install.packages("remotes")   # if needed
remotes::install_github("keyuan/ccube")
library(ccube)


base_dir <- "~/HMF/version1/somatic/DR-252-update1"  # <-- change
out_dir  <- file.path("~/clonality/Ccube")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#!/usr/bin/env Rscript
# run_ccube_from_purple.R
#
# Ccube SNV CCF estimation using HMF/Purple outputs.
# - Reads tumour+normal somatic VCF (selects tumour sample)
# - Extracts ref/alt counts robustly from many AD formats
# - Maps SNVs to Purple CN segments (major/minor/total CN)
# - Runs ccube::RunCcubePipeline per patient
# - Writes per-patient and combined outputs
#
# Requirements:
#   Bioconductor: VariantAnnotation, GenomicRanges, IRanges
#   CRAN: dplyr, stringr, purrr, readr, tidyr
#   ccube: GitHub keyuan/ccube (installed automatically if missing)

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(tidyr)
  library(GenomicRanges)
  library(IRanges)
  library(VariantAnnotation)
})



# ccube settings
numOfClusterPool <- 1:6
numOfRepeat <- 1
set.seed(123)


# ----------------------------
# Small utilities
# ----------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

find_one <- function(dir, patterns) {
  for (pat in patterns) {
    hits <- list.files(dir, pattern = pat, full.names = TRUE)
    if (length(hits) > 0) return(hits[1])
  }
  NA_character_
}

# ----------------------------
# Purple purity
# ----------------------------
read_purple_purity <- function(purple_dir) {
  cand <- list.files(purple_dir, full.names = TRUE)
  pri <- cand[str_detect(basename(cand), regex("purity|purple\\.purity|tumou?r\\.purity|purity\\.tsv|purity\\.csv", ignore_case = TRUE))]
  cand2 <- c(pri, cand)
  
  for (f in cand2) {
    if (file.info(f)$isdir) next
    
    # try TSV then CSV
    x <- tryCatch(suppressWarnings(readr::read_delim(f, delim = "\t", n_max = 5, show_col_types = FALSE)),
                  error = function(e) NULL)
    if (is.null(x) || ncol(x) == 0) {
      x <- tryCatch(suppressWarnings(readr::read_csv(f, n_max = 5, show_col_types = FALSE)),
                    error = function(e) NULL)
    }
    if (is.null(x) || ncol(x) == 0) next
    
    pcols <- which(str_detect(names(x), regex("purity", ignore_case = TRUE)))
    if (length(pcols) == 0) next
    
    vals <- suppressWarnings(as.numeric(unlist(x[1, pcols], use.names = FALSE)))
    vals <- vals[is.finite(vals)]
    if (length(vals) > 0 && vals[1] > 0 && vals[1] < 1) return(vals[1])
  }
  
  stop("Could not find a purity value in: ", purple_dir)
}

# ----------------------------
# Purple CN segments
# ----------------------------
read_purple_segments <- function(seg_file) {
  df <- readr::read_tsv(seg_file, show_col_types = FALSE) %>%
    transmute(
      chr = gsub("^chr", "", as.character(.data$chromosome %||% .data$chr %||% .data$Chromosome)),
      start = as.integer(.data$start %||% .data$Start),
      end   = as.integer(.data$end %||% .data$End),
      total_cn = suppressWarnings(as.numeric(.data$copyNumber %||% .data$totalCN %||% .data$total_cn)),
      major_cn = suppressWarnings(as.numeric(.data$majorAlleleCopyNumber %||% .data$majorCN %||% .data$major_cn)),
      minor_cn = suppressWarnings(as.numeric(.data$minorAlleleCopyNumber %||% .data$minorCN %||% .data$minor_cn))
    ) %>%
    filter(is.finite(start), is.finite(end), start < end)
  
  df <- df %>%
    mutate(total_cn = if_else(is.finite(total_cn), total_cn,
                              if_else(is.finite(major_cn) & is.finite(minor_cn), major_cn + minor_cn, NA_real_))) %>%
    filter(is.finite(total_cn), is.finite(major_cn), is.finite(minor_cn))
  
  GRanges(
    seqnames = df$chr,
    ranges   = IRanges(df$start, df$end),
    major_cn = df$major_cn,
    minor_cn = df$minor_cn,
    total_cn = df$total_cn
  )
}

# ----------------------------
# VCF helpers: select tumour sample
# ----------------------------
choose_tumour_index <- function(sample_names, tumour_sample = NULL) {
  if (!is.null(tumour_sample)) {
    t_idx <- match(tumour_sample, sample_names)
    if (!is.finite(t_idx)) stop("Requested tumour_sample not in VCF. Available: ", paste(sample_names, collapse = ", "))
    return(t_idx)
  }
  if (length(sample_names) == 1) return(1)
  
  # HMF/Purple typical: <pid>R (normal) and <pid>T (tumour)
  t_idx <- which(grepl("T$", sample_names))
  if (length(t_idx) != 1) t_idx <- which(grepl("tumou?r|tumor", sample_names, ignore.case = TRUE))
  if (length(t_idx) != 1) {
    not_r <- which(!grepl("R$", sample_names))
    if (length(not_r) == 1) t_idx <- not_r
  }
  if (length(t_idx) != 1) {
    stop("Cannot uniquely determine tumour sample. Samples: ", paste(sample_names, collapse = ", "),
         ". Provide tumour_sample= explicitly.")
  }
  t_idx
}

# ----------------------------
# Robust AD extraction (supports many VariantAnnotation representations)
# ----------------------------
extract_AD_counts <- function(AD, t_idx, sample_names) {
  ns <- length(sample_names)
  
  # 3D array [variant, sample, allele]
  if (!is.null(dim(AD)) && length(dim(AD)) == 3) {
    ref_counts <- as.integer(AD[, t_idx, 1])
    var_counts <- as.integer(AD[, t_idx, 2])
    return(list(ref = ref_counts, var = var_counts))
  }
  
  if (!is.matrix(AD)) stop("AD is not matrix/array. class=", paste(class(AD), collapse = ","))
  
  # character matrix with "ref,alt"
  if (is.character(AD)) {
    spl <- stringr::str_split(AD[, t_idx], ",", simplify = TRUE)
    ref_counts <- suppressWarnings(as.integer(spl[, 1]))
    var_counts <- suppressWarnings(as.integer(spl[, 2]))
    return(list(ref = ref_counts, var = var_counts))
  }
  
  # numeric/integer matrix
  if (is.numeric(AD) || is.integer(AD)) {
    # Common: AD is 2 columns (ref,alt) even if VCF has 2 samples
    if (ncol(AD) == 2) {
      ref_counts <- as.integer(AD[, 1])
      var_counts <- as.integer(AD[, 2])
      return(list(ref = ref_counts, var = var_counts))
    }
    # Another common: 2 columns per sample (ref,alt) => 2*ns columns
    if (ncol(AD) == 2 * ns) {
      ref_col <- 2 * (t_idx - 1) + 1
      alt_col <- 2 * (t_idx - 1) + 2
      ref_counts <- as.integer(AD[, ref_col])
      var_counts <- as.integer(AD[, alt_col])
      return(list(ref = ref_counts, var = var_counts))
    }
    stop("AD numeric matrix format not recognized. ncol(AD)=", ncol(AD), " nsamples=", ns)
  }
  
  # list-mode matrix (each cell is scalar/vec)
  if (typeof(AD) == "list") {
    if (ncol(AD) == 2) {
      ref_counts <- as.integer(vapply(AD[, 1], function(x) as.integer(x)[1], integer(1)))
      var_counts <- as.integer(vapply(AD[, 2], function(x) as.integer(x)[1], integer(1)))
      return(list(ref = ref_counts, var = var_counts))
    }
    if (ncol(AD) == 2 * ns) {
      ref_col <- 2 * (t_idx - 1) + 1
      alt_col <- 2 * (t_idx - 1) + 2
      ref_counts <- as.integer(vapply(AD[, ref_col], function(x) as.integer(x)[1], integer(1)))
      var_counts <- as.integer(vapply(AD[, alt_col], function(x) as.integer(x)[1], integer(1)))
      return(list(ref = ref_counts, var = var_counts))
    }
    stop("AD list-matrix format not recognized. ncol(AD)=", ncol(AD), " nsamples=", ns)
  }
  
  stop("AD matrix format not recognized (neither character, numeric/integer, nor list). class=",
       paste(class(AD), collapse = ","), " typeof=", typeof(AD))
}

# ----------------------------
# Read SNVs from somatic VCF and extract counts (tumour)
# ----------------------------
read_snv_counts_from_vcf <- function(vcf_file, tumour_sample = NULL) {
  vcf <- VariantAnnotation::readVcf(vcf_file)
  
  rr  <- rowRanges(vcf)
  chr <- gsub("^chr", "", as.character(GenomicRanges::seqnames(rr)))
  pos <- GenomicRanges::start(rr)
  ref <- as.character(VariantAnnotation::ref(vcf))
  alt <- vapply(VariantAnnotation::alt(vcf), function(x) as.character(x[1]), character(1))
  
  gt <- VariantAnnotation::geno(vcf)
  
  sample_names <- NULL
  if (!is.null(gt$GT)) sample_names <- colnames(gt$GT)
  if (is.null(sample_names)) {
    # fall back: try any geno field
    any_field <- gt[[1]]
    if (!is.null(any_field)) sample_names <- colnames(any_field)
  }
  if (is.null(sample_names) || length(sample_names) < 1) {
    stop("Could not determine sample names from VCF geno fields: ", vcf_file)
  }
  
  t_idx <- choose_tumour_index(sample_names, tumour_sample)
  
  n <- length(pos)
  ref_counts <- rep(NA_integer_, n)
  var_counts <- rep(NA_integer_, n)
  total_counts <- rep(NA_integer_, n)
  
  if (!is.null(gt$AD)) {
    out <- extract_AD_counts(gt$AD, t_idx, sample_names)
    ref_counts <- out$ref
    var_counts <- out$var
    total_counts <- ref_counts + var_counts
  } else {
    # Fallback: DP + AF (less ideal)
    if (is.null(gt$DP)) stop("VCF lacks AD and DP; cannot compute read counts: ", vcf_file)
    total_counts <- suppressWarnings(as.integer(gt$DP[, t_idx]))
    
    af <- NULL
    if (!is.null(gt$AF)) af <- suppressWarnings(as.numeric(gt$AF[, t_idx]))
    if (is.null(af)) stop("VCF lacks AD and AF; cannot infer var_counts reliably: ", vcf_file)
    
    var_counts <- as.integer(round(total_counts * af))
    ref_counts <- total_counts - var_counts
  }
  
  tibble(
    chr = chr,
    pos = as.integer(pos),
    ref = ref,
    alt = alt,
    mutation_id = paste(chr, pos, ref, alt, sep = "_"),
    ref_counts = as.integer(ref_counts),
    var_counts = as.integer(var_counts),
    total_counts = as.integer(total_counts)
  ) %>%
    filter(is.finite(ref_counts), is.finite(var_counts), is.finite(total_counts), total_counts > 0)
}

# ----------------------------
# Map SNVs to CN segments
# ----------------------------
annotate_with_cn <- function(snv_df, seg_gr) {
  snv_gr <- GRanges(seqnames = snv_df$chr, ranges = IRanges(snv_df$pos, snv_df$pos))
  hits <- findOverlaps(snv_gr, seg_gr, select = "first")
  
  out <- snv_df %>%
    mutate(
      major_cn = mcols(seg_gr)$major_cn[hits],
      minor_cn = mcols(seg_gr)$minor_cn[hits],
      total_cn = mcols(seg_gr)$total_cn[hits],
      normal_cn = 2
    ) %>%
    filter(is.finite(major_cn), is.finite(minor_cn), is.finite(total_cn))
  
  out
}

# ----------------------------
# Run Ccube for one patient
# ----------------------------
run_ccube_one_patient <- function(patient_dir) {
  pid <- basename(patient_dir)
  
  purple_dir <- file.path(patient_dir, "purple")
  if (!dir.exists(purple_dir)) return(NULL)
  
  # segment file patterns (adjust if needed)
  seg_file <- find_one(
    purple_dir,
    patterns = c(
      "\\.purple\\.cnv\\.somatic\\.tsv$",
      "\\.purple\\.cnv\\.tsv$",
      "cnv\\.tsv$",
      "somatic\\.cnv\\.tsv$",
      "segments\\.tsv$"
    )
  )
  if (is.na(seg_file)) stop("No CN segment file found in: ", purple_dir)
  
  # somatic VCF patterns (adjust if needed)
  vcf_file <- find_one(
    purple_dir,
    patterns = c(
      "\\.purple\\.somatic\\.vcf\\.gz$",
      "\\.somatic\\.vcf\\.gz$",
      "somatic\\.vcf\\.gz$",
      "\\.vcf\\.gz$"
    )
  )
  if (is.na(vcf_file)) stop("No somatic VCF (.vcf.gz) found in: ", purple_dir)
  
  purity <- read_purple_purity(purple_dir)
  seg_gr <- read_purple_segments(seg_file)
  
  # NOTE: tumour_sample can be forced if needed, e.g. paste0(pid, "T")
  snv <- read_snv_counts_from_vcf(vcf_file, tumour_sample = NULL) %>%
    annotate_with_cn(seg_gr)
  
  if (nrow(snv) < 50) {
    warning("Too few SNVs after filtering for ", pid, " (", nrow(snv), "). Skipping.")
    return(NULL)
  }
  
  ssm <- snv %>%
    transmute(
      mutation_id = mutation_id,
      ref_counts  = as.integer(ref_counts),
      var_counts  = as.integer(var_counts),
      total_counts = as.integer(total_counts),
      minor_cn    = as.numeric(minor_cn),
      major_cn    = as.numeric(major_cn),
      total_cn    = as.numeric(total_cn),
      purity      = as.numeric(purity),
      normal_cn   = 2
    )
  
  res <- ccube::RunCcubePipeline(
    ssm = ssm,
    numOfClusterPool = numOfClusterPool,
    numOfRepeat = numOfRepeat,
    runAnalysis = TRUE,
    runQC = TRUE
  )
  
  # Try to keep the best converged model if Ccube provides QC/model selection info
  if (!is.null(res$QC)) {
    qc <- res$QC
    if ("converged" %in% names(qc)) {
      if (any(qc$converged)) {
        best <- qc %>% filter(converged) %>% arrange(BIC %||% AIC %||% nLL) %>% slice(1)
        message("Selected converged model: ", paste(names(best), best, sep="=", collapse=", "))
      } else {
        message("WARNING: no converged models reported by QC; keeping default output.")
      }
    }
  }
  
  
  out <- res$ssm %>%
    mutate(patient_id = pid) %>%
    dplyr::select(patient_id, mutation_id, ccube_ccf, ccube_mult, everything())
  
  out_file <- file.path(out_dir, paste0(pid, ".ccube.tsv.gz"))
  readr::write_tsv(out, out_file)
  message("Wrote: ", out_file, " (", nrow(out), " SNVs)")
  
  out
}

# ----------------------------
# Run all patients
# ----------------------------
patient_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
patient_dirs <- patient_dirs[dir.exists(file.path(patient_dirs, "purple"))]
message("Found ", length(patient_dirs), " patient dirs with purple/")

all_ccube <- purrr::map_dfr(
  patient_dirs,
  ~tryCatch(run_ccube_one_patient(.x),
            error = function(e) {
              message("FAILED ", basename(.x), ": ", e$message)
              NULL
            })
)

out_all <- file.path(out_dir, "ccube_all_patients.tsv.gz")
readr::write_tsv(all_ccube, out_all)
message("DONE. Combined output: ", out_all)

# ----------------------------
# Optional: quick summary
# ----------------------------
if (nrow(all_ccube) > 0) {
  message("Summary ccube_ccf: ")
  print(summary(all_ccube$ccube_ccf))
}

