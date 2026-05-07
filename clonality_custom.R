#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(GenomicRanges)
  library(IRanges)
})

# ---------------------------
# USER SETTINGS
# ---------------------------
base_dir     <- "~/HMF/somatic/DR-252-update1"
gene_bed_tsv <- "~/refs/gene_coords_hg19.tsv"  # hg19: gene, chr, start, end
out_tsv      <- "drivers_clonality_all.tsv"

tol_integer <- 0.15

# ---------------------------
# Helpers
# ---------------------------
find_one <- function(dir, patterns, recursive = TRUE) {
  for (pat in patterns) {
    f <- list.files(dir, pattern = pat, full.names = TRUE, recursive = recursive, ignore.case = TRUE)
    if (length(f) > 0) return(f[[1]])
  }
  NA_character_
}

read_delim_safe <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  # PURPLE sometimes uses .tsv, sometimes .csv; auto-detect
  tryCatch({
    if (grepl("\\.csv$", path, ignore.case = TRUE)) {
      readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
    } else {
      readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
    }
  }, error = function(e) NULL)
}

weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

# ---------------------------
# PURPLE purity/ploidy
# ---------------------------
read_purple_purity <- function(purple_dir) {
  f <- find_one(
    purple_dir,
    patterns = c("\\.purple\\.purity\\.(csv|tsv)$", "purity\\.(csv|tsv)$"),
    recursive = TRUE
  )
  df <- read_delim_safe(f)
  if (is.null(df) || nrow(df) == 0) return(tibble(purity=NA_real_, ploidy=NA_real_, wgd=NA, minPurity=NA_real_, maxPurity=NA_real_))
  
  tibble(
    purity = suppressWarnings(as.numeric(df$purity[1])),
    ploidy = suppressWarnings(as.numeric(df$ploidy[1])),
    wgd    = df$wholeGenomeDuplication[1],
    minPurity = suppressWarnings(as.numeric(df$minPurity[1])),
    maxPurity = suppressWarnings(as.numeric(df$maxPurity[1]))
  )
}

# ---------------------------
# Driver catalog (your schema)
# ---------------------------
read_driver_catalog <- function(driver_file) {
  drv <- read_delim_safe(driver_file)
  if (is.null(drv) || nrow(drv) == 0) return(NULL)
  
  needed <- c("chromosome","gene","driver","category","driverLikelihood",
              "missense","nonsense","splice","inframe","frameshift",
              "biallelic","minCopyNumber","maxCopyNumber")
  for (x in needed) if (!x %in% names(drv)) drv[[x]] <- NA
  
  drv %>%
    transmute(
      chr = str_replace(as.character(chromosome), "^chr", ""),
      gene = as.character(gene),
      driver = as.logical(driver),
      category = as.character(category),
      driverLikelihood = suppressWarnings(as.numeric(driverLikelihood)),
      missense = suppressWarnings(as.numeric(missense)),
      nonsense = suppressWarnings(as.numeric(nonsense)),
      splice = suppressWarnings(as.numeric(splice)),
      inframe = suppressWarnings(as.numeric(inframe)),
      frameshift = suppressWarnings(as.numeric(frameshift)),
      biallelic = suppressWarnings(as.numeric(biallelic)),
      minCopyNumber = suppressWarnings(as.numeric(minCopyNumber)),
      maxCopyNumber = suppressWarnings(as.numeric(maxCopyNumber))
    )
}

get_driver_gene_sets <- function(drv) {
  snv_flag <- (coalesce(drv$missense,0) + coalesce(drv$nonsense,0) +
                 coalesce(drv$splice,0) + coalesce(drv$inframe,0) +
                 coalesce(drv$frameshift,0)) > 0
  
  cnv_flag <- is.finite(drv$minCopyNumber) | is.finite(drv$maxCopyNumber) |
    str_detect(coalesce(drv$category,""), regex("amp|del|loss|gain|cna|cnv", ignore_case=TRUE))
  
  list(
    snv_genes = unique(na.omit(drv$gene[snv_flag])),
    cnv_genes = unique(na.omit(drv$gene[cnv_flag]))
  )
}

# ---------------------------
# VCF parsing (gene + PURPLE fields)
# ---------------------------
parse_csq_gene <- function(vcf) {
  csq <- info(vcf)$CSQ
  n <- length(csq)
  gene <- rep(NA_character_, n)
  if (is.null(csq)) return(gene)
  
  desc <- tryCatch(info(header(vcf))["CSQ","Description"], error=function(e) NA_character_)
  if (!is.character(desc) || is.na(desc)) return(gene)
  
  fmt <- sub(".*Format: ", "", desc)
  fields <- strsplit(fmt, "\\|")[[1]]
  
  idx <- match("SYMBOL", fields)
  if (is.na(idx)) idx <- match("Gene", fields)
  if (is.na(idx)) idx <- match("HGNC", fields)
  if (is.na(idx)) return(gene)
  
  for (i in seq_len(n)) {
    if (is.na(csq[i]) || !nzchar(csq[i])) next
    first_ann <- strsplit(csq[i], ",", fixed = TRUE)[[1]][1]
    parts <- strsplit(first_ann, "\\|")[[1]]
    if (length(parts) >= idx && nzchar(parts[[idx]])) gene[i] <- parts[[idx]]
  }
  gene
}

get_info_numeric <- function(vcf, candidates) {
  for (k in candidates) if (k %in% names(info(vcf))) return(suppressWarnings(as.numeric(info(vcf)[[k]])))
  rep(NA_real_, nrow(vcf))
}

read_purple_vcf_variants <- function(vcf_file, genome = "hg19") {
  vcf <- readVcf(vcf_file, genome = genome)
  rr  <- rowRanges(vcf)
  
  tibble(
    chr = str_replace(as.character(seqnames(rr)), "^chr", ""),
    pos = start(rr),
    ref = as.character(ref(vcf)),
    alt = vapply(alt(vcf), function(a) as.character(a)[1], character(1)),
    gene = parse_csq_gene(vcf),
    purple_af  = get_info_numeric(vcf, c("PURPLE_AF","purple_af","AF","VAF")),
    purple_cn  = get_info_numeric(vcf, c("PURPLE_CN","purple_cn","CN","copyNumber")),
    purple_vcn = get_info_numeric(vcf, c("PURPLE_VCN","purple_vcn","VCN","variantCopyNumber","varCopyNumber"))
  )
}

# ---------------------------
# SNV CCF: prefer PURPLE_VCN; else purity-aware fallback from VAF+CN+purity
# ---------------------------
ccf_from_purple <- function(purple_vcn, purple_cn, purple_af, purity, normal_cn = 2) {
  # Returns list(ccf, m_hat, method, implied_vcn)
  cn_ok <- is.finite(purple_cn) && purple_cn > 0
  
  # 1) Primary: PURPLE VCN (already purity/CN modelled)
  if (is.finite(purple_vcn) && cn_ok) {
    m_hat <- round(purple_vcn)
    m_hat <- pmax(1, pmin(m_hat, ceiling(purple_cn)))
    ccf <- purple_vcn / m_hat
    ccf <- pmin(1, pmax(0, ccf))
    return(list(ccf=ccf, m_hat=m_hat, method="vcn", implied_vcn=purple_vcn))
  }
  
  # 2) Fallback: implied VCN from VAF + CN + purity
  p_ok  <- is.finite(purity) && purity > 0 && purity < 1
  af_ok <- is.finite(purple_af) && purple_af >= 0 && purple_af <= 1
  
  if (p_ok && af_ok && cn_ok) {
    implied_vcn <- purple_af * (purity * purple_cn + (1 - purity) * normal_cn) / purity
    
    # multiplicity guess: choose integer 1..ceil(CN) that yields ccf<=1 and is closest to 1
    ms <- seq_len(max(1L, ceiling(purple_cn)))
    ccf_cands <- implied_vcn / ms
    ccf_cands[ccf_cands <= 0] <- NA_real_
    ccf_cands[ccf_cands > 1.2] <- NA_real_  # discard very implausible
    # pick ccf closest to 1 (clonal) if any; else use m=1
    if (all(!is.finite(ccf_cands))) {
      m_hat <- 1
      ccf <- pmin(1, pmax(0, implied_vcn))
    } else {
      m_hat <- ms[which.min(abs(ccf_cands - 1))]
      ccf <- ccf_cands[which.min(abs(ccf_cands - 1))]
      ccf <- pmin(1, pmax(0, ccf))
    }
    return(list(ccf=ccf, m_hat=m_hat, method="vaf_purity_cn", implied_vcn=implied_vcn))
  }
  
  list(ccf=NA_real_, m_hat=NA_real_, method="missing_fields", implied_vcn=NA_real_)
}

# ---------------------------
# CNV clonality: ploidy-aware baseline + best integer allele state search
# ---------------------------
baseline_alleles_from_ploidy <- function(ploidy) {
  # crude but better than hard-coded 1/1 when WGD/aneuploid
  if (!is.finite(ploidy) || ploidy <= 0) return(c(1L, 1L))
  t <- max(2L, as.integer(round(ploidy)))
  # balanced split baseline
  minor <- floor(t / 2)
  major <- t - minor
  c(major, minor)
}

best_integer_state <- function(obs_major, obs_minor, obs_baf = NA_real_) {
  # choose (M,m) integer state that best matches observed allele CN and (optionally) baf
  if (!is.finite(obs_major) || !is.finite(obs_minor)) return(list(M=NA_integer_, m=NA_integer_, score=Inf))
  
  t_hat <- max(0L, as.integer(round(obs_major + obs_minor)))
  # allow +-1 total CN flexibility
  totals <- unique(pmax(0L, c(t_hat-1L, t_hat, t_hat+1L)))
  
  cand <- expand.grid(total = totals, m = 0:max(totals))
  cand <- cand %>%
    mutate(M = total - m) %>%
    filter(M >= m, M >= 0, m >= 0)
  
  if (nrow(cand) == 0) return(list(M=NA_integer_, m=NA_integer_, score=Inf))
  
  cand <- cand %>%
    mutate(
      score_cn = abs(obs_major - M) + abs(obs_minor - m),
      baf_exp = ifelse((M + m) > 0, m / (M + m), 0.5),
      score_baf = ifelse(is.finite(obs_baf), abs(obs_baf - baf_exp), 0),
      score = score_cn + 0.5 * score_baf
    ) %>%
    arrange(score)
  
  list(M = cand$M[1], m = cand$m[1], score = cand$score[1])
}

cnv_fraction_from_alleles <- function(obs_major, obs_minor, base_major, base_minor, ev_major, ev_minor) {
  sols <- c()
  if (is.finite(obs_major) && (ev_major - base_major) != 0) sols <- c(sols, (obs_major - base_major) / (ev_major - base_major))
  if (is.finite(obs_minor) && (ev_minor - base_minor) != 0) sols <- c(sols, (obs_minor - base_minor) / (ev_minor - base_minor))
  if (length(sols) == 0) return(list(f=NA_real_, method="unsolved"))
  f <- mean(sols, na.rm = TRUE)
  f <- max(0, min(1, f))
  list(f=f, method="mixture_2state")
}

# ---------------------------
# Gene coordinates (hg19)
# ---------------------------
genes_df <- read_tsv(gene_bed_tsv, show_col_types = FALSE) %>%
  transmute(
    gene = as.character(gene),
    chr  = str_replace(as.character(chr), "^chr", ""),
    start = as.integer(start),
    end   = as.integer(end)
  ) %>%
  filter(!is.na(gene), !is.na(chr), is.finite(start), is.finite(end), end >= start)

gene_gr_all <- GRanges(
  seqnames = genes_df$chr,
  ranges   = IRanges(start = genes_df$start, end = genes_df$end),
  gene     = genes_df$gene
)

# ---------------------------
# Per-sample processing
# ---------------------------
process_one <- function(pid_dir) {
  pid <- basename(pid_dir)
  purple_dir <- file.path(pid_dir, "purple")
  
  # purity/ploidy
  pp <- read_purple_purity(purple_dir)
  purity <- pp$purity
  ploidy <- pp$ploidy
  wgd    <- pp$wgd
  
  base <- baseline_alleles_from_ploidy(ploidy)
  base_major <- base[1]; base_minor <- base[2]
  
  # driver catalog
  driver_file <- find_one(purple_dir, patterns = c("driver.*catalog.*\\.(tsv|csv)$", ".*driver.*\\.(tsv|csv)$"))
  drv <- read_driver_catalog(driver_file)
  if (is.null(drv)) return(NULL)
  gs <- get_driver_gene_sets(drv)
  
  # --------------- SNV/INDEL driver CCF from VCF ---------------
  vcf_file <- find_one(purple_dir, patterns = c("somatic.*\\.vcf\\.gz$", "somatic.*\\.vcf$", "\\.vcf\\.gz$", "\\.vcf$"))
  snv_out <- NULL
  if (length(gs$snv_genes) > 0 && !is.na(vcf_file)) {
    v <- tryCatch(read_purple_vcf_variants(vcf_file, genome = "hg19"), error = function(e) NULL)
    
    if (!is.null(v) && nrow(v) > 0) {
      snv_out <- v %>%
        filter(!is.na(gene), gene %in% gs$snv_genes) %>%
        rowwise() %>%
        mutate(
          patient_id = pid,
          event_type = "SNV",
          purity = purity,
          ploidy = ploidy,
          wgd = wgd,
          tmp = list(ccf_from_purple(purple_vcn, purple_cn, purple_af, purity)),
          clonality = tmp$ccf,
          multiplicity_hat = tmp$m_hat,
          implied_vcn = tmp$implied_vcn,
          clonality_method = tmp$method,
          qc_low_purity = is.finite(purity) && purity < 0.2,
          clonality_simple = case_when(
            is.na(clonality) ~ NA_character_,
            clonality >= 0.9 ~ "clonal_like",
            TRUE ~ "subclonal_like"
          )
        ) %>%
        ungroup() %>%
        select(patient_id, event_type, gene, chr, pos, ref, alt,
               purity, ploidy, wgd,
               purple_af, purple_cn, purple_vcn, implied_vcn,
               multiplicity_hat, clonality, clonality_method, qc_low_purity, clonality_simple)
    }
  }
  
  # --------------- CNV driver clonality from CNV segments ---------------
  cnv_file <- find_one(
    purple_dir,
    patterns = c("\\.purple\\.cnv\\.somatic\\.(tsv|csv)$",
                 "\\.purple\\.cnv\\.(tsv|csv)$",
                 "cnv.*\\.(tsv|csv)$")
  )
  
  cnv_out <- NULL
  if (length(gs$cnv_genes) > 0 && !is.na(cnv_file)) {
    cnv <- read_delim_safe(cnv_file)
    if (!is.null(cnv) && nrow(cnv) > 0) {
      
      gene_gr <- gene_gr_all[mcols(gene_gr_all)$gene %in% gs$cnv_genes]
      
      cnv2 <- cnv %>%
        transmute(
          chr   = str_replace(as.character(chromosome), "^chr", ""),
          start = as.integer(start),
          end   = as.integer(end),
          total_cn = suppressWarnings(as.numeric(copyNumber)),
          baf      = suppressWarnings(as.numeric(baf)),
          observedBAF = suppressWarnings(as.numeric(observedBAF)),
          bafCount = suppressWarnings(as.numeric(bafCount)),
          depthWindowCount = suppressWarnings(as.numeric(depthWindowCount)),
          major_cn = suppressWarnings(as.numeric(majorAlleleCopyNumber)),
          minor_cn = suppressWarnings(as.numeric(minorAlleleCopyNumber))
        ) %>%
        filter(!is.na(chr), is.finite(start), is.finite(end), end >= start)
      
      if (nrow(cnv2) > 0 && length(gene_gr) > 0) {
        
        seg_gr <- GRanges(
          seqnames = cnv2$chr,
          ranges   = IRanges(start = cnv2$start, end = cnv2$end),
          total_cn = cnv2$total_cn,
          major_cn = cnv2$major_cn,
          minor_cn = cnv2$minor_cn,
          baf      = cnv2$baf,
          observedBAF = cnv2$observedBAF,
          bafCount = cnv2$bafCount,
          depthWindowCount = cnv2$depthWindowCount
        )
        
        hits <- findOverlaps(gene_gr, seg_gr, ignore.strand = TRUE)
        if (length(hits) > 0) {
          q <- queryHits(hits); s <- subjectHits(hits)
          ov <- pintersect(ranges(gene_gr)[q], ranges(seg_gr)[s])
          wlen <- width(ov)
          
          # confidence weight: overlap length × support (softly)
          bafCount <- mcols(seg_gr)$bafCount[s]
          depthCnt <- mcols(seg_gr)$depthWindowCount[s]
          conf <- sqrt(pmax(1, bafCount)) * sqrt(pmax(1, depthCnt))
          w <- wlen * conf
          
          df <- tibble(
            gene = mcols(gene_gr)$gene[q],
            w = w,
            wlen = wlen,
            total_cn = mcols(seg_gr)$total_cn[s],
            major_cn = mcols(seg_gr)$major_cn[s],
            minor_cn = mcols(seg_gr)$minor_cn[s],
            baf = mcols(seg_gr)$baf[s],
            observedBAF = mcols(seg_gr)$observedBAF[s]
          ) %>%
            group_by(gene) %>%
            summarise(
              total_cn = weighted_mean(total_cn, w),
              major_cn = weighted_mean(major_cn, w),
              minor_cn = weighted_mean(minor_cn, w),
              baf      = weighted_mean(baf, w),
              observedBAF = weighted_mean(observedBAF, w),
              overlap_bp = sum(wlen, na.rm = TRUE),
              .groups = "drop"
            ) %>%
            rowwise() %>%
            mutate(
              patient_id = pid,
              event_type = "CNV",
              purity = purity,
              ploidy = ploidy,
              wgd = wgd,
              baseline_major = base_major,
              baseline_minor = base_minor,
              integerish = (abs(major_cn - round(major_cn)) <= tol_integer) &
                (abs(minor_cn - round(minor_cn)) <= tol_integer),
              
              best = list(best_integer_state(major_cn, minor_cn, obs_baf = baf)),
              event_major_int = best$M,
              event_minor_int = best$m,
              int_fit_score = best$score,
              
              tmp = list(cnv_fraction_from_alleles(major_cn, minor_cn,
                                                   base_major, base_minor,
                                                   event_major_int, event_minor_int)),
              clonality = tmp$f,
              clonality_method = tmp$method,
              
              clonality_simple = case_when(
                is.na(clonality) ~ NA_character_,
                clonality >= 0.9 ~ "clonal_like",
                clonality <= 0.1 ~ "absent_or_baseline",
                TRUE ~ "subclonal_like"
              ),
              state_simple = case_when(
                is.na(total_cn) ~ NA_character_,
                total_cn < 0.5 ~ "HomDel_like",
                total_cn < 1.5 ~ "HetLoss_like",
                total_cn < 2.5 ~ "Neutral_like",
                total_cn < 3.5 ~ "Gain_like",
                TRUE ~ "Amp_like"
              )
            ) %>%
            ungroup() %>%
            select(patient_id, event_type, gene,
                   purity, ploidy, wgd,
                   total_cn, major_cn, minor_cn, baf, observedBAF,
                   baseline_major, baseline_minor,
                   event_major_int, event_minor_int, int_fit_score,
                   overlap_bp, integerish,
                   clonality, clonality_method, clonality_simple, state_simple)
          
          cnv_out <- df
        }
      }
    }
  }
  
  bind_rows(snv_out, cnv_out)
}

# ---------------------------
# Run cohort
# ---------------------------
patient_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

all <- purrr::map(patient_dirs, function(d) {
  pid <- basename(d)
  message("▶ Processing ", pid, " ...")
  
  res <- tryCatch(
    process_one(d),
    error = function(e) {
      message("✖ ERROR in ", pid, ": ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(res) || nrow(res) == 0) {
    message("  ↳ no driver events found")
  } else {
    message("  ↳ ", nrow(res), " driver events")
  }
  
  res
}) %>%
  purrr::compact() %>%
  bind_rows()


if (nrow(all) == 0) stop("No results produced. Check file patterns and inputs.")

write_tsv(all, out_tsv)
message("Wrote: ", out_tsv, " (rows: ", nrow(all), ")")
print(all %>% count(event_type, sort = TRUE))



