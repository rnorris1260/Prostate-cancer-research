
library(dplyr)
library(readr)
library(vroom)
library(VariantAnnotation)
library(GenomicRanges)
library(IRanges)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(forcats)
library(purrr)

## import res_all_pc_spec - relationship between driver mutations (prostate cancer specific) and ssgsea

setwd("~/HMF/subclonality/")
res_all_pc_spec <- read.csv("res_all_pc_spec.csv")

## Driver importance (from driver vs ssGSEA results) ----
compute_driver_importance <- function(res_all_pc_spec,
                                      alpha   = 0.015,
                                      min_sig = 0) {
  res_all %>%
    group_by(driver) %>%
    summarise(
      n_sig       = sum(adj.P.Val <= alpha, na.rm = TRUE),
      max_log10p  = suppressWarnings(max(-log10(adj.P.Val[adj.P.Val <= alpha]),
                                         na.rm = TRUE)),
      mean_abs_fc = mean(abs(logFC[adj.P.Val <= alpha]), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      across(c(n_sig, max_log10p, mean_abs_fc),
             ~ ifelse(is.finite(.x), .x, 0)),
      importance = n_sig + max_log10p + mean_abs_fc
    ) %>%
    filter(n_sig >= min_sig) %>%
    arrange(desc(importance))
}

## CNV reader for PURPLE ----
read_purple_cnv <- function(cnv_file) {
  cnv <- vroom::vroom(cnv_file, show_col_types = FALSE)
  
  pick <- function(cands, required = TRUE) {
    hit <- intersect(cands, names(cnv))
    if (!length(hit)) {
      if (required) stop("Missing required CNV column: one of ", paste(cands, collapse = ", "))
      return(NA_character_)
    }
    hit[1]
  }
  
  chr_col   <- pick(c("chromosome","chr","chrom","seqnames"))
  start_col <- pick(c("start","segmentStart","segment_start","start_position","minRegionStart"))
  end_col   <- pick(c("end","segmentEnd","segment_end","end_position","minRegionEnd"))
  maj_col   <- pick(c("majorAlleleCopyNumber","majorCN"), required = FALSE)
  min_col   <- pick(c("minorAlleleCopyNumber","minorCN"), required = FALSE)
  tot_col   <- pick(c("copyNumber","tumorCopyNumber","totalCN","CNt","maxCopyNumber"), required = FALSE)
  
  out <- cnv %>%
    transmute(
      chr_raw = .data[[chr_col]],
      start   = as.numeric(.data[[start_col]]),
      end     = as.numeric(.data[[end_col]]),
      majorCN = if (!is.na(maj_col)) as.numeric(.data[[maj_col]]) else NA_real_,
      minorCN = if (!is.na(min_col)) as.numeric(.data[[min_col]]) else NA_real_,
      totalCN = dplyr::coalesce(
        if (!is.na(tot_col)) as.numeric(.data[[tot_col]]) else NA_real_,
        if (!is.na(maj_col) && !is.na(min_col))
          as.numeric(.data[[maj_col]]) + as.numeric(.data[[min_col]]) else NA_real_
      )
    ) %>%
    mutate(
      chr = gsub("^chr", "", as.character(chr_raw)),
      totalCN = ifelse(is.finite(totalCN) & totalCN > 0, totalCN, 2),
      majorCN = ifelse(is.finite(majorCN), majorCN, NA_real_),
      minorCN = ifelse(is.finite(minorCN), minorCN, NA_real_)
    ) %>%
    dplyr::select(chr, start, end, majorCN, minorCN, totalCN)
  
  out
}

## Purity ----
read_purity <- function(purity_file) {
  p <- vroom::vroom(purity_file, show_col_types = FALSE) %>%
    dplyr::select(dplyr::matches("purity", ignore.case = TRUE)) %>%
    dplyr::slice(1) %>%
    unlist(use.names = FALSE) %>%
    as.numeric()
  p <- p[is.finite(p)][1]
  if (!is.finite(p) || p <= 0 || p >= 1) stop("Purity not found or invalid in: ", purity_file)
  p
}

## VCF picker (handles .vcf / .vcf.gz) ----
pick_vcf_path <- function(purple_dir, patient_id) {
  candidates <- c(
    file.path(purple_dir, paste0(patient_id, ".purple.somatic.vcf")),
    file.path(purple_dir, paste0(patient_id, ".purple.somatic.vcf.gz"))
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    stop("No VCF found for ", patient_id,
         " (looked for .vcf and .vcf.gz in: ", purple_dir, ")")
  }
  hit[1]
}

## Simple time score (0 = early clonal, 1 = late/subclonal) ----
add_time_score <- function(timed_df) {
  stopifnot(all(c("P_EARLY","P_LATE","P_SUBCLONAL") %in% names(timed_df)))
  timed_df %>%
    mutate(
      time_score = 0.0 * P_EARLY +
        0.5 * P_LATE +
        1.0 * P_SUBCLONAL
    )
}

## check
#drv_importance <- compute_driver_importance(res_all_pc_spec)

read_vcf_snvs <- function(vcf_file, min_dp = 20, keep_all = FALSE) {
  library(VariantAnnotation)
  library(tibble)
  library(dplyr)
  
  vcf <- readVcf(vcf_file)
  
  gr  <- rowRanges(vcf)
  chr <- as.character(seqnames(gr))
  pos <- start(gr)
  
  n_var <- length(vcf)
  
  ## helper: safely unlist anything into a numeric vector
  safe_unlist_num <- function(x) {
    if (is.null(x)) return(numeric(0))
    as.numeric(as.vector(unlist(x)))
  }
  
  ## ---- DP (depth) ----
  dp <- rep(NA_integer_, n_var)
  
  if ("DP" %in% names(geno(vcf))) {
    dp_raw <- geno(vcf)$DP
    if (is.matrix(dp_raw) || is.array(dp_raw)) {
      # variants x samples; take first sample
      dp <- suppressWarnings(as.integer(dp_raw[, 1]))
    } else {
      dp_vec <- safe_unlist_num(dp_raw)
      if (length(dp_vec) == n_var) {
        dp <- suppressWarnings(as.integer(dp_vec))
      } else if (length(dp_vec) > 0) {
        warning("DP length (", length(dp_vec), ") != n_var (", n_var,
                "); leaving some DP as NA")
        dp[seq_len(min(length(dp_vec), n_var))] <- suppressWarnings(as.integer(dp_vec[seq_len(min(length(dp_vec), n_var))]))
      }
    }
  } else if ("DP" %in% names(info(vcf))) {
    dp_raw <- info(vcf)$DP
    dp_vec <- safe_unlist_num(dp_raw)
    if (length(dp_vec) == n_var) {
      dp <- suppressWarnings(as.integer(dp_vec))
    } else if (length(dp_vec) > 0) {
      warning("INFO/DP length (", length(dp_vec), ") != n_var (", n_var,
              "); leaving some DP as NA")
      dp[seq_len(min(length(dp_vec), n_var))] <- suppressWarnings(as.integer(dp_vec[seq_len(min(length(dp_vec), n_var))]))
    }
  }
  
  ## ---- VAF from PURPLE_AF only ----
  vaf <- rep(NA_real_, n_var)
  
  if ("PURPLE_AF" %in% names(info(vcf))) {
    purple_raw <- info(vcf)$PURPLE_AF
    purple_vec <- safe_unlist_num(purple_raw)
    if (length(purple_vec) == n_var) {
      vaf <- purple_vec
    } else if (length(purple_vec) > 0) {
      warning("PURPLE_AF length (", length(purple_vec), ") != n_var (", n_var,
              "); filling up to available length")
      vaf[seq_len(min(length(purple_vec), n_var))] <- purple_vec[seq_len(min(length(purple_vec), n_var))]
    }
  } else {
    warning("No PURPLE_AF in INFO; VAF will be NA")
  }
  
  ## ---- Gene symbol from CSQ (VEP) ----
  gene <- rep(NA_character_, n_var)
  info_hdr <- info(header(vcf))
  
  if ("CSQ" %in% rownames(info_hdr)) {
    csq_desc <- info_hdr["CSQ", "Description"]
    fmt_str  <- sub(".*[Ff]ormat: *", "", csq_desc)
    csq_fields <- strsplit(fmt_str, "\\|")[[1]]
    csq_fields <- trimws(csq_fields)
    
    sym_idx <- which(csq_fields %in% c("SYMBOL","Gene","GENE","HGNC_ID"))[1]
    
    if (!is.na(sym_idx)) {
      csq <- info(vcf)$CSQ
      csq <- as.character(csq)
      csq_first <- vapply(strsplit(csq, ","), `[`, character(1), 1)
      
      gene <- vapply(strsplit(csq_first, "\\|"), function(x) {
        if (length(x) >= sym_idx) x[sym_idx] else NA_character_
      }, character(1))
    } else {
      warning("Could not find SYMBOL/Gene field in CSQ description; gene will be NA")
    }
  } else {
    warning("No CSQ field in INFO; gene will be NA")
  }
  
  ## ---- Build tibble ----
  snvs_raw <- tibble(
    chr = gsub("^chr", "", chr),
    pos = pos,
    ref = as.character(ref(vcf)),
    alt = vapply(alt(vcf), function(a) as.character(a)[1], character(1)),
    dp  = dp,
    vaf = vaf,
    gene = gene
  )
  
  passed <- with(snvs_raw,
                 !is.na(dp) & dp >= min_dp &
                   !is.na(vaf) & vaf > 0)
  
  snvs_raw$passed_filters <- passed
  
  if (keep_all) {
    snvs_raw
  } else {
    dplyr::filter(snvs_raw, passed_filters)
  }
}


vcf_test <- readVcf(vcf_file_test)

str(geno(vcf_test)$DP)
str(info(vcf_test)$PURPLE_AF)




map_snvs_to_cnv <- function(snvs, cnv_df) {
  if (!nrow(snvs) || !nrow(cnv_df)) return(tibble())
  
  gr_snvs <- GRanges(
    seqnames = snvs$chr,
    ranges   = IRanges(snvs$pos, snvs$pos)
  )
  gr_cnv <- GRanges(
    seqnames = cnv_df$chr,
    ranges   = IRanges(cnv_df$start, cnv_df$end)
  )
  
  hits <- findOverlaps(gr_snvs, gr_cnv)
  if (!length(hits)) {
    warning("No overlaps between SNVs and CNV segments")
    return(tibble())
  }
  
  sn_idx  <- queryHits(hits)
  cnv_idx <- subjectHits(hits)
  
  snvs_mapped <- snvs[sn_idx, , drop = FALSE]
  cnv_mapped  <- cnv_df[cnv_idx, , drop = FALSE]
  
  snvs_mapped %>%
    mutate(
      seg_id  = cnv_idx,
      majorCN = cnv_mapped$majorCN,
      minorCN = cnv_mapped$minorCN,
      totalCN = cnv_mapped$totalCN
    )
}


probabilistic_timing <- function(mapped,
                                 purity,
                                 jitter    = c(-0.05,-0.03,0,0.03,0.05),
                                 gain_rule = c("major_ge2","total_ge3")) {
  gain_rule <- match.arg(gain_rule)
  
  if (!nrow(mapped)) return(tibble())
  
  mapped <- mapped %>%
    mutate(
      totalCN_eff = ifelse(is.finite(totalCN) & totalCN > 0, totalCN, 2),
      ccf = ifelse(
        !is.na(vaf) & !is.na(purity) & purity > 0,
        pmin(1,
             vaf * (2 * (1 - purity) + totalCN_eff * purity) / purity),
        NA_real_
      ),
      is_gain = case_when(
        gain_rule == "major_ge2" & !is.na(majorCN) & majorCN >= 2 ~ TRUE,
        gain_rule == "total_ge3" & totalCN_eff >= 3               ~ TRUE,
        TRUE                                                      ~ FALSE
      )
    )
  
  # sanity: if ccf is still NA everywhere, you know vaf/purity are broken
  if (all(is.na(mapped$ccf))) {
    warning("All CCF values are NA. Check vaf and purity.")
  }
  
  # heuristic timing probabilities based on ccf and gain
  mapped <- mapped %>%
    mutate(
      P_SUBCLONAL = ifelse(
        is.na(ccf), NA_real_,
        pmin(1, pmax(0, 1 - ccf))
      ),
      P_CLONAL    = ifelse(
        is.na(ccf), NA_real_,
        1 - P_SUBCLONAL
      ),
      P_EARLY     = ifelse(
        is.na(P_CLONAL),
        NA_real_,
        ifelse(is_gain, P_CLONAL, P_CLONAL * 0.7)
      ),
      P_LATE      = ifelse(
        is.na(P_CLONAL) | is.na(P_EARLY),
        NA_real_,
        pmax(0, P_CLONAL - P_EARLY)
      )
    ) %>%
    mutate(
      P_sum       = P_EARLY + P_LATE + P_SUBCLONAL,
      P_EARLY     = ifelse(is.na(P_sum) | P_sum == 0, NA_real_, P_EARLY     / P_sum),
      P_LATE      = ifelse(is.na(P_sum) | P_sum == 0, NA_real_, P_LATE      / P_sum),
      P_SUBCLONAL = ifelse(is.na(P_sum) | P_sum == 0, NA_real_, P_SUBCLONAL / P_sum)
    ) %>%
    dplyr::select(-P_sum)
  
  mapped
}



read_purple_sv <- function(purple_dir, patient_id) {
  # Guess a PURPLE SV file; rename if needed
  sv_file <- file.path(purple_dir, paste0(patient_id, ".purple.sv.vcf.gz"))
  if (!file.exists(sv_file)) return(tibble())
  
  vcf <- readVcf(sv_file)
  gr  <- rowRanges(vcf)
  
  tibble(
    chr1 = gsub("^chr", "", as.character(seqnames(gr))),
    pos1 = start(gr),
    chr2 = gsub("^chr", "", as.character(seqnames(gr))),
    pos2 = end(gr),
    sv_id = names(gr)
  )
}

## assign SV timing from flanking CN segments
time_svs_from_cnv <- function(svs, cnv_df, seg_timing) {
  if (!nrow(svs) || !nrow(seg_timing)) return(tibble())
  
  gr_cnv <- GRanges(
    seqnames = cnv_df$chr,
    ranges   = IRanges(cnv_df$start, cnv_df$end)
  )
  
  # breakpoint 1
  gr1  <- GRanges(seqnames = svs$chr1, ranges = IRanges(svs$pos1, svs$pos1))
  hit1 <- findOverlaps(gr1, gr_cnv)
  seg1 <- rep(NA_integer_, nrow(svs))
  seg1[queryHits(hit1)] <- subjectHits(hit1)
  
  # breakpoint 2
  gr2  <- GRanges(seqnames = svs$chr2, ranges = IRanges(svs$pos2, svs$pos2))
  hit2 <- findOverlaps(gr2, gr_cnv)
  seg2 <- rep(NA_integer_, nrow(svs))
  seg2[queryHits(hit2)] <- subjectHits(hit2)
  
  # look up timing for segments
  s1 <- seg_timing$timing_label[match(seg1, seg_timing$seg_id)]
  s2 <- seg_timing$timing_label[match(seg2, seg_timing$seg_id)]
  t1 <- seg_timing$time_score[match(seg1, seg_timing$seg_id)]
  t2 <- seg_timing$time_score[match(seg2, seg_timing$seg_id)]
  
  tibble(
    sv_id      = svs$sv_id,
    chr1       = svs$chr1,
    pos1       = svs$pos1,
    chr2       = svs$chr2,
    pos2       = svs$pos2,
    seg1       = seg1,
    seg2       = seg2,
    sv_time_score = rowMeans(cbind(t1, t2), na.rm = TRUE),
    sv_timing_label = case_when(
      !is.na(s1) & !is.na(s2) & s1 == s2 ~ s1,
      !is.na(t1) & !is.na(t2) & sv_time_score < 0.33 ~ "EARLY",
      !is.na(t1) & !is.na(t2) & sv_time_score > 0.66 ~ "SUBCLONAL",
      TRUE ~ "LATE"
    )
  )
}


compute_segment_timing <- function(timed_snvs, cnv_df) {
  if (!nrow(timed_snvs)) return(tibble())
  
  timed_snvs <- add_time_score(timed_snvs)
  
  seg_timing <- timed_snvs %>%
    group_by(seg_id) %>%
    summarise(
      n_snvs      = n(),
      mean_P_EARLY = mean(P_EARLY, na.rm = TRUE),
      mean_P_LATE  = mean(P_LATE, na.rm = TRUE),
      mean_P_SUB   = mean(P_SUBCLONAL, na.rm = TRUE),
      time_score   = mean(time_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      timing_label = case_when(
        mean_P_EARLY >= mean_P_LATE & mean_P_EARLY >= mean_P_SUB ~ "EARLY",
        mean_P_LATE  >= mean_P_EARLY & mean_P_LATE  >= mean_P_SUB ~ "LATE",
        mean_P_SUB   >= mean_P_EARLY & mean_P_SUB   >= mean_P_LATE ~ "SUBCLONAL",
        TRUE ~ NA_character_
      )
    )
  
  seg_timing <- seg_timing %>%
    mutate(
      seg_id = seg_id,
      pseudo_time = rank(time_score, ties.method = "average") / max(rank(time_score))
    )
  
  # attach CN coordinates
  seg_timing %>%
    left_join(
      cnv_df %>% mutate(seg_id = row_number()),
      by = "seg_id"
    )
}



build_pseudotime_events <- function(sample_id,
                                    timed_snvs,
                                    seg_timing) {
  # SNV events
  events_snv <- timed_snvs %>%
    add_time_score() %>%
    mutate(
      event_type = "SNV",
      label = if ("gene" %in% names(.)) {
        gene
      } else {
        paste0("chr", chr, ":", pos)
      },
      pseudo_time = rank(time_score, ties.method = "average") /
        max(rank(time_score))
    ) %>%
    dplyr::select(sample_id,
           event_type, chr, pos, label,
           time_score, pseudo_time,
           P_EARLY, P_LATE, P_SUBCLONAL,
           dplyr::everything())
  
  # CN segment events
  events_cnv <- seg_timing %>%
    mutate(
      sample_id   = sample_id,
      event_type  = "CN_segment",
      label       = paste0(chr, ":", start, "-", end)
    ) %>%
    # seg_timing already has time_score + pseudo_time
    dplyr::select(sample_id,
           event_type, chr, start, end, label,
           time_score, pseudo_time,
           timing_label,
           dplyr::everything())
  
  bind_rows(
    events_snv,
    events_cnv
  )
}





run_pseudotime_for_patient_full <- function(
    patient_id,
    base_dir,
    driver_imp = NULL,
    min_dp            = 20,
    min_muts_per_seg  = 5,
    gain_rule         = c("major_ge2","total_ge3"),
    jitter            = c(-0.05,-0.03,0,0.03,0.05),
    save_timed_dir    = NULL,
    save_events_dir   = NULL
) {
  gain_rule  <- match.arg(gain_rule)
  purple_dir <- file.path(base_dir, patient_id, "purple")
  
  cnv_file    <- file.path(purple_dir, paste0(patient_id, ".purple.cnv.somatic.tsv"))
  purity_file <- file.path(purple_dir, paste0(patient_id, ".purple.purity.tsv"))
  vcf_file    <- pick_vcf_path(purple_dir, patient_id)
  
  if (!file.exists(cnv_file))    stop("Missing CNV file: ", cnv_file)
  if (!file.exists(purity_file)) stop("Missing purity file: ", purity_file)
  
  cnv    <- read_purple_cnv(cnv_file)
  purity <- read_purity(purity_file)
  
  snvs <- read_vcf_snvs(vcf_file, min_dp = min_dp)
  if (!nrow(snvs)) stop("No SNVs after filtering for ", patient_id)
  
  mapped <- map_snvs_to_cnv(snvs, cnv)
  timed  <- probabilistic_timing(mapped,
                                 purity    = purity,
                                 jitter    = jitter,
                                 gain_rule = gain_rule) %>%
    mutate(sample_id = patient_id)
  
  # safety: ensure gene column exists (for labels)
  if (!"gene" %in% names(timed)) {
    timed$gene <- NA_character_
  }
  
  if (!nrow(timed)) stop("No timed mutations for ", patient_id)
  
  # CN segment timing derived from SNVs
  seg_timing <- compute_segment_timing(timed, cnv)
  
  # Combined pseudotime events (SNV + CN segments only)
  events <- build_pseudotime_events(patient_id, timed, seg_timing)
  
  # Save if requested
  if (!is.null(save_timed_dir)) {
    dir.create(save_timed_dir, showWarnings = FALSE, recursive = TRUE)
    readr::write_rds(
      timed,
      file.path(save_timed_dir, paste0(patient_id, "_timed_snvs.rds")),
      compress = "gz"
    )
  }
  if (!is.null(save_events_dir)) {
    dir.create(save_events_dir, showWarnings = FALSE, recursive = TRUE)
    readr::write_rds(
      events,
      file.path(save_events_dir, paste0(patient_id, "_events_full.rds")),
      compress = "gz"
    )
  }
  
  list(
    timed_snvs  = timed,
    seg_timing  = seg_timing,
    events_full = events
  )
}




### Cohort runner

run_pseudotime_cohort_full <- function(
    patients,
    base_dir      = ".",
    res_all_path  = "~/HMF/Rws/res_all.csv",
    out_events    = "pseudo_events",
    out_timed     = "timed_snvs",
    min_dp        = 20,
    gain_rule     = c("major_ge2","total_ge3"),
    overwrite     = FALSE
) {
  dir.create(out_events, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_timed,  showWarnings = FALSE, recursive = TRUE)
  
  res_all    <- readr::read_csv(res_all_path, show_col_types = FALSE)
  driver_imp <- compute_driver_importance(res_all, alpha = 0.015)
  
  gain_rule <- match.arg(gain_rule)
  index <- tibble()
  
  for (pid in patients) {
    message("▶ ", pid)
    
    ev_file <- file.path(out_events, paste0(pid, "_events_full.rds"))
    if (!overwrite && file.exists(ev_file)) {
      message("  already exists; adding to index")
      ev <- readr::read_rds(ev_file)
      if (nrow(ev)) {
        index <- bind_rows(index,
                           ev %>%
                             summarise(
                               sample_id   = pid,
                               n_events    = n(),
                               n_snv       = sum(event_type == "SNV"),
                               n_cnv_seg   = sum(event_type == "CN_segment"),
                               min_time    = min(pseudo_time, na.rm = TRUE),
                               max_time    = max(pseudo_time, na.rm = TRUE)
                             ))
      }
      next
    }
    
    res <- tryCatch(
      run_pseudotime_for_patient_full(
        patient_id     = pid,
        base_dir       = base_dir,
        driver_imp     = driver_imp,
        min_dp         = min_dp,
        gain_rule      = gain_rule,
        save_timed_dir = out_timed,
        save_events_dir= out_events
      ),
      error = function(e) {
        message("  ✖ FAILED: ", e$message)
        return(NULL)
      }
    )
    
    if (!is.null(res)) {
      ev <- res$events_full
      if (nrow(ev)) {
        index <- bind_rows(index,
                           ev %>%
                             summarise(
                               sample_id   = pid,
                               n_events    = n(),
                               n_snv       = sum(event_type == "SNV"),
                               n_cnv_seg   = sum(event_type == "CN_segment"),
                               min_time    = min(pseudo_time, na.rm = TRUE),
                               max_time    = max(pseudo_time, na.rm = TRUE)
                             ))
      }
    }
  }
  
  index
}


allPats <- list.files("~/HMF/isofox/data_isofox/")

##### RUN PIPELINE #####

res_test <- run_pseudotime_for_patient_full(
  patient_id = allPats[2],
  base_dir = "~/HMF/somatic/DR-252-update1",
  save_timed_dir = "timed_snvs_full",
  save_events_dir= "pseudo_events_full"
)

head(res$events_full)

snvs   <- read_vcf_snvs(vcf_file_2, min_dp = 10)
mapped <- map_snvs_to_cnv(snvs, cnv)
timed  <- probabilistic_timing(mapped, purity = purity, jitter = jitter, gain_rule = gain_rule)

res_pat2 <- run_pseudotime_for_patient_full(
        patient_id = allPats[2],
        base_dir = "~/HMF/somatic/DR-252-update1",
        driver_imp = drv_importance,   # still available if you want to use it somewhere
        min_dp            = 20,
        min_muts_per_seg  = 5,
        gain_rule         = c("major_ge2","total_ge3"),
        jitter            = c(-0.05,-0.03,0,0.03,0.05),
        save_timed_dir = "timed_snvs_full",
        save_events_dir= "pseudo_events_full"
    )


#patients   <- list.files("~/HMF/isofox/data_isofox/")
patients <- allPats
base_dir   <- "~/HMF/somatic/DR-252-update1"

index_full <- run_pseudotime_cohort_full(
  patients      = patients,
  base_dir      = base_dir,
  res_all_path  = "~/HMF/subclonality/res_all_pc_spec.csv",
  out_events    = "pseudo_events",
  out_timed     = "timed_snvs_full",
  min_dp        = 20,
  gain_rule     = "major_ge2"
)



#### annotate ####

library(dplyr)
library(vroom)
library(GenomicRanges)
library(IRanges)
library(tibble)

BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("org.Hs.eg.db")


library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(dplyr)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(GenomeInfoDb)  # for seqlevelsStyle
library(dplyr)

library(GenomicRanges)
library(GenomeInfoDb)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(dplyr)

library(dplyr)
library(readr)
library(tibble)

# you already have this:
# read_vcf_snvs()
# assign_gene_by_coords()
# gene_regions_hg19_pc

library(dplyr)

annotate_sample_timed_snvs <- function(sample_id,
                                       base_dir,
                                       timed_snvs,      # already loaded
                                       gene_regions,
                                       min_dp = 10,
                                       padding = 0L) {
  # 1) VCF path (adjust extension if needed)
  vcf_file <- file.path(
    base_dir,
    sample_id,
    "purple",
    paste0(sample_id, ".purple.somatic.vcf")
  )
  
  # 2) Read SNVs and assign genes
  snvs_all <- read_vcf_snvs(vcf_file, min_dp = min_dp, keep_all = TRUE)
  
  snvs_gene <- assign_gene_by_coords(
    snvs        = snvs_all,
    gene_regions = gene_regions,
    padding     = padding
  )
  # snvs_gene has: gene, chr, pos, ref, alt, dp, vaf, passed_filters
  
  # Normalise types and chr naming
  snvs_gene_clean <- snvs_gene %>%
    mutate(
      chr = gsub("^chr", "", as.character(chr)),
      pos = as.integer(pos),
      ref = as.character(ref),
      alt = as.character(alt)
    )
  
  # 3) Subset timed_snvs for this sample and normalise
  timed_sample <- timed_snvs %>%
    filter(sample_id == !!sample_id) %>%
    mutate(
      chr = gsub("^chr", "", as.character(chr)),
      pos = as.integer(pos),
      ref = as.character(ref),
      alt = as.character(alt)
    )
  
  # 4) Join timing to gene annotation by chr+pos+ref+alt
  timed_annot <- timed_sample %>%
    left_join(
      snvs_gene_clean %>% dplyr::select(chr, pos, ref, alt, gene),
      by = c("chr", "pos", "ref", "alt")
    )
  
  # If timed_sample already had a gene column, coalesce and clean
  if ("gene.x" %in% names(timed_annot) && "gene.y" %in% names(timed_annot)) {
    timed_annot <- timed_annot %>%
      mutate(gene = coalesce(gene.x, gene.y)) %>%
      dplyr::select(-gene.x, -gene.y)
  } else if ("gene.y" %in% names(timed_annot)) {
    timed_annot <- timed_annot %>%
      rename(gene = gene.y)
  }
  
  timed_annot
}


timed_ar_sample <- annotate_sample_timed_snvs(
  sample_id   = "ACTN01020030T",
  base_dir    = "~/HMF/somatic/DR-252-update1",
  timed_snvs  = snvs_test,
  gene_regions = gene_regions_hg19_pc,
  min_dp      = 10,
  padding     = 0L
)

# sanity checks
table(is.na(timed_ar_sample$gene))
timed_ar_sample %>% filter(gene == "FOXA1") %>% head()
timed_ar_sample %>% filter(gene == "AR") %>% head()


timed_ar_sample <- timed_ar_sample %>%
  mutate(gene = coalesce(gene.x, gene.y)) %>%
  dplyr::select(-gene.x, -gene.y)

# Check AR mutations in this sample
timed_ar_sample %>%
  filter(gene == "AR")


timed_file    <- file.path(subclonality_dir, "timed_snvs_full", paste0(patient_id, "_timed_snvs.rds"))
cnv_gene_file <- file.path(purple_dir, paste0(patient_id, ".purple.cnv.gene.tsv"))

txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
timed_annot_hg19 <- annotate_timed_snvs_all_genes(timed_snvs, txdb = txdb_hg19)
table(is.na(timed_annot_hg19$gene_symbol))


vcf_file_test    <-  "~/HMF/somatic/DR-252-update1/ACTN01020030T/purple/ACTN01020030T.purple.somatic.vcf.gz"
snvs_test <- read_vcf_snvs(vcf_file_test, min_dp = 10)

# txdb_hg38 <- TxDb.Hsapiens.UCSC.hg38.knownGene
# timed_annot_hg38 <- annotate_timed_snvs_all_genes(timed_snvs, txdb = txdb_hg38)
# table(is.na(timed_annot_hg38$gene_symbol))

## hg19 gave more genes sp we'll use that

timed_annot <- timed_annot_hg19   # or hg19, whichever worked


patient_id <- allPats[2]
base_dir   <- "~/HMF/somatic/DR-252-update1"
purple_dir <- file.path(base_dir, patient_id, "purple")
subclonality_dir <- "~/HMF/subclonality/"


timed_snvs    <- readr::read_rds(timed_file)
genes_of_interest <- c("RB1","AR","ATM","KDM6A","KMT2D","APC")

annotate_snv_timing <- function(timed_snvs) {
  timed_snvs %>%
    mutate(
      timing_label = case_when(
        P_EARLY     >= P_LATE & P_EARLY     >= P_SUBCLONAL ~ "EARLY",
        P_LATE      >= P_EARLY & P_LATE     >= P_SUBCLONAL ~ "LATE",
        P_SUBCLONAL >= P_EARLY & P_SUBCLONAL >= P_LATE     ~ "SUBCLONAL",
        TRUE ~ NA_character_
      )
    )
}



##### GENES OF INTEREST #####

genes_of_interest <- c("RB1","AR","ATM","KDM6A","KMT2D","APC")

gene_timing_sample <- timed_labeled %>%
  filter(gene %in% genes_of_interest) %>%
  group_by(gene) %>%
  summarise(
    n_snv       = n(),
    n_early     = sum(timing_label == "EARLY", na.rm = TRUE),
    n_late      = sum(timing_label == "LATE", na.rm = TRUE),
    n_subclonal = sum(timing_label == "SUBCLONAL", na.rm = TRUE),
    gene_timing = case_when(
      n_early     > 0 ~ "EARLY",
      n_late      > 0 ~ "LATE",
      n_subclonal > 0 ~ "SUBCLONAL",
      TRUE            ~ "No SNV"
    ),
    .groups = "drop"
  )



# sanity checks
head(timed_snvs_annot$gene)
table(is.na(timed_snvs_annot$gene))
intersect(c("RB1","AR","ATM","KDM6A","KMT2D","APC"), timed_snvs_annot$gene)

## import drivers_pc
drivers_pc <- read.csv("~/HMF/drivers/drivers_pc.csv")

### sanity checks
drv_this <- drivers_pc %>%
  filter(patient_id == patient_id,
         gene == "AR")

drv_this



snvs_test <- read_vcf_snvs(vcf_file, min_dp = 20)


ar_chr <- "X"   # or "chrX", depending on your vcf 'chr' column
ar_start <- 66700000
ar_end   <- 67000000

snvs_AR <- snvs_debug %>%
  dplyr::filter(chr == ar_chr,
                pos >= ar_start, pos <= ar_end)

snvs_AR




############ add gene coordinates independently

library(tibble)

gene_regions_hg19_pc <- tribble(
  ~gene,       ~chr,   ~start,     ~end,
  
  # Core prostate cancer drivers
  "AR",        "X",    66763865,   66950461,
  "RB1",       "13",   48303751,   48935716,
  "PTEN",      "10",   87863119,   87971930,
  "TP53",      "17",   7565097,    7590856,
  
  # Homologous recombination / DDR
  "BRCA2",     "13",   32889617,   32973805,
  "BRCA1",     "17",   41196312,   41322266,
  "ATM",       "11",   108093558,  108239829,
  "CHEK2",     "22",   29083712,   29100047,
  
  # chr13 region
  "RNASEH2B",  "13",   51482371,   51490694,
  
  # Prostate cancer–relevant drivers
  "CDK12",     "17",   37677136,   37807366,
  "SPOP",      "17",   49603748,   49627980,
  "CHD1",      "5",    97738871,   97899801,
  "FOXA1",     "14",   38061669,   38079329,
  "NKX3_1",    "8",    23434336,   23437333,
  
  # Common oncogenes (amplicons)
  "MYC",       "8",    128748315,  128753680,
  "CCND1",     "11",   68887542,   68990020,
  "EGFR",      "7",    55086714,   55275823,
  
  # WNT pathway
  "APC",       "5",    112707498,  112846239,
  
  # **Epigenetic regulators you requested**
  "KDM6A",     "X",    45042558,   45148778,   # UTX
  "KMT2D",     "12",   49412701,   49634355    # MLL2
)

gene_regions_hg19_pc <- gene_regions_hg19_pc %>%
  mutate(start = start - 5000,
         end   = end + 5000)


library(dplyr)

library(dplyr)

assign_gene_by_coords <- function(snvs, gene_regions, padding = 0L) {
  # ensure expected columns exist in gene_regions
  stopifnot(all(c("gene", "chr", "start", "end") %in% names(gene_regions)))
  
  # drop any existing 'gene' column in snvs to avoid conflicts
  snvs_clean <- snvs %>% dplyr::select(-any_of("gene"))
  
  gene_regions %>%
    mutate(
      start = start - padding,
      end   = end   + padding
    ) %>%
    inner_join(
      snvs_clean,
      by = "chr",
      relationship = "many-to-many"  # expected: many SNVs per chr per gene
    ) %>%
    filter(pos >= start, pos <= end) %>%
    # keep a tidy set of columns
    dplyr::select(
      gene, chr, pos,
      ref, alt,
      dp, vaf,
      passed_filters
    )
}


snvs_test <- read_vcf_snvs(vcf_file_test, min_dp = 10)

snvs_gene <- assign_gene_by_coords(snvs_test, gene_regions_hg19_pc, padding = 0L)

snvs_gene %>% count(gene)
snvs_gene %>% filter(gene == "AR") %>% head()



#### let's go back to timed_labeled

annotate_snv_timing <- function(timed_snvs) {
  timed_snvs %>%
    mutate(
      timing_label = case_when(
        P_EARLY     >= P_LATE & P_EARLY     >= P_SUBCLONAL ~ "EARLY",
        P_LATE      >= P_EARLY & P_LATE     >= P_SUBCLONAL ~ "LATE",
        P_SUBCLONAL >= P_EARLY & P_SUBCLONAL >= P_LATE     ~ "SUBCLONAL",
        TRUE ~ NA_character_
      )
    )
}

timed_labeled <- annotate_snv_timing(timed_annot)
table(timed_labeled$timing_label, useNA = "ifany")


