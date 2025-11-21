### alts_vs_ssgsea - prostate specific
pc_drivers <- build_pc_driver_dataframe(pats, base_dir, driver_genes = driver_genes)

# prostate tumour supressors
TSG_pc <- c("PTEN","RB1","TP53","BRCA2","ATM","CDK12","ARID1A","CHEK2","PALB2","FANCA","RAD51C","RAD51D","BAP1")
# prostae oncogenes
ONC_pc <- c("AR","MYC","CCND1","MDM2","MET","BRAF","PIK3CA","PIK3CB","AKT1")

pc_drivers <- pc_drivers %>%
  mutate(driver_class = case_when(
    gene %in% TSG_pc ~ "TSG",
    gene %in% ONC_pc ~ "Oncogene",
    TRUE ~ "Other"
  ))


# PURPLE’s likelihood scoring assumes pan-cancer priors:
#   e.g., RB1 deep deletion is rare pan-cancer but common in CRPC → PURPLE underweights it.
# 
# We fix that.
# 
# ✔ Prostate-specific CN rules
# TSGs:
#   
#   HomDel (minCN < 0.5) → strong driver (score = 4)
# 
# HetDel + LOH (minCN ~1) → moderate–strong (score 3)
# 
# LOH only → weak (score = 1)
# 
# No CN event → 0
# 
# Oncogenes:
#   
#   High-level amp (totalCN ≥ 6) → strong (4)
# 
# Amp (totalCN 4–6) → moderate (3)
# 
# Gain (totalCN 3) → weak (1)

drivers_pc <- pc_drivers %>%
  left_join(cn_gene_sub, by = c("sampleId","gene")) %>%
  mutate(
    prostate_CN_score = case_when(
      driver_class == "TSG" & minCopyNumber < 0.5 ~ 4,
      driver_class == "TSG" & minCopyNumber < 1.5 ~ 3,
      driver_class == "TSG" & minorCN == 0        ~ 1,
      driver_class == "Oncogene" & totalCN >= 6   ~ 4,
      driver_class == "Oncogene" & totalCN >= 4   ~ 3,
      driver_class == "Oncogene" & totalCN == 3   ~ 1,
      TRUE                                        ~ 0
    )
  )



### wrapper function
# Required packages
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(VariantAnnotation)

prostate_genes = c( "AR","FOXA1","SPOP","ERG","CCND1","CDK4","CDK6","CDKN2B","MDM2",
                    "RB1","TP53","KDM5A","KDM6A","KMT2C","KMT2D","PARP1","MUTYH",
                    "ERCC2","ERCC5","ATM","ATR","PRKDC","CDK12","BAP1","BARD1","BLM",
                    "BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","FANCA","PALB2","RAD50",
                    "RAD51","RAD51B","RAD51C","RAD51D","MLH1","MSH2","MSH6","BRAF",
                    "EP300","B2M","MYC","ARID1A","MET","PIK3CA","PIK3CB","PTEN",
                    "AKT1","CTNNB1","APC"
                      )


build_prostate_driver_catalog_simple <- function(
    samples,
    purple_dir,
    prostate_genes = c(
      "AR","FOXA1","SPOP","ERG","CCND1","CDK4","CDK6","CDKN2B","MDM2",
      "RB1","TP53","KDM5A","KDM6A","KMT2C","KMT2D","PARP1","MUTYH",
      "ERCC2","ERCC5","ATM","ATR","PRKDC","CDK12","BAP1","BARD1","BLM",
      "BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","FANCA","PALB2","RAD50",
      "RAD51","RAD51B","RAD51C","RAD51D","MLH1","MSH2","MSH6","BRAF",
      "EP300","B2M","MYC","ARID1A","MET","PIK3CA","PIK3CB","PTEN",
      "AKT1","CTNNB1","APC"
    )
) {
  library(dplyr)
  library(readr)
  library(purrr)
  
  ## -----------------------------
  ## 1. Driver catalogs (all char → then convert)
  ## -----------------------------
  driver_list <- lapply(samples, function(sid) {
    f <- file.path(purple_dir, sid, "purple",
                   paste0(sid, ".purple.driver.catalog.somatic.tsv"))
    if (!file.exists(f)) {
      warning("Driver file not found for sample: ", sid)
      return(NULL)
    }
    
    df <- read_tsv(
      f,
      col_types = cols(.default = col_character()),
      na = c("NA","."),
      show_col_types = FALSE
    )
    df$sampleId <- sid
    df
  })
  
  drivers_all_raw <- bind_rows(driver_list)
  
  # let readr guess types *once* on the combined table
  drivers_all <- type_convert(drivers_all_raw, na = c("NA","."),
                              guess_integer = TRUE)
  
  ## -----------------------------
  ## 2. CNV gene tables
  ## -----------------------------
  cn_gene_list <- lapply(samples, function(sid) {
    f <- file.path(purple_dir, sid, "purple",
                   paste0(sid, ".purple.cnv.gene.tsv"))
    if (!file.exists(f)) {
      warning("CNV gene file not found for sample: ", sid)
      return(NULL)
    }
    
    df <- read_tsv(
      f,
      col_types = cols(.default = col_character()),
      na = c("NA","."),
      show_col_types = FALSE
    )
    df$sampleId <- sid
    df
  })
  
  cn_gene_raw <- bind_rows(cn_gene_list)
  cn_gene <- type_convert(cn_gene_raw, na = c("NA","."), guess_integer = TRUE)
  
  cn_gene_sub <- cn_gene %>%
    dplyr::select(sampleId, gene, minCopyNumber, maxCopyNumber)
  
  ## -----------------------------
  ## 3. Prostate filter + role assignment
  ## -----------------------------
  drivers_pc <- drivers_all %>%
    filter(gene %in% prostate_genes)
  
  TSG_pc <- c("PTEN","RB1","TP53","BRCA2","ATM","CDK12","ARID1A",
              "CHEK2","PALB2","FANCA","RAD51C","RAD51D","BAP1",
              "MSH2","MSH6","MLH1")
  ONC_pc <- c("AR","MYC","CCND1","MDM2","MET","BRAF",
              "PIK3CA","PIK3CB","AKT1")
  
  drivers_pc <- drivers_pc %>%
    mutate(
      driver_class = case_when(
        gene %in% TSG_pc ~ "TSG",
        gene %in% ONC_pc ~ "Oncogene",
        TRUE            ~ "Other"
      )
    ) %>%
    left_join(cn_gene_sub, by = c("sampleId","gene"))
  
  drivers_pc <- drivers_pc[, c(1:19)]
  names(drivers_pc)[c(16,17)] <- c("minCopyNumber", "maxCopyNumber") 
  
  ## -----------------------------
  ## 4. Prostate-specific CN score
  ## -----------------------------
  drivers_pc <- drivers_pc %>%
    mutate(
      prostate_CN_score = case_when(
        driver_class == "TSG" & !is.na(minCopyNumber) & minCopyNumber < 0.5 ~ 4,
        driver_class == "TSG" & !is.na(minCopyNumber) & minCopyNumber < 1.5 ~ 3,
        driver_class == "Oncogene" & !is.na(maxCopyNumber) & maxCopyNumber >= 6 ~ 4,
        driver_class == "Oncogene" & !is.na(maxCopyNumber) & maxCopyNumber >= 4 ~ 3,
        driver_class == "Oncogene" & !is.na(maxCopyNumber) & maxCopyNumber == 3 ~ 1,
        TRUE ~ 0
      )
    )
  
  ## -----------------------------
  ## 5. Prostate-specific SNV score from driver catalog flags
  ## -----------------------------
  # ensure flag columns exist; if not, create as 0
  flag_cols <- c("missense","nonsense","truncating","splice","frameshift")
  for (cc in flag_cols) {
    if (!cc %in% names(drivers_pc)) {
      drivers_pc[[cc]] <- 0
    }
  }
  
  drivers_pc <- drivers_pc %>%
    mutate(
      missense    = as.numeric(missense),
      nonsense    = as.numeric(nonsense),
      truncating  = as.numeric(truncating),
      splice      = as.numeric(splice),
      frameshift  = as.numeric(frameshift),
      missense    = ifelse(is.na(missense),   0, missense),
      nonsense    = ifelse(is.na(nonsense),   0, nonsense),
      truncating  = ifelse(is.na(truncating), 0, truncating),
      splice      = ifelse(is.na(splice),     0, splice),
      frameshift  = ifelse(is.na(frameshift), 0, frameshift),
      # high-impact LoF flag
      lof_flag = (nonsense + truncating + splice + frameshift) > 0
    ) %>%
    mutate(
      prostate_SNV_score = case_when(
        gene == "SPOP" & missense > 0          ~ 4, # canonical prostate driver
        lof_flag                                 ~ 4,
        missense > 0                             ~ 2,
        TRUE                                     ~ 0
      )
    )
  
  ## -----------------------------
  ## 6. Combined prostate driver score
  ## -----------------------------
  drivers_pc <- drivers_pc %>%
    mutate(
      prostate_driver_score = prostate_CN_score + prostate_SNV_score,
      prostate_driver_strength = case_when(
        prostate_driver_score >= 6 ~ "Strong",
        prostate_driver_score >= 3 ~ "Moderate",
        TRUE                       ~ "Weak"
      )
    )
  
  list(
    drivers_raw = drivers_all,   # all original PURPLE drivers
    cn_gene     = cn_gene,       # full CNV gene table
    cn_gene_sub = cn_gene_sub,   # CN features used for scoring
    drivers_pc  = drivers_pc     # prostate-specific driver catalog with scores
  )
}



samples   <- allPats  # or c("pat1","pat2",...)
purple_dir <- "~/HMF/somatic/DR-252-update1"

driver_cat <- build_prostate_driver_catalog(samples, purple_dir)

drivers_pc  <- res$drivers_pc
cn_gene_sub <- res$cn_gene_sub
mut_gene    <- res$mut_gene


drivers_pc_mod_strong <- drivers_pc[drivers_pc$prostate_driver_strength != "Weak", ]  

res_pc_spec2 <- analyze_driver_ssgsea(scores, metadata, drivers_pc_mod_strong,
              covars = c("tumorPurity"),
              use_spline_for_tumorPurity = TRUE, spline_df = 3)

plot_resall_bubbles(res_pc_spec$all, alpha = 0.01)
