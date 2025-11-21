library(msigdbr)
library(dplyr)

brca_core <- msigdbr(species = "Homo sapiens") %>%
  filter(gs_name %in% c(
    "PUJANA_BRCA2_PCC_NETWORK",
    "PUJANA_BRCA1_PCC_NETWORK",
    "HALLMARK_DNA_REPAIR",
    "KEGG_HOMOLOGOUS_RECOMBINATION",
    "GOBP_HOMOLOGOUS_RECOMBINATION"
  ))

library(msigdbr)
library(dplyr)

pathway_names <- c(
  # PI3K/AKT
  "REACTOME_PI3K_AKT_ACTIVATION",
  "KEGG_PI3K_AKT_SIGNALING_PATHWAY",
  "REACTOME_SIGNALING_BY_AKT",
  
  # SRC / FAK / Integrin
  "REACTOME_SIGNALING_BY_SRC_FAMILY_KINASES",
  "REACTOME_FOCAL_ADHESION",
  "KEGG_FOCAL_ADHESION",
  "PID_SRC_PATHWAY",
  "PID_FAK_PATHWAY",
  
  # MAPK / ERK
  "REACTOME_MAPK1_ERK2_ACTIVATION",
  "KEGG_MAPK_SIGNALING_PATHWAY",
  "PID_MAPK_PATHWAY",
  "PID_ERK_PATHWAY",
  
  # ATR–ATF2 & checkpoint
  "REACTOME_ATR_PATHWAY",
  "REACTOME_ACTIVATION_OF_ATF2_BY_ATR",
  "GOBP_ATR_SIGNALING",
  "PID_ATR_PATHWAY"
)

msig_df <- msigdbr(species = "Homo sapiens") %>%
  filter(gs_name %in% pathway_names)

pws <- split(msig_df$gene_symbol, msig_df$gs_name)

DNAvsRNA <- function(gene,
                     pws,
                     pats    = NULL,      # optional: user-supplied mutant patients
                     use_cnv = FALSE,     # TRUE to also include CNV deletions
                     logFC_cutoff = 1.5,
                     fdr_cutoff   = 0.01,
                     top_n        = 10,
                     minSize      = 1,
                     maxSize      = 5000) {
  
  ## 0. Basic checks
  if (length(gene) != 1) {
    stop("DNAvsRNA() only accepts a single gene. Please pass one gene symbol.")
  }
  if (missing(pws) || is.null(pws)) {
    stop("Please provide a named list of pathways in 'pws'.")
  }
  if (!is.list(pws) || is.null(names(pws))) {
    stop("pws must be a *named list* of pathways (each is a character vector of genes).")
  }
  
  ## 1. Choose mutant patients
  # If pats is NULL -> derive from mutations
  # If pats is provided -> use that directly
  if (is.null(pats)) {
    mut_pats <- get_mut_pats(gene)
    message("Using get_mut_pats(): ", length(mut_pats), " mutant patients for ", gene)
  } else {
    mut_pats <- as.character(pats)
    message("Using user-supplied pats: ", length(mut_pats), " mutant patients")
  }
  
  ## 2. Optionally add CNV-deleted patients
  if (use_cnv) {
    cnv_df   <- getCNV(gene, "DEL")          # your existing helper
    cnv_pats <- unique(cnv_df$patientID)
    message("Adding ", length(cnv_pats), " CNV-del patients")
    
    mut_pats <- unique(c(mut_pats, cnv_pats))
  }
  
  ## 3. Define universe of patients and non-mutant group
  # Assuming counts has all the samples you care about
  all_pats <- colnames(counts)
  
  # Keep only IDs that actually exist in counts
  mut_pats <- intersect(mut_pats, all_pats)
  if (length(mut_pats) == 0) {
    stop("No mutant patients found in counts after matching IDs.")
  }
  
  NON_mut_pats <- setdiff(all_pats, mut_pats)
  if (length(NON_mut_pats) == 0) {
    stop("No non-mutant patients left (NON_mut_pats is empty).")
  }
  
  message("Performing DE for ",
          length(mut_pats), " mutant vs ",
          length(NON_mut_pats), " non-mutant patients")
  
  
  mut_cnv_pats_df    <- counts %>% dplyr::select(all_of(unlist(mut_pats)))
  NONmut_cnv_pats_df <- counts %>% dplyr::select(all_of(unlist(NON_mut_pats)))
  count_df           <- cbind(mut_cnv_pats_df, NONmut_cnv_pats_df)
  
  group  <- factor(c(rep("T", length(mut_pats)),
                     rep("C", length(NON_mut_pats))))
  design <- model.matrix(~ batch + group)  # batch assumed global
  
  dge <- edgeR::DGEList(counts = count_df)
  keep <- edgeR::filterByExpr(dge, group = group)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- edgeR::calcNormFactors(dge)
  
  message("Estimating dispersion")
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design)
  qlf <- edgeR::glmQLFTest(fit, coef = "groupT")
  
  DE  <- edgeR::topTags(qlf, n = nrow(dge$counts))
  res <- DE$table
  
  ## 2. Volcano plot
  res$gene   <- rownames(res)
  res$logFDR <- -log10(res$FDR)
  
  res$Significant <- "NS"
  res$Significant[res$logFC >  logFC_cutoff & res$FDR < fdr_cutoff] <- "Up"
  res$Significant[res$logFC < -logFC_cutoff & res$FDR < fdr_cutoff] <- "Down"
  
  message("Differential analysis complete (",
          sum(res$Significant == "Up"), " up, ",
          sum(res$Significant == "Down"), " down )")
  
  sig_hits <- res[res$Significant != "NS", ]
  
  Volcano <- ggplot2::ggplot(res, ggplot2::aes(x = logFC, y = logFDR)) +
    ggplot2::geom_point(ggplot2::aes(color = Significant), alpha = 0.6) +
    ggplot2::geom_hline(yintercept = -log10(fdr_cutoff),
                        linetype = "dashed", color = "black") +
    ggplot2::geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff),
                        linetype = "dotted", color = "black") +
    ggplot2::scale_color_manual(
      values = c(
        "Up"   = "#b2182b",
        "Down" = "#2166ac",
        "NS"   = "#bdbdbd"
      )
    ) +
    ggrepel::geom_text_repel(
      data = sig_hits,
      ggplot2::aes(label = gene),
      max.overlaps = 30,
      size = 3,
      box.padding = 0.3,
      point.padding = 0.2,
      segment.color = "grey50"
    ) +
    ggplot2::labs(
      title = paste("Differential expression for:", gene, "mut vs non-mut"),
      x = "Log2 Fold Change",
      y = "-Log10 FDR"
    ) +
    ggplot2::theme_minimal()
  
  ## 3. Ranked list for GSEA
  set.seed(123)
  gene_list <- res$logFC
  names(gene_list) <- rownames(res)  # assume symbols
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  gene_list_names <- names(gene_list)
  
  ## 4. DEBUG: pathway sizes & overlaps (this is the code you were asking about)
  path_sizes    <- sapply(pws, length)
  path_overlap  <- sapply(pws, function(g) length(intersect(gene_list_names, g)))
  
  message("Pathway sizes (total genes per set):")
  print(path_sizes)
  
  message("Overlap with ranked gene list (per pathway):")
  print(path_overlap)
  
  ## 5. Run fgsea
  fgseaRes <- fgsea::fgseaMultilevel(
    pathways = pws,
    stats    = gene_list,
    minSize  = minSize,
    maxSize  = maxSize
  )
  
  fgseaRes_df <- as.data.frame(fgseaRes)
  
  if (nrow(fgseaRes_df) == 0) {
    warning("fgsea returned no pathways (likely zero overlap or minSize too high).")
    return(list(
      gene      = gene,
      DE_res    = res,
      Vol       = Volcano,
      GSEA      = fgseaRes_df,
      GSEA_plot = NULL,
      pathway_sizes   = path_sizes,
      pathway_overlap = path_overlap
    ))
  }
  
  top_n <- min(top_n, nrow(fgseaRes_df))
  topPathways_df <- fgseaRes_df[order(fgseaRes_df$padj), ][1:top_n, ]
  
  ## 6. Clear significance vs NS in GSEA plot
  topPathways_df$signif <- ifelse(topPathways_df$padj < fdr_cutoff,
                                  "Significant", "Not significant")
  topPathways_df$direction <- ifelse(topPathways_df$NES >= 0,
                                     "Up", "Down")
  
  topPathways_df$cat <- paste(topPathways_df$direction,
                              topPathways_df$signif,
                              sep = "_")
  topPathways_df$cat <- factor(
    topPathways_df$cat,
    levels = c(
      "Up_Significant",
      "Up_Not significant",
      "Down_Significant",
      "Down_Not significant"
    )
  )
  
  message("Pre-ranked GSEA complete; showing top ", top_n, " pathways.")
  
  path_enrich_plot <- ggplot2::ggplot(
    topPathways_df,
    ggplot2::aes(x = reorder(pathway, NES), y = NES, fill = cat)
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = "Pathway",
      y = "Normalized Enrichment Score (NES)",
      title = paste("Top enriched pathways for", gene, "DE")
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Up_Significant"       = "#b2182b",
        "Up_Not significant"   = "#f4a582",
        "Down_Significant"     = "#2166ac",
        "Down_Not significant" = "#92c5de"
      ),
      labels = c(
        "Up_Significant"       = "Up (NES>0, padj<0.05)",
        "Up_Not significant"   = "Up (NES>0, NS)",
        "Down_Significant"     = "Down (NES<0, padj<0.05)",
        "Down_Not significant" = "Down (NES<0, NS)"
      )
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.title       = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank()
    )
  
  ## 7. Return everything, including debug info
  list(
    gene            = gene,
    DE_res          = res,
    Vol             = Volcano,
    GSEA            = topPathways_df,
    GSEA_plot       = path_enrich_plot,
    pathway_sizes   = path_sizes,
    pathway_overlap = path_overlap
  )
}



BRCA2_DE_gsea <- DNAvsRNA("BRCA2", pws = pws)
BRCA2_DE_gsea_87 <- DNAvsRNA("BRCA2", pws = pws, pats = l)


BRCA2_DE_gsea$GSEA_plot
