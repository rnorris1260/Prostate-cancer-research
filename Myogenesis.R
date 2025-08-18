## heatmap of myogenesis genes down-regulated in patients with BRCA2 mutations across whole cohort 

BRCA2_DE <- BRCA2_allres_qlf$DE_res
ranks_df <- data.frame("gene" = BRCA2_DE$gene, "logFC" = BRCA2_DE$logFC, "FDR" = BRCA2_DE$FDR)

BRCA2_myog_ranks <- path_genes_vol(ranks_df, path, "blue")

BRCA2_myog <- subset(BRCA2_myog_ranks, BRCA2_myog_ranks$in_path == TRUE)
BRCA2_myog_sig <- subset(BRCA2_myog, BRCA2_myog$logFC < -1.5 & BRCA2_myog$FDR < 0.01)

counts_myog_sig <- subset(counts, rownames(counts) %in% BRCA2_myog_sig$gene)

counts_myog_sig_mat <- as.matrix(counts_myog_sig)

# Log-transform (variance stabilizing)
log_mat <- log2(counts_myog_sig_mat + 1)

library(pheatmap)
# Plot heatmap
pm_BRCA2 <- pheatmap(
        log_mat,
        scale = "row", # z-score per gene
        clustering_distance_rows = "euclidean",
        clustering_distance_cols = "euclidean",
        clustering_method = "complete",
        color = colorRampPalette(c("blue", "white", "red"))(100)
)

# Cut into 2 clusters
clusters <- cutree(pm_BRCA2$tree_col, k = 2)

# Get samples from cluster 1
samples_cluster2 <- names(clusters[clusters == 2])

# lets pull out the metadata for these samples
metadata_myog_up <- subset(metadata, metadata$sampleId %in% samples_cluster2)
Myog_down_pats <- setdiff(pats, samples_cluster2)
metadata_myog_down <- subset(metadata, metadata$sampleId %in% Myog_down_pats)

write.csv(metadata_myog_up, "metadata_myog_up.csv")
write.csv(metadata_myog_down, "metadata_myog_down.csv")

# Specify samples to highlight
highlight_samples <- mut_pats_rna

# Create annotation dataframe
annotation_col <- data.frame(
  Highlight = ifelse(colnames(log_mat) %in% highlight_samples, "Yes", "No")
)
rownames(annotation_col) <- colnames(log_mat)

# Annotation colors
ann_colors <- list(Highlight = c(Yes = "gold", No = "white"))

# Plot heatmap with highlighted samples
pheatmap(pm_BRCA2 <- pheatmap(
  log_mat,
  scale = "row", # z-score per gene
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  annotation_col = annotation_col, annotation_colors = ann_colors)
)


### lets do a differential analysis between samples in cluster 1 vs cluster 2

### need to find out what mutations the patients with up-regulated myogenesis have
