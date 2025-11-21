library(dplyr)
library(ggplot2)
library(edgeR)
library(ggrepel)
library(fgsea)
library(msigdbr)

### Driver ###

get_mut_pats <- function(gene){
  setwd("~/HMF/somatic/DR-252-update1")
  mut_pats <- list()
  for (pat in pats){
    driver <- read.delim(paste(pat,"/purple/",pat,".purple.driver.catalog.somatic.tsv",sep=""), sep="\t")
    genes <- driver$gene
    if (gene %in% genes == T){
      #print(pat)
      mut_pats <- append(mut_pats, pat)
    }
  }
  return(mut_pats)
}


### CNV ###

getCNV <- function(gene, CNV_arg){
  cnv_df <- data.frame(patientID = character(), CN = numeric())
  for (i in 1:length(pats)){
    #    print(pats[i])
    gene_value = gene
    cnv <- read.delim(paste("~/HMF/somatic/DR-252-update1/",pats[i], "/purple/", pats[i], ".purple.cnv.gene.tsv", sep=""))
    cnv_gene <- subset(cnv, cnv$gene == gene_value)
    CN = cnv_gene$maxCopyNumber
    #    print(CN)
    cnv_df[i, ] <- c(pats[i], CN)
  }
  if(CNV_arg == "DEL"){
    CNVs <- subset(cnv_df, cnv_df$CN < 1.1)
  }
  else if(CNV_arg == "AMP"){
    CNVs <- subset(cnv_df, cnv_df$CN > 6)
  }
  else if(CNV_arg == "ALL"){CNVs <- subset(cnv_df, cnv_df$CN < 1.1 | cnv_df$CN > 6)}
  return(CNVs)
}


### RNA ###

# load rna-seq count matrix
counts <- read.csv("~/HMF/gene_counts_from_fittedfragments/HMF_rna_counts.csv")
rownames(counts) <- counts$X
counts <- counts[-1]

pats <- colnames(counts)


### metadata ###

# Should have columns: sample, group, batch
metadata <- read.csv("~/HMF/metadata/metadata.csv", sep = "\t")

# Ensure order of metadata matches count data columns
metadata <- metadata[metadata$sampleId %in% colnames(counts), ]
metadata_purity <- data.frame("patientID" = metadata$sampleId, "batch" = metadata$tumorPurity)
metadata_purity$purity_scaled <- scale(metadata_purity$batch, center = TRUE, scale = TRUE)
#metadata_purity$batch <- round(metadata_purity$batch *10)
# batch
batch <- metadata_purity$batch


# Human Hallmark gene sets
m_df <- msigdbr(species = "Homo sapiens", collection = "H")  # H = Hallmark
pathways <- split(m_df$gene_symbol, m_df$gs_name)




### Main function

DNAvsRNA <- function(gene){
  
  if (length(gene) > 1){
    mut_pats <- list()
    for (i in 1:length(gene)){
      mut_pats_i <- get_mut_pats(gene[i])
      print(paste("found", length(mut_pats_i), "patients for", gene[i], sep=" "))
      mut_pats <- append(mut_pats, mut_pats_i)
      mut_pats <- unique(mut_pats)
    }}
    else {mut_pats <- get_mut_pats(gene)}
  
  print(paste(length(mut_pats),"mutations loaded", sep =" "))
  
  # CNV <- getCNV(gene, cnv_arg)
  # CNV_pats <- CNV$patientID
  # print(paste(nrow(CNV),"CNV loaded", sep =" "))
  #muts_CNV <- unique(c(muts, CNV_pats))
  
  NON_mut_pats <- setdiff(pats, mut_pats)
  
  ### differential expression ###
  # create count df for mut vs control
  print(paste("performing differential expression analysis for", length(mut_pats), "patients", sep=" "))      
  mut_cnv_pats_df <- counts %>% dplyr::select(all_of(unlist(mut_pats)))
  NONmut_cnv_pats_df <- counts %>% dplyr::select(all_of(unlist(NON_mut_pats)))
  count_df <- cbind(mut_cnv_pats_df, NONmut_cnv_pats_df)
  
  # create group
  group <- factor(c(rep("T",length(mut_pats)), rep("C",length(NON_mut_pats))))
  
  # design
  design <- model.matrix(~ batch + group)
  
  dge <- DGEList(counts = count_df)
  
  # Filter lowly expressed genes
  keep <- filterByExpr(dge, group = group)
  dge <- dge[keep,,keep.lib.sizes = FALSE]
  
  # Normalize
  dge <- calcNormFactors(dge)
  
  # Estimate dispersion with batch + group design
  print("estimating dispersion")
  dge <- estimateDisp(dge, design)
  
  # Fit GLM
  fit <- glmQLFit(dge, design)
  
  # Run QL on the group coefficient (mutation vs control)
  qlf <- glmQLFTest(fit, coef = "groupT")
  
  # Extract DE genes
  DE <- topTags(qlf, n = nrow(dge$counts))
  res <- DE$table
  
  ### volcano plot
  res$gene <- rownames(res)
  # Add columns for plotting
  res$logFDR <- -log10(res$FDR)
  
  # Define significance thresholds
  logFC_cutoff <- 1.5
  pval_cutoff <- 0.01
  
  # Add a column for significance and direction
  res$Significant <- "NS"  # Not significant
  res$Significant[res$logFC >  logFC_cutoff & res$FDR < pval_cutoff] <- "Up"
  res$Significant[res$logFC < -logFC_cutoff & res$FDR < pval_cutoff] <- "Down"
  
  print(" differential analysis complete")
  
  # volcano plot with gene name labels
  
  if (length(gene) > 1){
    # Collapse into a single string
    gene_title <- paste(gene, collapse = ", ")}
  else{gene_title <- gene}
  
  sig_hits <- res[res$Significant != "NS", ]

  Volcano <- ggplot(res, aes(x = logFC, y = logFDR)) +
    geom_point(aes(color = Significant), alpha = 0.6) +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-1.5, 1.5), linetype = "dotted", color = "black") +
    scale_color_manual(values = c("Down" = "blue", "Up" = "red", "NS" = "gray")) +
    geom_text_repel(
      data = sig_hits,  # label only significant points
      aes(label = gene),                        # replace 'gene' with your label column
      max.overlaps = 15,
      size = 3,
      box.padding = 0.3,
      point.padding = 0.2,
      segment.color = 'grey50'
    ) +
    labs(
     title = paste("Differential expression for:", gene_title, "mut vs non-mut", sep=" "),
      x = "Log2 Fold Change",
      y = "-Log10 FDR"
    ) +
    theme_minimal()
  
  
  # Create ranked list (ENTREZ IDs recommended for many gene sets)
  set.seed(123)
  gene_list <- res$logFC
  names(gene_list) <- rownames(res)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  fgseaRes <- fgseaMultilevel(pathways = pathways,
                    stats = gene_list,
                    minSize = 15,
                    maxSize = 500)
  
  # View top results
  # head(fgseaRes[order(padj), ])
  
  # Assuming 'fgseaRes' is your fgsea result data.table
  topPathways <- fgseaRes[order(padj)][1:10]
  topPathways_df <- as.data.frame(topPathways)
  
  print("pre-ranked GSEA complete")  
  
  path_enrich_plot <- ggplot(topPathways, aes(x = reorder(pathway, NES), y = NES)) +
    geom_col(aes(fill = padj < 0.05)) +
    coord_flip() +
    labs(x = "Pathway", y = "Normalized Enrichment Score (NES)",
         title = paste("Top Enriched Pathways", gene, "DE", sep= " ")) +
    scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "gray")) +
    theme_minimal()
  
  
  return(list(gene = gene, DE_res = res, Vol = Volcano, GSEA = topPathways_df, GSEA_plot = path_enrich_plot))
  print("All analysis complete and results returned")
}


BRCA2_allres_qlf <- DNAvsRNA("BRCA2")

PALB2_allres <- DNAvsRNA("PALB2")

RAD51_allres_qlf <- DNAvsRNA(gene = c("RAD51", "RAD51B", "RAD51C", "RAD51D"))

BRCA2_PALB2_RAD51 <- DNAvsRNA(gene = c("BRCA2", "PALB2", "RAD51", "RAD51B", "RAD51C", "RAD51D"))



##########################################################################################################################


BRCA2_GSEA <- BRCA2_allres_qlf$GSEA
BRCA2_myog_LE <- unlist(BRCA2_GSEA[1,8])


#### GSEA from  SU2C data via cbioportal

##### you can do a very basic DE analysis on cbioportal

SU2C_2015_BRCA2 <- read.delim("~/SU2C/prad_su2c_2015/BRCA2_mut_vs_nonmut_DE/table.tsv")
sig <- subset(SU2C_2015_BRCA2, SU2C_2015_BRCA2$q.Value < 0.01)

head(SU2C_2015_BRCA2)

### convert log ratio to fold change
SU2C_2015_BRCA2$Log2.FC = 2^SU2C_2015_BRCA2$Log2.Ratio

# Create ranked list (ENTREZ IDs recommended for many gene sets)
set.seed(123)

gene_list <- SU2C_2015_BRCA2$Log2.FC
names(gene_list) <- SU2C_2015_BRCA2$Gene
gene_list <- sort(gene_list, decreasing = TRUE)

fgseaRes <- fgseaMultilevel(pathways = pathways,
                            stats = gene_list,
                            minSize = 15,
                            maxSize = 500)

topPathways <- fgseaRes[order(padj)][1:5]
topPathways_df <- as.data.frame(topPathways)

ggplot(topPathways, aes(x = reorder(pathway, NES), y = NES)) +
  geom_col(aes(fill = padj < 0.05)) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score (NES)",
       title = paste("Top Enriched Pathways", "BRCA2 SU2C 2015", "DE", sep= " ")) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "gray")) +
  theme_minimal()

