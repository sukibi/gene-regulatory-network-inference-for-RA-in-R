setwd("/Users/suki/Desktop/M2 ST4Health/ML for Biological networks/raproject")
rm(list = ls())
options(stringsAsFactors = FALSE)


if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_pkgs <- c("GEOquery", "limma", "WGCNA", "sva", "GENIE3", "minet")
cran_pkgs <- c("tidyverse", "igraph", "data.table")

for (p in c(bioc_pkgs, cran_pkgs)) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p %in% bioc_pkgs) BiocManager::install(p) else install.packages(p)
  }
  library(p, character.only = TRUE)
}
allowWGCNAThreads()

preprocess_geo <- function(gse) {
  expr <- exprs(gse)
  fdata <- fData(gse)
  sym_col <- grep("Symbol", colnames(fdata), ignore.case = TRUE, value = TRUE)[1]
  geneSymbols <- sapply(strsplit(as.character(fdata[[sym_col]]), " /// "), `[`, 1)
  
  keep <- !is.na(geneSymbols) & geneSymbols != ""
  expr <- expr[keep, ]
  expr_gene <- limma::avereps(expr, ID = geneSymbols[keep])
  v <- apply(expr_gene, 1, var)
  return(expr_gene[v > quantile(v, 0.5), ])
}


gse74143 <- getGEO("GSE74143", GSEMatrix = TRUE)[[1]]
gse93272 <- getGEO("GSE93272", GSEMatrix = TRUE)[[1]]

expr74143 <- preprocess_geo(gse74143)
expr93272 <- preprocess_geo(gse93272)

# Parse Phenotypes
p93272 <- pData(gse93272)
pheno93272 <- data.frame(
  sample = rownames(p93272),
  condition = ifelse(grepl("HC", p93272$characteristics_ch1), "Control", "RA"),
  row.names = rownames(p93272)
)

ra_93272_samples <- rownames(pheno93272)[pheno93272$condition == "RA"]
hc_93272_samples <- rownames(pheno93272)[pheno93272$condition == "Control"]

message("GSE93272: ", length(ra_93272_samples), " RA + ", length(hc_93272_samples), " HC")
message("GSE74143: ", ncol(expr74143), " RA samples")

# hub discovery
datExpr <- t(expr93272)
powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
softPower <- sft$fitIndices$Power[which(sft$fitIndices$SFT.R.sq >= 0.85)[1]]
if (is.na(softPower)) softPower <- 6

net <- blockwiseModules(datExpr, power = softPower, TOMType = "unsigned", 
                        minModuleSize = 30, numericLabels = TRUE, verbose = 0)

# Identify Module correlated with RA
condition_bin <- as.numeric(factor(pheno93272$condition, levels=c("Control", "RA"))) - 1
moduleTraitCor <- cor(net$MEs, condition_bin, use = "p")
best_mod_idx <- which.max(abs(moduleTraitCor))
module_label <- colnames(net$MEs)[best_mod_idx]

# Extract Top 200 Hubs
datKME <- signedKME(datExpr, net$MEs)
kME_col <- paste0("kME", gsub("ME", "", module_label))
if(!kME_col %in% colnames(datKME)) kME_col <- paste0("k", module_label)

hub_values <- datKME[[kME_col]]
names(hub_values) <- rownames(datKME)
top_hubs <- names(sort(abs(hub_values), decreasing = TRUE))[1:200]

message("Selected ", length(top_hubs), " hub genes from module ", module_label)


# Batch correction

common_genes <- intersect(intersect(rownames(expr93272), rownames(expr74143)), top_hubs)
message("Common genes across datasets: ", length(common_genes))

# --- CORRECT RA SAMPLES ONLY ---
ra_93272_expr <- expr93272[common_genes, ra_93272_samples]
ra_74143_expr <- expr74143[common_genes, ]

ra_merged <- cbind(ra_93272_expr, ra_74143_expr)
ra_batch <- c(rep("GSE93272", ncol(ra_93272_expr)), rep("GSE74143", ncol(ra_74143_expr)))

message("Correcting batch effects in RA samples...")
ra_corrected <- ComBat(dat = ra_merged, batch = ra_batch, mod = NULL)

# --- HEALTHY CONTROLS (No batch correction needed - single dataset) ---
hc_expr <- expr93272[common_genes, hc_93272_samples]

message("RA mega-cohort: ", ncol(ra_corrected), " samples")
message("HC cohort: ", ncol(hc_expr), " samples")


genie_ra_mega <- GENIE3(as.matrix(ra_corrected))
links_genie_ra <- getLinkList(genie_ra_mega)
write.table(links_genie_ra, "Links_RA_MegaCohort_GENIE3.txt", sep="\t", row.names=F, quote=F)

message("Running ARACNe on RA Mega-Cohort...")
mim_ra <- build.mim(t(ra_corrected), estimator = "spearman")
net_aracne_ra <- aracne(mim_ra, eps = 0.1)
links_aracne_ra <- as.data.frame(as.table(net_aracne_ra))
colnames(links_aracne_ra) <- c("Regulator", "Target", "Weight")
links_aracne_ra <- links_aracne_ra[links_aracne_ra$Weight > 0, ]
write.table(links_aracne_ra, "Links_RA_MegaCohort_ARACNe.txt", sep="\t", row.names=F, quote=F)


genie_hc <- GENIE3(as.matrix(hc_expr))
links_genie_hc <- getLinkList(genie_hc)
write.table(links_genie_hc, "Links_HC_GENIE3.txt", sep="\t", row.names=F, quote=F)

message("Running ARACNe on Healthy Controls...")
mim_hc <- build.mim(t(hc_expr), estimator = "spearman")
net_aracne_hc <- aracne(mim_hc, eps = 0.1)
links_aracne_hc <- as.data.frame(as.table(net_aracne_hc))
colnames(links_aracne_hc) <- c("Regulator", "Target", "Weight")
links_aracne_hc <- links_aracne_hc[links_aracne_hc$Weight > 0, ]
write.table(links_aracne_hc, "Links_HC_ARACNe.txt", sep="\t", row.names=F, quote=F)

#export for mIIC
# RA Mega-Cohort
write.csv(as.data.frame(t(ra_corrected)), 
          "MIIC_RA_MegaCohort_Input.csv", 
          row.names = FALSE)
message("Exported MIIC_RA_MegaCohort_Input.csv: ", ncol(ra_corrected), 
        " samples x ", nrow(ra_corrected), " genes")

# Healthy Controls
write.csv(as.data.frame(t(hc_expr)), 
          "MIIC_HC_Input.csv", 
          row.names = FALSE)
message("Exported MIIC_HC_Input.csv: ", ncol(hc_expr), 
        " samples x ", nrow(hc_expr), " genes")


# verifying batch correction optional block

# PCA before correction
ra_merged_uncorrected <- cbind(ra_93272_expr, ra_74143_expr)
pca_before <- prcomp(t(ra_merged_uncorrected), scale. = TRUE)
pca_before_df <- data.frame(
  PC1 = pca_before$x[,1],
  PC2 = pca_before$x[,2],
  Dataset = ra_batch
)

# PCA after correction
pca_after <- prcomp(t(ra_corrected), scale. = TRUE)
pca_after_df <- data.frame(
  PC1 = pca_after$x[,1],
  PC2 = pca_after$x[,2],
  Dataset = ra_batch
)

# Save plots
pdf("Batch_Correction_Verification.pdf", width = 12, height = 5)
par(mfrow = c(1, 2))
plot(pca_before_df$PC1, pca_before_df$PC2, 
     col = ifelse(pca_before_df$Dataset == "GSE93272", "red", "blue"),
     pch = 19, main = "Before Batch Correction",
     xlab = "PC1", ylab = "PC2")
legend("topright", legend = c("GSE93272", "GSE74143"), 
       col = c("red", "blue"), pch = 19)

plot(pca_after_df$PC1, pca_after_df$PC2, 
     col = ifelse(pca_after_df$Dataset == "GSE93272", "red", "blue"),
     pch = 19, main = "After Batch Correction",
     xlab = "PC1", ylab = "PC2")
legend("topright", legend = c("GSE93272", "GSE74143"), 
       col = c("red", "blue"), pch = 19)
dev.off()

message("Batch correction verification plot saved")



g_ra <- graph_from_data_frame(links_genie_ra[1:300, 1:2], directed=TRUE)
g_hc <- graph_from_data_frame(links_genie_hc[1:300, 1:2], directed=TRUE)

g_disease_specific <- difference(g_ra, g_hc)
message("Disease-specific edges (RA but not HC): ", ecount(g_disease_specific))

top_regulators <- sort(degree(g_disease_specific, mode = "out"), decreasing = TRUE)
message("\nTop 10 Preliminary Master Regulators:")
print(head(top_regulators, 10))


saveRDS(list(
  common_genes = common_genes,
  ra_corrected_colnames = colnames(ra_corrected),
  hc_colnames = colnames(hc_expr),
  batch_info = ra_batch
), "megacohort_metadata.rds")


message("Run rscriptpart4_2bconsensus for consensus integration")
message("\nKey advantage: ", ncol(ra_corrected), " RA samples vs ", 
        length(ra_93272_samples), " in original approach")