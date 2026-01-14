

setwd("/Users/suki/Desktop/M2 ST4Health/ML for Biological networks/raproject")
rm(list = ls())
options(stringsAsFactors = FALSE)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_pkgs <- c("GEOquery", "limma", "WGCNA", "GENIE3", "minet")
cran_pkgs <- c("data.table", "tidyverse", "igraph")

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
  geneSymbols <- sapply(strsplit(as.character(fdata$`Gene Symbol`), " /// "), `[`, 1)
  geneSymbols <- trimws(geneSymbols)
  keep <- !is.na(geneSymbols) & geneSymbols != ""
  expr <- expr[keep, ]
  expr_gene <- limma::avereps(expr, ID = geneSymbols[keep])
  # Filter for common variance (top 50%)
  v <- apply(expr_gene, 1, var)
  return(expr_gene[v > quantile(v, 0.5), ])
}


gse74143 <- getGEO("GSE74143", GSEMatrix = TRUE)[[1]]
gse93272 <- getGEO("GSE93272", GSEMatrix = TRUE)[[1]]

expr74143 <- preprocess_geo(gse74143)
expr93272 <- preprocess_geo(gse93272)

message("GSE93272 dimensions: ", nrow(expr93272), " genes x ", ncol(expr93272), " samples")
message("GSE74143 dimensions: ", nrow(expr74143), " genes x ", ncol(expr74143), " samples")


# GSE93272 has both RA and HC; GSE74143 is all RA
p93272 <- pData(gse93272)
pheno93272 <- data.frame(
  sample = rownames(p93272),
  condition = ifelse(grepl("HC", p93272$characteristics_ch1), "Control", "RA"),
  row.names = rownames(p93272)
)

ra_93272_samples <- rownames(pheno93272)[pheno93272$condition == "RA"]
hc_93272_samples <- rownames(pheno93272)[pheno93272$condition == "Control"]

message("GSE93272 RA samples: ", length(ra_93272_samples))
message("GSE93272 HC samples: ", length(hc_93272_samples))
message("GSE74143 RA samples: ", ncol(expr74143))

# Check if we have enough samples for analysis
if(length(ra_93272_samples) < 10 | length(hc_93272_samples) < 10) {
  warning("Sample sizes may be too small for robust network inference!")
}

message("\nRunning WGCNA for module detection...")
run_wgcna <- function(expr) {
  datExpr <- t(expr)
  powers <- 1:20
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
  softPower <- sft$fitIndices$Power[which(sft$fitIndices$SFT.R.sq >= 0.85)[1]]
  if (is.na(softPower)) softPower <- 6
  message("Using soft threshold power: ", softPower)
  net <- blockwiseModules(datExpr, power = softPower, TOMType = "unsigned", 
                          minModuleSize = 30, numericLabels = TRUE, verbose = 0)
  return(list(net = net, datExpr = datExpr))
}

wgcna_discovery <- run_wgcna(expr93272)

condition_bin <- as.numeric(factor(pheno93272$condition, levels=c("Control", "RA"))) - 1
moduleTraitCor <- cor(wgcna_discovery$net$MEs, condition_bin, use = "p")
best_module <- which.max(abs(moduleTraitCor))
module_label <- colnames(wgcna_discovery$net$MEs)[best_module]
message("\nTarget Module: ", module_label)
message("Module-Disease Correlation: ", round(moduleTraitCor[best_module], 3))


message("\nIdentifying hub genes by module connectivity...")
datKME <- signedKME(wgcna_discovery$datExpr, wgcna_discovery$net$MEs)

# Get the correct kME column name
kME_target_col <- paste0("kME", gsub("ME", "", module_label))

if (!kME_target_col %in% colnames(datKME)) {
  # Try alternative naming
  kME_target_col <- paste0("k", module_label)
}

if (!kME_target_col %in% colnames(datKME)) {
  stop("Cannot find kME column. Available columns: ", paste(colnames(datKME), collapse=", "))
}

message("Using connectivity column: ", kME_target_col)

# Extract top 200 hubs by absolute connectivity
hub_values <- datKME[[kME_target_col]]
names(hub_values) <- rownames(datKME)
top_hubs <- names(sort(abs(hub_values), decreasing = TRUE))[1:200]

message("Selected ", length(top_hubs), " hub genes")


valid_hubs <- intersect(top_hubs, rownames(expr74143))
message("Valid hubs across both datasets: ", length(valid_hubs))

if(length(valid_hubs) < 50) {
  warning("Only ", length(valid_hubs), " hubs overlap between datasets. Consider relaxing filters.")
}

# Save for Script 2
saveRDS(valid_hubs, "valid_hubs.rds")
saveRDS(list(ra = ra_93272_samples, hc = hc_93272_samples), "sample_ids.rds")

#genie3

message("GENIE3 for RA samples (GSE93272)...")
links_RA_1 <- getLinkList(GENIE3(expr93272[valid_hubs, ra_93272_samples]))
write.table(links_RA_1, "Links_GSE93272_RA.txt", sep="\t", row.names=F, quote=F)

message("GENIE3 for Healthy samples (GSE93272)...")
links_HC <- getLinkList(GENIE3(expr93272[valid_hubs, hc_93272_samples]))
write.table(links_HC, "Links_GSE93272_Healthy.txt", sep="\t", row.names=F, quote=F)

message("GENIE3 for RA validation (GSE74143)...")
links_RA_2 <- getLinkList(GENIE3(expr74143[valid_hubs, ]))
write.table(links_RA_2, "Links_GSE74143_RA_Validation.txt", sep="\t", row.names=F, quote=F)


#aracne inference

run_aracne <- function(expr_matrix, hubs, sample_ids) {
  data_subset <- t(expr_matrix[hubs, sample_ids])
  mim <- build.mim(data_subset, estimator = "spearman")
  net_aracne <- aracne(mim, eps = 0.1)
  links <- as.data.frame(as.table(net_aracne))
  colnames(links) <- c("Regulator", "Target", "Weight")
  links <- links[links$Weight > 0, ]
  links <- links[order(links$Weight, decreasing = TRUE), ]
  return(links)
}

message("ARACNe for RA samples (GSE93272)...")
links_ARACNE_RA1 <- run_aracne(expr93272, valid_hubs, ra_93272_samples)
write.table(links_ARACNE_RA1, "ARACNe_GSE93272_RA.txt", sep="\t", row.names=F, quote=F)

message("ARACNe for Healthy samples (GSE93272)...")
links_ARACNE_HC <- run_aracne(expr93272, valid_hubs, hc_93272_samples)
write.table(links_ARACNE_HC, "ARACNe_GSE93272_HC.txt", sep="\t", row.names=F, quote=F)

message("ARACNe for RA validation (GSE74143)...")
links_ARACNE_RA2 <- run_aracne(expr74143, valid_hubs, colnames(expr74143))
write.table(links_ARACNE_RA2, "ARACNe_GSE74143_RA_Validation.txt", sep="\t", row.names=F, quote=F)


# Export 1: RA samples (Discovery)
miic_RA_discovery <- as.data.frame(t(expr93272[valid_hubs, ra_93272_samples]))
write.csv(miic_RA_discovery, "MIIC_GSE93272_kme_RA.csv", row.names = FALSE)
message("Exported MIIC_GSE93272_kme_RA.csv: ", nrow(miic_RA_discovery), " samples x ", 
        ncol(miic_RA_discovery), " genes")

# Export 2: Healthy samples (Discovery)
miic_HC_discovery <- as.data.frame(t(expr93272[valid_hubs, hc_93272_samples]))
write.csv(miic_HC_discovery, "MIIC_GSE93272_kme_Healthy.csv", row.names = FALSE)
message("Exported MIIC_GSE93272_kme_Healthy.csv: ", nrow(miic_HC_discovery), " samples x ", 
        ncol(miic_HC_discovery), " genes")

# Export 3: RA samples (Validation)
miic_RA_validation <- as.data.frame(t(expr74143[valid_hubs, ]))
write.csv(miic_RA_validation, "MIIC_GSE74143_kme_RA.csv", row.names = FALSE)
message("Exported MIIC_GSE74143_kme_RA.csv: ", nrow(miic_RA_validation), " samples x ", 
        ncol(miic_RA_validation), " genes")



# Load top edges for quick comparison
g_RA <- graph_from_data_frame(links_RA_1[1:300, 1:2], directed=TRUE)
g_HC <- graph_from_data_frame(links_HC[1:300, 1:2], directed=TRUE)
g_Val <- graph_from_data_frame(links_RA_2[1:300, 1:2], directed=TRUE)

# Disease-specific links (in RA but NOT in Healthy)
g_disease_specific <- difference(g_RA, g_HC)
message("Disease-specific edges: ", ecount(g_disease_specific))

# Validated links (appear in BOTH RA datasets)
g_robust_RA <- intersection(g_RA, g_Val)
message("Cross-validated edges: ", ecount(g_robust_RA))

# Preliminary master regulators
top_regulators <- sort(degree(g_disease_specific, mode = "out"), decreasing = TRUE)
message("\nTop 10 Preliminary Master Regulators (GENIE3 only):")
print(head(top_regulators, 10))



message("Run Script rscriptpart2_1bconsensus for  consensus integration")