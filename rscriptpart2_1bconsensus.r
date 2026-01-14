
setwd("/Users/suki/Desktop/M2 ST4Health/ML for Biological networks/raproject")
rm(list = ls())
options(stringsAsFactors = FALSE)


library(tidyverse)
library(igraph)

if (!requireNamespace("clusterProfiler", quietly=TRUE)) {
  BiocManager::install("clusterProfiler")
}
library(clusterProfiler)
library(org.Hs.eg.db)

valid_hubs <- readRDS("valid_hubs.rds")
sample_info <- readRDS("sample_ids.rds")

message("Working with ", length(valid_hubs), " validated hub genes")
message("GSE93272 RA samples: ", length(sample_info$ra))
message("GSE93272 HC samples: ", length(sample_info$hc))


# Discovery RA
genie_ra_disc <- read.table("Links_GSE93272_RA.txt", header=T)
# Force column names regardless of what they are
colnames(genie_ra_disc) <- c("Regulator", "Target", "Genie_Weight_RA")
genie_ra_disc <- genie_ra_disc %>% filter(Regulator %in% valid_hubs, Target %in% valid_hubs)

# Discovery Healthy
genie_hc_disc <- read.table("Links_GSE93272_Healthy.txt", header=T)
colnames(genie_hc_disc) <- c("Regulator", "Target", "Genie_Weight_HC")
genie_hc_disc <- genie_hc_disc %>% filter(Regulator %in% valid_hubs, Target %in% valid_hubs)

# Validation RA
genie_ra_val <- read.table("Links_GSE74143_RA_Validation.txt", header=T)
colnames(genie_ra_val) <- c("Regulator", "Target", "Genie_Weight_Val")
genie_ra_val <- genie_ra_val %>% filter(Regulator %in% valid_hubs, Target %in% valid_hubs)

message("GENIE3 RA discovery edges: ", nrow(genie_ra_disc))
message("GENIE3 HC discovery edges: ", nrow(genie_hc_disc))
message("GENIE3 RA validation edges: ", nrow(genie_ra_val))


# Discovery RA
aracne_ra_disc <- read.table("ARACNe_GSE93272_RA.txt", header=T)
colnames(aracne_ra_disc) <- c("Regulator", "Target", "Aracne_MI_RA")
aracne_ra_disc <- aracne_ra_disc %>% filter(Regulator %in% valid_hubs, Target %in% valid_hubs)

# Discovery Healthy
aracne_hc_disc <- read.table("ARACNe_GSE93272_HC.txt", header=T)
colnames(aracne_hc_disc) <- c("Regulator", "Target", "Aracne_MI_HC")
aracne_hc_disc <- aracne_hc_disc %>% filter(Regulator %in% valid_hubs, Target %in% valid_hubs)

# Validation RA
aracne_ra_val <- read.table("ARACNe_GSE74143_RA_Validation.txt", header=T)
colnames(aracne_ra_val) <- c("Regulator", "Target", "Aracne_MI_Val")
aracne_ra_val <- aracne_ra_val %>% filter(Regulator %in% valid_hubs, Target %in% valid_hubs)

message("ARACNe RA discovery edges: ", nrow(aracne_ra_disc))
message("ARACNe HC discovery edges: ", nrow(aracne_hc_disc))
message("ARACNe RA validation edges: ", nrow(aracne_ra_val))


# Function to process MIIC output
process_miic <- function(filepath, edge_label) {
  if(!file.exists(filepath)) {
    warning("MIIC file not found: ", filepath)
    return(NULL)
  }
  
  miic_raw <- read.table(filepath, header=T, sep="\t")
  
  # Diagnostic check (first run only)
  if(edge_label == "RA_Discovery") {
 
    message("First 5 edges from RA Discovery:")
    print(head(miic_raw[, c("x", "y", "ort_inferred", "type")], 5))
    message("ort_inferred encoding: positive = x->y, negative = y->x")
   
  }
  
  miic_processed <- miic_raw %>%
    filter(type == "P") %>%  # Only prospective/confirmed edges
    mutate(
      # Adjust direction based on ort_inferred
      # Positive values: x ->y
      # Negative values: y->x
      Final_Regulator = ifelse(ort_inferred > 0, x, y),
      Final_Target = ifelse(ort_inferred > 0, y, x),
      MIIC_Confidence = info_shifted_bits,
      Edge_Label = edge_label
    )
  
  # Use dplyr::select explicitly
  miic_processed <- miic_processed %>%
    dplyr::select(Regulator = Final_Regulator, 
                  Target = Final_Target, 
                  MIIC_Confidence,
                  is_causal,
                  Edge_Label) %>%
    filter(Regulator %in% valid_hubs, Target %in% valid_hubs)
  
  return(miic_processed)
}

# Load MIIC results
miic_ra_disc <- process_miic("Final_MIIC_GSE93272_kme_RA_edgesList.miic.summary.txt", "RA_Discovery")
miic_hc_disc <- process_miic("Final_MIIC_GSE93272_kme_Healthy_edgesList.miic.summary.txt", "HC_Discovery")
miic_ra_val <- process_miic("Final_MIIC_GSE74143_kme_RA_edgesList.miic.summary.txt", "RA_Validation")

if(!is.null(miic_ra_disc)) message("MIIC RA discovery edges: ", nrow(miic_ra_disc))
if(!is.null(miic_hc_disc)) message("MIIC HC discovery edges: ", nrow(miic_hc_disc))
if(!is.null(miic_ra_val)) message("MIIC RA validation edges: ", nrow(miic_ra_val))


# Merge RA discovery networks
ra_discovery_edges <- genie_ra_disc %>%
  full_join(aracne_ra_disc, by = c("Regulator", "Target"))

if(!is.null(miic_ra_disc)) {
  # FIX: Use dplyr::select explicitly to avoid igraph conflict
  ra_discovery_edges <- ra_discovery_edges %>%
    full_join(miic_ra_disc %>% dplyr::select(Regulator, Target, MIIC_Confidence, is_causal), 
              by = c("Regulator", "Target"))
}
# Create HC edge identifiers for filtering
hc_edges_genie <- genie_hc_disc %>%
  filter(Genie_Weight_HC > quantile(Genie_Weight_HC, 0.75)) %>%  # Top 25% of HC edges
  mutate(edge_id = paste(Regulator, Target, sep = "->"))

hc_edges_aracne <- aracne_hc_disc %>%
  filter(Aracne_MI_HC > quantile(Aracne_MI_HC, 0.75)) %>%
  mutate(edge_id = paste(Regulator, Target, sep = "->"))

hc_edges_miic <- NULL
if(!is.null(miic_hc_disc)) {
  hc_edges_miic <- miic_hc_disc %>%
    mutate(edge_id = paste(Regulator, Target, sep = "->"))
}

# Filter out edges strong in healthy controls
disease_specific_edges <- ra_discovery_edges %>%
  mutate(edge_id = paste(Regulator, Target, sep = "->")) %>%
  filter(
    !edge_id %in% hc_edges_genie$edge_id,
    !edge_id %in% hc_edges_aracne$edge_id
  )

if(!is.null(hc_edges_miic)) {
  disease_specific_edges <- disease_specific_edges %>%
    filter(!edge_id %in% hc_edges_miic$edge_id)
}

# Count method agreement
disease_specific_edges <- disease_specific_edges %>%
  rowwise() %>%
  mutate(
    Methods_Count = sum(!is.na(Genie_Weight_RA), 
                        !is.na(Aracne_MI_RA), 
                        !is.na(MIIC_Confidence))
  ) %>%
  ungroup() %>%
  filter(Methods_Count >= 2) %>%  # Require at least 2 methods
  arrange(desc(Methods_Count), desc(Genie_Weight_RA))

message("Disease-specific edges (2+ methods, absent in HC): ", nrow(disease_specific_edges))


# Check which disease-specific edges replicate in validation dataset
validated_edges <- disease_specific_edges %>%
  mutate(
    Validated_GENIE = paste(Regulator, Target, sep = "->") %in% 
      paste(genie_ra_val$Regulator, genie_ra_val$Target, sep = "->"),
    Validated_ARACNE = paste(Regulator, Target, sep = "->") %in% 
      paste(aracne_ra_val$Regulator, aracne_ra_val$Target, sep = "->")
  )

if(!is.null(miic_ra_val)) {
  validated_edges <- validated_edges %>%
    mutate(
      Validated_MIIC = paste(Regulator, Target, sep = "->") %in% 
        paste(miic_ra_val$Regulator, miic_ra_val$Target, sep = "->")
    ) %>%
    rowwise() %>%
    mutate(
      Validation_Score = sum(Validated_GENIE, Validated_ARACNE, Validated_MIIC, na.rm=TRUE)
    ) %>%
    ungroup()
} else {
  validated_edges <- validated_edges %>%
    rowwise() %>%
    mutate(
      Validation_Score = sum(Validated_GENIE, Validated_ARACNE, na.rm=TRUE)
    ) %>%
    ungroup()
}

# Keep only edges validated in at least one method
final_consensus <- validated_edges %>%
  filter(Validation_Score >= 1) %>%
  arrange(desc(Validation_Score), desc(Methods_Count), desc(Genie_Weight_RA))

message("Validated disease-specific edges: ", nrow(final_consensus))




g_final <- graph_from_data_frame(
  final_consensus %>% dplyr::select(Regulator, Target), 
  directed = TRUE
)

# Calculate hub metrics
hub_out_degree <- degree(g_final, mode = "out")
hub_in_degree <- degree(g_final, mode = "in")
hub_betweenness <- betweenness(g_final, directed = TRUE)

# Compile master regulator table
master_regulators <- data.frame(
  Gene = V(g_final)$name,
  OutDegree = hub_out_degree[V(g_final)$name],
  InDegree = hub_in_degree[V(g_final)$name],
  Betweenness = hub_betweenness[V(g_final)$name]
) %>%
  mutate(
    Total_Degree = OutDegree + InDegree,
    # Normalize betweenness to prevent it from overwhelming the score
    Regulator_Score = (OutDegree * 2) + (hub_betweenness / max(hub_betweenness + 1) * 10)
  ) %>%
  arrange(desc(OutDegree), desc(Betweenness)) %>%
  mutate(Rank = row_number())


write.csv(final_consensus, "RA_Final_Consensus_Network.csv", row.names=F)
write.csv(master_regulators, "RA_Master_Regulators_Final.csv", row.names=F)
write.csv(disease_specific_edges, "RA_Disease_Specific_Edges_All.csv", row.names=F)

message("Exported:")
message("  - RA_Final_Consensus_Network.csv (", nrow(final_consensus), " edges)")
message("  - RA_Master_Regulators_Final.csv (", nrow(master_regulators), " genes)")
message("  - RA_Disease_Specific_Edges_All.csv (", nrow(disease_specific_edges), " edges)")

message("\n TOP 20 MASTER REGULATORS")
print(master_regulators %>% 
        dplyr::select(Rank, Gene, OutDegree, InDegree, Betweenness, Regulator_Score) %>%
        head(20))


library(ggplot2)

# Plot Top 20 Master Regulators
top_20_plot <- ggplot(head(master_regulators, 20), aes(x = reorder(Gene, OutDegree), y = OutDegree, fill = OutDegree)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "skyblue", high = "royalblue4") +
  labs(title = "Top 20 Master Regulators in RA",
       subtitle = "Ranked by Out-Degree (Number of Target Genes)",
       x = "Gene Symbol", y = "Out-Degree") +
  theme_minimal()

ggsave("RA_Top20_MasterRegulators_Barplot.pdf", plot = top_20_plot, width = 8, height = 6)


top_genes <- master_regulators$Gene[1:min(50, nrow(master_regulators))]
valid_genes <- top_genes[top_genes %in% keys(org.Hs.eg.db, keytype = "SYMBOL")]

message("Analyzing ", length(valid_genes), " genes for GO enrichment...")

if(length(valid_genes) >= 5) {
  tryCatch({
    # Biological Process
    ego_bp <- enrichGO(
      gene = valid_genes, 
      OrgDb = org.Hs.eg.db, 
      keyType = "SYMBOL", 
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
    
    if(nrow(ego_bp@result) > 0) {
      message("Found ", nrow(ego_bp@result), " enriched BP terms")
      write.csv(ego_bp@result, "RA_MasterRegulators_GO_BP.csv", row.names=F)
      
      pdf("RA_MasterRegulators_GO_BP_Dotplot.pdf", width=10, height=8)
      print(dotplot(ego_bp, showCategory=20))
      dev.off()
      message("Saved GO BP enrichment plot")
    } else {
      message("No significant GO BP terms found")
    }
    
    # Molecular Function
    ego_mf <- enrichGO(
      gene = valid_genes, 
      OrgDb = org.Hs.eg.db, 
      keyType = "SYMBOL", 
      ont = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05
    )
    
    if(nrow(ego_mf@result) > 0) {
      write.csv(ego_mf@result, "RA_MasterRegulators_GO_MF.csv", row.names=F)
    }
    
  }, error = function(e) {
    message("GO enrichment error: ", e$message)
  })
} else {
  message("Insufficient valid gene symbols for enrichment (need >= 5)")
}


message("Total nodes: ", vcount(g_final))
message("Total edges: ", ecount(g_final))
message("Network density: ", round(edge_density(g_final), 4))
message("Average out-degree: ", round(mean(hub_out_degree), 2))
message("Genes with OutDegree >= 5: ", sum(hub_out_degree >= 5))
message("Genes with OutDegree >= 10: ", sum(hub_out_degree >= 10))

# Identify hub categories
tier1_hubs <- master_regulators %>% filter(OutDegree >= 10) %>% pull(Gene)
tier2_hubs <- master_regulators %>% filter(OutDegree >= 5, OutDegree < 10) %>% pull(Gene)

message("\nTier 1 Master Regulators (≥10 targets): ", length(tier1_hubs))
if(length(tier1_hubs) > 0) print(tier1_hubs)

message("\nTier 2 Master Regulators (5-9 targets): ", length(tier2_hubs))
if(length(tier2_hubs) > 0) print(tier2_hubs)

