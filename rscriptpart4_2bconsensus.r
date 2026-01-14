

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


meta <- readRDS("megacohort_metadata.rds")
common_genes <- meta$common_genes

message("Working with ", length(common_genes), " genes across mega-cohort")
message("RA mega-cohort samples: ", length(meta$ra_corrected_colnames))
message("HC samples: ", length(meta$hc_colnames))

# RA Mega-Cohort
genie_ra_mega <- read.table("Links_RA_MegaCohort_GENIE3.txt", header=T)
colnames(genie_ra_mega) <- c("Regulator", "Target", "Genie_Weight_RA")
genie_ra_mega <- genie_ra_mega %>%
  filter(Regulator %in% common_genes, Target %in% common_genes)

# Healthy Controls
genie_hc <- read.table("Links_HC_GENIE3.txt", header=T)
colnames(genie_hc) <- c("Regulator", "Target", "Genie_Weight_HC")
genie_hc <- genie_hc %>%
  filter(Regulator %in% common_genes, Target %in% common_genes)

message("GENIE3 RA mega-cohort edges: ", nrow(genie_ra_mega))
message("GENIE3 HC edges: ", nrow(genie_hc))



# RA Mega-Cohort
aracne_ra_mega <- read.table("Links_RA_MegaCohort_ARACNe.txt", header=T)
colnames(aracne_ra_mega) <- c("Regulator", "Target", "Aracne_MI_RA")
aracne_ra_mega <- aracne_ra_mega %>%
  filter(Regulator %in% common_genes, Target %in% common_genes)

# Healthy Controls
aracne_hc <- read.table("Links_HC_ARACNe.txt", header=T)
colnames(aracne_hc) <- c("Regulator", "Target", "Aracne_MI_HC")
aracne_hc <- aracne_hc %>%
  filter(Regulator %in% common_genes, Target %in% common_genes)

message("ARACNe RA mega-cohort edges: ", nrow(aracne_ra_mega))
message("ARACNe HC edges: ", nrow(aracne_hc))


# Function to process MIIC output
process_miic <- function(filepath, edge_label) {
  if(!file.exists(filepath)) {
    warning("MIIC file not found: ", filepath, " - Skipping MIIC integration")
    return(NULL)
  }
  
  miic_raw <- read.table(filepath, header=T, sep="\t")
  
  # Diagnostic check (first run only)
  if(edge_label == "RA_MegaCohort") {
    message("\n=== MIIC DIRECTION DIAGNOSTIC ===")
    message("First 5 edges from RA Mega-Cohort:")
    print(head(miic_raw[, c("x", "y", "ort_inferred", "type")], 5))
  }
  
  miic_processed <- miic_raw %>%
    filter(type == "P") %>%  # Only prospective/confirmed edges
    mutate(
      
      Final_Regulator = ifelse(ort_inferred > 0, x, y),
      Final_Target = ifelse(ort_inferred > 0, y, x),
      MIIC_Confidence = info_shifted_bits,
      Edge_Label = edge_label
    )
  
 
  miic_processed <- miic_processed %>%
    dplyr::select(Regulator = Final_Regulator, 
                  Target = Final_Target, 
                  MIIC_Confidence,
                  is_causal,
                  Edge_Label) %>%
    filter(Regulator %in% common_genes, Target %in% common_genes)
  
  return(miic_processed)
}

# Load MIIC results
miic_ra_mega <- process_miic("Final_MIIC_RA_MegaCohort_edgesList.miic.summary.txt", "RA_MegaCohort")
miic_hc <- process_miic("Final_MIIC_HC_Input_edgesList.miic.summary.txt", "HC")

if(!is.null(miic_ra_mega)) message("MIIC RA mega-cohort edges: ", nrow(miic_ra_mega))
if(!is.null(miic_hc)) message("MIIC HC edges: ", nrow(miic_hc))


# Start with GENIE3 as base (typically has most edges)
ra_mega_edges <- genie_ra_mega %>%
  full_join(aracne_ra_mega, by = c("Regulator", "Target"))

# Add MIIC if available
if(!is.null(miic_ra_mega)) {
  ra_mega_edges <- ra_mega_edges %>%
    full_join(miic_ra_mega %>% dplyr::select(Regulator, Target, MIIC_Confidence, is_causal), 
              by = c("Regulator", "Target"))
}

# Count method agreement
ra_mega_edges <- ra_mega_edges %>%
  rowwise() %>%
  mutate(
    Methods_Count = sum(!is.na(Genie_Weight_RA), 
                        !is.na(Aracne_MI_RA), 
                        !is.na(MIIC_Confidence))
  ) %>%
  ungroup()

# Filter: Require at least 2 methods agree
ra_consensus <- ra_mega_edges %>%
  filter(Methods_Count >= 2) %>%
  arrange(desc(Methods_Count), desc(Genie_Weight_RA))

message("RA mega-cohort consensus edges (2+ methods): ", nrow(ra_consensus))


# We want to remove edges that are STRONG in healthy controls
hc_edges_genie <- genie_hc %>%
  filter(Genie_Weight_HC > quantile(Genie_Weight_HC, 0.75)) %>%  # Top 25% HC edges
  mutate(edge_id = paste(Regulator, Target, sep = "->"))

hc_edges_aracne <- aracne_hc %>%
  filter(Aracne_MI_HC > quantile(Aracne_MI_HC, 0.75)) %>%
  mutate(edge_id = paste(Regulator, Target, sep = "->"))

hc_edges_miic <- NULL
if(!is.null(miic_hc)) {
  hc_edges_miic <- miic_hc %>%
    mutate(edge_id = paste(Regulator, Target, sep = "->"))
}

# Filter out edges strong in healthy controls
disease_specific_edges <- ra_consensus %>%
  mutate(edge_id = paste(Regulator, Target, sep = "->")) %>%
  filter(
    !edge_id %in% hc_edges_genie$edge_id,
    !edge_id %in% hc_edges_aracne$edge_id
  )

if(!is.null(hc_edges_miic)) {
  disease_specific_edges <- disease_specific_edges %>%
    filter(!edge_id %in% hc_edges_miic$edge_id)
}

# Tag edges with their disease-specificity status
disease_specific_edges <- disease_specific_edges %>%
  mutate(
    Disease_Specific = TRUE,
    HC_Evidence = "Absent in HC (top 25%)"
  )

message("Disease-specific edges (present in RA, weak/absent in HC): ", 
        nrow(disease_specific_edges))

# Also create a table showing which edges ARE present in both (for comparison)
shared_edges <- ra_consensus %>%
  mutate(edge_id = paste(Regulator, Target, sep = "->")) %>%
  filter(
    edge_id %in% hc_edges_genie$edge_id | edge_id %in% hc_edges_aracne$edge_id
  ) %>%
  mutate(
    Disease_Specific = FALSE,
    HC_Evidence = "Also strong in HC"
  )

message("Shared edges (present in both RA and HC): ", nrow(shared_edges))


all_edges_annotated <- bind_rows(disease_specific_edges, shared_edges) %>%
  arrange(desc(Disease_Specific), desc(Methods_Count), desc(Genie_Weight_RA))


# Tier 1: Disease-Specific Master Regulators
g_disease <- graph_from_data_frame(
  disease_specific_edges %>% dplyr::select(Regulator, Target), 
  directed = TRUE
)

disease_out_degree <- degree(g_disease, mode = "out")
disease_in_degree <- degree(g_disease, mode = "in")
disease_betweenness <- betweenness(g_disease, directed = TRUE)

master_regulators_disease <- data.frame(
  Gene = names(disease_out_degree),
  OutDegree_Disease = disease_out_degree,
  InDegree_Disease = disease_in_degree,
  Betweenness_Disease = disease_betweenness,
  stringsAsFactors = FALSE
) %>%
  arrange(desc(OutDegree_Disease), desc(Betweenness_Disease))

# Tier 2: Overall RA Network Hubs (including shared edges)
g_all_ra <- graph_from_data_frame(
  ra_consensus %>% dplyr::select(Regulator, Target), 
  directed = TRUE
)

all_out_degree <- degree(g_all_ra, mode = "out")
all_betweenness <- betweenness(g_all_ra, directed = TRUE)

master_regulators_all <- data.frame(
  Gene = names(all_out_degree),
  OutDegree_Total = all_out_degree,
  Betweenness_Total = all_betweenness,
  stringsAsFactors = FALSE
)

# Combine tiers
master_regulators_final <- master_regulators_disease %>%
  full_join(master_regulators_all, by = "Gene") %>%
  replace_na(list(
    OutDegree_Disease = 0, 
    InDegree_Disease = 0,
    Betweenness_Disease = 0,
    OutDegree_Total = 0,
    Betweenness_Total = 0
  )) %>%
  mutate(
    # Priority score: Disease-specific regulation is weighted higher
    Priority_Score = OutDegree_Disease * 3 + 
      (Betweenness_Disease / max(Betweenness_Disease + 1)) * 5 +
      OutDegree_Total * 1,
    Tier = case_when(
      OutDegree_Disease >= 10 ~ "Tier 1: High-Impact Disease Regulator",
      OutDegree_Disease >= 5 ~ "Tier 2: Moderate Disease Regulator",
      OutDegree_Disease >= 2 ~ "Tier 3: Minor Disease Regulator",
      OutDegree_Total >= 10 ~ "Tier 4: General RA Network Hub",
      TRUE ~ "Tier 5: Supporting Node"
    )
  ) %>%
  arrange(desc(Priority_Score)) %>%
  mutate(Rank = row_number())


write.csv(disease_specific_edges, 
          "RA_MegaCohort_Disease_Specific_Network.csv", 
          row.names=F)

write.csv(all_edges_annotated, 
          "RA_MegaCohort_All_Edges_Annotated.csv", 
          row.names=F)

write.csv(master_regulators_final, 
          "RA_MegaCohort_Master_Regulators_Final.csv", 
          row.names=F)

message("Exported:")
message("  - RA_MegaCohort_Disease_Specific_Network.csv (", 
        nrow(disease_specific_edges), " edges)")
message("  - RA_MegaCohort_All_Edges_Annotated.csv (", 
        nrow(all_edges_annotated), " edges)")
message("  - RA_MegaCohort_Master_Regulators_Final.csv (", 
        nrow(master_regulators_final), " genes)")


top20 <- master_regulators_final %>%
  filter(OutDegree_Disease > 0) %>%
  dplyr::select(Rank, Gene, Tier, OutDegree_Disease, OutDegree_Total, 
                Betweenness_Disease, Priority_Score) %>%
  head(20)

print(top20)

# Count by tier
tier_summary <- master_regulators_final %>%
  filter(OutDegree_Disease > 0) %>%
  count(Tier) %>%
  arrange(Tier)

message("\n MASTER REGULATOR TIER DISTRIBUTION ")
print(tier_summary)



# Analyze top disease-specific regulators
top_disease_genes <- master_regulators_final %>%
  filter(OutDegree_Disease >= 2) %>%
  pull(Gene)

message("Analyzing ", length(top_disease_genes), " disease-specific regulators for GO enrichment...")

if(length(top_disease_genes) >= 5) {
  # Check which genes are valid
  valid_genes <- top_disease_genes[top_disease_genes %in% keys(org.Hs.eg.db, keytype = "SYMBOL")]
  
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
        write.csv(ego_bp@result, "Mega_RA_MegaCohort_MasterRegulators_GO_BP.csv", row.names=F)
        
        pdf("RA_MegaCohort_GO_BP_Dotplot.pdf", width=10, height=8)
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
        write.csv(ego_mf@result, "Mega_RA_MegaCohort_MasterRegulators_GO_MF.csv", row.names=F)
        message("Saved GO MF enrichment results")
      }
      
      # KEGG Pathways
      # Convert to Entrez IDs for KEGG
      gene_entrez <- bitr(valid_genes, 
                          fromType = "SYMBOL",
                          toType = "ENTREZID", 
                          OrgDb = org.Hs.eg.db)
      
      if(nrow(gene_entrez) >= 5) {
        kegg <- enrichKEGG(
          gene = gene_entrez$ENTREZID,
          organism = "hsa",
          pvalueCutoff = 0.05
        )
        
        if(nrow(kegg@result) > 0) {
          write.csv(kegg@result, "Mega_RA_MegaCohort_MasterRegulators_KEGG.csv", row.names=F)
          message("Saved KEGG pathway enrichment")
        }
      }
      
    }, error = function(e) {
      message("GO enrichment error: ", e$message)
    })
  } else {
    message("Only ", length(valid_genes), " valid gene symbols found (need >= 5)")
  }
} else {
  message("Insufficient disease-specific regulators for enrichment (need >= 5)")
}

message("Disease-Specific Network:")
message("  Nodes: ", vcount(g_disease))
message("  Edges: ", ecount(g_disease))
message("  Density: ", round(edge_density(g_disease), 4))
message("  Avg out-degree: ", round(mean(disease_out_degree), 2))

message("\nOverall RA Network:")
message("  Nodes: ", vcount(g_all_ra))
message("  Edges: ", ecount(g_all_ra))
message("  Density: ", round(edge_density(g_all_ra), 4))

# Method contribution analysis
method_stats <- all_edges_annotated %>%
  summarise(
    GENIE3_only = sum(!is.na(Genie_Weight_RA) & is.na(Aracne_MI_RA) & is.na(MIIC_Confidence)),
    ARACNe_only = sum(is.na(Genie_Weight_RA) & !is.na(Aracne_MI_RA) & is.na(MIIC_Confidence)),
    MIIC_only = sum(is.na(Genie_Weight_RA) & is.na(Aracne_MI_RA) & !is.na(MIIC_Confidence)),
    Two_methods = sum(Methods_Count == 2),
    All_three = sum(Methods_Count == 3)
  )


print(method_stats)

message("\nCREATING VISUALIZATION FILES")

# Node attributes for disease-specific network
node_attrs <- master_regulators_final %>%
  filter(Gene %in% V(g_disease)$name) %>%
  dplyr::select(Gene, OutDegree_Disease, Betweenness_Disease, Tier, Priority_Score)

write.csv(node_attrs, "Mega_RA_Disease_Network_Nodes.csv", row.names=F)

# Edge list with attributes
edge_attrs <- disease_specific_edges %>%
  dplyr::select(Regulator, Target, Genie_Weight_RA, Aracne_MI_RA, 
                MIIC_Confidence, Methods_Count, is_causal)

write.csv(edge_attrs, "Mega_RA_Disease_Network_Edges.csv", row.names=F)

message("Exported Cytoscape-ready files:")
message("  - Mega_RA_Disease_Network_Nodes.csv")
message("  - Mega_RA_Disease_Network_Edges.csv")


message("\nKey findings:")
message("  Disease-specific edges: ", nrow(disease_specific_edges))
message("  Tier 1 master regulators: ", 
        sum(master_regulators_final$Tier == "Tier 1: High-Impact Disease Regulator"))
message("  Tier 2 master regulators: ", 
        sum(master_regulators_final$Tier == "Tier 2: Moderate Disease Regulator"))
message("\nRecommendations for next steps:")
