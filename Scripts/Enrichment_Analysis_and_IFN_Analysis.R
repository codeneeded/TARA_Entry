#Load Required Libraries
library(Nebulosa)
library(Seurat)
library(SingleCellExperiment)
library(viridis)
library(scater)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(hdf5r)
library(dsb)
library(data.table)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(SeuratWrappers)
library(Azimuth)
library(ggrepel)
library(patchwork)
library(scCustomize)
library(reticulate)
library(circlize)
library(ComplexHeatmap)
library(readxl)
library(Polychrome)
library(EnhancedVolcano)
library(ggpubr)
library(enrichR)
library(openxlsx)

##set path to load data


setwd('C:/Users/axi313/Documents/TARA_Entry/WNN')

load.path <- "C:/Users/axi313/Documents/TARA_Entry/saved_R_data/"

load(paste0(load.path,'Seuratv5_WNN_Complete.RData'))

### Identify the colours used in each cluster

palette36_colors <- palette36.colors(36)  # Generates 36 colors

# Get the number of clusters in your dataset
num_clusters <- length(levels(seurat_isotype$wsnn_res.0.8))

# Generate the Polychrome color palette for the number of clusters
polychrome_colors <- palette36.colors(num_clusters)

# Name the colors by the clusters
names(polychrome_colors) <- levels(seurat_isotype$wsnn_res.0.8)

# View the hex codes for each cluster
print(polychrome_colors)


#### Genes of interest - EXTRA
setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/IFN')

DefaultAssay(seurat_isotype) <-'RNA'

### Subset
# Subset the Seurat object to only include cells with Condition = "HEI"
seurat_isotype_HEI <- subset(seurat_isotype, subset = Condition == "HEI")

# Subset the Seurat object to only include cells with Condition = "HEI" or "HEU"
seurat_isotype_HEI_HEU <- subset(seurat_isotype, subset = Condition %in% c("HEI", "HEU"))

#### GENE LISTS #############

# Primarily Type I IFN-Stimulated Genes
type1_genes <- c("MX1", "OAS1", "ISG15", "IFITM3", "EIF2AK2", "IFIT1", 
                 "IRF9", "BST2", "RSAD2", "IFI6", "ZC3HAV1", "TRIM22", 
                 "DDX58", "IFI27", "SAMHD1", "APOBEC3G", "USP18", "IFNAR1", 
                 "IFIT3", "IFIT5","IFNAR1","IFNAR2")

# Genes Responsive to Both Type I and Type II IFNs
type1_and_type2_genes <- c("STAT1", "STAT2", "IRF1", "IRF7", "SOCS1", 
                           "SOCS3", "CXCL10", "SP100", "ADAR", 
                           "TRIM21", "GBP1", "PML", "HLA-E", "XAF1")

# Genes Responsive to Type I and Type III IFNs (with Tissue-Specific Expression)
type1_and_type3_genes <- c("ISG15", "IFITM3", "OAS1", "MX1", "DDX58", 
                           "IFI6", "IFI16", "DHX58", "ZBP1", "RSAD2", 
                           "HERC5")

sting_pathway_genes <- c("TMEM173", "CGAS", "TBK1", "IRF3", "IFI16", 
                         "DDX41", "NFKB1", 
                         "IKBKB", "CHUK", "NFKBIA", "TMEM173", "TBKBP1", 
                          "ULK1", "TRAF6", "TRAF3", "TRAF2", 
                         "TRAFD1")

### Violin plots for all genes

### Function for saving vln plots
# Define the function
create_violin_plots <- function(seurat_obj, gene_list, split_by) {
  # Loop through each gene in the gene list
  for (gene in gene_list) {
    # Create the violin plot with significance testing
    p <- VlnPlot_scCustom(seurat_obj, features = gene, split.by = split_by, pt.size = 0.1) + 
      stat_compare_means(label = "p.signif", method = "wilcox.test",hide.ns = T) +
      labs(title = paste("Gene:", gene))
    
    # Save the plot to the working directory with the gene name in the file name
    ggsave(filename = paste0(gene, "_violin_plot.png"), plot = p, 
           width = 11, height = 7, dpi = 300, bg='white')

  }
}
### HIGHVLvsLOWVL


setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/IFN/VLN/HighVLvsLowVL/Type1')
create_violin_plots(seurat_isotype_HEI, type1_genes, "Viral_Load_Category")

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/IFN/VLN/HighVLvsLowVL/Type1_Type2')
create_violin_plots(seurat_isotype_HEI, type1_and_type2_genes, "Viral_Load_Category")

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/IFN/VLN/HighVLvsLowVL/Type1_Type3')
create_violin_plots(seurat_isotype_HEI, type1_and_type3_genes, "Viral_Load_Category")

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/IFN/VLN/HighVLvsLowVL/STING')
create_violin_plots(seurat_isotype_HEI, sting_pathway_genes, "Viral_Load_Category")


### HEI vs HEU


setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/IFN/VLN/HEIvsHEU/Type1')
create_violin_plots(seurat_isotype_HEI_HEU, type1_genes, "Condition")

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/IFN/VLN/HEIvsHEU/Type1_Type2')
create_violin_plots(seurat_isotype_HEI_HEU, type1_and_type2_genes, "Condition")

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/IFN/VLN/HEIvsHEU/Type1_Type3')
create_violin_plots(seurat_isotype_HEI_HEU, type1_and_type3_genes, "Condition")

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/IFN/VLN/HEIvsHEU/STING')
create_violin_plots(seurat_isotype_HEI_HEU, sting_pathway_genes, "Condition")


####### Differential Expression for Genelists ################

# Define the function with error handling
perform_diff_exp_by_cluster <- function(seurat_obj, split_by = "Viral_Load_Category", ident_high = "High", ident_low = "Low") {
  
  # Initialize a list to store results
  diff_exp_results <- list()
  
  # Loop through each cluster to perform differential expression analysis
  for (cluster in unique(Idents(seurat_obj))) {
    # Perform differential expression with error handling
    tryCatch({
      # Perform differential expression between high and low viral load for each cluster
      markers <- FindMarkers(seurat_obj, ident.1 = ident_high, ident.2 = ident_low, 
                             group.by = split_by, subset.ident = cluster, test.use = "MAST")
      
      # Add cluster and gene information to the results
      markers$cluster <- cluster
      markers$gene <- rownames(markers)
      
      # Store results in the list
      diff_exp_results[[cluster]] <- markers
      
    }, error = function(e) {
      # Print the error message but continue the loop
      message(paste("Error in cluster", cluster, ":", e$message))
    })
  }
  
  # Combine all results into a single data frame
  all_results <- do.call(rbind, diff_exp_results)
  
  # Filter significant results (adjusted p-value < 0.05)
  significant_results <- all_results %>%
    filter(p_val_adj < 0.05)
  
  # Return the combined and filtered data frame
  return(significant_results)
}

# Define the function to filter by gene list
filter_by_genelist <- function(diff_exp_results, gene_list) {
  # Filter results to include only genes in the gene list
  filtered_results <- diff_exp_results %>%
    filter(gene %in% gene_list)
  filtered_results$cluster <- factor(filtered_results$cluster, levels = as.character(sort(unique(as.numeric(filtered_results$cluster)))))
  filtered_results <- filtered_results %>% filter(!is.na(cluster))
  return(filtered_results)
}

# Define the function to plot and save
plot_and_save_DGE <- function(data, title = "DGE High vs Low VL Type1 IFN", filename = "DGE_plot.png") {
  
  # Generate the plot
  p <- ggplot(data, aes(x = cluster, y = gene, color = avg_log2FC, size = abs(avg_log2FC))) +
    geom_point() +
    scale_color_gradient2(low = "#008080", mid = "white", high = "#FF69B4", midpoint = 0) +
    labs(title = title,
         x = "Cluster", y = "Gene", color = "Log2 Fold Change") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )+
    guides(size = "none") 
  
  # Save the plot to the working directory
  ggsave(filename = filename, plot = p, width = 8, height = 6, dpi = 300, bg='white')
  
  # Print confirmation
  message(paste("Plot saved as", filename, "in the working directory."))
}

# High Vs Low VL
setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/IFN/Dotplot/Viral_Load')

VL_DGE_ALL <- perform_diff_exp_by_cluster(seurat_isotype_HEI)

### Type 1
type1_DGE <- filter_by_genelist(significant_diff_exp_results,type1_genes)
plot_and_save_DGE(type1_DGE, title = "DGE of High vs Low VL (Type 1 IFN genelist)", filename = "DGE_Type1_IFN_plot.png")

### Type 2
type1_and_type2_DGE <- filter_by_genelist(significant_diff_exp_results,type1_and_type2_genes)
plot_and_save_DGE(type1_and_type2_DGE, title = "DGE High vs Low VL (Type 1 and Type 2 IFN genelist)", filename = "DGE_Type1_Type2_IFN_plot.png")

### Type 3
type1_and_type3 <- filter_by_genelist(significant_diff_exp_results,type1_and_type3_genes)
plot_and_save_DGE(type1_and_type3, title = "DGE High vs Low VL (Type 1 and Type 3 IFN genelist)", filename = "DGE_Type1_Type3_IFN_plot.png")

### STING
sting_pathway_DGE <- filter_by_genelist(significant_diff_exp_results,sting_pathway_genes)
plot_and_save_DGE(sting_pathway_DGE, title = "DGE High vs Low VL (STING Pathway genelist)", filename = "DGE_STING_IFN_plot.png")


# HEI vs HEU
setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/IFN/Dotplot/HIV_Status')

HIV_DGE_ALL <- perform_diff_exp_by_cluster(seurat_isotype_HEI_HEU,split_by = "Condition", ident_high = "HEI", ident_low = "HEU")

### Type 1
type1_DGE <- filter_by_genelist(HIV_DGE_ALL,type1_genes)
plot_and_save_DGE(type1_DGE, title = "DGE HEI vs HEU (Type 1 IFN genelist)", filename = "DGE_Type1_IFN_plot.png")

### Type 2
type1_and_type2_DGE <- filter_by_genelist(HIV_DGE_ALL,type1_and_type2_genes)
plot_and_save_DGE(type1_and_type2_DGE, title = "DGE HEI vs HEU (Type 1 and Type 2 IFN genelist)", filename = "DGE_Type1_Type2_IFN_plot.png")

### Type 3
type1_and_type3 <- filter_by_genelist(HIV_DGE_ALL,type1_and_type3_genes)
plot_and_save_DGE(type1_and_type3, title = "DGE HEI vs HEU (Type 1 and Type 3 IFN genelist)", filename = "DGE_Type1_Type3_IFN_plot.png")

### STING
sting_pathway_DGE <- filter_by_genelist(HIV_DGE_ALL,sting_pathway_genes)
plot_and_save_DGE(sting_pathway_DGE, title = "DGE HEI vs HEU (STING Pathway genelist)", filename = "DGE_STING_IFN_plot.png")

###################### Lesley Fas Genes + Proteins #######################################
setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Manuscript/Dotplot')

genelist_1 <- c("FAS","HLA-A","HLA-B","HLA-C","HLA-E","HLA-G")
genelist_2 <- c("FASLG","KIR2DL1", "KIR2DL3", "KLRB1", "KLRD1", "KLRG1", "KLRK1")
protlist_1 <- c("FAS","HLA-A","HLA-E")
protlist_2 <- c("KIR2DL1", "KIR2DL3", "KLRB1", "KLRD1", "KLRG1", "KLRK1")



DefaultAssay(seurat_isotype_HEI) <- "RNA"
VL_DGE_ALL <- perform_diff_exp_by_cluster(seurat_isotype_HEI)
DefaultAssay(seurat_isotype_HEI) <- "ADT"
VL_DPE_ALL <- perform_diff_exp_by_cluster(seurat_isotype_HEI)

DefaultAssay(seurat_isotype_HEI_HEU) <- "RNA"
HIV_DGE_ALL <- perform_diff_exp_by_cluster(seurat_isotype_HEI_HEU, split_by = "Condition", ident_high = "HEI", ident_low = "HEU")
DefaultAssay(seurat_isotype_HEI_HEU) <- "ADT"
HIV_DPE_ALL <- perform_diff_exp_by_cluster(seurat_isotype_HEI_HEU, split_by = "Condition", ident_high = "HEI", ident_low = "HEU")
########### SAVE ###########
save(VL_DGE_ALL, file = paste0(load.path, "VL_DGE_ALL.RData"))
save(VL_DPE_ALL, file = paste0(load.path, "VL_DPE_ALL.RData"))
save(VL_DGE_ALL, file = paste0(load.path, "HIV_DGE_ALL.RData"))
save(VL_DPE_ALL, file = paste0(load.path, "HIV_DPE_ALL.RData"))

load.path.2 <- 'C:/Users/axi313/Documents/TARA_Entry/WNN/Manuscript/DE/'
write.csv(VL_DGE_ALL, file = paste0(load.path.2, "VL_DGE_ALL.csv"), row.names = FALSE)
write.csv(VL_DPE_ALL, file = paste0(load.path.2, "VL_DPE_ALL.csv"), row.names = FALSE)
write.csv(HIV_DGE_ALL, file = paste0(load.path.2, "HIV_DGE_ALL.csv"), row.names = FALSE)
write.csv(HIV_DPE_ALL, file = paste0(load.path.2, "HIV_DPE_ALL.csv"), row.names = FALSE)

############## LOAD ##########################
load(paste0(load.path,'VL_DGE_ALL.RData'))
load(paste0(load.path,'VL_DPE_ALL.RData'))

###################################

DGE_1 <- filter_by_genelist(VL_DGE_ALL,genelist_1)
DGE_2 <- filter_by_genelist(VL_DGE_ALL,genelist_2)
DPE_1 <- filter_by_genelist(VL_DPE_ALL,protlist_1)
DPE_2 <- filter_by_genelist(VL_DPE_ALL,protlist_2)

### Filter to needed clusters only
clusters_to_keep_1 <- c(0, 2, 3, 6, 8, 11, 21)
clusters_to_keep_2 <- c(4, 15, 16, 19)

# Subset the data to only include these clusters
DGE_1_filtered <- DGE_1[DGE_1$cluster %in% clusters_to_keep_1, ]
DPE_1_filtered <- DPE_1[DPE_1$cluster %in% clusters_to_keep_1, ]
DGE_2_filtered <- DGE_2[DGE_2$cluster %in% clusters_to_keep_2, ]
DPE_2_filtered <- DPE_2[DPE_2$cluster %in% clusters_to_keep_2, ]

#### Plots

ggplot(DGE_1_filtered, aes(x = cluster, y = gene, color = avg_log2FC, size = abs(avg_log2FC))) +
  geom_point() +
  scale_color_gradient2(low = "#008080", mid = "white", high = "#FF69B4", midpoint = 0) +
  labs(title= "DGE High vs Low VL Load CD4 Clusters",
       x = "Cluster", y = "Gene", color = "Log2 Fold Change") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )+
  guides(size = "none") 


# Save the plot to the working directory
ggsave("DGE_CD4_Dotplot.png", width = 5, height = 3, dpi = 300, bg='white')

ggplot(DPE_1_filtered, aes(x = cluster, y = gene, color = avg_log2FC, size = abs(avg_log2FC))) +
  geom_point() +
  scale_color_gradient2(low = "#008080", mid = "white", high = "#FF69B4", midpoint = 0) +
  labs(title= "DPE High vs Low VL Load CD4 Clusters",
       x = "Cluster", y = "Gene", color = "Log2 Fold Change") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )+
  guides(size = "none") 


# Save the plot to the working directory
ggsave("DPE_CD4_Dotplot.png", width = 5, height = 3, dpi = 300, bg='white')

###
ggplot(DGE_2_filtered, aes(x = cluster, y = gene, color = avg_log2FC, size = abs(avg_log2FC))) +
  geom_point() +
  scale_color_gradient2(low = "#008080", mid = "white", high = "#FF69B4", midpoint = 0) +
  labs(title= "DGE High vs Low VL Load NK Clusters",
       x = "Cluster", y = "Gene", color = "Log2 Fold Change") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )+
  guides(size = "none") 


# Save the plot to the working directory
ggsave("DGE_NK_Dotplot.png", width = 5, height = 3, dpi = 300, bg='white')

ggplot(DPE_2_filtered, aes(x = cluster, y = gene, color = avg_log2FC, size = abs(avg_log2FC))) +
  geom_point() +
  scale_color_gradient2(low = "#008080", mid = "white", high = "#FF69B4", midpoint = 0) +
  labs(title= "DPE High vs Low VL Load NK Clusters",
       x = "Cluster", y = "Gene", color = "Log2 Fold Change") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )+
  guides(size = "none") 


# Save the plot to the working directory
ggsave("DPE_NK_Dotplot.png", width = 5, height = 3, dpi = 300, bg='white')



################################# ENRICHMENT ANALYSIS ###################################################
# Get unique clusters
clusters <- unique(VL_DGE_ALL$cluster)
# Define databases for enrichment analysis
databases <- c("TRRUST_Transcription_Factors_2019", "ChEA_2022","TRANSFAC_and_JASPAR_PWMs","KEGG_2021_Human", 
               "WikiPathways_2024_Human","GO_Biological_Process_2023","MSigDB_Hallmark_2020",
               "Panther_2016","Reactome_2022","BioPlanet_2019")
# Initialize a list to store results
enrichment_results <- list()

#### Viral Load

### mRNA
# Initialize a list to store results
enrichment_results <- list()
# Perform enrichment for each cluster
for (i in 0:29) {
  # Get genes for the current cluster
  gene_list <- VL_DGE_ALL$gene[VL_DGE_ALL$cluster == i]
  
  # Run enrichment analysis
  enrichment_results[[as.character(i)]] <- enrichR::enrichr(genes = gene_list, databases = databases)
}

# Define the output directory
output_dir <- "C:/Users/axi313/Documents/TARA_Entry/WNN/Manuscript/Enrichment_Viral_Load/mRNA/Data/"

if (!dir.exists(output_dir)) dir.create(output_dir)

# Create separate Excel files for each cluster
for (cluster in names(enrichment_results)) {
  # Create a new workbook for this cluster
  wb <- createWorkbook()
  
  for (db_name in names(enrichment_results[[cluster]])) {
    # Get the filtered results for this database
    db_results <- enrichment_results[[cluster]][[db_name]]
    
    # Create a valid sheet name (truncate to 31 characters)
    sheet_name <- substr(db_name, 1, 31)
    
    # Add the results to a new sheet
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, db_results)
  }
  
  # Define the file name for the current cluster
  output_file <- paste0(output_dir, "Cluster_", cluster, "_Enrichment_Results.xlsx")
  
  # Save the workbook
  saveWorkbook(wb, output_file, overwrite = TRUE)
}

#### ENRICHMENT PLOTS - Viral Load

# Define the output directory
output_dir <- "C:/Users/axi313/Documents/TARA_Entry/WNN/Manuscript/Enrichment_Viral_Load/mRNA/Plots/"

# Define transcription factor and pathway databases
tf_databases <- c("TRRUST_Transcription_Factors_2019", "ChEA_2022", "TRANSFAC_and_JASPAR_PWMs")
pathway_databases <- c("KEGG_2021_Human", "WikiPathways_2024_Human", "GO_Biological_Process_2023", 
                       "MSigDB_Hallmark_2020", "Panther_2016", "Reactome_2022", "BioPlanet_2019")

# Iterate over each cluster in the filtered enrichment results
for (cluster_to_plot in names(enrichment_results)) {
  
  # Initialize a list to store the top pathways for all databases
  top_pathways_list <- list()
  
  # Iterate over each database for the cluster
  for (db_name in names(enrichment_results[[cluster_to_plot]])) {
    # Extract the results for the database
    db_results <- enrichment_results[[cluster_to_plot]][[db_name]]
    
    # Filter for significant pathways (adjusted p-value < 0.05)
    significant_results <- db_results %>%
      filter(`Adjusted.P.value` < 0.05)
    
    # Proceed only if there are significant results
    if (nrow(significant_results) > 0) {
      # Select the top 10 pathways by adjusted p-value
      top_10 <- significant_results %>%
        arrange(`Combined.Score`) %>%
        slice_head(n = 10)
      
      # Add a column for the database name
      top_10 <- top_10 %>%
        mutate(Database = db_name)
      
      # Append to the list
      top_pathways_list[[db_name]] <- top_10
    }
  }
  
  # Combine all results into a single data frame
  top_pathways_df <- bind_rows(top_pathways_list)
  
  # Filter into transcription factors and pathway subsets
  tf_df <- top_pathways_df %>%
    filter(Database %in% tf_databases) %>%
    arrange(desc(Combined.Score)) %>%
    slice(1:20)
  
  pathway_df <- top_pathways_df %>%
    filter(Database %in% pathway_databases) %>%
    arrange(desc(Combined.Score)) %>%
    slice(1:20)
  
  # Function to create and save plots
  create_and_save_plot <- function(data, title_suffix, filename_suffix) {
    if (nrow(data) > 0) {
      p <- ggplot(data, aes(x = reorder(Term, `Combined.Score`), y = Combined.Score, fill = Database)) +
        geom_bar(stat = "identity", position = "dodge") +
        coord_flip() +
        labs(
          title = paste("Top", title_suffix, "for Cluster", cluster_to_plot),
          x = "Pathway",
          y = "Enrichment Score"
        ) +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 12, color = "black"),  # Y-axis tick labels
          axis.text.x = element_text(size = 12),        # X-axis tick labels
          axis.title.y = element_text(size = 14),       # Y-axis title
          axis.title.x = element_text(size = 14),       # X-axis title
          plot.title = element_text(size = 16, hjust = 0.5)  # Plot title
        )
      
      ggsave(
        filename = paste0(output_dir, "Cluster_", cluster_to_plot, "_", filename_suffix, ".png"),
        plot = p,
        width = 15,
        height = 12,
        dpi = 300,
        bg = "white"
      )
      
      cat("Saved", title_suffix, "plot for Cluster", cluster_to_plot, "\n")
    } else {
      message(paste("No significant", title_suffix, "found for cluster", cluster_to_plot, "across the databases."))
    }
  }
  
  # Create and save plots
  create_and_save_plot(tf_df, "Transcription Factors", "Transcription_Factors")
  create_and_save_plot(pathway_df, "Pathways", "Pathways")
}



### HIV Status

### mRNA
enrichment_results <- list()

# Perform enrichment for each cluster
for (i in 0:29) {
  # Get genes for the current cluster
  gene_list <- HIV_DGE_ALL$gene[HIV_DGE_ALL$cluster == i]
  
  # Run enrichment analysis
  enrichment_results[[as.character(i)]] <- enrichR::enrichr(genes = gene_list, databases = databases)
}

# Define the output directory
output_dir <- "C:/Users/axi313/Documents/TARA_Entry/WNN/Manuscript/Enrichment_HIV_Status/mRNA/Data/"

if (!dir.exists(output_dir)) dir.create(output_dir)

# Create separate Excel files for each cluster
for (cluster in names(enrichment_results)) {
  # Create a new workbook for this cluster
  wb <- createWorkbook()
  
  for (db_name in names(enrichment_results[[cluster]])) {
    # Get the filtered results for this database
    db_results <- enrichment_results[[cluster]][[db_name]]
    
    # Create a valid sheet name (truncate to 31 characters)
    sheet_name <- substr(db_name, 1, 31)
    
    # Add the results to a new sheet
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, db_results)
  }
  
  # Define the file name for the current cluster
  output_file <- paste0(output_dir, "Cluster_", cluster, "_Enrichment_Results.xlsx")
  
  # Save the workbook
  saveWorkbook(wb, output_file, overwrite = TRUE)
}


#### ENRICHMENT PLOTS - HIV Status

# Define the output directory
output_dir <- "C:/Users/axi313/Documents/TARA_Entry/WNN/Manuscript/Enrichment_HIV_Status/mRNA/Plots/"

# Define transcription factor and pathway databases
tf_databases <- c("TRRUST_Transcription_Factors_2019", "ChEA_2022", "TRANSFAC_and_JASPAR_PWMs")
pathway_databases <- c("KEGG_2021_Human", "WikiPathways_2024_Human", "GO_Biological_Process_2023", 
                       "MSigDB_Hallmark_2020", "Panther_2016", "Reactome_2022", "BioPlanet_2019")

# Iterate over each cluster in the filtered enrichment results
for (cluster_to_plot in names(enrichment_results)) {
  
  # Initialize a list to store the top pathways for all databases
  top_pathways_list <- list()
  
  # Iterate over each database for the cluster
  for (db_name in names(enrichment_results[[cluster_to_plot]])) {
    # Extract the results for the database
    db_results <- enrichment_results[[cluster_to_plot]][[db_name]]
    
    # Filter for significant pathways (adjusted p-value < 0.05)
    significant_results <- db_results %>%
      filter(`Adjusted.P.value` < 0.05)
    
    # Proceed only if there are significant results
    if (nrow(significant_results) > 0) {
      # Select the top 10 pathways by adjusted p-value
      top_10 <- significant_results %>%
        arrange(`Combined.Score`) %>%
        slice_head(n = 10)
      
      # Add a column for the database name
      top_10 <- top_10 %>%
        mutate(Database = db_name)
      
      # Append to the list
      top_pathways_list[[db_name]] <- top_10
    }
  }
  
  # Combine all results into a single data frame
  top_pathways_df <- bind_rows(top_pathways_list)
  
  # Filter into transcription factors and pathway subsets
  tf_df <- top_pathways_df %>%
    filter(Database %in% tf_databases) %>%
    arrange(desc(Combined.Score)) %>%
    slice(1:20)
  
  pathway_df <- top_pathways_df %>%
    filter(Database %in% pathway_databases) %>%
    arrange(desc(Combined.Score)) %>%
    slice(1:20)
  
  # Function to create and save plots
  create_and_save_plot <- function(data, title_suffix, filename_suffix) {
    if (nrow(data) > 0) {
      p <- ggplot(data, aes(x = reorder(Term, `Combined.Score`), y = Combined.Score, fill = Database)) +
        geom_bar(stat = "identity", position = "dodge") +
        coord_flip() +
        labs(
          title = paste("Top", title_suffix, "for Cluster", cluster_to_plot),
          x = "Pathway",
          y = "Enrichment Score"
        ) +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 12, color = "black"),  # Y-axis tick labels
          axis.text.x = element_text(size = 12),        # X-axis tick labels
          axis.title.y = element_text(size = 14),       # Y-axis title
          axis.title.x = element_text(size = 14),       # X-axis title
          plot.title = element_text(size = 16, hjust = 0.5)  # Plot title
        )
      
      ggsave(
        filename = paste0(output_dir, "Cluster_", cluster_to_plot, "_", filename_suffix, ".png"),
        plot = p,
        width = 15,
        height = 12,
        dpi = 300,
        bg = "white"
      )
      
      cat("Saved", title_suffix, "plot for Cluster", cluster_to_plot, "\n")
    } else {
      message(paste("No significant", title_suffix, "found for cluster", cluster_to_plot, "across the databases."))
    }
  }
  
  # Create and save plots
  create_and_save_plot(tf_df, "Transcription Factors", "Transcription_Factors")
  create_and_save_plot(pathway_df, "Pathways", "Pathways")
}









############ TYPE 1 IFN ####################






#### FILTER FOR TYPE 1 IFN
# Combine all interferon-related gene lists
all_ifn_genes <- unique(c(type1_genes, type1_and_type2_genes, type1_and_type3_genes))

# Filter the enrichment results
filtered_enrichment_results <- list()



for (cluster in names(enrichment_results)) {
  # Initialize a list to store filtered databases for this cluster
  filtered_cluster_results <- list()
  
  # Iterate over each database in the enrichment results for this cluster
  for (db_name in names(enrichment_results[[cluster]])) {
    # Extract the data frame for this database
    db_results <- enrichment_results[[cluster]][[db_name]]
    
    # Ensure db_results is not empty before proceeding
    if (nrow(db_results) == 0) next
    
    # Filter rows where the 'Genes' column includes any of the interferon-related genes
    db_results_filtered <- db_results[sapply(db_results$Genes, function(genes) {
      any(all_ifn_genes %in% unlist(strsplit(as.character(genes), ";")))
    }), ]
    
    # Store filtered results for this database if any rows are left
    if (nrow(db_results_filtered) > 0) {
      filtered_cluster_results[[db_name]] <- db_results_filtered
    }
  }
  
  # Store filtered results for this cluster if any databases have results
  if (length(filtered_cluster_results) > 0) {
    filtered_enrichment_results[[cluster]] <- filtered_cluster_results
  }
}

# Define the output directory
output_dir <- "Enrichment_Results_IFN/"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Create separate Excel files for each cluster
for (cluster in names(enrichment_results)) {
  # Create a new workbook for this cluster
  wb <- createWorkbook()
  
  for (db_name in names(enrichment_results[[cluster]])) {
    # Get the filtered results for this database
    db_results <- enrichment_results[[cluster]][[db_name]]
    
    # Create a valid sheet name (truncate to 31 characters)
    sheet_name <- substr(db_name, 1, 31)
    
    # Add the results to a new sheet
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, db_results)
  }
  
  # Define the file name for the current cluster
  output_file <- paste0(output_dir, "Cluster_", cluster, "_Enrichment_Results.xlsx")
  
  # Save the workbook
  saveWorkbook(wb, output_file, overwrite = TRUE)
}


#### ENRICHMENT PLOTS #############
# Define the output directory
output_dir <- "Enrichment_Plots_IFN/"
if (!dir.exists(output_dir)) dir.create(output_dir)


# Iterate over each cluster in the filtered enrichment results
for (cluster_to_plot in names(filtered_enrichment_results)) {
  
  # Initialize a list to store the top pathways for all databases
  top_pathways_list <- list()
  
  # Iterate over each database for the cluster
  for (db_name in names(filtered_enrichment_results[[cluster_to_plot]])) {
    # Extract the results for the database
    db_results <- filtered_enrichment_results[[cluster_to_plot]][[db_name]]
    
    # Filter for significant pathways (adjusted p-value < 0.05)
    significant_results <- db_results %>%
      filter(`Adjusted.P.value` < 0.05)
    
    # Proceed only if there are significant results
    if (nrow(significant_results) > 0) {
      # Select the top 10 pathways by adjusted p-value
      top_10 <- significant_results %>%
        arrange(`Adjusted.P.value`) %>%
        slice_head(n = 10)
      
      # Add a column for the database name
      top_10 <- top_10 %>%
        mutate(Database = db_name)
      
      # Append to the list
      top_pathways_list[[db_name]] <- top_10
    }
  }
  
  # Combine all results into a single data frame
  top_pathways_df <- bind_rows(top_pathways_list)
  
  # Check if there are any pathways to plot
  if (nrow(top_pathways_df) > 0) {
    # Plot the results
    p <- ggplot(top_pathways_df, aes(x = reorder(Term, -`Adjusted.P.value`), y = -log10(`Adjusted.P.value`), fill = Database)) +
      geom_bar(stat = "identity", position = "dodge") +
      coord_flip() +
      labs(
        title = paste("Top Pathways for Cluster", cluster_to_plot),
        x = "Pathway",
        y = "-log10(Adjusted P-value)"
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 12, color = "black"),  # Y-axis tick labels
        axis.text.x = element_text(size = 12),        # X-axis tick labels
        axis.title.y = element_text(size = 14),       # Y-axis title
        axis.title.x = element_text(size = 14),       # X-axis title
        plot.title = element_text(size = 16, hjust = 0.5)  # Plot title
      )
    
    # Save the plot
    ggsave(
      filename = paste0(output_dir, "Cluster_", cluster_to_plot, "_Top_Pathways.png"),
      plot = p,
      width = 12,
      height = 16,
      dpi = 300,
      bg = "white"  # Set background to white
    )
    
    cat("Saved plot for Cluster", cluster_to_plot, "\n")
  } else {
    message(paste("No significant pathways found for cluster", cluster_to_plot, "across the databases."))
  }
}

