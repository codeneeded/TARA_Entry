###
#Load Required Libraries
library(Nebulosa)
library(Seurat)
library(SingleCellExperiment)
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
##set path to load data


setwd('C:/Users/ammas/Documents/TARA_Entry/WNN')

load.path <- "C:/Users/ammas/Documents/TARA_Entry/saved_R_data/"

load(paste0(load.path,'Seuratv5_WNN_Complete.RData'))

##### Subset NK Cells #######

setwd('C:/Users/ammas/Documents/TARA_Entry/WNN/NK_Only')
NK <- subset(seurat_isotype, idents=c(4,15,16,19))

DimPlot(NK, reduction='wnn.umap',label=T,label.size = 4.5,repel=T, cols='polychrome', split.by = 'Condition')
ggsave('NK_WNN_ClusterPlot_res0.8.png',dpi=500, width = 10)

### Feature Plots

### RNA
setwd('C:/Users/ammas/Documents/TARA_Entry/WNN/NK_Only/Feature_Plot/RNA')

DefaultAssay(NK) <- 'RNA'

rna.features <-  c("NCAM1", "FCGR3A", "B3GAT1", "KLRC1", "KLRC2", "KLRK1", "NCR1", 
                   "NCR2", "NCR3", "PRF1", "GZMB", "HLA-DRB5", "HLA-DRB1", "FCER1G", "LAG3"
)
# Feature Plots


for (i in rna.features) {
  fea.pl <- FeaturePlot(NK, reduction = 'wnn.umap', features = i)
  ggsave(paste0(i,'_Featureplot.png'),dpi=500, width = 8, fea.pl)
}
### ADT
setwd('C:/Users/ammas/Documents/TARA_Entry/WNN/NK_Only/Feature_Plot/ADT')

DefaultAssay(NK) <- 'ADT'

adt.features <-  c( "LAG3",
                   "NCAM1", "FCGR3A", "B3GAT1","KLRK1", "NCR1"
)
# Feature Plots


for (i in adt.features) {
  fea.pl <- FeaturePlot(NK, reduction = 'wnn.umap', features = i)
  ggsave(paste0(i,'_Featureplot.png'),dpi=500, width = 8, fea.pl)
}



############## NK PLOTS #########################
library(pheatmap)

seurat_isotype <- ScaleData(seurat_isotype, assay = "ADT", vars.to.regress = NULL)
clusters_of_interest <- c("1","8", "10")
seurat_subset <- subset(seurat_isotype, idents = clusters_of_interest)
seurat_subset <- RenameIdents(object = seurat_subset, `1` = "CD8_1")
seurat_subset <- RenameIdents(object = seurat_subset, `8` = "CD8_2")
seurat_subset <- RenameIdents(object = seurat_subset, `10` = "CTL")

# Calculate average expression
avg_expression <- AverageExpression(seurat_subset, assays = 'RNA', verbose = FALSE)

# Calculate standard deviation for each gene
gene_variability <- apply(avg_expression$RNA, 1, sd) # Assuming average expression data is in avg_expression$RNA

# Sort genes by variability
sorted_genes <- sort(gene_variability, decreasing = TRUE)

# Exclude ribosomal proteins (genes starting with RPL and RPS)
non_ribosomal_genes <- !grepl("^RPL|^RPS", names(sorted_genes))

# Get names of the top 15 most variable non-ribosomal genes
top_genes <- names(sorted_genes[non_ribosomal_genes])
top15_non_ribosomal_genes <- top_genes[1:15]

### Chemokines ####
# List of specific genes of interest
### Plotting for senescence and Exhaution
s_genes <- c("HGF", "MMP2", "MMP9", "CCL2", "CCL20", "CXCL1", "CXCL8", "CXCL10", "IL10", "IL13",
             "IL15", "IL18", "IL1A", "IL1B", "IL2", "IL32", "IL6", "IL7", "SPP1", "TNF", "VEGFA", 
             "IL6ST", "SERPINE1","MKI67", "ATM", "NFKB1", "BCL2", "ACTB", "TGFB1", "E2F1", "ATR",
             "NFATC1", "BCL2L1", "GAPDH", "IFNG", "CDKN2A", "H2AFX", "RELA", "BCL2L11", "TUBA1A", "LMNB1",
             "CDKN1A", "STAT3", "BAD", "PUM1", "GLB1", "TP53", "PRDM1", "TERT", "CAV1")

e_genes <- c("PDCD1", "LAG3", "TIGIT", "HAVCR2", "CTLA4")
common_genes_s <- intersect(s_genes, rownames(avg_expression$RNA))
common_genes_e <- intersect(e_genes, rownames(avg_expression$RNA))

# Subset the average expression data for the selected genes
avg_expression_subset_s <- as.matrix(avg_expression$RNA[common_genes_s, ])
avg_expression_subset_s <- avg_expression_subset_s[apply(avg_expression_subset_s, 1, var) != 0, ]
avg_expression_subset_e <- as.matrix(avg_expression$RNA[common_genes_e, ])
avg_expression_subset_e <- avg_expression_subset_e[apply(avg_expression_subset_e, 1, var) != 0, ]

setwd('C:/Users/ammas/Documents/TARA_Entry/WNN/Exhaution_Markers')

# Plot the heatmap
pheatmap(avg_expression_subset_s,
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",
         scale = "row",
         treeheight_row = 0, # Hides row dendrogram
         treeheight_col = 0,
         main= "Senescence Comparisons",
         filename = "Senescence_Comparisons_Heatmap.png",  # Specify the file name and format
         width = 10,  # Width of the output file (in inches)
         height = 8)  # Scale rows to make patterns more visible

# Plot the heatmap
pheatmap(avg_expression_subset_e,
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_cols = FALSE,
         clustering_distance_rows = "euclidean",
         scale = "row",
         treeheight_row = 0, # Hides row dendrogram
         treeheight_col = 0,
         main= "Exhaution Comparisons",
         filename = "Exhaution_Comparisons_Heatmap.png",  # Specify the file name and format
         width = 6,  # Width of the output file (in inches)
         height = 4) # Scale rows to make patterns more visible

