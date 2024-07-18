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
library(scRepertoire)

##set path to load data


setwd('C:/Users/axi313/Documents/TARA_Entry/WNN')

load.path <- "C:/Users/axi313/Documents/TARA_Entry/saved_R_data/"

load(paste0(load.path,'Seuratv5_WNN_labled.RData'))

### Fix CTLGrouping metadata

metadata <- seurat_isotype@meta.data
metadata$CTLGrouping[which(str_detect(metadata$cells, "SAAH29"))] <- "high"
metadata$CTLGrouping[which(str_detect(metadata$cells, "SATY021"))] <- "low"


seurat_isotype@meta.data <- metadata

seurat_HEI <- subset(seurat_isotype, subset = Condition == 'HEI' )

DimPlot(seurat_HEI, label = T, split.by = 'CTLGrouping', reduction = 'wnn.umap')
DimPlot(seurat_HEI, label = T, split.by = 'Viral_Load_Category', reduction = 'wnn.umap')

#### VDJ ####

#Filenames

in.path <- "C:/Users/axi313/Documents/TARA_Entry/Raw_Data/"

# Function to get all folder names within a specified path
get_folder_names <- function(in.path) {
  # List all directories within the specified path without recursion
  folder_names <- list.dirs(in.path, full.names = FALSE, recursive = FALSE)
  
  # Filter out the root path itself if present
  folder_names <- folder_names[folder_names != ""]
  
  return(folder_names)
}

f_names <- get_folder_names(in.path)



colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                            "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))


for (name in f_names) {
  # Construct the file paths
  t_file <- paste0(in.path, name, '/per_sample_outs/TCR/filtered_contig_annotations.csv')
  b_file <- paste0(in.path, name, '/per_sample_outs/BCR/filtered_contig_annotations.csv')
  
  # Read the files
  t_data <- read.csv(t_file)
  b_data <- read.csv(b_file)
  
  # Create dynamically named variables
  assign(paste0(name, ".TCR"), t_data)
  assign(paste0(name, ".BCR"), b_data)
}

# Create Contig list

f_names.TCR <- paste(f_names, ".TCR", sep="")
f_names.BCR <- paste(f_names, ".BCR", sep="")
contig_list.TCR <- as.list(mget(f_names.TCR))
contig_list.BCR <- as.list(mget(f_names.BCR))


#Combine For downstream Analysis

combined.TCR <- combineTCR(contig_list.TCR,samples = f_names)

combined.BCR <- combineBCR(contig_list.BCR,samples = f_names)

############ Merge Seurat #####################
# Access the cell barcodes (Assuming they are in the column names of the data slot)
barcodes <- rownames(seurat_isotype[[]])

# Use gsub to modify the barcodes, removing everything before and including the fourth '_'
modified_barcodes <- gsub(".*_.*_.*_.*_(.*)", "\\1", barcodes)
modified_barcodes <- paste0(seurat_isotype$orig.ident, "_", modified_barcodes)

# Assign the modified barcodes back to the Seurat object
seurat_isotype <- RenameCells(seurat_isotype, new.names = modified_barcodes)


seurat.tcr <- combineExpression(combined.TCR, 
                                seurat_isotype, 
                                cloneCall="strict",
                                group.by = 'sample',
                                cloneSize = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded =
                                                500),
                                proportion = FALSE)

seurat.bcr <- combineExpression(combined.BCR, 
                                seurat_isotype, 
                                cloneCall="strict",
                                group.by = 'sample',
                                cloneSize = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded =
                                                500),
                                proportion = FALSE)

###################################################### Hyper Expansion Plots ##################################

cols <- c("Hyperexpanded (100 < X <= 500)"="#F0F921","Large (20 < X <= 100)" = "#F69441",
          "Medium (5 < X <= 20)" = "#CA4778" ,"Small (1 < X <= 5)" = "#7D06A5",
          "Single (0 < X <= 1)" = "#0D0887")

## first ordering the Clone Size as a factor, this prevents the coloring from being in alphabetical order. 
slot(seurat.tcr, "meta.data")$cloneSize <- factor(slot(seurat.tcr, "meta.data")$cloneSize, 
                                                  levels = c("Hyperexpanded (100 < X <= 500)", 
                                                             "Large (20 < X <= 100)", 
                                                             "Medium (5 < X <= 20)", 
                                                             "Small (1 < X <= 5)", 
                                                             "Single (0 < X <= 1)", NA))
slot(seurat.bcr, "meta.data")$cloneSize <- factor(slot(seurat.bcr, "meta.data")$cloneSize, 
                                                  levels = c("Hyperexpanded (100 < X <= 500)", 
                                                             "Large (20 < X <= 100)", 
                                                             "Medium (5 < X <= 20)", 
                                                             "Small (1 < X <= 5)", 
                                                             "Single (0 < X <= 1)", NA))


### TCR
setwd("C:/Users/axi313/Documents/TARA_Entry/WNN/VDJ/Seurat_Plots")

tcr.condition <- DimPlot(seurat.tcr, group.by = "cloneSize", split.by = 'Condition', reduction = 'wnn.umap') +
  scale_color_manual(values=cols) + ggtitle('Expanded TCR Clone Size')

tcr.VL <- DimPlot(seurat.tcr, group.by = "cloneSize", split.by = 'CTLGrouping', reduction = 'wnn.umap') +
  scale_color_manual(values=cols) + ggtitle('Expanded TCR Clone Size')

tcr.CTL <- DimPlot(seurat.tcr, group.by = "cloneSize", split.by = 'Viral_Load_Category', reduction = 'wnn.umap') +
  scale_color_manual(values=cols) + ggtitle('Expanded TCR Clone Size')

ggsave('TCR_Seurat_Condition.png', width = 15, dpi = 500, tcr.condition)
ggsave('TCR_Seurat_VL.png', width = 12, dpi = 500, tcr.VL)
ggsave('TCR_Seurat_CTL.png', width = 12, dpi = 500, tcr.CTL)


DimPlot(subset(seurat.tcr, ident='20'), group.by = "cloneSize", reduction = 'wnn.umap') +
  scale_color_manual(values=cols) + ggtitle('Expanded TCR Clone Size')



DimPlot(seurat_isotype, label=T, reduction='wnn.umap')
DimPlot(subset(seurat_isotype, ident='singleton'), label=T, reduction='wnn.umap')

### BCR
bcr.condition <- DimPlot(seurat.bcr, group.by = "cloneSize", split.by = 'condition', reduction = 'wnn.umap') +
  scale_color_manual(values=cols) + ggtitle('Expanded BCR Clone Size')
bcr.stim <- DimPlot(seurat.bcr, group.by = "cloneSize", split.by = 'stim', reduction = 'wnn.umap') +
  scale_color_manual(values=cols) + ggtitle('Expanded BCR Clone Size')
#bcr.sample <- DimPlot(seurat.bcr, group.by = "cloneSize", split.by = 'sample', reduction = 'wnn.umap') +
#  scale_color_manual(values=cols) + ggtitle('Expanded BCR Clone Size')

ggsave('BCR_Seurat_Condition.png', width = 12, dpi = 500, bcr.condition)
ggsave('BCR_Seurat_Stim.png', width = 12, dpi = 500, bcr.stim)
#ggsave('BCR_Seurat_Sample.png', width = 18, dpi = 500, bcr.sample)

### Recluster

seurat_isotype <- FindClusters(seurat_isotype, graph.name = "wsnn", 
                               group.singletons=F,
                               algorithm = 3, 
                               resolution = 0.8,
                               random.seed = 1990)

save(seurat_isotype, file=paste0(load.path,"Seuratv5_WNN_Complete_2.RData"))
