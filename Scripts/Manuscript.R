###
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

##### Interferon Stimulated Genes ###########
setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Manuscript/ISG')

type1_IFN_vector <- c("MX1", "OAS1", "ISG15", "IFITM3", "EIF2AK2", "IFIT1", 
                      "IRF9", "BST2", "RSAD2", "IFI6", "ZC3HAV1", "TRIM22", 
                      "DDX58", "IFI27", "SAMHD1", "APOBEC3G")

IFN_gamma_vector <- c("GBP1", "CXCL10", "IRF9", "SOCS1", "MX1")


DefaultAssay(seurat_isotype) <- 'RNA'

### Feature Plots


pal <- viridis(n = 10, option = "A") # Magma Colour Scheme

### Type 1 IFN stim genes

## All
fea_al_custom <- FeaturePlot_scCustom(seurat_isotype, reduction = 'wnn.umap', features = type1_IFN_vector, colors_use = pal,keep.scale = 'feature',order=TRUE)
ggsave('Type1_ISG_Featureplot_ALL_RNA_magma.png',dpi=500, width = 13,height = 12, fea_al_custom)

## All split by VL
fea_al_custom_VL <- FeaturePlot_scCustom(seurat_isotype, reduction = 'wnn.umap', 
                                         split.by='Viral_Load_Category',features = type1_IFN_vector, colors_use = pal,keep.scale = 'feature',order=TRUE)
ggsave('Type_1_ISG_Featureplot_ALL_RNA_Viral_Load_magma.png',dpi=500, width = 12,height = 48,fea_al_custom_VL)

### Clustered Dot Plot
clust_dot_plot_IFN <- Clustered_DotPlot(seurat_object = seurat_isotype, features = type1_IFN_vector, split.by = 'Viral_Load_Category',  cluster_ident = FALSE,
                                      colors_use_idents='polychrome'
)

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Manuscript/ISG/Violin_Plot')
for (i in type1_IFN_vector) {
  
  vln.pl <-VlnPlot_scCustom(seurat_isotype,features = i ,split.by = 'Viral_Load_Category')
  ggsave(paste0(i,'_VLNplot_VL_HighvsLow.png'),dpi=500, width = 12, vln.pl)
}
#### CD4 Analysis

### Subset CD4 Cells
#CD4 Clusters, : 

CD4_sub <- subset(seurat_isotype,idents=c(0,2,3,6,8,11,21),subset = Condition == 'HEI')
DefaultAssay(CD4_sub) <-'RNA'
# Define the CD4 clusters of interest
cd4_clusters <- c("0", "2", "3", "6", "8", "11", "21")

# Extract the hex codes for the specified CD4 clusters (without names)
CD4_colors <- unname(polychrome_colors[cd4_clusters])

# View the resulting vector of hex codes
print(CD4_colors)



### Annotation
setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Manuscript/CD4_Annotation')

# Find variable features in the CD4_sub dataset
CD4_sub <- FindVariableFeatures(CD4_sub)

# Extract the top 30 most variable features
top50_features <- head(VariableFeatures(CD4_sub), 50)

# View the top 30 most variable features
print(top50_features)

png("CD4_Clust_dotplot_RNA_Top50.png",width=8,height=7,units="in",res=1200)

Clustered_DotPlot(seurat_object = CD4_sub, features = top50_features,  cluster_ident = FALSE,
                  colors_use_idents= CD4_colors
)
dev.off()
# Vector of CD4+ T-cell subset marker genes
cd4_marker_genes <- c(
  
  # Th1 markers
  "TBX21", "IFNG", "STAT1", "CXCR3",  # Loww CCR6
  
  # TH1/TH17
  "CXCR3", "CCR6",
  
  # Th2 markers
  "GATA3", "IL4", "IL5", "IL13", "CCR4", # Low CXCR3
  
  # Treg markers
  "FOXP3", "IL2RA", "CTLA4", "IKZF2", "TGFB1",
  
  # Th17 markers
  "RORC", "IL17A", "IL17F", "IL22", "CCR6",  # Low CXCR3
  
  # Tfh markers
  "BCL6", "CXCR5", "ICOS", "PDCD1", "IL21",
  
  # Th9 markers
  "IL9", "IRF4", "STAT6",
  
  # Th22 markers
  "IL22", "AHR",
  
  # Th0 (Naive) markers
  "CCR7", "SELL", "LEF1", "TCF7", 
  
  # Memory CD4+ T cell markers
  "IL7R", "CCR7", "CD27", 'CD45RO',
  
  # Additional activation or differentiation markers
  "CD69", "CD40LG", "KLRG1", "HLA-DRA"
)
DefaultAssay(CD4_sub) <- 'RNA'


png("CD4_Clust_dotplot_RNA_CD4markers.png",width=8,height=7,units="in",res=1200)

Clustered_DotPlot(seurat_object = CD4_sub, features = cd4_marker_genes,  cluster_ident = FALSE,
                  colors_use_idents= CD4_colors
)
dev.off()

DefaultAssay(CD4_sub) <- 'ADT'

ADT.features <- CD4_sub@assays$ADT@counts@Dimnames[[1]]

png("CD4_Clust_dotplot_ADT_all.png",width=10,height=18,units="in",res=1200)

Clustered_DotPlot(seurat_object = CD4_sub, features = ADT.features,  cluster_ident = FALSE,
                  colors_use_idents= CD4_colors
)
dev.off()

## 11 - > Treg

###   RNA 

DefaultAssay(CD4_sub) <- 'RNA'

rna.features <- c(
  "IRF4", "NCR3LG1", "PDCD1", "CXCR5", "IL7R", "IFNGR2", "ICOS", "FOXP3", "TGFB1", 
  "CCR7", "TCF7", "HLA-DRA", "HLA-E", "HLA-DPB1", "HLA-DQB1", "HLA-DRB1", 
  "IL2RA", "KLRG1", "IFNG", "STAT1",  "RORC", "NECTIN2", "TNFRSF1A", "MKI67", "PVR", "GATA3", "HLA-C", "IL7", 
   "IFNGR1", "LEF1", "CCR4", "MICA", "CD40LG", "CXCR3", "CD69", 
   "BCL6", "FAS", "CD27", "HLA-B", "TNFRSF1B", "SELL", "PRF1", "TBX21", 
  "AHR", "IL15", "MICB", "STAT6", "CCR6", "HLA-A", "CTLA4", "IKZF2"
)

png("CD4_Clust_dotplot_RNA_Flowsom_Ligand_Compliment.png",width=10,height=14,units="in",res=1200)
Clustered_DotPlot(seurat_object = CD4_sub, features = rna.features,  cluster_ident = FALSE,
                  colors_use_idents= CD4_colors
)
dev.off()

png("CD4_Clust_dotplot_RNA_Flowsom_Ligand_Compliment_byVL.png",width=13,height=14,units="in",res=1200)
Clustered_DotPlot(seurat_object = CD4_sub, features = rna.features,  cluster_ident = FALSE,
                  split.by = 'Viral_Load_Category'
)
dev.off()

# VLN Plots
setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Manuscript/CD4_Annotation/Violin/RNA')

for (i in rna.features) {
  
  vln.pl <-VlnPlot_scCustom(CD4_sub,features = i ,split.by = 'Viral_Load_Category')
  ggsave(paste0(i,'_CD4_VLNplot_VL_HighvsLow.png'),dpi=500, width = 10, vln.pl)
}

# Feature Plots
setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Manuscript/CD4_Annotation/Feature/RNA')

for (i in rna.features) {
  pal <- viridis(n = 10, option = "G")
  fea.pl <- FeaturePlot_scCustom(CD4_sub, reduction = 'wnn.umap', features = i,split.by = 'Viral_Load_Category', 
                                 colors_use = pal,order=TRUE)
  ggsave(paste0(i,'_Featureplot_HEIvsHEU_Mako.png'),dpi=500, width = 10, fea.pl)

}

### ADT
DefaultAssay(CD4_sub) <- 'ADT'

# VLN Plots
setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Manuscript/CD4_Annotation/Violin/ADT')

for (i in ADT.features) {
  
  vln.pl <-VlnPlot_scCustom(CD4_sub,features = i ,split.by = 'Viral_Load_Category')
  ggsave(paste0(i,'_CD4_VLNplotADT_VL_HighvsLow.png'),dpi=500, width = 10, vln.pl)
}

# Feature Plots
setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Manuscript/CD4_Annotation/Feature/ADT')

for (i in ADT.features) {
  pal <- viridis(n = 10, option = "G")
  fea.pl <- FeaturePlot_scCustom(CD4_sub, reduction = 'wnn.umap', features = i,split.by = 'Viral_Load_Category', 
                                 colors_use = pal,order=TRUE)
  ggsave(paste0(i,'_FeatureplotADT_HEIvsHEU_Mako.png'),dpi=500, width = 10, fea.pl)
  
}

DefaultAssay(CD4_sub) <- 'RNA'

#### DGE ###
setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Manuscript/CD4_DGE')


#Identify clusters in CD4_sub
clusters <- unique(Idents(CD4_sub))

# Loop through each cluster to perform DGE
dge_results <- list()

for (cluster in clusters) {
  
  # Subset data for the current cluster
  cluster_data <- subset(CD4_sub, idents = cluster)
  
  # Perform DGE analysis for high vs low viral load within the cluster
  dge <- FindMarkers(
    cluster_data, 
    ident.1 = "High", 
    ident.2 = "Low", 
    group.by = "Viral_Load_Category", 
    logfc.threshold = 0.25,
    test.use = 'MAST'
  )
  
  # Save DGE results for the cluster
  dge_file <- paste0("DGE_results_cluster_", cluster, ".csv")
  write.csv(dge, dge_file)
  
  # Append to results list
  dge_results[[as.character(cluster)]] <- dge
}

### VOlcano
for (cluster in names(dge_results)) {
  
  # Extract DGE results for the cluster
  dge <- dge_results[[cluster]]
  
  # Sort to find top up and downregulated genes
  top_genes <- dge %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 10) %>%
    bind_rows(dge %>% arrange(avg_log2FC) %>% slice_head(n = 10))
  
  # Generate volcano plot
  plot <- EnhancedVolcano(
    dge,
    lab = rownames(dge),
    x = "avg_log2FC",
    y = "p_val_adj",
    title = paste("Volcano Plot of Cluster", cluster),
    selectLab = rownames(top_genes),
    subtitle = bquote(italic('High VL vs Low VL')),
    pCutoff = 0.05,
    FCcutoff = 0.25,
    labSize = 3.0,
    pointSize = 1.5,  # Smaller point size
    drawConnectors = TRUE,  # Add lines connecting labels
    widthConnectors = 0.5,  # Width of connector lines
    colConnectors = "black",  # Color of connector lines
    legendPosition = 'top',
    legendLabSize = 10,
    legendIconSize = 3.0
  )
  
  # Save the plot to the working directory
  plot_file <- paste0("Volcano_Plot_cluster_", cluster, ".png")
  ggsave(plot_file, plot = plot, width = 8, height = 6)
}

### NK Clusters ###########

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Manuscript/NK_Annotation')
NK <- subset(seurat_isotype,idents=c(4,15,16,19),subset = Condition == 'HEI')

# Define the NK clusters of interest
NK_clusters <- c("4", "15", "16", "19")

# Extract the hex codes for the specified NK clusters (without names)
NK_colors <- unname(polychrome_colors[NK_clusters])

# View the resulting vector of hex codes
print(NK_colors)


### RNA
rna.features <- c(
  "B3GAT1", "CD226", "CD69", "CX3CR1", "FASLG", "FCER1G", "FCGR3A", "GNLY", 
  "GZMB", "HLA-DRB1", "HLA-DRB5", "IFNG", "IL2RB", "IL7R", "ITGA2", "KIR2DL1", 
  "KIR3DL1", "KLRB1", "KLRC1", "KLRC2", "KLRD1", "KLRG1", "KLRK1", "LAG3", 
  "MKI67", "NCAM1", "NCR1", "NCR2", "NCR3", "NKG7", "PRF1", "SIGLEC7", 
  "TIGIT", "TNF", "TNFRSF18", "XCL1", "XCL2"
)

DimPlot_scCustom(NK, reduction='wnn.umap',label=F, repel=T, colors_use =NK_colors)
ggsave('NK_WNN_ClusterPlot_res0.8_HEI.png',dpi=500)
DimPlot_scCustom(NK, reduction='wnn.umap',label=F, repel=T, colors_use =NK_colors,split.by = 'Viral_Load_Category')
ggsave('NK_WNN_ClusterPlot_res0.8_splitbyVL.png',width=11, dpi=500)


# Feature Plots
setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Manuscript/NK_Annotation/Feature_Plot/RNA')

for (i in rna.features) {
  pal <- viridis(n = 10, option = "A")
  fea.pl <- FeaturePlot_scCustom(NK, reduction = 'wnn.umap', features = i,split.by = 'Viral_Load_Category', 
                                 colors_use = pal,order=TRUE)
  ggsave(paste0(i,'_Featureplot_HEIvsHEU_Magma.png'),dpi=500, width = 10, fea.pl)
}

# VLN Plots
setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Manuscript/NK_Annotation/Violin_Plot/RNA')

for (i in rna.features) {
  
  vln.pl <-VlnPlot_scCustom(NK,features = i ,split.by = 'Condition')
  ggsave(paste0(i,'_VLNplot_HEIvsHEU.png'),dpi=500, width = 10, vln.pl)
}

### Dot Plots

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Manuscript/NK_Annotation/Dot_Plot')

png("NK_Clust_dotplot_RNA.png",width=8,height=9,units="in",res=1200)

Clustered_DotPlot(seurat_object = NK, features = rna.features,  cluster_ident = FALSE,
                  colors_use_idents= NK_colors
)
dev.off()

png("NK_Clust_dotplot_RNA_VL.png",width=8,height=9,units="in",res=1200)

Clustered_DotPlot(seurat_object = NK, features = rna.features,  cluster_ident = FALSE,
                  split.by = 'Viral_Load_Category'
)
dev.off()
