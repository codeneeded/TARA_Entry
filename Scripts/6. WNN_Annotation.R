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


setwd('C:/Users/axi313/Documents/TARA_Entry/WNN')

load.path <- "C:/Users/axi313/Documents/TARA_Entry/saved_R_data/"

load(paste0(load.path,'Seuratv5_isotype_Assay3.RData'))

################################################# Integrate RNA #############################################
DefaultAssay(seurat_isotype) <- 'RNA'
DefaultAssay(seurat_isotype)


rna.list <- SplitObject(seurat_isotype, split.by = "orig.ident")

for (i in 1:length(rna.list)) {
  rna.list[[i]] <- NormalizeData(rna.list[[i]],assay = 'RNA')
  rna.list[[i]] <- FindVariableFeatures(rna.list[[i]], selection.method = "vst", 
                                        nfeatures = 2000,assay = 'RNA')
}

rna.anchors <- FindIntegrationAnchors(object.list = rna.list, dims = 1:30)

rna.integrated <- IntegrateData(anchorset = rna.anchors, dims = 1:30)

rna.integrated <- RenameAssays(rna.integrated, integrated = 'rna.integrated')

# Run the standard workflow for visualization and clustering
rna.integrated <- ScaleData(rna.integrated,assay = 'rna.integrated')
rna.integrated <- RunPCA(rna.integrated, npcs = 30, reduction.name = "pca.rna.integrated")
rna.integrated <- RunUMAP(rna.integrated, reduction = "pca.rna.integrated", dims = 1:30, reduction.name = "umap.rna.integrated")

save(rna.integrated, file=paste0(load.path,"Seuratv5_WNN_rna_integrated.RData"))
#load(paste0(load.path,'Seuratv5_WNN_rna_integrated.RData'))

################################################# Integrate ADT ###############################################
DefaultAssay(rna.integrated) <- 'ADT'
DefaultAssay(rna.integrated)
adt.list <- SplitObject(rna.integrated, split.by = "orig.ident")

# define proteins to use in clustering (non-isptype controls)
prots <- rownames(seurat_isotype@assays$ADT@data)
isotype_genes <- c('Mouse-IgG1', 'Mouse-IgG2a', 'Mouse-IgG2b', 'Rat-IgG2b', 'Rat-IgG1', 'Rat-IgG2a', 'Hamster-IgG')
prots <- setdiff(prots, isotype_genes)
for (i in 1:length(adt.list)) {
  VariableFeatures(adt.list[[i]]) = prots
}

features <- SelectIntegrationFeatures(object.list = adt.list)

adt.list <- lapply(X = adt.list, FUN = function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
})

anchors <- FindIntegrationAnchors(object.list = adt.list, reduction = "rpca", 
                                  dims = 1:30)

rna_adt_integrated <- IntegrateData(anchorset = anchors, dims = 1:30,new.assay.name = "integrated")
rna_adt_integrated <- RenameAssays(rna_adt_integrated, integrated = 'adt.integrated')

rna_adt_integrated <- ScaleData(rna_adt_integrated,assay = 'adt.integrated')
rna_adt_integrated <- RunPCA(rna_adt_integrated, npcs = 30, reduction.name = "pca.adt.integrated")
rna_adt_integrated <- RunUMAP(rna_adt_integrated, reduction = "pca.adt.integrated", dims = 1:30, reduction.name = "umap.adt.integrated")
save(rna_adt_integrated, file=paste0(load.path,"Seuratv5_WNN_rna_adt_integrated.RData"))

####################################### Seurat WNN default with PCA on dsb normalized protein ###############

seurat_isotype <- rna_adt_integrated


########### Check PCA variance
# SET rna features 
rna.features <- rownames(seurat_isotype@assays$rna.integrated)
seurat_isotype <- ScaleData(seurat_isotype,assay = 'rna.integrated')
seurat_isotype <- RunPCA(seurat_isotype, npcs = 30, reduction.name = "pca.rna.integrated",features=rna.features, assay = 'rna.integrated')
seurat_isotype <- RunUMAP(seurat_isotype, reduction = "pca.rna.integrated", dims = 1:30, reduction.name = "umap.rna.integrated")

### RNA

ElbowPlot(object = seurat_isotype, 
          ndims = 50, reduction = 'pca.rna.integrated')


# Determine percent of variation associated with each PC
pct <- seurat_isotype[["pca.rna.integrated"]]@stdev / sum(seurat_isotype[["pca.rna.integrated"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

# Minimum of the two calculation
pcs <- min(co1, co2)

pcs
# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

# Elbow plot to visualize 

ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

### PCs=20

### ADT
ElbowPlot(object = seurat_isotype, 
          ndims = 50, reduction = 'pca.adt.integrated')

# Determine percent of variation associated with each PC
pct <- seurat_isotype[["pca.adt.integrated"]]@stdev / sum(seurat_isotype[["pca.adt.integrated"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

# Minimum of the two calculation
pcs <- min(co1, co2)

pcs
# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

# Elbow plot to visualize 

ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()


### apca=14

# run WNN 
seurat_isotype = FindMultiModalNeighbors(
  seurat_isotype, reduction.list = list("pca.rna.integrated", "pca.adt.integrated"), 
  dims.list = list(1:20, 1:14)
)

# cluster 
seurat_isotype <-RunUMAP(seurat_isotype, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
seurat_isotype <- FindClusters(seurat_isotype, graph.name = "wsnn", 
                               algorithm = 3, 
                               resolution = 0.8,
                               group.singletons = F,
                               random.seed = 1990)

#Remove clusters with <50 cells
!Idents(seurat_isotype) %in% small_clusters
# Identify clusters with fewer than 20 cells
cluster_sizes <- table(Idents(seurat_isotype))
small_clusters <- names(cluster_sizes[cluster_sizes < 20])
large_clusters <- names(cluster_sizes[cluster_sizes > 20])



seurat_isotype <- subset(seurat_isotype, idents = large_clusters)


############################# Cluster Counts #####################################


# Extract the current levels
current_levels <- levels(droplevels(seurat_isotype$wsnn_res.0.8))

# Convert the levels to numeric and sort them
sorted_levels <- sort(as.numeric(current_levels))
sorted_levels <- c(sorted_levels,'singleton')


# Set the new sorted levels back to the factor
seurat_isotype$wsnn_res.0.8 <- factor(seurat_isotype$wsnn_res.0.8, levels = as.character(sorted_levels))

# Verify the new levels
levels(seurat_isotype$wsnn_res.0.8)


Idents(seurat_isotype) <- seurat_isotype$wsnn_res.0.8

# Extract the cluster information
cluster_counts <- table(Idents(seurat_isotype))

# Convert to a data frame for easier manipulation if needed
cluster_counts_df <- as.data.frame(cluster_counts)

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Cluster_Counts')

# Create a bar plot
# Assuming you have already created the plot and stored it in a variable
p <- ggplot(cluster_counts_df, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of Cells in Each Cluster", x = "Cluster", y = "Number of Cells") +
  theme_minimal()

# Save the plot to a file
ggsave("cluster_barplot.png", plot = p, width = 8, height = 6, bg = 'white')
# Save the data frame to a CSV file
write.csv(cluster_counts_df, "cluster_counts.csv", row.names = FALSE)

### Proportion Barplot

# Extract the metadata and cluster information
metadata <- seurat_isotype@meta.data

# Calculate the counts of each sample in each cluster
sample_cluster_counts <- metadata %>%
  group_by(orig.ident, wsnn_res.0.8) %>%
  summarise(count = n()) %>%
  ungroup()

# Calculate the total number of cells in each cluster
cluster_totals <- sample_cluster_counts %>%
  group_by(wsnn_res.0.8) %>%
  summarise(total = sum(count))

# Calculate the proportion of each sample in each cluster
sample_cluster_proportions <- sample_cluster_counts %>%
  left_join(cluster_totals, by = "wsnn_res.0.8") %>%
  mutate(proportion = count / total)

# Convert cluster and sample to factors for plotting
sample_cluster_proportions$wsnn_res.0.8 <- factor(sample_cluster_proportions$wsnn_res.0.8)
sample_cluster_proportions$orig.ident <- factor(sample_cluster_proportions$orig.ident, levels = unique(metadata$orig.ident))

# Create a proportion bar plot
proportion_plot <- ggplot(sample_cluster_proportions, aes(x = wsnn_res.0.8, y = proportion, fill = orig.ident)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Proportion of Each Sample Contributing to Each Cluster",
       x = "Cluster",
       y = "Proportion",
       fill = "Sample") +
  theme_minimal()

# Display the plot
print(proportion_plot)

# Save the plot to a file
ggsave("proportion_barplot.png", plot = proportion_plot, width = 10, height = 7, bg='white')

# Save the proportions data frame to a CSV file
write.csv(sample_cluster_proportions, "sample_cluster_proportions.csv", row.names = FALSE)

DimPlot(seurat_isotype, reduction='wnn.umap',label=T,label.size = 5,repel=T, cols='polychrome', split.by = 'Condition')

save(seurat_isotype, file=paste0(load.path,"Seuratv5_WNN_Complete.RData"))


########## Annotation #############

### Cluster Plots

seurat_isotype@meta.data$wsnn_res.0.8
setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Cluster_Plots')

DimPlot(seurat_isotype, reduction='wnn.umap',label=T,label.size = 5,repel=T, group.by = 'predicted.celltype.l1')
ggsave('WNN_ClusterPlot_Azimuth1.png',dpi=500, width = 8)

DimPlot(seurat_isotype, reduction='wnn.umap',label=T,label.size = 5,repel=T, group.by = 'predicted.celltype.l2')
ggsave('WNN_ClusterPlot_Azimuth2.png',dpi=500, width = 9.5)

DimPlot(seurat_isotype, reduction='wnn.umap',label=T,label.size = 5,repel=T, group.by = 'predicted.celltype.l3')
ggsave('WNN_ClusterPlot_Azimuth3.png',dpi=500, width = 11)

DimPlot(seurat_isotype, reduction='wnn.umap',label=T,label.size = 5,repel=T, cols='polychrome')
ggsave('WNN_ClusterPlot_res0.8.png',dpi=500, width = 8)

DimPlot(seurat_isotype, reduction='wnn.umap',label=T,label.size = 5,repel=T, split.by = 'Condition', cols='polychrome')
ggsave('WNN_ClusterPlot_res0.8_bycondition.png',dpi=500, width = 14)

DimPlot(seurat_isotype, reduction='wnn.umap',label=T,label.size = 5,repel=T, split.by = 'CTLGrouping', cols='polychrome')
ggsave('WNN_ClusterPlot_res0.8_byCTLGrouping.png',dpi=500, width = 14)

DimPlot(seurat_isotype, reduction='wnn.umap',label=T,label.size = 5,repel=T, split.by = 'Viral_Load_Category', cols='polychrome')
ggsave('WNN_ClusterPlot_res0.8_byViral_Load_Category.png',dpi=500, width = 14)

seurat_HEI <- subset(seurat_isotype, subset= Condition == 'HEI')

DimPlot(seurat_HEI, reduction='wnn.umap',label=T,label.size = 5,repel=T, split.by = 'CTLGrouping', cols='polychrome')
ggsave('HEI_WNN_ClusterPlot_res0.8_byCTLGrouping.png',dpi=500, width = 14)

DimPlot(seurat_HEI, reduction='wnn.umap',label=T,label.size = 5,repel=T, split.by = 'Viral_Load_Category', cols='polychrome')
ggsave('HEI_WNN_ClusterPlot_res0.8_byViral_Load_Category.png',dpi=500, width = 14)

