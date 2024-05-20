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
                               random.seed = 1990)

save(seurat_isotype, file=paste0(load.path,"Seuratv5_WNN_Complete.RData"))

####
DimPlot(seurat_isotype, reduction='wnn.umap',label=T,label.size = 5,repel=T)

############################# Cluster Counts #####################################
# Extract the cluster information
cluster_counts <- table(Idents(seurat_isotype))

# Convert to a data frame for easier manipulation if needed
cluster_counts_df <- as.data.frame(cluster_counts)

# Identify clusters with fewer than 10 cells
clusters_to_keep <- names(cluster_counts[cluster_counts >= 10])

seurat_isotype <- subset(seurat_isotype, idents = clusters_to_keep)

# Extract the current levels
current_levels <- levels(seurat_isotype$wsnn_res.0.8)

# Convert the levels to numeric and sort them
sorted_levels <- sort(as.numeric(current_levels))

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
sample_cluster_proportions$wsnn_res.0.8 <- factor(sample_cluster_proportions$wsnn_res.0.8, levels = sort(as.numeric(levels(seurat_isotype$wsnn_res.0.8))))
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

DimPlot(seurat_isotype, reduction='wnn.umap',label=T,label.size = 5,repel=T)
ggsave('WNN_ClusterPlot_res0.8.png',dpi=500, width = 8)

DimPlot(seurat_isotype, reduction='wnn.umap',label=T,label.size = 5,repel=T, split.by = 'Condition')
ggsave('WNN_ClusterPlot_res0.8_bycondition.png',dpi=500, width = 14)


#################################### Heatmaps ##############################################
# create multimodal heatmap 


# find marker genes for the joint clusters 

Idents(seurat_isotype) = "wsnn_res.0.8"
DefaultAssay(seurat_isotype)  = "RNA"
rnade = FindAllMarkers(seurat_isotype, features = rna.features, only.pos = T)
gene_plot = rnade %>% 
  dplyr::filter(avg_log2FC > 1 ) %>%  
  dplyr::group_by(cluster) %>% 
  dplyr::top_n(3) %$% gene %>% unique 


cite_data = GetAssayData(seurat_isotype,slot = 'data',assay = 'ADT') %>% t()
rna_subset = GetAssayData(seurat_isotype,assay = 'RNA',slot = 'data')[gene_plot, ] %>%
  as.data.frame() %>% 
  t() %>% 
  as.matrix()

# combine into dataframe 
d_r = cbind(seurat_isotype@meta.data, rna_subset) 
d_p = cbind(seurat_isotype@meta.data, cite_data) 
# calculate the median protein expression per cluster
dat_plot_r = d_r %>% 
  dplyr::group_by(wsnn_res.0.8) %>% 
  dplyr::summarize_at(.vars = c(gene_plot), .funs = median) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("wsnn_res.0.8") 

dat_plot_p = d_p %>% 
  dplyr::group_by(wsnn_res.0.8) %>% 
  dplyr::summarize_at(.vars = c(prots), .funs = median) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("wsnn_res.0.8") 

# protein heatmap 
# protein heatmap 
prot_col = circlize::colorRamp2(breaks = seq(-1,25, by = 1), 
                                colors = viridis::viridis(n = 27, option = "B"))
p1 = Heatmap(t(dat_plot_p)[prots, ], 
             name = "protein", 
             col = prot_col, 
             use_raster = T,
             row_names_gp = gpar(color = "black", fontsize = 5)
)
p1


# mRNA heatmap 
mrna = t(dat_plot_r)[gene_plot, ]
rna_col = circlize::colorRamp2(breaks = c(-2,-1,0,1,2), 
                               colors = colorspace::diverge_hsv(n = 5))

zero_variance_rows <- which(apply(mrna, 1, var) == 0)
mrna_cleaned <- mrna[-zero_variance_rows, ]


p2 = Heatmap(t(scale(t(mrna_cleaned))), 
             name = "mRNA", 
             col = rna_col,
             use_raster = T, 
             clustering_method_columns = 'average',
             column_names_gp = gpar(color = "black", fontsize = 7), 
             row_names_gp = gpar(color = "black", fontsize = 5))


p2

# heatmaps
setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Heatmap')

tiff("ADT_Heatmap_ALLPROT.tiff", width = 2500, height = 4000, res = 300) # Adjust width, height, and resolution as needed
draw(p1)
dev.off()
tiff("RNA_Heatmap_TopGene.tiff", width = 2000, height = 2000, res = 300) # Adjust width, height, and resolution as needed
draw(p2)
dev.off()

########################################## ADT Feature Plots and VLN Plots ###############################################
DefaultAssay(seurat_isotype) <- 'ADT'

features <- seurat_isotype@assays$ADT@counts@Dimnames[[1]]

# Violin Plot

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Violin_Plot/ADT')

for (i in features) {
  vln.pl <- VlnPlot(seurat_isotype, features = i, assay = 'ADT')
  ggsave(paste0(i,'_VLNplot.png'),dpi=500, width = 13, vln.pl)
  vln.pl.2 <- VlnPlot(seurat_isotype, features = i, assay = 'ADT', group.by = 'predicted.celltype.l2')
  ggsave(paste0(i,'_VLNplot_Azimuth.png'),dpi=500, width = 13, vln.pl.2)
}

# Feature Plots

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Feature_Plot/ADT')

for (i in features) {
  fea.pl <- FeaturePlot(seurat_isotype, reduction = 'wnn.umap', features = i)
  ggsave(paste0(i,'_Featureplot.png'),dpi=500, width = 8, fea.pl)
}

########################################## RNA Feature Plots and VLN Plots ###############################################

DefaultAssay(seurat_isotype) <- 'RNA'

features <- c('CD14','FCGR2B','SERPING1','CCR7','CD27','TCF7','CCL5','FCGR3A','PRF1','CD40LG','IRF8','TNFRSF4',
              'CD8A','TNFRSF9','XCL2','CD7','CD8B','NELL2','C1QBP','CD3E','ICOS','IGFBP2','IGFBP4','LDHA',
              'CCND3','MIR155HG','NR4A1','CTLA4','FOXP3','IL2RA','CD19','CD79A','IGHM','EBI3','HLA-DPA1',
              'HLA-DRB1','CTSW','KLRC1','TNFRSF18','CCR4','IRF4','MALAT1','IKZF2','TRDV1','TRGC2',
              'CD3D','CXCR3','GZMK','CCL2','HLA-DRA','SERPINA1','GNLY','NKG7','TIGIT','LTB','MAL','SELL',
              'CCL4L2','CD70','IFNG','IL2RB','KLRD1','TRBC1','HAVCR2','LGALS1','NCAM1','CD36','CD4','IFI30',
              'CXCL8','ITGAX','IL18BP','TNF','TRDV2','TRGV9','FABP5','MT-ND1','MT-ND5','CCL3','IL1B','TNFAIP2',
              'CD40','MS4A1','XCL1','HIST1H4C','LTA','MKI67')

# Violin Plot

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Violin_Plot/RNA')

for (i in features) {
  vln.pl <- VlnPlot(seurat_isotype, features = i, assay = 'RNA')
  ggsave(paste0(i,'_VLNplot.png'),dpi=500, width = 13, vln.pl)
  vln.pl.2 <- VlnPlot(seurat_isotype, features = i, assay = 'RNA', group.by = 'predicted.celltype.l2')
  ggsave(paste0(i,'_VLNplot_Azimuth.png'),dpi=500, width = 13, vln.pl.2)
}

# Feature Plots

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Feature_Plot/RNA')

for (i in features) {
  fea.pl <- FeaturePlot(seurat_isotype, reduction = 'wnn.umap', features = i)
  ggsave(paste0(i,'_Featureplot.png'),dpi=500, width = 8, fea.pl)
}



















###### EXTRA #################


#############################################################################################################

###### Cluster Distribution by donor
setwd('C:/Users/axi313/Documents/ECCITE_Seq_Katherine/WNN/Cluster_Proportions')

# Step 1: Extract relevant information
predicted_celltypes <- seurat_isotype_h@meta.data$predicted.celltype.l2
donor_ids <- seurat_isotype_h@meta.data$orig.ident

# Step 2: Calculate cell counts per cluster per donor
cluster_counts <- table(donor_ids, predicted_celltypes)

# Step 3: Create a table with raw numbers
raw_numbers_table <- as.data.frame.matrix(cluster_counts)

# Print or save the two tables
print("Raw Numbers:")
print(raw_numbers_table)

# Export raw_numbers_table as a CSV file
write.csv(raw_numbers_table, file = "raw_numbers_table.csv", row.names = TRUE)

# Create a stacked bar plot
raw <- seurat_isotype_h@meta.data %>% 
  ggplot(aes(x=orig.ident, fill=predicted.celltype.l2)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")


# Calculate percentage distribution
percentage_distribution <- prop.table(table(seurat_isotype_h@meta.data$orig.ident, seurat_isotype_h@meta.data$predicted.celltype.l2), margin = 1) * 100

# Convert the percentage distribution table to a data frame
percentage_df <- as.data.frame(percentage_distribution)

# Rename the columns for clarity
colnames(percentage_df) <- c("Donor ID", "Cluster", "Percentage")


# Plot the stacked bar plot with % values as fill
percent <- ggplot(percentage_df, aes(x = `Donor ID`, y = Percentage, fill = Cluster)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Percentage Distribution of Cell Clusters by Donor")

ggsave('Cluster_Raw_byDonor.png',dpi=500, width = 7, raw)
ggsave('Cluster_Prop_byDonor.png',dpi=500, width = 7, percent)

###########################
# load(paste0(load.path,"Seuratv5_WNN_Complete.RData"))

################################## Differential Expression ############################################
setwd('C:/Users/axi313/Documents/ECCITE_Seq_Katherine/WNN/Differential_Expression/Azimuth')


# Assuming 'seurat_isotype_h' is your Seurat object

# Define the clusters
clusters <- levels(as.factor(seurat_isotype_h@meta.data$predicted.celltype.l2))
clusters
# Create an empty list to store differential expression results
de_results <- list()

# Perform differential expression analysis for each cluster
for (cluster in clusters) {
  # Subset the data for the current cluster
  subset_data <- subset(seurat_isotype_h, subset= predicted.celltype.l2 == cluster)
  
  # Check if both conditions have at least 20 cells
  if (sum(subset_data$stim == "Med") >= 20 & sum(subset_data$stim == "HIV_peptide") >= 20) {
    # Perform differential expression analysis comparing HIV peptide stimulation to media stimulation
    de_result <- FindMarkers(subset_data, ident.1 = "HIV_peptide", ident.2 = "Med", group.by = 'stim', test.use = "MAST")
    
    # Add the result to the list
    de_results[[cluster]] <- de_result
  } else {
    # If the condition is not met, skip the cluster and print a message
    cat("Skipping cluster", cluster, "due to insufficient cells for differential expression analysis.\n")
  }
}


# Save the results as CSV files
for (result_name in names(de_results)) {
  write.csv(de_results[[result_name]], file = paste0(result_name, "_HIV_peptide_VS_Media__DGE_results_Azimuth.csv"), row.names = TRUE)
}
