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
library(reticulate)
library(clusterProfiler)
library(BiocParallel)
library(fgsea)
library(DOSE)
library(enrichplot)
library(org.Hs.eg.db)
library(pathview)
library(magrittr)
use_python("/path/to/your/python")
##set path to load data


setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Vinh')

load.path <- "C:/Users/axi313/Documents/TARA_Entry/saved_R_data/"

load(paste0(load.path,'Seuratv5_WNN_labled.RData'))

seurat_NK <- subset(seurat_isotype, idents = c( 6,8,17,18,31), subset = Timepoint== 'Entry')

#### NK Subset

### PCA
rna.features <- rownames(seurat_NK@assays$rna.integrated)
adt.features <- rownames(seurat_NK@assays$adt.integrated)

seurat_NK <- ScaleData(seurat_NK,assay = 'rna.integrated')
seurat_NK <- RunPCA(seurat_NK, npcs = 30, reduction.name = "pca.rna.integrated",features=rna.features, assay = 'rna.integrated')
seurat_NK <- RunUMAP(seurat_NK, reduction = "pca.rna.integrated", dims = 1:30, reduction.name = "umap.rna.integrated")
seurat_NK <- ScaleData(seurat_NK,assay = 'adt.integrated')
seurat_NK <- RunPCA(seurat_NK, npcs = 30, reduction.name = "pca.adt.integrated",features=adt.features, assay = 'adt.integrated')
seurat_NK <- RunUMAP(seurat_NK, reduction = "pca.adt.integrated", dims = 1:30, reduction.name = "umap.adt.integrated"
                     )

### WNN ####

# run WNN 
seurat_NK = FindMultiModalNeighbors(
  seurat_NK, reduction.list = list("pca.rna.integrated", "pca.adt.integrated"), 
  dims.list = list(1:20, 1:15)
)

# cluster 
seurat_NK <-RunUMAP(seurat_NK, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

DimPlot(seurat_NK, label=T, split.by = 'Condition')


seurat_NK <- subset(seurat_NK, subset = Condition == 'HEI' | Condition == 'HEU')

### Cluster Plots ####

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Vinh/Cluster_Plots')

DimPlot(seurat_NK, reduction='wnn.umap',repel=T, group.by = 'wsnn_res.0.8')
ggsave('NK_WNN_ClusterPlot.png',dpi=500, width = 9.5)
DimPlot(seurat_NK, reduction='wnn.umap',label=T,label.size = 5,repel=T, group.by = 'wsnn_res.0.8')
ggsave('NK_WNN_ClusterPlot_Labled.png',dpi=500, width = 9.5)

### DGE ####

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Vinh/DGE')
DefaultAssay(seurat_NK) <- 'RNA'
for (i in levels(Idents(seurat_NK))) {
  dge <- FindMarkers(seurat_NK,ident.1 = 'HEI', ident.2 = 'HEU', group.by = 'Condition',
                     subset.ident = i, test.use = 'MAST')
  write.csv(dge, paste0('Cluster_',i,'_DGE_HEIvsHEU.csv'), row.names = TRUE)
}

### DPE ####

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Vinh/DPE')
DefaultAssay(seurat_NK) <- 'ADT'
for (i in levels(Idents(seurat_NK))) {
  dge <- FindMarkers(seurat_NK,ident.1 = 'HEI', ident.2 = 'HEU', group.by = 'Condition',
                     subset.ident = i, test.use = 'MAST')
  write.csv(dge, paste0('Cluster_',i,'_DPE_HEIvsHEU.csv'), row.names = TRUE)
}

### Violin Plots ####

###
rna.features <- c("CD4", "CD8A", "IL2RG", "KLRC4", "KLRK1", "GZMM", "GZMK", "CD69",
              "IL2RA", "NCAM1", "NCR2", "NCR3", "GZMB", "IL2RB", "GNLY", "KLRC1", 
              "KLRC2", "KLRC3", "GZMH", "NCR1", "KLRD1", "KLRF1", "NKG7", "FCGR3A", 
              "GZMA", "PRF1", "TIGIT", "KLRB1", "KLRG1")

prot.features <- seurat_NK@assays$ADT@counts@Dimnames[[1]]

# RNA

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Vinh/Violin_Plot/RNA')

for (i in rna.features) {
  vln.pl <- VlnPlot(seurat_NK, features = i, assay = 'RNA', split.by = 'Condition')
  ggsave(paste0(i,'_VLNplot_splitbyCondition.png'),dpi=500, width = 13, vln.pl)
  vln.pl <- VlnPlot(seurat_NK, features = i, assay = 'RNA')
  ggsave(paste0(i,'_VLNplot.png'),dpi=500, width = 13, vln.pl)
}

# ADT

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Vinh/Violin_Plot/ADT')

for (i in prot.features) {
  vln.pl <- VlnPlot(seurat_NK, features = i, assay = 'ADT', split.by = 'Condition')
  ggsave(paste0(i,'_VLNplot_splitbyCondition.png'),dpi=500, width = 13, vln.pl)
  vln.pl <- VlnPlot(seurat_NK, features = i, assay = 'ADT')
  ggsave(paste0(i,'_VLNplot.png'),dpi=500, width = 13, vln.pl)
}

#### Feature Plots ####
DefaultAssay(seurat_NK) <- 'ADT'

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Vinh/Feature_Plot/ADT')


for (i in prot.features) {
  fea.pl <- FeaturePlot(seurat_NK, reduction = 'wnn.umap', features = i)
  ggsave(paste0(i,'_Featureplot.png'),dpi=500, width = 8, fea.pl)
}

DefaultAssay(seurat_NK) <- 'RNA'

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Vinh/Feature_Plot/RNA')

for (i in rna.features) {
  fea.pl <- FeaturePlot(seurat_NK, reduction = 'wnn.umap', features = i)
  ggsave(paste0(i,'_Featureplot.png'),dpi=500, width = 8, fea.pl)
}

### Heatmap ####

seurat_NK <- ScaleData(seurat_NK, assay = 'RNA')
seurat_NK <- ScaleData(seurat_NK, assay = 'ADT')

heatmap.rna <- DoHeatmap(seurat_NK, features = rna.features, size = 3, assay = 'RNA') + 
  scale_fill_gradientn(colors = c("blue", "white", "red"))

heatmap.prot <- DoHeatmap(seurat_NK, features = adt.features, size = 3, assay = 'ADT') + 
  scale_fill_gradientn(colors = c("blue", "white", "red"))

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Vinh/Heatmap')

ggsave('RNA_Heatmap_CellSplit.png',dpi=500, width = 13, heatmap.rna)

ggsave('ADT_Heatmap_CellSplit.png',dpi=500, width = 12, height = 15, heatmap.prot)

### Find All Markers ####
setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Vinh/All_Markers_Comp')

dpe_all <- FindAllMarkers(seurat_NK,assay = 'ADT',test.use = 'MAST')
write.csv(dpe_all, 'ADT_All_Markers.csv', row.names = TRUE)

dge_all <- FindAllMarkers(seurat_NK,assay = 'RNA',test.use = 'MAST')
write.csv(dge_all, 'RNA_All_Markers.csv', row.names = TRUE)



setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Vinh/Proportions_Plot')



### Proportions ####
setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Vinh/Proportions_Plot')

seurat_NK$wsnn_res.0.8

# Create a table showing the number of cells in each cluster per donor
cluster_counts <- seurat_NK@meta.data %>%
  group_by(orig.ident, wsnn_res.0.8) %>%
  summarise(count = n()) %>%
  ungroup()

# Print the table
print(cluster_counts)

# Calculate the proportion of cells in each cluster per donor
cluster_proportions <- cluster_counts %>%
  group_by(orig.ident) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# Print the proportion table
print(cluster_proportions)

# Plot the number of cells in each cluster per donor
count_plot <- ggplot(cluster_counts, aes(x = orig.ident, y = count, fill = as.factor(wsnn_res.0.8))) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Number of Cells in Each Cluster per Donor",
       x = "Donor",
       y = "Number of Cells",
       fill = "Cluster") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot the proportion of cells in each cluster per donor
proportion_plot <- ggplot(cluster_proportions, aes(x = orig.ident, y = proportion, fill = as.factor(wsnn_res.0.8))) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Proportion of Cells in Each Cluster per Donor",
       x = "Donor",
       y = "Proportion of Cells",
       fill = "Cluster") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Export tables to CSV files
write_csv(cluster_counts, "cluster_counts.csv")
write_csv(cluster_proportions, "cluster_proportions.csv")

# Print the plots
ggsave('Absolute_Counts_BySample_PerCluster.png',dpi=500, width = 13, bg='white', count_plot)

ggsave('Prop_Counts_BySample_PerCluster.png',dpi=500, width = 13, bg='white', proportion_plot)



#### Pathway Analysis ####
#### Function to format and create gene list for Pathway Analysis
create_gene_list <- function(group1,group2,group,ident, sub.ident = NULL) {
  Idents(seurat_NK) <- ident
  DefaultAssay(seurat_NK) <- 'RNA'
  dif_genes <- FindMarkers(seurat_isotype, ident.1 = group1, ident.2= group2, group.by =  group,
                           subset.ident = sub.ident ,min.pct=0.1, logfc.threshold=0.5, test.use = 'MAST')
  dif_genes <- subset(dif_genes,p_val_adj < 0.05)
  
  # we want the named vector with log2 fold change
  genelist <- setNames(dif_genes$avg_log2FC,rownames(dif_genes))
  # omit any NA values 
  genelist<-na.omit(genelist)
  # sort the list in decreasing order (required for clusterProfiler)
  genelist = sort(genelist, decreasing = TRUE)
  return(genelist)
}

create_gene_list_Entrez <- function(group1,group2,group,ident, sub.ident = NULL) {
  Idents(seurat_NK) <- ident
  DefaultAssay(seurat_NK) <- 'RNA'
  dif_genes <- FindMarkers(seurat_NK, ident.1 = group1, ident.2= group2, group.by =  group,
                           subset.ident = sub.ident ,min.pct=0.1, logfc.threshold=0.5, test.use = 'MAST')
  dif_genes <- subset(dif_genes,p_val_adj < 0.05)
  ## Ad columns with EntrezID
  x <- bitr(rownames(dif_genes), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db",drop = F)
  #Remove duplicate if present
  x <- subset(x, ENTREZID != "7795")
  dif_genes$SYMBOL <- rownames(dif_genes)
  dif_genes <- merge(dif_genes,x,by="SYMBOL")
  dif_genes <- na.omit(dif_genes) #Remove Missing Values
  genelist <- extract2(dif_genes, 'avg_log2FC') %>% set_names(dif_genes$ENTREZID)
  # sort the list in decreasing order (required for clusterProfiler)
  genelist = sort(genelist, decreasing = TRUE)
  return(genelist)
}



#### Create folders to save output ####
base_dir <- "C:/Users/axi313/Documents/TARA_Entry/WNN/Vinh/Pathway_Analysis"

Idents(seurat_NK)

####################### Working Dirs To Save Output ##########################


all.sample.names <-c('6','8','17','18','31')

for (runx in all.sample.names) {
  
  ### Move to Working Directory and prepare Genelists for gsea and e.gsea
  wd <- paste0('C:/Users/axi313/Documents/TARA_Entry/WNN/Vinh/Pathway_Analysis/',runx)
  setwd(wd)
  genelist <- create_gene_list(group1='HEI',group2='HEU',group = 'Condition',ident='wsnn_res.0.8' , sub.ident = runx)
  egenelist <- create_gene_list_Entrez(group1='HEI',group2='HEU',group = 'Condition',ident='wsnn_res.0.8' , sub.ident = runx)
  # Directory names
  directories <- c("GO", "KEGG", "KEGG_Pathview")
  
  # Create each directory if it doesn't already exist
  for (dir in directories) {
    if (!dir.exists(dir)) {
      dir.create(dir)
    } else {
      message(paste("Directory", dir, "already exists"))
    }
  }
  
  
  ############## GO Analysis ###########
  
  ### GO Enrichment (Overexpression Analysis)
  
  setwd(paste0(wd,'/GO'))
  
  ### GSEA Analysis
  go.gsea <- gseGO(geneList     = genelist,
                   OrgDb        = org.Hs.eg.db,
                   keyType = 'SYMBOL',
                   ont          = "All",
                   minGSSize    = 5,
                   maxGSSize    = 500,
                   pvalueCutoff = 0.05,
                   verbose      = TRUE)
  
  ### Results CSV
  
  df.gsea <- as.data.table(go.gsea)
  write.csv(df.gsea,paste0(runx,'_GO_GSEA_Enrichment.csv'))
  
  ### GO Dotplot
  
  dp.gsea <- dotplot(go.gsea, showCategory=30) + ggtitle("dotplot for GO GSEA")
  ggsave(paste0(runx,'_GO_Dotplot.png'),dpi=500, height =16, width = 13, dp.gsea)
  
  ### Gene-Concept Network
  
  go.gc.p1 <- cnetplot(go.gsea,color.params = list(foldChange = genelist))
  go.gc.p2 <- cnetplot(go.gsea, color.params = list(foldChange = genelist), circular = TRUE, colorEdge = TRUE) 
  go.gc.p3 <- cowplot::plot_grid(go.gc.p1, go.gc.p2, ncol=2, labels=LETTERS[1:2], rel_widths=c(.8, 1.2))
  ggsave(paste0(runx,'_GO_Cnetplot.png'),dpi=500, height =16, width = 21,bg='white', go.gc.p3)
  
  ### Heatplot
  hp.go.gsea <- heatplot(go.gsea, showCategory = 10, foldChange=genelist) 
  ggsave(paste0(runx,'_GO_Heatplot.png'),dpi=500, height =11, width = 28,bg='white', hp.go.gsea)
  
  ### Emap Plot
  go.gsea.tree <- pairwise_termsim(go.gsea)
  emap.go.gsea <- emapplot(go.gsea.tree)
  ggsave(paste0(runx,'_GO_Emmap.png'),dpi=500, height =11, width = 28,bg='white', emap.go.gsea)
  
  
  #### KEGG GSEA ####
  
  setwd(paste0(wd,'/KEGG'))
  
  kegg.gsea <- gseKEGG(
    geneList = egenelist,
    keyType = 'ncbi-geneid',
    organism = "hsa")
  
  ### Results CSV
  df.KEGG <- as.data.table(kegg.gsea)
  write.csv(df.KEGG,paste0(runx,'_KEGG_Enrichment.csv'))
  
  ### Dot Plot
  
  dp.kegg <- dotplot(kegg.gsea, showCategory=30) + ggtitle("dotplot for KEGG GSEA")
  ggsave(paste0(runx,'_KEGG_Dotplot.png'),dpi=500, height =16, width = 13, dp.kegg)
  
  ## Gene-Concept Network
  
  kegg.gsea.2 <- setReadable(kegg.gsea, 'org.Hs.eg.db', 'ENTREZID')
  kegg.p1 <- cnetplot(kegg.gsea.2,color.params = list(foldChange = genelist))
  kegg.p2 <- cnetplot(kegg.gsea.2, color.params = list(foldChange = genelist), circular = TRUE, colorEdge = TRUE) 
  kegg.p3 <- cowplot::plot_grid(kegg.p1, kegg.p2, ncol=2, labels=LETTERS[1:2], rel_widths=c(.8, 1.2))
  ggsave(paste0(runx,'_KEGG_Cnetplot.png'),dpi=500, height =16, width = 21,bg='white', kegg.p3)
  
  ### Heatplot
  
  hp.kegg.gsea <- heatplot(kegg.gsea.2, showCategory = 20, foldChange=genelist) 
  ggsave(paste0(runx,'_KEGG_Heatplot.png'),dpi=500, height =11, width = 28,bg='white', hp.kegg.gsea)
  
  ### Emap Plot
  go.kegg.tree <- pairwise_termsim(kegg.gsea.2)
  
  emap.kegg.gsea <- emapplot(go.kegg.tree)
  ggsave(paste0(runx,'_KEGG_Emmap.png'),dpi=500, height =11, width = 28,bg='white', emap.kegg.gsea)
  
  #### Pathway Overlay
  image_directory <- paste0(wd,'/KEGG_Pathview')
  setwd(image_directory)
  
  # Get Ids of pathways where adj Pval < 0.05
  
  KEGG.pathways <- df.KEGG[df.KEGG$p.adjust <0.05,]$ID
  
  for (KP in KEGG.pathways) {
    tryCatch({
      # Attempt to generate the plot with pathview
      pathview(
        gene.data = egenelist,
        pathway.id = KP,
        species = "hsa",
        limit = list(gene = max(abs(egenelist)), cpd = 1),
        new.signature = FALSE
      )
      dev.off()
    }, error = function(e) {
      # In case of an error, print the pathway ID and the error message
      message(sprintf("Error with pathway %s: %s", KP, e$message))
      # The loop will continue despite the error
    })
  }
  
  #Delete Extra Pathview files
  
  # List all files in the directory
  all_files <- list.files(image_directory, full.names = TRUE)
  
  # Identify the files that do not end with '.pathview.png'
  files_to_delete <- all_files[!grepl("\\.pathview\\.png$", all_files)]
  
  # Remove the identified files
  file.remove(files_to_delete)
  
}
