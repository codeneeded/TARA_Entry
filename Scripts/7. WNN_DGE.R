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


setwd('C:/Users/axi313/Documents/TARA_Entry/WNN')

load.path <- "C:/Users/axi313/Documents/TARA_Entry/saved_R_data/"

load(paste0(load.path,'Seuratv5_WNN_Complete.RData'))

#################################### Heatmaps ##############################################

# create multimodal heatmap 
rna.features <-  c('CD14','FCGR2B','SERPING1','CCR7','CD27','TCF7','CCL5','FCGR3A','PRF1','CD40LG','IRF8','TNFRSF4',
                              'CD8A','TNFRSF9','XCL2','CD7','CD8B','NELL2','C1QBP','CD3E','ICOS','IGFBP2','IGFBP4','LDHA',
                              'CCND3','MIR155HG','NR4A1','CTLA4','FOXP3','IL2RA','CD19','CD79A','IGHM','EBI3','HLA-DPA1',
                              'HLA-DRB1','CTSW','KLRC1','TNFRSF18','CCR4','IRF4','MALAT1','IKZF2','TRDV1','TRGC2',
                              'CD3D','CXCR3','GZMK','CCL2','HLA-DRA','SERPINA1','GNLY','NKG7','TIGIT','LTB','MAL','SELL',
                              'CCL4L2','CD70','IFNG','IL2RB','KLRD1','TRBC1','HAVCR2','LGALS1','NCAM1','CD36','CD4','IFI30',
                              'CXCL8','ITGAX','IL18BP','TNF','TRDV2','TRGV9','FABP5','MT-ND1','MT-ND5','CCL3','IL1B','TNFAIP2',
                              'CD40','MS4A1','XCL1','HIST1H4C','LTA','MKI67')

prots <- rownames(seurat_isotype@assays$ADT@data)

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



DoHeatmap(seurat_isotype, features = rna.features, slot='data') + 
  scale_fill_gradientn(colors = c("blue", "white", "red"))

# Save the heatmap as a file if needed
ggsave("RNA_Seurat_heatmap.png", plot = last_plot(), width = 10, height = 14)
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

########################## CD57 CD28 Status ##############################

seurat_isotype@meta.data$CD57_CD28_Status <- ifelse(
  FetchData(seurat_isotype, vars = c("B3GAT1", "CD28"))$B3GAT1 > 3 & FetchData(seurat_isotype, vars = c("B3GAT1", "CD28"))$CD28 > 3, "CD57+CD28+",
  ifelse(
    FetchData(seurat_isotype, vars = c("B3GAT1", "CD28"))$B3GAT1 < 3 & FetchData(seurat_isotype, vars = c("B3GAT1", "CD28"))$CD28 < 3, "CD57-CD28-",
    ifelse(
      FetchData(seurat_isotype, vars = c("B3GAT1", "CD28"))$B3GAT1 > 3 & FetchData(seurat_isotype, vars = c("B3GAT1", "CD28"))$CD28 < 3, "CD57+CD28-",
      ifelse(
        FetchData(seurat_isotype, vars = c("B3GAT1", "CD28"))$B3GAT1 < 3 & FetchData(seurat_isotype, vars = c("B3GAT1", "CD28"))$CD28 > 3, "CD57-CD28+",
        "Other"
      )
    )
  )
)

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/CD57_CD28_Status')

cd57 <- DimPlot(seurat_isotype, reduction = 'wnn.umap', split.by = 'CD57_CD28_Status',repel=T, cols='polychrome')
ggsave('CD57_CD28_Status.png',dpi=500, width = 12, cd57)

save(seurat_isotype, file=paste0(load.path,"Seuratv5_WNN_Complete.RData"))

################################### Differential Gene Expression #############################################

##### HEI Vs HEU ##############################
seurat_entry <- subset(seurat_isotype, subset = Timepoint =='Entry'& Condition == c('HEI','HEU'))


### RNA

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Differential_Expression/HEIVsHEU/RNA')
DefaultAssay(seurat_entry) <- 'RNA'
for (i in levels(Idents(seurat_entry))) {
  dge <- FindMarkers(seurat_entry,ident.1 = 'HEI', ident.2 = 'HEU', group.by = 'Condition',
                     subset.ident = i, test.use = 'MAST',min.pct = 0.1)
  write.csv(dge, paste0('Cluster_',i,'_DGE_HEI_Vs_HEU.csv'), row.names = TRUE)
}

### ADT

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Differential_Expression/HEIVsHEU/ADT')
DefaultAssay(seurat_entry) <- 'ADT'
for (i in levels(Idents(seurat_entry))) {
  dge <- FindMarkers(seurat_entry,ident.1 = 'HEI', ident.2 = 'HEU', group.by = 'Condition',
                     subset.ident = i, test.use = 'MAST',min.pct = 0.1)
  write.csv(dge, paste0('Cluster_',i,'_DGE_HEI_Vs_HEU.csv'), row.names = TRUE)
}




### Subset

seurat_viral <- subset(seurat_isotype, subset = Timepoint =='Entry'& Condition == 'HEI')

### RNA

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Differential_Expression/HighVL_Vs_LowVL/RNA')
DefaultAssay(seurat_viral) <- 'RNA'
for (i in levels(Idents(seurat_viral))) {
  dge <- FindMarkers(seurat_viral,ident.1 = 'High', ident.2 = 'Low', group.by = 'Viral_Load_Category',
                     subset.ident = i, test.use = 'MAST',min.pct = 0.1)
  write.csv(dge, paste0('Cluster_',i,'_DGE_HighVL_Vs_LowVL.csv'), row.names = TRUE)
}

### ADT

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Differential_Expression/HighVL_Vs_LowVL/ADT')
DefaultAssay(seurat_viral) <- 'ADT'
for (i in levels(Idents(seurat_viral))) {
  dge <- FindMarkers(seurat_viral,ident.1 = 'High', ident.2 = 'Low', group.by = 'Viral_Load_Category',
                     subset.ident = i, test.use = 'MAST',min.pct = 0.1)
  write.csv(dge, paste0('Cluster_',i,'_DGE_HighVL_Vs_LowVL.csv'), row.names = TRUE)
}



########## DGE CTL-expanding vs CTL non-expanding ######



setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Differential_Expression/CTL_expanding_Vs_CTL-non_expanding/RNA')

DefaultAssay(seurat_viral) <- 'RNA'

for (i in levels(Idents(seurat_viral))) {
  dge <- FindMarkers(seurat_viral,ident.1 = 'high', ident.2 = 'low', group.by = 'CTLGrouping',
                     subset.ident = i, test.use = 'MAST',min.pct = 0.1)
  write.csv(dge, paste0('Cluster_',i,'_DGE_High_Vs_Low_CTLGrouping.csv'), row.names = TRUE)
}

### ADT

setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Differential_Expression/CTL_expanding_Vs_CTL-non_expanding/ADT')
DefaultAssay(seurat_viral) <- 'ADT'
for (i in levels(Idents(seurat_viral))) {
  dge <- FindMarkers(seurat_viral,ident.1 = 'high', ident.2 = 'low', group.by = 'CTLGrouping',
                     subset.ident = i, test.use = 'MAST',min.pct = 0.1)
  write.csv(dge, paste0('Cluster_',i,'_DGE_High_Vs_Low_CTLGrouping.csv'), row.names = TRUE)
}
 
