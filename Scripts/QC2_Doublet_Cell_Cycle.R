###Citeseq Pipeline
#Cite seurat and ds packages
#Load Required Libraries
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
library(DoubletFinder)
library(ggplot2)
library(gridExtra)
library(SeuratWrappers)
library(Nebulosa)
library(tricycle)
library(org.Hs.eg.db)
library(Azimuth)
library(scDblFinder)

# Data Input, Creation of Merged Seurat Object

##set path to load data

setwd("C:/Users/axi313/Documents/TARA_Entry/Preprocessing/Cell_Cycle")
load.path <- "C:/Users/axi313/Documents/TARA_Entry/saved_R_data/"

# Load Data

load(paste0(load.path,'Seuratv5_filtered_seurat.RData'))
load(paste0(load.path,'Cell_Cycle_Genes.RData'))

###################################Complete Initialisation Steps ######################################################
DefaultAssay(filtered_seurat) <-'RNA'

### Azimuth Analysis
seurat_phase <- RunAzimuth(filtered_seurat, reference = "pbmcref")

### Basic Transformations
seurat_phase <- NormalizeData(seurat_phase)
seurat_phase <- FindVariableFeatures(seurat_phase, selection.method = "vst")
seurat_phase <- ScaleData(seurat_phase, features = rownames(seurat_phase))
seurat_phase <- RunPCA(seurat_phase, features = VariableFeatures(seurat_phase), ndims.print = 6:10, nfeatures.print = 10)
seurat_phase <- FindNeighbors(seurat_phase, dims = 1:30, reduction = "pca")
seurat_phase <- FindClusters(seurat_phase, resolution = 2, cluster.name = "unintegrated_clusters")
seurat_phase <- RunUMAP(seurat_phase, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
seurat_phase$predicted.celltype.l2
DimPlot(seurat_phase, group.by = 'predicted.celltype.l2', label = T)

############################################## EXPLORE CELL CYCLE ########################################################

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

#Evaluate effects of Cell Cycle on Data

# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m.genes, 
                                 s.features = s.genes,
                                 set.ident = TRUE)

RidgePlot(seurat_phase, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

seurat_phase <- RunPCA(seurat_phase, features = c(s.genes, g2m.genes))
DimPlot(seurat_phase)

# View cell cycle scores and phases assigned to cells 
png(file="Cell_Cycle_Phase_PerSample.png", width = 900, height = 600)
VlnPlot(seurat_phase, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", pt.size = 0.1)
dev.off()

# Plot the PCA colored by cell cycle phase
png(file="PCA_by_Cell_Cycle_Split.png", width = 900, height = 600)
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
dev.off()

png(file="PCA_by_Cell_Cycle.png", width = 900, height = 600)
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase")
dev.off()

#Check effect of Condition
png(file="PCA_by_Cell_Cycle_Split_by_phasae_groupby_Condition.png", width = 900, height = 600)
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Condition",
        split.by = "Phase")
dev.off()

################################################### Tricycle ########################################################

seurat_phase <- Runtricycle(object = seurat_phase, assay = 'RNA', slot = "data", reduction.name = "tricycleEmbedding", 
                            reduction.key = "tricycleEmbedding_", gname = NULL,
                            gname.type = "SYMBOL", species = "human")

### Visualisation

plot.df <- FetchData(object = seurat_phase, vars = c("tricycleEmbedding_1", "tricycleEmbedding_2", 
                                                     "tricyclePosition", "TOP2A"))
names(plot.df)[4] <- "Top2a"

###
p <- tricycle:::.plot_emb_circle_scale(emb.m = plot.df[, 1:2], color.value = plot.df$tricyclePosition, 
                                       color_by = "tricyclePosition", point.size = 1.5, point.alpha = 2) + theme_light()
legend <- circle_scale_legend(text.size = 5, alpha = 3)
p.1 <- plot_grid(p, legend, ncol = 2, rel_widths = c(1, 0.4)) + theme_light()

ggsave("Tricycle_Embedings.png",dpi = 400,p.1)

#PCA
plot.df <- FetchData(object = seurat_phase, vars = c("PC_1", "PC_2", 
                                                     "tricyclePosition", "TOP2A"))

p <- tricycle:::.plot_emb_circle_scale(emb.m = plot.df[, 1:2], color.value = plot.df$tricyclePosition,
                                       color_by = "tricyclePosition", point.size = 1.5, point.alpha = 0.8)+
  theme_bw(base_size = 14)
legend <- circle_scale_legend(text.size = 5, alpha = 0.8)
p.1 <- grid.arrange(p, legend, ncol = 2,heights=c(5), widths=c(4,2)) 


ggsave("Tricycle_Embedings_PCA.png",dpi = 500,p.1)


### By Condition

plot.df.hivpos <- FetchData(object = subset(seurat_phase,subset= `condition` == 'HIV+'), vars = c("PC_1", "PC_2", 
                                                                                                  "tricyclePosition", "TOP2A"))
plot.df.hivneg <- FetchData(object = subset(seurat_phase,subset= `condition` == 'HIV-'), vars = c("PC_1", "PC_2", 
                                                                                                  "tricyclePosition", "TOP2A"))
p.hp <- tricycle:::.plot_emb_circle_scale(emb.m = plot.df.hivpos[, 1:2], color.value = plot.df.hivpos$tricyclePosition,
                                          color_by = "tricyclePosition", point.size = 1.5, point.alpha = 0.8)+
  theme_bw(base_size = 14)
p.hn <- tricycle:::.plot_emb_circle_scale(emb.m = plot.df.hivneg[, 1:2], color.value = plot.df.hivneg$tricyclePosition,
                                          color_by = "tricyclePosition", point.size = 1.5, point.alpha = 0.8)+
  theme_bw(base_size = 14)
legend <- circle_scale_legend(text.size = 5, alpha = 0.8)
p.hp.f <- grid.arrange(p.hp, legend, ncol = 2,heights=c(5), widths=c(4,2), top = grid::textGrob('HIV+', gp=grid::gpar(fontsize=18)))
p.hn.f <- grid.arrange(p.hn, legend, ncol = 2,heights=c(5), widths=c(4,2), top = grid::textGrob('HIV-', gp=grid::gpar(fontsize=18)))
vsplotscondition <- grid.arrange(p.hp.f,p.hn.f)

ggsave("Tricycle_Embedings_PCA_byCondition.png",dpi = 500,vsplotscondition)

#####################################################################################################################
##As an alternative, we suggest regressing out the difference between the G2M and S phase scores. 
##This means that signals separating non-cycling cells and cycling cells will be maintained, 
##but differences in cell cycle phase among proliferating cells (which are often uninteresting), 
##will be regressed out of the data

############################################Doublet Finder###########################################################
### 
#save(seurat_phase, file = paste0(load.path,'SeuratV5_SeuratPhase.RData'))
#load(paste0(load.path,'SeuratV5_SeuratPhase.RData'))
# Split Seurat Object For Doublet Filtration

split_seurat <- SplitObject(seurat_phase, split.by = "orig.ident")

samples <-  levels(as.factor(seurat_phase$orig.ident))

setwd("C:/Users/axi313/Documents/TARA_Entry/Preprocessing/Doublets")

for (i in samples) {
  sce <- scDblFinder(GetAssayData(split_seurat[[i]], slot="counts"))
  split_seurat[[i]]$scDblFinder.score <- sce$scDblFinder.score
  split_seurat[[i]]$scDblFinder.class <- sce$scDblFinder.class
  plot <- DimPlot(split_seurat[[i]],group.by = 'scDblFinder.class') + ggtitle(paste0(i,' Doublets'))
  ggsave(paste0(i,'_doublets.png'),device = 'png', width = 6, 
         height = 6, dpi = 600,plot)
}

save(split_seurat, file = paste0(load.path,'SeuratV5_SplitSeurat_Preintegration_WithDoublets.RData'))

### Remove Doublets

for (i in samples) {
  split_seurat[[i]]<- subset(split_seurat[[i]],subset= `scDblFinder.class`== 'singlet')
}

save(split_seurat, file = paste0(load.path,'SeuratV5_SplitSeurat_Preintegration.RData'))
