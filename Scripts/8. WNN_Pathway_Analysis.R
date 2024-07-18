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
library(reticulate)
library(clusterProfiler)
library(BiocParallel)
library(fgsea)
library(DOSE)
library(enrichplot)
library(org.Hs.eg.db)
library(pathview)
library(magrittr)
##set path to load data


setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Pathway_Analysis')

load.path <- "C:/Users/axi313/Documents/TARA_Entry/saved_R_data/"

load(paste0(load.path,'Seuratv5_WNN_Complete.RData'))

#### Pathway Analysis HIgh Vs Low Viral####

seurat_viral <- subset(seurat_isotype, subset = Timepoint =='Entry'& Condition == 'HEI')


#### Function to format and create gene list for Pathway Analysis
create_gene_list <- function(x_path) {
  dif_genes <- read.csv(x_path)
  dif_genes <- subset(dif_genes,p_val_adj < 0.05)
  
  # we want the named vector with log2 fold change
  genelist <- setNames(dif_genes$avg_log2FC,dif_genes$X)
  # omit any NA values 
  genelist<-na.omit(genelist)
  # sort the list in decreasing order (required for clusterProfiler)
  genelist = sort(genelist, decreasing = TRUE)
  return(genelist)
}


create_gene_list_Entrez <- function(x_path) {
  dif_genes <- read.csv(x_path)
  dif_genes <- subset(dif_genes,p_val_adj < 0.05)
  ## Ad columns with EntrezID
  x <- bitr(dif_genes$X, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db",drop = F)
  #Remove duplicate if present
  x <- subset(x, ENTREZID != "7795")
  dif_genes$SYMBOL <- dif_genes$X
  dif_genes <- merge(dif_genes,x,by="SYMBOL")
  dif_genes <- na.omit(dif_genes) #Remove Missing Values
  genelist <- extract2(dif_genes, 'avg_log2FC') %>% set_names(dif_genes$ENTREZID)
  # sort the list in decreasing order (required for clusterProfiler)
  genelist = sort(genelist, decreasing = TRUE)
  return(genelist)
}



#### Create folders to save output ####
base_dir <- "C:/Users/axi313/Documents/TARA_Entry/WNN/Vinh/Pathway_Analysis"

Idents(seurat_viral)

####################### Working Dirs To Save Output ##########################


### Viral Grouping #######

all.sample.names <- levels(Idents(seurat_viral))[0:30]
viral_path <- 'C:/Users/axi313/Documents/TARA_Entry/WNN/Differential_Expression/HighVL_Vs_LowVL/RNA/'

for (runx in c('13','14','15','16','17','18','19','20')) {
  
  setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Pathway_Analysis/Viral_Load')
  dir.create(runx, showWarnings = FALSE)
  cat("Directory", runx, "created or already exists.\n")
  
  ### Move to Working Directory and prepare Genelists for gsea and e.gsea
  
  wd <- paste0('C:/Users/axi313/Documents/TARA_Entry/WNN/Pathway_Analysis/Viral_Load/',runx)
  setwd(wd)
  genelist <- create_gene_list(paste0(viral_path,'Cluster_',runx,'_DGE_HighVL_Vs_LowVL.csv'))
  egenelist <- create_gene_list_Entrez(paste0(viral_path,'Cluster_',runx,'_DGE_HighVL_Vs_LowVL.csv'))
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



### CTL Grouping ########

all.sample.names <- levels(Idents(seurat_viral))[0:30]
viral_path <- 'C:/Users/axi313/Documents/TARA_Entry/WNN/Differential_Expression/CTL_expanding_Vs_CTL-non_expanding/RNA/'

for (runx in c('14','15','16','17','18','19','20')) {
  
  setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Pathway_Analysis/CTL_grouping')
  dir.create(runx, showWarnings = FALSE)
  cat("Directory", runx, "created or already exists.\n")
  
  ### Move to Working Directory and prepare Genelists for gsea and e.gsea
  
  wd <- paste0('C:/Users/axi313/Documents/TARA_Entry/WNN/Pathway_Analysis/CTL_grouping/',runx)
  setwd(wd)
  genelist <- create_gene_list(paste0(viral_path,'Cluster_',runx,'_DGE_High_Vs_Low_CTLGrouping.csv'))
  egenelist <- create_gene_list_Entrez(paste0(viral_path,'Cluster_',runx,'_DGE_High_Vs_Low_CTLGrouping.csv'))
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

################### HEI vs HEU #####################################

seurat_entry <- subset(seurat_isotype, subset = Timepoint =='Entry'& Condition == c('HEI','HEU'))


### Viral Grouping

all.sample.names <- levels(Idents(seurat_entry))[0:30]
viral_path <- 'C:/Users/axi313/Documents/TARA_Entry/WNN/Differential_Expression/HEIvsHEU/RNA/'

for (runx in all.sample.names) {
  
  setwd('C:/Users/axi313/Documents/TARA_Entry/WNN/Pathway_Analysis/HEIvsHEU')
  dir.create(runx, showWarnings = FALSE)
  cat("Directory", runx, "created or already exists.\n")
  
  ### Move to Working Directory and prepare Genelists for gsea and e.gsea
  
  wd <- paste0('C:/Users/axi313/Documents/TARA_Entry/WNN/Pathway_Analysis/HEIvsHEU/',runx)
  setwd(wd)
  genelist <- create_gene_list(paste0(viral_path,'Cluster_',runx,'_DGE_HEI_Vs_HEU.csv'))
  egenelist <- create_gene_list_Entrez(paste0(viral_path,'Cluster_',runx,'_DGE_HEI_Vs_HEU.csv'))
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
