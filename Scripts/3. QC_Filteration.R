
# Quality Control Visualizations
##Ref -> https://www.singlecellcourse.org/scrna-seq-analysis-with-bioconductor.html
##Ref -> http://bioconductor.org/books/3.15/OSCA.intro/getting-scrna-seq-datasets.html
#Ref -> https://hbctraining.github.io/scRNA-seq_online/schedule/links-to-lessons.html

library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(hdf5r)
library(data.table)
library(ggplot2)
library(biomaRt)

############################################ Global Variables #######################################################
setwd("C:/Users/axi313/Documents/TARA_Entry/Preprocessing/Post_DSB_QC")
in.path <- "C:/Users/axi313/Documents/TARA_Entry/Raw_Data/"
out.path <- "C:/Users/axi313/Documents/TARA_Entry/saved_R_data/"


load(paste0(out.path,"Seuratv5_CITEseq_dsbnorm_merged_Seurat.RData"))

############################################## Cell Level QC #########################################################

#Post dsb QC, pre-final filteration

colnames(merged_seurat@meta.data)[4] <- "nCounts_ADT"
colnames(merged_seurat@meta.data)[5] <- "nFeatures_ADT"
colnames(merged_seurat@meta.data)

## Number of Cells per Sample
#The cell numbers can also vary by protocol, producing cell numbers that are much higher than what we loaded. 
#For example, during the inDrops protocol, the cellular barcodes are present in the hydrogels, which are encapsulated in the droplets 
#with a single cell and lysis/reaction mixture. While each hydrogel should have a single cellular barcode associated with it, 
#occasionally a hydrogel can have more than one cellular barcode. 
#Similarly, with the 10X protocol there is a chance of obtaining only a barcoded bead in the emulsion droplet (GEM) and no actual cell.

metadata <- merged_seurat@meta.data

png(file="Cells_per_sample.png", width = 900, height = 600)
metadata %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
dev.off()

png(file="Cells_per_Condition.png", width = 900, height = 600)
metadata %>% 
  ggplot(aes(x=Condition, fill=orig.ident)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
dev.off()

# QC Features rough look

feats.1 <- c("nCount_RNA", "nFeature_RNA","nCounts_ADT", "nFeatures_ADT",
             "percent_mito","percent_ribo","percent_hb","percent_plat")



png(file="Pre-QC_features_grouped.png", width = 900, height = 800)
VlnPlot(merged_seurat, group.by = "orig.ident", features = feats.1, pt.size = 0.1, ncol = 2) +
  NoLegend()
dev.off()

#UMI Count
#Should generally be above 500, that is the low end of what we expect.
#If UMI counts are between 500-1000 counts, it is usable but the cells probably should have been sequenced more deeply.

png(file="UMI_Count.png", width = 900, height = 600)
metadata %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
dev.off()

#nGenes
#For high quality data, the proportional histogram should contain a single large peak that represents cells that were encapsulated. 
#If we see a small shoulder to the left of the major peak (not present in our data), or a bimodal distribution of the cells, 
#that can indicate a couple of things. It might be that there are a set of cells that failed for some reason. 
#It could also be that there are biologically different types of cells (i.e. quiescent cell populations, less complex cells of interest), 
#and/or one type is much smaller than the other (i.e. cells with high counts may be cells that are larger in size).

png(file="nGenes.png", width = 900, height = 600)
metadata %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 250)
dev.off()

#Complexity Score
## Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
#We can evaluate each cell in terms of how complex the RNA species are by using a measure called the novelty score. 
#The novelty score is computed by taking the ratio of nGenes over nUMI. If there are many captured transcripts (high nUMI) 
#and a low number of genes detected in a cell, this likely means that you only captured a low number of genes and simply sequenced transcripts 
#from those lower number of genes over and over again. These low complexity (low novelty) cells could represent a specific cell type 
#(i.e. red blood cells which lack a typical transcriptome), or could be due to an artifact or contamination. 
#Generally, we expect the novelty score to be above 0.80 for good quality cells.

png(file="Complexity_Score.png", width = 900, height = 600)
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
dev.off()

#Mito Ratio
#This metric can identify whether there is a large amount of mitochondrial contamination from dead or dying cells. 
#We define poor quality samples for mitochondrial counts as cells which surpass the 0.2 mitochondrial ratio mark, 
#unless of course you are expecting this in your sample.
# Visualize the distribution of mitochondrial gene expression detected per cell

png(file="Mito_Ratio.png", width = 900, height = 600)
metadata %>% 
  ggplot(aes(color=orig.ident, x=percent_mito, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 15)
dev.off()

# Ribo Ratio

png(file="Ribo_Ratio.png", width = 900, height = 600)
metadata %>% 
  ggplot(aes(color=orig.ident, x=percent_ribo, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 5)
dev.off()

#Heme Ratio

png(file="Heme_Ratio.png", width = 900, height = 600)
metadata %>% 
  ggplot(aes(color=orig.ident, x=percent_hb, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 20)
dev.off()

#Platelet Ratio

png(file="Platlet_Ratio.png", width = 900, height = 600)
metadata %>% 
  ggplot(aes(color=orig.ident, x=percent_plat, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 2)
dev.off()

#################################################### FILTERING ######################################################
#Cells with at least 200 detected genes and genes need to be expressed in at least 3 cells.

merged_seurat <- JoinLayers(merged_seurat)
selected_c <- WhichCells(merged_seurat, expression = nFeature_RNA > 200)
selected_f <- rownames(merged_seurat)[Matrix::rowSums(merged_seurat@assays$RNA@layers$counts) > 3]


#Proteins
DefaultAssay(merged_seurat) <- 'ADT'
selected_p <- rownames(merged_seurat)[Matrix::rowSums(merged_seurat) > 3]


filtered_seurat <- subset(merged_seurat, features = c(selected_f,selected_p), cells = selected_c)

DefaultAssay(merged_seurat) <- 'RNA'
DefaultAssay(filtered_seurat) <- 'RNA'


# Filter out low quality cells using selected thresholds -> these will change with experiment
##Cell Level Filtering
#nUMI > 500
#nGene > 250
#log10GenesPerUMI > 0.8
#mitoRatio < 20%
#riboratio > 5%
#Heme <20%
#Platlet <2%


filtered_seurat <- subset(x = filtered_seurat, 
                          subset= (nCount_RNA >= 500) & 
                            (nFeature_RNA >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (percent_mito < 15) &
                            (percent_ribo > 5) &
                            (percent_hb < 20) &
                            (percent_plat < 2)
)


#Within our data we will have many genes with zero counts. 
#These genes can dramatically reduce the average expression for a cell and so we will remove them from our data. 

# Extract counts
counts <- LayerData(object = filtered_seurat, layer = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0


#filtering by prevalence. If a gene is only expressed in a handful of cells, it is not particularly meaningful 
#as it still brings down the averages for all other cells it is not expressed in. 
#For our data we choose to keep only genes which are expressed in 10 or more cells. 
#By using this filter, genes which have zero counts in all cells will effectively be removed.

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
#filtered_seurat@assays$RNA@layers$counts <- filtered_counts

# Reassign to filtered Seurat object
filtered_seurat[['RNA']] <- CreateAssayObject(filtered_counts)


# Cells Post QC
png(file="Post-QC_Cells_per_sample.png", width = 1100, height = 800)
filtered_seurat@meta.data %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
dev.off()

png(file="Post-QC_features_grouped.png", width = 900, height = 800)
VlnPlot(filtered_seurat, group.by = "orig.ident", features = feats.1, pt.size = 0.1, ncol = 2) +
  NoLegend()
dev.off()



f_metadata <- filtered_seurat@meta.data

### Gene based Filtering

par(mar = c(6, 8, 2, 1))
C <- filtered_seurat@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]

png(file="Post-QC_Genes_Expressed_per_cell.png", width = 900, height = 600)
par(mar = c(4, 8, 2, 1))
boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
dev.off()


# Check Cell Counts Pre and Post Filteration
table(merged_seurat$orig.ident)
table(filtered_seurat$orig.ident)

######################################## Sex Determination #######################################################
genes.file = 'C:/Users/axi313/Desktop/Analysis/Pediatric/Scripts/Gene_Files/genes.table.csv'

if (!file.exists(genes.file)) {
  suppressMessages(require(biomaRt))
  
  # initialize connection to mart, may take some time if the sites are
  # unresponsive.
  mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  
  # fetch chromosome info plus some other annotations
  genes.table <- try(biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                                   "description", "gene_biotype", "chromosome_name", "start_position"), mart = mart,
                                    useCache = F))
  
  if (!dir.exists("C:/Users/axi313/Desktop/Analysis/Pediatric/Scripts/Gene_Files/")) {
    dir.create("C:/Users/axi313/Desktop/Analysis/Pediatric/Scripts/Gene_Files/")
  }
  if (is.data.frame(genes.table)) {
    write.csv(genes.table, file = genes.file)
  }
  
  if (!file.exists(genes.file)) {
    download.file("https://raw.githubusercontent.com/NBISweden/workshop-scRNAseq/master/labs/misc/genes.table.csv",
                  destfile = "data/results/genes.table.csv")
    genes.table = read.csv(genes.file)
  }
  
} else {
  genes.table = read.csv(genes.file)
}
Idents(filtered_seurat) <- 'orig.ident'
genes.table <- genes.table[genes.table$external_gene_name %in% rownames(filtered_seurat),]
chrY.gene = genes.table$external_gene_name[genes.table$chromosome_name == "Y"]

filtered_seurat$pct_chrY = colSums(filtered_seurat@assays$RNA@counts[chrY.gene, ])/colSums(filtered_seurat@assays$RNA@counts)
FeatureScatter(filtered_seurat, feature1 = "XIST", feature2 = "pct_chrY")

png(file="Chromosome_Sex_Determination.png", width = 900, height = 600)
VlnPlot(filtered_seurat, features = c("XIST", "pct_chrY"))
dev.off()
levels(as.factor(filtered_seurat$orig.ident))

### Add gender to Meadata
# Extract the metadata
metadata <- filtered_seurat@meta.data

# Create a new column 'Sex' and set default values to "Female"
metadata$Sex <- "Female"

# Define the entries that should be labeled as "Male"
male_entries <- c("CE037_entry", "CP003_entry", "CP013_entry", "CP042_entry")

# Update the 'Sex' column for the specified 'orig.ident' values
metadata$Sex[metadata$orig.ident %in% male_entries] <- "Male"

# Assign the updated metadata back to the Seurat object
filtered_seurat@meta.data <- metadata

save(filtered_seurat, file="Seuratv5_filtered_seurat.RData")

