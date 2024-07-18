###Citeseq Pipeline
#Cite seurat and ds packages
#Load Required Libraries
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(hdf5r)
library(dsb)
library(data.table)
library(ggplot2)

###############################################Path + GLobal Variables ##################################################
##set path to load data
setwd("C:/Users/axi313/Documents/TARA_Entry/Preprocessing")
in.path <- "C:/Users/axi313/Documents/TARA_Entry/Raw_Data/"
out.path <- "C:/Users/axi313/Documents/TARA_Entry/saved_R_data/"
r.str <- "/multi/GEX_FB/raw_feature_bc_matrix.h5"
f.str <- "/per_sample_outs/GEX_FB/sample_filtered_feature_bc_matrix.h5"

### 
##Load 10x data as String names
# Function to get all folder names within a specified path
get_folder_names <- function(in.path) {
  # List all directories within the specified path without recursion
  folder_names <- list.dirs(in.path, full.names = FALSE, recursive = FALSE)
  
  # Filter out the root path itself if present
  folder_names <- folder_names[folder_names != ""]
  
  return(folder_names)
}
list.dirs(in.path, full.names = FALSE, recursive = FALSE)

# Get folder names
f_names <- get_folder_names(in.path)


####################################### Load Data + Basic Pre-processing ##############################################################
for (i in f_names) {
  #Read in Raw and filtered data from 10x output
  raw <- Read10X_h5(paste0(in.path,i,r.str))
  cells<- Read10X_h5(paste0(in.path,i,f.str))
  raw$`Antibody Capture`@Dimnames[[1]] <-  gsub(pattern = "_TotalSeqC", replacement = "", 
                                                x = raw$`Antibody Capture`@Dimnames[[1]])
  cells$`Antibody Capture`@Dimnames[[1]] <-  gsub(pattern = "_TotalSeqC", replacement = "", 
                                                  x = cells$`Antibody Capture`@Dimnames[[1]])
  # define a vector of cell-containing barcodes and remove them from unfiltered data 
  stained_cells <- colnames(cells$`Gene Expression`)
  background <- setdiff(colnames(raw$`Gene Expression`), stained_cells)
  
  # split the data into separate matrices per assay 
  prot <- raw$`Antibody Capture`
  rna <- raw$`Gene Expression`
  
  # create metadata of droplet QC stats used in standard scRNAseq processing
  rna.size <- log10(Matrix::colSums(rna))
  prot.size <- log10(Matrix::colSums(prot))
  nCount_RNA <- Matrix::colSums(rna) #Molecules per cell
  nCount_ADT <- Matrix::colSums(prot) #Molecules per cell
  nFeature_RNA <- Matrix::colSums(rna > 0) #genes per cell
  nFeature_ADT <- Matrix::colSums(prot > 0) #proteins per cell
  mtgene <- grep(pattern = "^MT-", rownames(rna), value = TRUE)
  mt.prop <- Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna)
  md <- as.data.frame(cbind( nCount_RNA, nFeature_RNA, nCount_ADT, nFeature_ADT, rna.size, prot.size,mt.prop))
  
  # add indicator for barcodes Cell Ranger called as cells
  md$drop.class <- ifelse(rownames(md) %in% stained_cells, 'cell', 'background')
  
  # remove barcodes with no evidence of capture in the experiment
  md <- md[md$rna.size > 0 & md$prot.size > 0, ]
  
  ##FIX PROTEIN NAMES
  prot@Dimnames[[1]][prot@Dimnames[[1]] == "human(TCR)"] <- "TCR-AB"
  prot@Dimnames[[1]][prot@Dimnames[[1]] == "human(TCR) "] <- "TCR-vA7.2"
  prot@Dimnames[[1]][prot@Dimnames[[1]] == "Human-TCR "] <- "TCR-vA7.2"
  prot@Dimnames[[1]][prot@Dimnames[[1]] == "human(TCR).1"] <- "TCR-vD2"
  prot@Dimnames[[1]][prot@Dimnames[[1]] == "Human-TCR"] <- "TCR-vD2"
  prot@Dimnames[[1]][prot@Dimnames[[1]] == "PTPRC"] <- "CD45RA"
  prot@Dimnames[[1]][prot@Dimnames[[1]] == "PTPRC.1"] <- "CD45RO"
  prot@Dimnames[[1]][prot@Dimnames[[1]] == "PTPRC.2"] <- "CD45"
  prot@Dimnames[[1]][prot@Dimnames[[1]] == "human(Ig)"] <- "Human-Ig"
  #Assign to final output
  assign(paste0(i,'.md'),md)
  assign(paste0(i,'.prot'),prot)
  assign(paste0(i,'.rna'),rna)
}


###  Droplet Settings

for (i in f_names) {
  #Set metadata object and protein object to prevent constant eval parse calling  
  md <- eval(parse(text= paste0(i,'.md')))
  prot <- eval(parse(text = paste0(i,'.prot')))
  rna <- eval(parse(text = paste0(i,'.rna')))
  
  #Output Plot for detected genes vs the protein library size for cells vs background drops
  png(paste0('Droplet_Thresholds/',i,'_genevsprotlibsize.png'),width = 800, height = 600)
  p <- ggplot(md, aes(x = log10(nFeature_RNA), y = prot.size )) + 
    theme_bw() + 
    geom_bin2d(bins = 300) + 
    scale_fill_viridis_c(option = "C") + 
    facet_wrap(~drop.class)
  print(p)
  dev.off()
}


############################################### Droplet Thresholds for DSB Norm###################################################

# Define the function to process each entry
process_entry <- function(entry_name) {
  # Dynamically create the variable name for the data frame
  md_var_name <- paste0(entry_name, ".md")
  
  # Get the data frame from the variable name
  md <- get(md_var_name)
  
  # Filter the rows based on the given conditions
  filtered_md <- md[md$prot.size > 0.5 & md$prot.size < 4 & md$rna.size < 2.5, ]
  
  # Extract the row names of the filtered rows
  background_drops <- rownames(filtered_md)
  
  # Dynamically assign the result to the new variable
  assign(paste0(entry_name, ".background_drops"), background_drops, envir = .GlobalEnv)
}


# Loop through each entry and process it
for (entry_name in f_names) {
  process_entry(entry_name)
}


###################################Initial QC + DSB Normalization###################################################


for (i in f_names) {
  #Set metadata object and protein object to prevent constant eval parse calling  
  md <- eval(parse(text= paste0(i,'.md')))
  prot <- eval(parse(text = paste0(i,'.prot')))
  rna <- eval(parse(text = paste0(i,'.rna')))
  background_drops <- eval(parse(text = paste0(i,'.background_drops')))
  
  background.adt.mtx = as.matrix(prot[ , background_drops])
  
  cellmd = md[md$drop.class == 'cell', ]  #Define Cell Metadata on only cells not background  
  
  # filter drops with + / - 3 median absolute deviations from the median library size
  rna.mult = (3*mad(cellmd$rna.size))
  prot.mult = (3*mad(cellmd$prot.size))
  rna.lower = median(cellmd$rna.size) - rna.mult
  rna.upper = median(cellmd$rna.size) + rna.mult
  prot.lower = median(cellmd$prot.size) - prot.mult
  prot.upper = median(cellmd$prot.size) + prot.mult
  
  # filter rows based on droplet qualty control metrics
  qc_cells = rownames(
    cellmd[cellmd$prot.size > prot.lower & 
             cellmd$prot.size < prot.upper & 
             cellmd$rna.size > rna.lower & 
             cellmd$rna.size < rna.upper & 
             cellmd$mt.prop < 0.25, ]
  )
  
  # Output thresholds for quality control metrics as in any standard scRNAseq analysis
  png(paste0(i,'_qc_thresholds.png'),width = 800, height = 600)
  plot_aes = list(theme_bw(), geom_point(shape = 21 , stroke = 0, size = 0.7), scale_fill_viridis_c(option = "C"))
  p1 = ggplot(cellmd, aes(x = rna.size )) + geom_histogram(bins = 50) + theme_bw() + xlab("log10 RNA library size")
  p2 = ggplot(cellmd, aes(x = mt.prop)) + geom_histogram(bins = 50) + theme_bw() + xlab("mitochondrial read proportion")
  p3 = ggplot(cellmd, aes(x = log10(nFeature_RNA), y = rna.size, fill = mt.prop )) + plot_aes
  p4 = ggplot(cellmd, aes(x = nFeature_RNA, y = prot.size, fill = mt.prop )) + plot_aes
  print (p1+p2+p3+p4)
  dev.off()
  cell.adt.raw = as.matrix(prot[ , qc_cells])
  cell.rna.raw = rna[ ,qc_cells]
  cellmd = cellmd[qc_cells, ]
  
  #Proteins without Staining
  pm = sort(apply(cell.adt.raw, 1, max))
  pm2 = apply(background.adt.mtx, 1, max)
  head(pm2)
  
  #Assign it to your final output
  assign(paste0(i,'.cell.adt.raw'),cell.adt.raw)
  assign(paste0(i,'.cell.rna.raw'),cell.rna.raw)
  assign(paste0(i,'.background.adt.mtx'),background.adt.mtx)
  assign(paste0(i,'.cellmd'),cellmd)
  assign(paste0(i,'.pm'),pm)
}


# Check if you need to remove proteins without staining
#https://www.rdocumentation.org/packages/dsb/versions/0.3.0
#prot.expres.total <- rbindlist(adt.list)

#In this case we do not

### DSB Normalisation
#Set isotype control

isotype.controls <- c('Mouse-IgG1', 'Mouse-IgG2a','Mouse-IgG2b', 'Rat-IgG2b'
                      ,'Rat-IgG1', 'Rat-IgG2a', 'Hamster-IgG')


#normalize protein data for the cell containing droplets with the dsb method. 

for (i in f_names) {
  cell.adt.raw <- eval(parse(text = paste0(i,'.cell.adt.raw')))
  background.adt.mtx <- eval(parse(text = paste0(i,'.background.adt.mtx')))
  # normalize and denoise with dsb with 
  cells.dsb.norm = DSBNormalizeProtein(
    cell_protein_matrix = cell.adt.raw, 
    empty_drop_matrix = background.adt.mtx, 
    denoise.counts = TRUE, 
    use.isotype.control = TRUE, 
    isotype.control.name.vec = isotype.controls
  )
  cells.dsb.norm <- Matrix(as.matrix(cells.dsb.norm),sparse=TRUE)
  assign(paste0(i,'.cells.dsb.norm'),cells.dsb.norm)
}


########################################Create + Merge Seurat Object (norm dsb+RNA)########################################

# Create Seurat Object

for (i in f_names) {
  cellmd <- eval(parse(text= paste0(i,'.cellmd')))
  cell.adt.raw <- eval(parse(text= paste0(i,'.cell.adt.raw')))
  cells.dsb.norm <- eval(parse(text= paste0(i,'.cells.dsb.norm')))
  cell.rna.raw <- eval(parse(text= paste0(i,'.cell.rna.raw')))
  
  # integrating with Seurat
  stopifnot(isTRUE(all.equal(rownames(cellmd), colnames(cell.adt.raw))))
  stopifnot(isTRUE(all.equal(rownames(cellmd), colnames(cell.rna.raw))))
  
  # create Seurat object note: min.cells is a gene filter, not a cell filter
  s = Seurat::CreateSeuratObject(counts = cell.rna.raw, 
                                 meta.data = cellmd,
                                 assay = "RNA", 
                                 min.cells = 20)
  
  # add dsb normalized matrix "dsb_norm_prot" to the "CITE" assay data slot
  s[["ADT"]] <- CreateAssayObject(data = cells.dsb.norm)
  s$orig.ident <- i
  #Assign
  assign(paste0(i,'.seurat'),s)
}


#Merge Seurat Objects

merged_seurat <- merge(
  x = CE021_entry.seurat, 
  y = c(CE025_entry.seurat, CE037_entry.seurat, CP002_entry.seurat, CP003_entry.seurat, CP006_12m.seurat, CP006_entry.seurat, 
    CP011_entry.seurat, CP013_entry.seurat, CP016_entry.seurat, CP017_entry.seurat, CP018_entry.seurat, CP022_entry.seurat, 
    CP042_entry.seurat, CS005_entry.seurat, CS015_entry.seurat, SAAH29_entry.seurat, SATY021_entry.seurat), 
  add.cell.id = f_names)


###################################################### Edit Seurat Metadata ###############################################


### Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

### Compute percent mitochondrial genes per cell
merged_seurat <- PercentageFeatureSet(merged_seurat, pattern = "^MT-", col.name = "percent_mito")

### Compute percentage of ribosomal genes per cell
merged_seurat <- PercentageFeatureSet(merged_seurat, pattern = "^RP[SL]", col.name = "percent_ribo")

# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
merged_seurat <- PercentageFeatureSet(merged_seurat, "^HB[^(P)]", col.name = "percent_hb")

#Same for platelets
merged_seurat <- PercentageFeatureSet(merged_seurat, "PECAM1|PF4", col.name = "percent_plat")


## Create a metadata Object from seurat to work on

### Create metadata dataframe
metadata <- merged_seurat@meta.data

### Add cell IDs to metadata
metadata$cells <- rownames(metadata)

### Create Condition Column
metadata$Condition <- NA
metadata$Condition[which(str_detect(metadata$cells, "CE"))] <- "HEU"
metadata$Condition[which(str_detect(metadata$cells, "CP"))] <- "HEI"
metadata$Condition[which(str_detect(metadata$cells, "CS"))] <- "HUU"
metadata$Condition[which(str_detect(metadata$cells, "SA"))] <- "HEI"


### Create Timepoint Column
metadata$Timepoint <- NA
metadata$Timepoint[which(str_detect(metadata$cells, "_entry"))] <- "Entry"
metadata$Timepoint[which(str_detect(metadata$cells, "_12m"))] <- "12m"


### Create Cohort Column
metadata$Cohort <- NA
metadata$Cohort <- ifelse(str_detect(metadata$cells, "SA"), "EARTH", "TARA")

### Add Viral Load Column
metadata$Viral_Load <- NA
metadata$Viral_Load[which(str_detect(metadata$cells, "CE021_entry"))] <- "0"
metadata$Viral_Load[which(str_detect(metadata$cells, "CE025_entry"))] <- "0"
metadata$Viral_Load[which(str_detect(metadata$cells, "CE037_entry"))] <- "0"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP002_entry"))] <- "4284389"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP003_entry"))] <- "656769"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP006_12m"))] <- "73"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP006_entry"))] <- "10000000"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP011_entry"))] <- "36965"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP013_entry"))] <- "3434"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP016_entry"))] <- "1978332"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP017_entry"))] <- "3167384"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP018_entry"))] <- "176970"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP022_entry"))] <- "5075764"
metadata$Viral_Load[which(str_detect(metadata$cells, "CP042_entry"))] <- "6113"
metadata$Viral_Load[which(str_detect(metadata$cells, "CS005_entry"))] <- "0"
metadata$Viral_Load[which(str_detect(metadata$cells, "CS015_entry"))] <- "0"
metadata$Viral_Load[which(str_detect(metadata$cells, "SAAH29_entry"))] <- "24769"
metadata$Viral_Load[which(str_detect(metadata$cells, "SATY021_entry"))] <- "488"

# Classify into high and low viral load
# Define the threshold for high and low viral load
threshold <- 100000

# Create a new column 'Viral_Load_Category' based on the threshold
metadata$Viral_Load_Category <- ifelse(as.numeric(as.character(metadata$Viral_Load)) > threshold, "High", "Low")

# Split into CTL Expressing Groups

##### Adding CTL-metadata ####
# Create a mapping of orig.ident to the new grouping based on the image
group_mapping <- list(
  "CP006_entry" = "low",
  "CP006_12m" = "low",
  "CP011_entry" = "high",
  "CP013_entry" = "low",
  "CP016_entry" = "high",
  "CP017_entry" = "high",
  "CP022_entry" = "high",
  "CP042_entry" = "low",
  "SAAH29_entry" = "high",
  "SATY021_entry" = "low",
  "CP003_entry" = "high",
  "CP002_entry" = "high",
  "CP018_entry" = "low"
)

# Function to assign the group based on orig.ident
assign_group <- function(orig.ident) {
  if (orig.ident %in% names(group_mapping)) {
    return(group_mapping[[orig.ident]])
  } else {
    return("HEU/HUU")
  }
}

# Create the new column 'Grouping'
metadata$CTLGrouping <- sapply(metadata$orig.ident, assign_group)

### After verifying created metadata is correct, add metadata back to seurat object
merged_seurat@meta.data <- metadata


# Create .RData object to load at any time
save(merged_seurat, file=paste0(out.path,"Seuratv5_CITEseq_dsbnorm_merged_Seurat.RData"))

########################################################################################################################