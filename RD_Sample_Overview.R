# R Packages
# ==========
library(dplyr)
library(ggplot2)
library(reshape)
library(ggpubr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(stringr)
library(tidyr)
library(RColorBrewer)
library(ggfortify)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Change directory
setwd("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/04_Summary/450k/modified_InfiniumPurify/Enmix_quantile_with_QCinfo_Harman_All/03_Sample_Overview/")

# get the 450k annotation data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k <- as.data.frame(ann450k)
head(ann450k)

# Define CpGs of interest
CpG_islands <- ann450k %>% filter(Relation_to_Island == "Island") %>% data.frame() %>% .$Name
OpenSea <- ann450k %>% filter(Relation_to_Island == "OpenSea") %>% data.frame() %>% .$Name
S_Shelf <- ann450k %>% filter(Relation_to_Island == "S_Shelf") %>% data.frame() %>% .$Name
S_Shore <- ann450k %>% filter(Relation_to_Island == "S_Shore") %>% data.frame() %>% .$Name
N_Shelf <- ann450k %>% filter(Relation_to_Island == "N_Shelf") %>% data.frame() %>% .$Name
N_Shore <- ann450k %>% filter(Relation_to_Island == "N_Shore") %>% data.frame() %>% .$Name

# Read in methylation data
bval <- readRDS("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/01_InfiniumPurify/450k/modified_InfiniumPurify/Enmix_quantile_with_QCinfo_Harman_All/beta_combat_purified.rds")
mval <- readRDS("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/01_InfiniumPurify/450k/modified_InfiniumPurify/Enmix_quantile_with_QCinfo_Harman_All/mval_combat_purified.rds")

dim(bval)
dim(mval)

# Read in sample annotations
Sample_Annotations <- read.table(file = "/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/ANN/20231031_sample_annotations.tsv", sep = "\t", header = TRUE, na.strings = c("","NA")) %>%
  subset(TYPE!="TFL") %>%
  subset(TIME_POINT!="T2") %>%
  filter(SAMPLE_ID %in% colnames(bval))

# Read in sample sheet with batch information
targets <- read.csv(file = "/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/SampleSheet/SampleSheet.csv", skip = 1, header = TRUE, na.strings = c("","NA"))
targets$ID <- gsub(" ", "", paste(targets$Sample_Name,".",targets$Sentrix_Position))

# Remove repeated samples to obtain correct batch information
rep_rm_samples <- c("LY_FL_535_T1.R08C01", "LY_FL_498_T1.R06C01", "LY_FL_158_T1.R01C01",
                    "LY_FL_479_T1.R02C01", "LY_FL_488_T1.R06C01", "LY_FL_523_T1.R01C01",
                    "LY_FL_524_T1.R02C01", "LY_FL_525_T1.R03C01", "LY_FL_527_T1.R01C01",
                    "LY_FL_529_T1.R01C01", "LY_FL_159_T1.R08C01", "LY_FL_536_T1.R01C01",
                    "HCT116_DKO_methylated.R08C01", "HCT116_DKO_methylated.R04C01")

targets <- targets[!targets$ID %in% rep_rm_samples,]
targets$Sample_Name <- str_replace(targets$Sample_Name, "LY_FL_159_T1_rep", "LY_FL_159_T1")
targets <- targets %>% dplyr::select(Sample_Name, Batch)
colnames(targets)[1]  <- "SAMPLE_ID"
targets <- targets %>% filter(SAMPLE_ID %in% colnames(bval))

# Read in the Benign_LN sample annotation
Benign_LN_sample_ann <- read.table(file = "/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/Benign_LN/SampleSheet/Benign_LN_sample_annotations.tsv",
                                   sep = "\t",
                                   header = TRUE,
                                   na.strings=c("","NA")) %>%
  dplyr::select(SAMPLE_ID, TYPE)

Benign_LN_sample_ann <- Benign_LN_sample_ann %>% filter(SAMPLE_ID %in% colnames(bval))
Benign_LN_sample_ann$Batch <- "Benign_LN"

# Read in the NIH sample annotation
NIH_sample_ann <- read.table(file = "/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/NIH/GSE237299/SampleSheet/GSE237299_sample_annotations.txt",
                             sep = "\t",
                             header = TRUE,
                             na.strings=c("","NA")) %>%
  dplyr::select(SAMPLE_ID, TYPE)

NIH_sample_ann$SAMPLE_ID <- gsub('-', '_', NIH_sample_ann$SAMPLE_ID)
NIH_sample_ann <- NIH_sample_ann %>% filter(SAMPLE_ID %in% colnames(bval))
NIH_sample_ann$Batch <- "GSE237299"

# Read in the EGA sample annotation
EGA_sample_ann <- read.table(file = "/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/EGA/SampleSheet/Sample_Annotations_EGAD00010001974.txt",
                             sep = "\t",
                             header = TRUE,
                             na.strings=c("","NA")) %>%
  dplyr::select(SAMPLE_ID, TYPE)

EGA_sample_ann$SAMPLE_ID <- gsub('-', '_', EGA_sample_ann$SAMPLE_ID)
EGA_sample_ann$SAMPLE_ID <- gsub('\\.', '_', EGA_sample_ann$SAMPLE_ID)
EGA_sample_ann$Batch <- "EGAD00010001974"
EGA_sample_ann <- EGA_sample_ann %>% filter(SAMPLE_ID %in% colnames(bval))

# Read in blueprint sample annotation
Blueprint_sample_ann <- read.csv(file = "/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/blueprint/Sample_Annotation_EGAS00001001196.txt",
                                 sep = "\t",
                                 header = TRUE,
                                 na.strings=c("","NA")) %>%
#  dplyr::select(Sample_Name, Source_Name, CELL_TYPE, TISSUE_TYPE, DONOR_SEX)
  dplyr::select(Sample_Name, Source_Name, CELL_TYPE, TISSUE_TYPE)

#colnames(Blueprint_sample_ann) <- c("SAMPLE_ID", "LY_FL_ID", "TYPE", "SITE_BIOPSY", "SEX")
colnames(Blueprint_sample_ann) <- c("SAMPLE_ID", "LY_FL_ID", "TYPE", "SITE_BIOPSY")
Blueprint_sample_ann$Batch <- 'Blueprint'
Blueprint_sample_ann$SAMPLE_ID <- gsub('-', '_', Blueprint_sample_ann$SAMPLE_ID)
Blueprint_sample_ann$TYPE <- gsub('-', '_', Blueprint_sample_ann$TYPE)
# Blueprint_sample_ann$SEX <- gsub('Male', 'M', Blueprint_sample_ann$SEX)
# Blueprint_sample_ann$SEX <- gsub('Female', 'F', Blueprint_sample_ann$SEX)
Blueprint_sample_ann <- Blueprint_sample_ann %>% filter(SAMPLE_ID %in% colnames(bval))

# Read in the GEO sample annotation
GEO_sample_ann <- read.table(file = "/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/GEO/GEO_Annotation.tsv",
                             sep = "\t",
                             header = TRUE,
                             na.strings=c("","NA")) %>%
  dplyr::select(SAMPLE_ID, TYPE)
GEO_sample_ann$Batch <- 'GSE255869'
GEO_sample_ann <- GEO_sample_ann %>% filter(SAMPLE_ID %in% colnames(bval))

# Read in clinical data
clinical_data <- read.csv(file = "/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/CLN/20240116_clinical_data.csv", 
                          header = TRUE, na.strings = c("","NA")) %>%
  dplyr::select(SAMPLE_ID, AGE_AT_DIAGNOSIS, AGE_CAT, SEX, ANN_ARBOR_STAGE, LDH_ELEVATED, HEMOGLOBIN_LESS_120, MORE_4_NODAL_SITES, FLIPI_BINARY, GRADE, T_14_18, PRIM_TX_CAT, CODE_PFS, CODE_TRANSF, CODE_OS, PFS, TTT, OS)

clinical_data <- clinical_data[rowSums(is.na(clinical_data)) < ncol(clinical_data), ]

# Left outer join sample sheet and sample annotation
df1 <- merge(x = Sample_Annotations,
             y = targets,
             by.x = "SAMPLE_ID",
             by.y = "SAMPLE_ID",
             all.x = TRUE)

# df2 <- merge(x = df1,
#              y = clinical_data,
#              by.x = "LY_FL_ID",
#              by.y = "LY_FL_ID",
#              all.x = TRUE)

df2 <- merge(x = df1,
             y = clinical_data,
             by.x = "SAMPLE_ID",
             by.y = "SAMPLE_ID",
             all.x = TRUE)

Sample_Annotations <- df2

dim(Sample_Annotations)
remove(df1, df2)

# Combine sample annotations
Sample_Annotations <- bind_rows(Sample_Annotations, Benign_LN_sample_ann, NIH_sample_ann, EGA_sample_ann, GEO_sample_ann, Blueprint_sample_ann)
dim(Sample_Annotations)
table(Sample_Annotations$TYPE)
table(Sample_Annotations$Batch)

row.names(Sample_Annotations) <- Sample_Annotations$SAMPLE_ID

# subset of merged methylation data
bval <- bval[, Sample_Annotations$SAMPLE_ID]
mval <- mval[, Sample_Annotations$SAMPLE_ID]
dim(bval)
dim(mval)

gc()

# Remove Normal_GCB
Sample_Annotations <- Sample_Annotations %>% filter(TYPE != "Normal_GCB") %>% data.frame()

# Remove Blueprint samples before clustering
# Sample_Annotations <- Sample_Annotations %>% filter(TYPE == "FL" | TYPE == "FL_A" | TYPE == "FL_B" | TYPE == "DLBCL" | TYPE == "DLBCL_EGA" | TYPE == "DLBCL_EN" | TYPE == "DLBCL_GC" | TYPE == "DLBCL_nonGC" | TYPE == "t_DLBCL" | TYPE == "Benign_LN" | TYPE == "LN_NN_A" | TYPE == "LN_NN_B" | TYPE == "RLN") %>% data.frame()

# Calculate Samples by Type per Cluster
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'FL_A'] <- 'FL'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'FL_B'] <- 'FL'

Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'DLBCL_EGA'] <- 'DLBCL'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'DLBCL_EN'] <- 'DLBCL'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'DLBCL_GC'] <- 'DLBCL'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'DLBCL_nonGC'] <- 'DLBCL'
# Sample_Annotations$TYPE[Sample_Annotations$TYPE == 't_DLBCL'] <- 'DLBCL'

Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'Benign_LN'] <- 'Normal'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'LN_NN_A'] <- 'Normal'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'LN_NN_B'] <- 'Normal'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'RLN'] <- 'Normal'

# Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'Normal_GCB'] <- 'gcBC'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'bm_PC'] <- 'bm_PC'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'gcBC'] <- 'gcBC'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'memBC'] <- 'memBC'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'naiBC'] <- 'naiBC'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 't_PC'] <- 't_PC'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 't_naiBC'] <- 't_naiBC'

table(Sample_Annotations$TYPE)
table(Sample_Annotations$Batch)

# subset of merged methylation data
bval <- bval[, Sample_Annotations$SAMPLE_ID]
mval <- mval[, Sample_Annotations$SAMPLE_ID]
dim(bval)
dim(mval)

# Write TYPE and Batch Information
df1 <- Sample_Annotations[c("SAMPLE_ID", "TYPE","Batch")]
write.table(df1, file="TYPE_INFO.txt", sep="\t", quote = FALSE, row.names = FALSE)

gc()

# Define the clinical columns
clinical_cols <- c("AGE_AT_DIAGNOSIS", "AGE_CAT", "SEX", 
                   "ANN_ARBOR_STAGE", "LDH_ELEVATED", "HEMOGLOBIN_LESS_120", 
                   "MORE_4_NODAL_SITES", "FLIPI_BINARY", "GRADE", 
                   "T_14_18", "PRIM_TX_CAT", "CODE_PFS", 
                   "CODE_TRANSF", "CODE_OS", "PFS", 
                   "TTT", "OS")

# Create the new "Clinical" column
Sample_Annotations$Clinical <- ifelse(rowSums(is.na(Sample_Annotations[clinical_cols])) == length(clinical_cols), "Not Profiled", "Profiled")

# Subset Sample_Annotations
Sample_Annotations_Sub <- Sample_Annotations %>% select(SAMPLE_ID, INSTITUTION, TYPE, Batch, Clinical)
dim(Sample_Annotations_Sub)

table(Sample_Annotations_Sub$INSTITUTION)

#Sample_Annotations_Sub <- Sample_Annotations_Sub %>% mutate(INSTITUTION = coalesce(INSTITUTION, Batch))

# Add Repository Column
Sample_Annotations_Sub <- Sample_Annotations_Sub %>% mutate(Repository = ifelse(is.na(INSTITUTION), Batch, NA))

Sample_Annotations_Sub <- Sample_Annotations_Sub %>% 
  mutate(Repository = str_replace(Repository, regex("Benign_LN", ignore_case = TRUE), "GSE237299"))

# Load RNA-seq data
# =================
# FL RNA-Seq, n=282
# DLBCL RNA-Seq, n=8

# mRNA <- read.table("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/RNA_Seq/2022-10-25_df_counts_adj.txt", sep = " ", header = TRUE)
passing.RNAseq.QC <- read.table("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/RNA_Seq/metadata_passed_290.txt", sep = "\t", header = TRUE) %>% subset(phenotype == "FL" | phenotype == "DLBCL")

table(passing.RNAseq.QC$phenotype)
length(unique(passing.RNAseq.QC$sample_id))

# mRNA <- mRNA[,passing.RNAseq.QC$id]
# dim(mRNA)
# mRNA[1:5, 1:5]

# Rename column names
colnames(passing.RNAseq.QC)[which(names(passing.RNAseq.QC) == "sample_id")] <- "SAMPLE_ID"
colnames(passing.RNAseq.QC)[which(names(passing.RNAseq.QC) == "phenotype")] <- "TYPE"

passing.RNAseq.QC <- passing.RNAseq.QC %>% select(SAMPLE_ID, TYPE)
dim(passing.RNAseq.QC)

# Extract only the part of SAMPLE_ID that starts with "LY_"
passing.RNAseq.QC <- passing.RNAseq.QC %>%
  mutate(SAMPLE_ID = str_extract(SAMPLE_ID, "LY_.*"))

# Add RNA_Seq column based on matching SAMPLE_ID
Sample_Annotations_Sub <- Sample_Annotations_Sub %>%
  mutate(RNA_Seq = ifelse(SAMPLE_ID %in% passing.RNAseq.QC$SAMPLE_ID, "Profiled", "Not Profiled"))

table(Sample_Annotations_Sub$RNA_Seq)

# Load SNPs data
# ==============
ClusterAIC <- read.table("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/ClusterAIC/2023-05-06_14_GMM_Cluster_Labels_flexmix_clusters.txt", sep = "\t", header = TRUE) %>% select(SAMPLE_ID, ClusterAIC)

# Add SNPs column based on matching SAMPLE_ID
Sample_Annotations_Sub <- Sample_Annotations_Sub %>%
  mutate(SNPs = ifelse(SAMPLE_ID %in% ClusterAIC$SAMPLE_ID, "Profiled", "Not Profiled"))

table(Sample_Annotations_Sub$SNPs)

# Load epiTypes
# =============
FL_epitypes <- read.table("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/02_Clustering/450k/modified_InfiniumPurify/Enmix_quantile_with_QCinfo_Harman_All/FL_epitypes.tsv", sep = "\t", header = TRUE)

Sample_Annotations_Sub <- merge(x = Sample_Annotations_Sub,
                                y = FL_epitypes,
                                by.x = "SAMPLE_ID",
                                by.y = "SAMPLE_ID",
                                all.x = TRUE)

# Write TYPE and Batch Information
df2 <- Sample_Annotations_Sub[c("SAMPLE_ID", "TYPE", "epitypes")]
write.table(df2, file="Epitypes_Info.txt", sep="\t", quote = FALSE, row.names = FALSE)

gc()

# Replace NAs with Blueprint
Sample_Annotations_Sub$epitypes[is.na(Sample_Annotations_Sub$epitypes)] <- "Blueprint"
table(Sample_Annotations_Sub$epitypes)

# Add Methylation column based on matching SAMPLE_ID
Sample_Annotations_Sub <- Sample_Annotations_Sub %>%
  mutate(Methylation = ifelse(SAMPLE_ID %in% FL_epitypes$SAMPLE_ID, "Profiled", "Not Profiled"))

table(Sample_Annotations_Sub$Methylation)

Sample_Annotations_Sub$Methylation[Sample_Annotations_Sub$Methylation == "Not Profiled"] <- "Profiled"

table(Sample_Annotations_Sub$Methylation)

# Rename column names
colnames(Sample_Annotations_Sub)[which(names(Sample_Annotations_Sub) == "INSTITUTION")] <- "Institute"
colnames(Sample_Annotations_Sub)[which(names(Sample_Annotations_Sub) == "TYPE")] <- "Type"

# Reorder columns using select
Sample_Annotations_Sub <- Sample_Annotations_Sub %>%
  select(SAMPLE_ID, Type, epitypes, Batch, Institute, Repository, RNA_Seq, SNPs, Methylation, Clinical)

# Ensure column names are correct
colnames(Sample_Annotations_Sub)

# Convert `epitype` to a factor with a defined order
Sample_Annotations_Sub$epitypes <- factor(Sample_Annotations_Sub$epitypes, 
                                         levels = c("C1", "C2", "C3", "UC", "Blueprint"))

# Define Type sorting order based on epitype
Sample_Annotations_Sub$Type <- factor(Sample_Annotations_Sub$Type, levels = c("FL", "t_DLBCL", "DLBCL", "Normal", "naiBC", "t_naiBC", "gcBC", "t_PC", "memBC", "bm_PC", 
                                                                              setdiff(unique(Sample_Annotations_Sub$Type), c("FL", "t_DLBCL", "DLBCL", "Normal", "naiBC", "t_naiBC", "gcBC", "t_PC", "memBC", "bm_PC"))))

# Sort first by epitype, then by Type (custom order within Blueprint)
Sample_Annotations_Sub <- Sample_Annotations_Sub[order(Sample_Annotations_Sub$epitypes, 
                                                       Sample_Annotations_Sub$Type), ]

# Sort first by epitype, then by Type within C1 & C2
Sample_Annotations_Sub <- Sample_Annotations_Sub[order(Sample_Annotations_Sub$epitypes, 
                                                       Sample_Annotations_Sub$Type), ]

# Adjusted colors for Epitypes
epitypes_colors <- c(
  "C1" = "#4D9900",
  "C2" = "#B20000",
  "C3" = "#5A00B0",
  "UC" = "#595959",
  "Blueprint" = "#A9C8E4"
)

# Adjusted colors for Types
type_colors <- c(
  "Normal" = "#4A89F3",
  "FL" = "#70C76F",
  "t_DLBCL" = "#FF8E42",
  "DLBCL" = "#D95842",
  "naiBC" = "#9C7746",
  "t_naiBC" = "#C7A6A6",
  "gcBC" = "#E8C21D",
  "t_PC" = "#6FBFBA",
  "memBC" = "#F58DA1",
  "bm_PC" = "#5B3D89"
)

# Batch colors for better contrast
batch_colors <- c(
  "Benign_LN" = "#9B30FF",           
  "Blueprint" = "#2F4F4F",           
  "EGAD00010001974" = "#f5d007",     
  "Genome-Quebec" = "#F4A300",       
  "GSE237299" = "#9ACD32",           
  "GSE255869" = "#1E8B8D",           
  "PM-Genomics-Centre" = "#8B4513",  
  "PM-OICR-TGL" = "#2a69d4"          
)

# Institute colors for better contrast
institute_colors <- c(
  "AARHUS" = "#D9534F",     
  "BCCA" = "#000000",       
  "Blueprint" = "#2F4F4F",   
  "E2408" = "#FF6347",       
  "EGAD00010001974" = "#f5d007", 
  "GSE237299" = "#9ACD32",   
  "GSE255869" = "#1E8B8D",   
  "JGH" = "#800080",        
  "KINGSTON" = "#A9C8E4",    
  "OSLO" = "#800000",        
  "UHN" = "#808000"          
)

# Repository colors
repository_colors <- c(
  "Blueprint" = "#2F4F4F",           
  "EGAD00010001974" = "#f5d007",
  "GSE237299" = "#9ACD32",
  "GSE255869" = "#1E8B8D"
)

methylation_colors <- c("Profiled" = "#008080", "Not Profiled" = "#D0EFFF")
snps_colors <- c("Profiled" = "#800080", "Not Profiled" = "#F8DFFF")
# rna_seq_colors <- c("Profiled" = "#FF4500", "Not Profiled" = "#FFE5D0")
clinical_colors <- c("Profiled" = "#4C5C67", "Not Profiled" = "#D1D8DC")

# Create COLUMN annotations with Epitype and Type
col_anno <- HeatmapAnnotation(
  Epitype = Sample_Annotations_Sub$epitype,
  Type = Sample_Annotations_Sub$Type,
  Methylation = Sample_Annotations_Sub$Methylation,
  Targeted_DNA_Seq = Sample_Annotations_Sub$SNPs,
  # RNA_Seq = Sample_Annotations_Sub$RNA_Seq,
  Clinical_Info = Sample_Annotations_Sub$Clinical,
  Institute = Sample_Annotations_Sub$Institute,
  Repository = Sample_Annotations_Sub$Repository,
  Batch = Sample_Annotations_Sub$Batch,
  col = list(
    Epitype = epitypes_colors,
    Type = type_colors,
    Methylation = methylation_colors,
    Targeted_DNA_Seq = snps_colors,
    # RNA_Seq = rna_seq_colors,
    Clinical_Info = clinical_colors,
    Institute = institute_colors,
    Repository = repository_colors,
    Batch = batch_colors
  ),
  na_col = "white",  
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 12, fontface = "bold"),  # Make legend titles bold
    # labels_gp = gpar(fontsize = 10, fontface = "bold"), # Make legend labels bold
    title_gp = gpar(fontsize = 12),  
    labels_gp = gpar(fontsize = 10),
    ncol = 1,  
    legend_height = unit(10, "cm")  
  ),
  annotation_name_side = "right"
)

# Extract sorted sample IDs
sorted_sample_ids <- Sample_Annotations_Sub$SAMPLE_ID

# Create a dummy matrix with the sorted order
dummy_matrix <- matrix(0, nrow = 1, ncol = length(sorted_sample_ids))  
colnames(dummy_matrix) <- sorted_sample_ids

# Generate heatmap with annotations
ht <- Heatmap(
  dummy_matrix,
  name = "Annotations",
  col = c("0" = "white"),
  show_heatmap_legend = FALSE,
  top_annotation = col_anno,
  column_split = Sample_Annotations_Sub$epitype,
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),  # Make split titles bold
  show_column_names = FALSE,
  cluster_columns = FALSE,
  cluster_rows = FALSE
)

# Save to PDF
pdf("Sample_Overview.pdf", width = 12, height = 4.2)
draw(ht, annotation_legend_side = "bottom", heatmap_legend_side = "bottom", 
     padding = unit(c(5, 2, 2, 2), "mm"))
dev.off()
