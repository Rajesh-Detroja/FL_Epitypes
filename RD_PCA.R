# R Packages
# ==========
library(dplyr)
library(ggplot2)
library(reshape)
library(ggpubr)
library(factoextra)
library(readxl)
library(FactoMineR)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(stringr)
library(tidyr)
library(RColorBrewer)
library(cluster)
library(pheatmap)
library(GGally)
library(ggfortify)

# Change directory
setwd("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/05_PCA/450k/modified_InfiniumPurify/Enmix_quantile_with_QCinfo_Harman_All/01_After_Tumor_Purity/")

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
  dplyr::select(Sample_Name, Source_Name, CELL_TYPE, TISSUE_TYPE, DONOR_SEX)

colnames(Blueprint_sample_ann) <- c("SAMPLE_ID", "LY_FL_ID", "TYPE", "SITE_BIOPSY", "SEX")
Blueprint_sample_ann$Batch <- 'Blueprint'
Blueprint_sample_ann$SAMPLE_ID <- gsub('-', '_', Blueprint_sample_ann$SAMPLE_ID)
Blueprint_sample_ann$TYPE <- gsub('-', '_', Blueprint_sample_ann$TYPE)
Blueprint_sample_ann$SEX <- gsub('Male', 'M', Blueprint_sample_ann$SEX)
Blueprint_sample_ann$SEX <- gsub('Female', 'F', Blueprint_sample_ann$SEX)
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
  dplyr::select(LY_FL_ID, ANN_ARBOR_STAGE, SEX, T_14_18, PRIM_TX_CAT, CODE_PFS, CODE_TRANSF, CODE_OS, PFS, TTT, OS)

# Left outer join sample sheet and sample annotation
df1 <- merge(x = Sample_Annotations,
             y = targets,
             by.x = "SAMPLE_ID",
             by.y = "SAMPLE_ID",
             all.x = TRUE)

df2 <- merge(x = df1,
             y = clinical_data,
             by.x = "LY_FL_ID",
             by.y = "LY_FL_ID",
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

# Create list of groups
naiBC_Annotations <- Sample_Annotations %>% filter(TYPE == "naiBC") %>% data.frame()
gcBC_Annotations <- Sample_Annotations %>% filter(TYPE == "gcBC") %>% data.frame()
memBC_Annotations <- Sample_Annotations %>% filter(TYPE == "memBC") %>% data.frame()

# Subset methylation data group
naiBC_bval <- bval[, naiBC_Annotations$SAMPLE_ID]
gcBC_bval <- bval[, gcBC_Annotations$SAMPLE_ID]
memBC_bval <- bval[, memBC_Annotations$SAMPLE_ID]

naiBC_mval <- mval[, naiBC_Annotations$SAMPLE_ID]
gcBC_mval <- mval[, gcBC_Annotations$SAMPLE_ID]
memBC_mval <- mval[, memBC_Annotations$SAMPLE_ID]

dim(naiBC_bval)
dim(naiBC_mval)

dim(gcBC_bval)
dim(gcBC_mval)

dim(memBC_bval)
dim(memBC_mval)

# Calculate Mean Methylation for Each CpG in Naive, GC, and Memory B Cells
naiBC_means <- rowMeans(naiBC_bval)
gcBC_means <- rowMeans(gcBC_bval)
memBC_means <- rowMeans(memBC_bval)

length(naiBC_means)
length(gcBC_means)
length(memBC_means)

# Identify CpG Probes with Similar Methylation Profiles
threshold <- 0.1

# Identify CpGs with similar methylation profiles to naive B cells
gcBC_naiBC <- abs(gcBC_means - naiBC_means) <= threshold
table(gcBC_naiBC)

memBC_naiBC <- abs(memBC_means - naiBC_means) <= threshold
table(memBC_naiBC)

# Combine the results to find CpGs similar to naive in both GC and memory B cells
naiBC_gcBC_memBC <- gcBC_naiBC & memBC_naiBC
table(naiBC_gcBC_memBC)

# CpG probes to keep (those that are NOT similar across all B cell types)
keep <- !naiBC_gcBC_memBC
remove <- naiBC_gcBC_memBC
table(keep)
table(remove)

# Filter Out Similar CpG Probes
bval_keep <- bval[keep, ]
mval_keep <- mval[keep, ]

bval_remove <- bval[remove, ]
mval_remove <- mval[remove, ]

dim(bval_keep)
dim(mval_keep)

dim(bval_remove)
dim(mval_remove)

gc()

# Remove Blueprint samples before clustering
Sample_Annotations <- Sample_Annotations %>% filter(TYPE == "FL" | TYPE == "FL_A" | TYPE == "FL_B" | TYPE == "DLBCL" | TYPE == "DLBCL_EGA" | TYPE == "DLBCL_EN" | TYPE == "DLBCL_GC" | TYPE == "DLBCL_nonGC" | TYPE == "t_DLBCL" | TYPE == "Benign_LN" | TYPE == "LN_NN_A" | TYPE == "LN_NN_B" | TYPE == "RLN") %>% data.frame()

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

Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'Normal_GCB'] <- 'gcBC'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'bm_PC'] <- 'bm_PC'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'gcBC'] <- 'gcBC'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'memBC'] <- 'memBC'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'naiBC'] <- 'naiBC'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 't_PC'] <- 't_PC'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 't_naiBC'] <- 't_naiBC'

table(Sample_Annotations$TYPE)
table(Sample_Annotations$Batch)

# Convert type to a factor with a defined order
Sample_Annotations$TYPE <- factor(Sample_Annotations$TYPE, levels = c("Normal", "FL", "t_DLBCL", "DLBCL"))

table(Sample_Annotations$TYPE)
table(Sample_Annotations$Batch)

# subset of merged methylation data
bval <- bval[, Sample_Annotations$SAMPLE_ID]
mval <- mval[, Sample_Annotations$SAMPLE_ID]
dim(bval)
dim(mval)

bval_keep <- bval_keep[, Sample_Annotations$SAMPLE_ID]
mval_keep <- mval_keep[, Sample_Annotations$SAMPLE_ID]
dim(bval_keep)
dim(mval_keep)

bval_remove <- bval_remove[, Sample_Annotations$SAMPLE_ID]
mval_remove <- mval_remove[, Sample_Annotations$SAMPLE_ID]
dim(bval_remove)
dim(mval_remove)

# Write TYPE and Batch Information
df1 <- Sample_Annotations[c("SAMPLE_ID", "TYPE","Batch")]
write.table(df1, file="TYPE_INFO.txt", sep="\t", quote = FALSE, row.names = FALSE)

gc()

# Probes with Low Variance
# Very strict: Cut-off near 0.05 SD—useful for clustering when you want to focus on highly variable probes.
# Moderate: Cut-off around 0.02 SD—removes probes with extremely low variance, while retaining most informative probes.
# Loose: Cut-off near 0.01 SD—retains almost all probes but removes the most non-variable ones.

# Define variability threshold
variability_threshold <- 0.05

# Compute standard deviation for each CpG in total, keep, and remove datasets
probe_var_total <- apply(bval, 1, sd, na.rm = TRUE)
probe_var_keep <- apply(bval_keep, 1, sd, na.rm = TRUE)
probe_var_remove <- apply(bval_remove, 1, sd, na.rm = TRUE)

# Get the total count of CpGs in each category and format with commas
num_total <- formatC(length(probe_var_total), format = "d", big.mark = ",")
num_keep <- formatC(length(probe_var_keep), format = "d", big.mark = ",")
num_remove <- formatC(length(probe_var_remove), format = "d", big.mark = ",")

# Plot histograms with formatted CpG count in the title
pdf("hist_total.pdf")
hist(probe_var_total, breaks = 100, main = paste("Distribution of CpG Probe SDs -", num_total, "CpGs"), xlab = "Standard Deviation")
abline(v = variability_threshold, col = "red")  
dev.off()

pdf("hist_keep.pdf")
hist(probe_var_keep, breaks = 100, main = paste("Distribution of CpG Probe SDs -", num_keep, "CpGs"), xlab = "Standard Deviation")
abline(v = variability_threshold, col = "red")  
dev.off()

pdf("hist_remove.pdf")
hist(probe_var_remove, breaks = 100, main = paste("Distribution of CpG Probe SDs -", num_remove, "CpGs"), xlab = "Standard Deviation")
abline(v = variability_threshold, col = "red")  
dev.off()

# Subset Methylation Data
bval <- bval[keep, ]
mval <- mval[keep, ]

dim(bval)
dim(mval)

remove(bval_keep, mval_keep, bval_remove, mval_remove)
gc()

# subset matrix based on most variable probes
set.seed(1234)
#probe_num = "10000"
probe_num = nrow(bval)

# subset random probes
bval_sub <- bval[sample(nrow(bval), probe_num), ]
mval_sub <- mval[intersect(rownames(bval_sub), rownames(mval)), ]

dim(bval_sub)
bval_sub[1:5, 1:5]

dim(mval_sub)
mval_sub[1:5, 1:5]

gc()

# Preparing input data
bval_sub <- data.matrix(bval_sub)
class(bval_sub)
dim(bval_sub)
bval_sub[1:5, 1:5]

mval_sub <- data.matrix(mval_sub)
class(mval_sub)
dim(mval_sub)
mval_sub[1:5, 1:5]

# # Remove unclassifiable samples before clustering analysis # n = 39 (FL: 33, DLBCL: 6)
# unclassified_samples <- c("2411_04_01BD",
#                           "DL13_D1108",
#                           "DL9_D1096",
#                           "GSM8081924",
#                           "GSM8081960",
#                           "HP_Z001",
#                           "LY_FL_001_T1",
#                           "LY_FL_013_T1",
#                           "LY_FL_020_T1",
#                           "LY_FL_062_T1",
#                           "LY_FL_095_T1",
#                           "LY_FL_1092_T1",
#                           "LY_FL_1093_T1",
#                           "LY_FL_1101_T1",
#                           "LY_FL_1112_T1",
#                           "LY_FL_1152_T1",
#                           "LY_FL_1160_T1",
#                           "LY_FL_1161_T1",
#                           "LY_FL_1174_T1",
#                           "LY_FL_127_T1",
#                           "LY_FL_136_T1",
#                           "LY_FL_143_T1",
#                           "LY_FL_152_T1",
#                           "LY_FL_158_T1",
#                           "LY_FL_172_T1",
#                           "LY_FL_210_T1",
#                           "LY_FL_237_T1",
#                           "LY_FL_249_T1",
#                           "LY_FL_260_T1",
#                           "LY_FL_264_T1",
#                           "LY_FL_295_T1",
#                           "LY_FL_307_T1",
#                           "LY_FL_464_T1",
#                           "LY_FL_488_T1",
#                           "LY_FL_524_T1",
#                           "LY_FL_609_T1",
#                           "GSM8081947",
#                           "HP_W678",
#                           "LY_FL_602_T1")
# 
# length(unclassified_samples)
# 
# # Remove unclassified samples from clustering input
# bval_sub <- bval_sub[, !colnames(bval_sub) %in% unclassified_samples]
# bval_sub[1:5, 1:5]
# dim(bval_sub)
# 
# mval_sub <- mval_sub[, !colnames(mval_sub) %in% unclassified_samples]
# mval_sub[1:5, 1:5]
# dim(mval_sub)
# 
# gc()

# bval - All CpG
bval_sub_t <- t(bval_sub)
bval_sub_t <- merge(bval_sub_t, Sample_Annotations, by=0)
names(bval_sub_t)[1] <- ""
bval_sub_t <- data.frame(bval_sub_t[,-1], row.names=bval_sub_t[,1])
bval_sub_t_df <- bval_sub_t[1:probe_num]

# mval - All CpG
mval_sub_t <- t(mval_sub)
mval_sub_t <- merge(mval_sub_t, Sample_Annotations, by=0)
names(mval_sub_t)[1] <- ""
mval_sub_t <- data.frame(mval_sub_t[,-1], row.names=mval_sub_t[,1])
mval_sub_t_df <- mval_sub_t[1:probe_num]

# bval - CpG Islands
bval_sub_CpG_islands <- bval_sub[intersect(CpG_islands, row.names(bval_sub)), ]
bval_sub_CpG_islands_t <- t(bval_sub_CpG_islands)
bval_sub_CpG_islands_t <- merge(bval_sub_CpG_islands_t, Sample_Annotations, by=0)
names(bval_sub_CpG_islands_t)[1] <- ""
bval_sub_CpG_islands_t <- data.frame(bval_sub_CpG_islands_t[,-1], row.names=bval_sub_CpG_islands_t[,1])
bval_sub_CpG_islands_t_df <- bval_sub_CpG_islands_t[1:ncol(t(bval_sub_CpG_islands))]

# mval - CpG Islands
mval_sub_CpG_islands <- mval_sub[intersect(CpG_islands, row.names(mval_sub)), ]
mval_sub_CpG_islands_t <- t(mval_sub_CpG_islands)
mval_sub_CpG_islands_t <- merge(mval_sub_CpG_islands_t, Sample_Annotations, by=0)
names(mval_sub_CpG_islands_t)[1] <- ""
mval_sub_CpG_islands_t <- data.frame(mval_sub_CpG_islands_t[,-1], row.names=mval_sub_CpG_islands_t[,1])
mval_sub_CpG_islands_t_df <- mval_sub_CpG_islands_t[1:ncol(t(mval_sub_CpG_islands))]

# bval - OpenSea
bval_sub_OpenSea <- bval_sub[intersect(OpenSea, row.names(bval_sub)), ]
bval_sub_OpenSea_t <- t(bval_sub_OpenSea)
bval_sub_OpenSea_t <- merge(bval_sub_OpenSea_t, Sample_Annotations, by=0)
names(bval_sub_OpenSea_t)[1] <- ""
bval_sub_OpenSea_t <- data.frame(bval_sub_OpenSea_t[,-1], row.names=bval_sub_OpenSea_t[,1])
bval_sub_OpenSea_t_df <- bval_sub_OpenSea_t[1:ncol(t(bval_sub_OpenSea))]

# mval - OpenSea
mval_sub_OpenSea <- mval_sub[intersect(OpenSea, row.names(mval_sub)), ]
mval_sub_OpenSea_t <- t(mval_sub_OpenSea)
mval_sub_OpenSea_t <- merge(mval_sub_OpenSea_t, Sample_Annotations, by=0)
names(mval_sub_OpenSea_t)[1] <- ""
mval_sub_OpenSea_t <- data.frame(mval_sub_OpenSea_t[,-1], row.names=mval_sub_OpenSea_t[,1])
mval_sub_OpenSea_t_df <- mval_sub_OpenSea_t[1:ncol(t(mval_sub_OpenSea))]

# bval - S_Shelf
bval_sub_S_Shelf <- bval_sub[intersect(S_Shelf, row.names(bval_sub)), ]
bval_sub_S_Shelf_t <- t(bval_sub_S_Shelf)
bval_sub_S_Shelf_t <- merge(bval_sub_S_Shelf_t, Sample_Annotations, by=0)
names(bval_sub_S_Shelf_t)[1] <- ""
bval_sub_S_Shelf_t <- data.frame(bval_sub_S_Shelf_t[,-1], row.names=bval_sub_S_Shelf_t[,1])
bval_sub_S_Shelf_t_df <- bval_sub_S_Shelf_t[1:ncol(t(bval_sub_S_Shelf))]

# mval - S_Shelf
mval_sub_S_Shelf <- mval_sub[intersect(S_Shelf, row.names(mval_sub)), ]
mval_sub_S_Shelf_t <- t(mval_sub_S_Shelf)
mval_sub_S_Shelf_t <- merge(mval_sub_S_Shelf_t, Sample_Annotations, by=0)
names(mval_sub_S_Shelf_t)[1] <- ""
mval_sub_S_Shelf_t <- data.frame(mval_sub_S_Shelf_t[,-1], row.names=mval_sub_S_Shelf_t[,1])
mval_sub_S_Shelf_t_df <- mval_sub_S_Shelf_t[1:ncol(t(mval_sub_S_Shelf))]

# bval - S_Shore
bval_sub_S_Shore <- bval_sub[intersect(S_Shore, row.names(bval_sub)), ]
bval_sub_S_Shore_t <- t(bval_sub_S_Shore)
bval_sub_S_Shore_t <- merge(bval_sub_S_Shore_t, Sample_Annotations, by=0)
names(bval_sub_S_Shore_t)[1] <- ""
bval_sub_S_Shore_t <- data.frame(bval_sub_S_Shore_t[,-1], row.names=bval_sub_S_Shore_t[,1])
bval_sub_S_Shore_t_df <- bval_sub_S_Shore_t[1:ncol(t(bval_sub_S_Shore))]

# mval - S_Shore
mval_sub_S_Shore <- mval_sub[intersect(S_Shore, row.names(mval_sub)), ]
mval_sub_S_Shore_t <- t(mval_sub_S_Shore)
mval_sub_S_Shore_t <- merge(mval_sub_S_Shore_t, Sample_Annotations, by=0)
names(mval_sub_S_Shore_t)[1] <- ""
mval_sub_S_Shore_t <- data.frame(mval_sub_S_Shore_t[,-1], row.names=mval_sub_S_Shore_t[,1])
mval_sub_S_Shore_t_df <- mval_sub_S_Shore_t[1:ncol(t(mval_sub_S_Shore))]

# bval - N_Shelf
bval_sub_N_Shelf <- bval_sub[intersect(N_Shelf, row.names(bval_sub)), ]
bval_sub_N_Shelf_t <- t(bval_sub_N_Shelf)
bval_sub_N_Shelf_t <- merge(bval_sub_N_Shelf_t, Sample_Annotations, by=0)
names(bval_sub_N_Shelf_t)[1] <- ""
bval_sub_N_Shelf_t <- data.frame(bval_sub_N_Shelf_t[,-1], row.names=bval_sub_N_Shelf_t[,1])
bval_sub_N_Shelf_t_df <- bval_sub_N_Shelf_t[1:ncol(t(bval_sub_N_Shelf))]

# mval - N_Shelf
mval_sub_N_Shelf <- mval_sub[intersect(N_Shelf, row.names(mval_sub)), ]
mval_sub_N_Shelf_t <- t(mval_sub_N_Shelf)
mval_sub_N_Shelf_t <- merge(mval_sub_N_Shelf_t, Sample_Annotations, by=0)
names(mval_sub_N_Shelf_t)[1] <- ""
mval_sub_N_Shelf_t <- data.frame(mval_sub_N_Shelf_t[,-1], row.names=mval_sub_N_Shelf_t[,1])
mval_sub_N_Shelf_t_df <- mval_sub_N_Shelf_t[1:ncol(t(mval_sub_N_Shelf))]

# bval - N_Shore
bval_sub_N_Shore <- bval_sub[intersect(N_Shore, row.names(bval_sub)), ]
bval_sub_N_Shore_t <- t(bval_sub_N_Shore)
bval_sub_N_Shore_t <- merge(bval_sub_N_Shore_t, Sample_Annotations, by=0)
names(bval_sub_N_Shore_t)[1] <- ""
bval_sub_N_Shore_t <- data.frame(bval_sub_N_Shore_t[,-1], row.names=bval_sub_N_Shore_t[,1])
bval_sub_N_Shore_t_df <- bval_sub_N_Shore_t[1:ncol(t(bval_sub_N_Shore))]

# mval - N_Shore
mval_sub_N_Shore <- mval_sub[intersect(N_Shore, row.names(mval_sub)), ]
mval_sub_N_Shore_t <- t(mval_sub_N_Shore)
mval_sub_N_Shore_t <- merge(mval_sub_N_Shore_t, Sample_Annotations, by=0)
names(mval_sub_N_Shore_t)[1] <- ""
mval_sub_N_Shore_t <- data.frame(mval_sub_N_Shore_t[,-1], row.names=mval_sub_N_Shore_t[,1])
mval_sub_N_Shore_t_df <- mval_sub_N_Shore_t[1:ncol(t(mval_sub_N_Shore))]

## TYPE
## ====
#type_colors <- c(Normal="#7618dc", FL="#6dcc04", DLBCL="#fc1417", t_DLBCL= "#6E0A0A", na.value="#c9c9c9")
#type_colors <- c(Benign_LN="#5086F3", bm_PC="#7f8992", DLBCL="#D02060", DLBCL_EGA="#fc1417", DLBCL_EN="#f26ded", DLBCL_GC="#D87087", DLBCL_nonGC="#9A3B8A", t_DLBCL= "#6E0A0A", Normal_GCB="#2E2EAA", FL="#6dcc04", FL_A="#02bdaa", FL_B="#ebc602", gcBC="#0303a8", LN_NN_A="#FF7F00", LN_NN_B="#B15928", memBC="#b5b80b", naiBC="#2480b9", RLN="#7618dc", t_naiBC="#6e988c", t_PC="#383838", na.value="#c9c9c9")
#type_colors <- c(Normal="#7618dc", FL="#6dcc04", DLBCL="#fc1417", bm_PC="#7f8992", DLBCL_GC="#f26ded", DLBCL_nonGC="#D87087", t_DLBCL= "#6E0A0A", Normal_GCB="#2E2EAA", gcBC="#0303a8", memBC="#b5b80b", naiBC="#2480b9", t_naiBC="#6e988c", t_PC="#383838", na.value="#c9c9c9")
#type_colors <- c(Normal="#7618dc", FL="#6dcc04", DLBCL="#fc1417", bm_PC="#7f8992", t_DLBCL= "#6E0A0A", Normal_GCB="#2E2EAA", gcBC="#0303a8", memBC="#b5b80b", naiBC="#2480b9", t_naiBC="#6e988c", t_PC="#383838", na.value="#c9c9c9")
#type_colors <- c(Benign_LN="#5086F3", DLBCL="#A01641", DLBCL_EGA="#fc1417", DLBCL_EN="#f26ded", FL="#6dcc04", FL_A="#02bdaa", FL_B="#ebc602", LN_NN_A="#FF7F00", LN_NN_B="#B15928", RLN="#7618dc", na.value="#c9c9c9")
#type_colors <- c(Normal="#7618dc", FL="#6dcc04", DLBCL="#fc1417", bm_PC="#7f8992", gcBC="#0303a8", memBC="#b5b80b", naiBC="#2480b9",  t_naiBC="#6e988c", t_PC="#383838", na.value="#c9c9c9")

# Define custom colors for improved visibility
type_colors <- c(
  "Normal" = "#5086F3",
  "FL" = "#59A14F",
  "t_DLBCL" = "#595959", # "#F28E2B"
  "DLBCL" = "#E15759"
  # "naiBC" = "#8C6D31",
  # "t_naiBC" = "#D1A1A6",
  # "gcBC" = "#E4B321",
  # "t_PC" = "#76B7B2",
  # "memBC" = "#FF9DA7",
  # "bm_PC" = "#6A4D8C"
)

# Beta-Value Clustering - TYPE - All CpGs
t <- bval_sub_t %>% drop_na(TYPE)
df <- bval_sub_t_df[intersect(row.names(t), row.names(bval_sub_t_df)), ]
pca_res <- prcomp(df, scale. = TRUE)

pdf("bval-pca-plot-type-All-CpGs.pdf", width = 6.2, height = 6)
autoplot(pca_res,
         data = t,
         colour = 'TYPE',
         frame = TRUE
) +
  theme_bw() +
  theme(title = element_text(size=16, face='bold'),
        panel.border = element_blank(),
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        axis.line = element_line(),
        axis.line.x.top = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size=14),
        legend.position="bottom",
        legend.title = element_text(size=14)) +
  ggtitle(paste("Beta value -", prettyNum(probe_num, big.mark=",", scientific=FALSE), "CpGs")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = type_colors) +
  scale_fill_manual(values = type_colors) +
  geom_point(aes(colour = factor(TYPE)), size = 1.5) +
  coord_cartesian(clip = 'off') +
  theme(plot.margin = margin(0.4, 0.8, 0.4, 0.4, "cm")) +
  labs(fill='Type') +
  labs(color='Type') +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
dev.off()


# M-Value Clustering - Type - All CpGs
t <- mval_sub_t %>% drop_na(TYPE)
df <- mval_sub_t_df[intersect(row.names(t), row.names(mval_sub_t_df)), ]
pca_res <- prcomp(df, scale. = TRUE)

pdf("mval-pca-plot-type-All-CpGs.pdf", width = 6.2, height = 6)
autoplot(pca_res,
         data = t,
         colour = 'TYPE',
         frame = TRUE
) +
  theme_bw() +
  theme(title = element_text(size=16, face='bold'),
        panel.border = element_blank(),
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        axis.line = element_line(),
        axis.line.x.top = element_blank(),
        axis.line.y.right = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text = element_text(size=14),
        legend.position="bottom",
        legend.title = element_text(size=14)) + 
  ggtitle(paste("M value -", prettyNum(probe_num, big.mark=",", scientific=FALSE), "CpGs")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = type_colors) +
  scale_fill_manual(values = type_colors) +
  geom_point(aes(colour = factor(TYPE)), size = 1.5) +
  coord_cartesian(clip = 'off') +
  theme(plot.margin = margin(0.4, 0.8, 0.4, 0.4, "cm")) +
  labs(fill='Type') +
  labs(color='Type') +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
dev.off()


# # Beta-Value Clustering - Type - CpG Islands
# t <- bval_sub_CpG_islands_t %>% drop_na(TYPE)
# df <- bval_sub_CpG_islands_t_df[intersect(row.names(t), row.names(bval_sub_CpG_islands_t_df)), ]
# pca_res <- prcomp(df, scale. = TRUE)
# 
# pdf("bval-pca-plot-type-CpG-Islands.pdf")
# autoplot(pca_res,
#          data = t,
#          colour = 'TYPE',
#          frame = TRUE
# ) +
#   theme_bw() +
#   theme(title = element_text(size=16, face='bold'),
#         axis.title.x = element_text(size = 16),
#         axis.text.x = element_text(size = 16),
#         axis.title.y = element_text(size = 16),
#         axis.text.y = element_text(size = 16),
#         legend.text = element_text(size=14),
#         legend.position="bottom",
#         legend.title = element_text(size=14)) + 
#   ggtitle(paste("Beta value - CpG Islands", prettyNum(nrow(bval_sub_CpG_islands), big.mark=",", scientific=FALSE), "/", prettyNum(probe_num, big.mark=",", scientific=FALSE), "CpGs")) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_color_manual(values = type_colors) +
#   scale_fill_manual(values = type_colors) +
#   geom_point(aes(colour = factor(TYPE)), size = 1.5) +
#   coord_cartesian(clip = 'off') +
#   theme(plot.margin = margin(0.4, 0.8, 0.4, 0.4, "cm")) +
#   labs(fill='Type') +
#   labs(color='Type') +
#   guides(color = guide_legend(nrow = 2, byrow = TRUE))
# dev.off()
# 
# 
# # M-Value Clustering - Type - CpG Islands
# t <- mval_sub_CpG_islands_t %>% drop_na(TYPE)
# df <- mval_sub_CpG_islands_t_df[intersect(row.names(t), row.names(mval_sub_CpG_islands_t_df)), ]
# pca_res <- prcomp(df, scale. = TRUE)
# 
# pdf("mval-pca-plot-type-CpG-Islands.pdf")
# autoplot(pca_res,
#          data = t,
#          colour = 'TYPE',
#          frame = TRUE
# ) +
#   theme_bw() +
#   theme(title = element_text(size=16, face='bold'),
#         axis.title.x = element_text(size = 16),
#         axis.text.x = element_text(size = 16),
#         axis.title.y = element_text(size = 16),
#         axis.text.y = element_text(size = 16),
#         legend.text = element_text(size=14),
#         legend.position="bottom",
#         legend.title = element_text(size=14)) + 
#   ggtitle(paste("M value - CpG Islands", prettyNum(nrow(mval_sub_CpG_islands), big.mark=",", scientific=FALSE), "/", prettyNum(probe_num, big.mark=",", scientific=FALSE), "CpGs")) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_color_manual(values = type_colors) +
#   scale_fill_manual(values = type_colors) +
#   geom_point(aes(colour = factor(TYPE)), size = 1.5) +
#   coord_cartesian(clip = 'off') +
#   theme(plot.margin = margin(0.4, 0.8, 0.4, 0.4, "cm")) +
#   labs(fill='Type') +
#   labs(color='Type') +
#   guides(color = guide_legend(nrow = 2, byrow = TRUE))
# dev.off()
# 
# 
# # Beta-Value Clustering - Type - OpenSea
# t <- bval_sub_OpenSea_t %>% drop_na(TYPE)
# df <- bval_sub_OpenSea_t_df[intersect(row.names(t), row.names(bval_sub_OpenSea_t_df)), ]
# pca_res <- prcomp(df, scale. = TRUE)
# 
# pdf("bval-pca-plot-type-OpenSea.pdf")
# autoplot(pca_res,
#          data = t,
#          colour = 'TYPE',
#          frame = TRUE
# ) +
#   theme_bw() +
#   theme(title = element_text(size=16, face='bold'),
#         axis.title.x = element_text(size = 16),
#         axis.text.x = element_text(size = 16),
#         axis.title.y = element_text(size = 16),
#         axis.text.y = element_text(size = 16),
#         legend.text = element_text(size=14),
#         legend.position="bottom",
#         legend.title = element_text(size=14)) + 
#   ggtitle(paste("Beta value - OpenSea", prettyNum(nrow(bval_sub_OpenSea), big.mark=",", scientific=FALSE), "/", prettyNum(probe_num, big.mark=",", scientific=FALSE), "CpGs")) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_color_manual(values = type_colors) +
#   scale_fill_manual(values = type_colors) +
#   geom_point(aes(colour = factor(TYPE)), size = 1.5) +
#   coord_cartesian(clip = 'off') +
#   theme(plot.margin = margin(0.4, 0.8, 0.4, 0.4, "cm")) +
#   labs(fill='Type') +
#   labs(color='Type') +
#   guides(color = guide_legend(nrow = 2, byrow = TRUE))
# dev.off()
# 
# 
# # M-Value Clustering - Type - OpenSea
# t <- mval_sub_OpenSea_t %>% drop_na(TYPE)
# df <- mval_sub_OpenSea_t_df[intersect(row.names(t), row.names(mval_sub_OpenSea_t_df)), ]
# pca_res <- prcomp(df, scale. = TRUE)
# 
# pdf("mval-pca-plot-type-OpenSea.pdf")
# autoplot(pca_res,
#          data = t,
#          colour = 'TYPE',
#          frame = TRUE
# ) +
#   theme_bw() +
#   theme(title = element_text(size=16, face='bold'),
#         axis.title.x = element_text(size = 16),
#         axis.text.x = element_text(size = 16),
#         axis.title.y = element_text(size = 16),
#         axis.text.y = element_text(size = 16),
#         legend.text = element_text(size=14),
#         legend.position="bottom",
#         legend.title = element_text(size=14)) + 
#   ggtitle(paste("M value - OpenSea", prettyNum(nrow(mval_sub_OpenSea), big.mark=",", scientific=FALSE), "/", prettyNum(probe_num, big.mark=",", scientific=FALSE), "CpGs")) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_color_manual(values = type_colors) +
#   scale_fill_manual(values = type_colors) +
#   geom_point(aes(colour = factor(TYPE)), size = 1.5) +
#   coord_cartesian(clip = 'off') +
#   theme(plot.margin = margin(0.4, 0.8, 0.4, 0.4, "cm")) +
#   labs(fill='Type') +
#   labs(color='Type') +
#   guides(color = guide_legend(nrow = 2, byrow = TRUE))
# dev.off()
# 
# 
# # Beta-Value Clustering - Type - S_Shelf
# t <- bval_sub_S_Shelf_t %>% drop_na(TYPE)
# df <- bval_sub_S_Shelf_t_df[intersect(row.names(t), row.names(bval_sub_S_Shelf_t_df)), ]
# pca_res <- prcomp(df, scale. = TRUE)
# 
# pdf("bval-pca-plot-type-S_Shelf.pdf")
# autoplot(pca_res,
#          data = t,
#          colour = 'TYPE',
#          frame = TRUE
# ) +
#   theme_bw() +
#   theme(title = element_text(size=16, face='bold'),
#         axis.title.x = element_text(size = 16),
#         axis.text.x = element_text(size = 16),
#         axis.title.y = element_text(size = 16),
#         axis.text.y = element_text(size = 16),
#         legend.text = element_text(size=14),
#         legend.position="bottom",
#         legend.title = element_text(size=14)) + 
#   ggtitle(paste("Beta value - S_Shelf", prettyNum(nrow(bval_sub_S_Shelf), big.mark=",", scientific=FALSE), "/", prettyNum(probe_num, big.mark=",", scientific=FALSE), "CpGs")) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_color_manual(values = type_colors) +
#   scale_fill_manual(values = type_colors) +
#   geom_point(aes(colour = factor(TYPE)), size = 1.5) +
#   coord_cartesian(clip = 'off') +
#   theme(plot.margin = margin(0.4, 0.8, 0.4, 0.4, "cm")) +
#   labs(fill='Type') +
#   labs(color='Type') +
#   guides(color = guide_legend(nrow = 2, byrow = TRUE))
# dev.off()
# 
# # M-Value Clustering - Type - S_Shelf
# t <- mval_sub_S_Shelf_t %>% drop_na(TYPE)
# df <- mval_sub_S_Shelf_t_df[intersect(row.names(t), row.names(mval_sub_S_Shelf_t_df)), ]
# pca_res <- prcomp(df, scale. = TRUE)
# 
# pdf("mval-pca-plot-type-S_Shelf.pdf")
# autoplot(pca_res,
#          data = t,
#          colour = 'TYPE',
#          frame = TRUE
# ) +
#   theme_bw() +
#   theme(title = element_text(size=16, face='bold'),
#         axis.title.x = element_text(size = 16),
#         axis.text.x = element_text(size = 16),
#         axis.title.y = element_text(size = 16),
#         axis.text.y = element_text(size = 16),
#         legend.text = element_text(size=14),
#         legend.position="bottom",
#         legend.title = element_text(size=14)) + 
#   ggtitle(paste("M value - S_Shelf", prettyNum(nrow(mval_sub_S_Shelf), big.mark=",", scientific=FALSE), "/", prettyNum(probe_num, big.mark=",", scientific=FALSE), "CpGs")) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_color_manual(values = type_colors) +
#   scale_fill_manual(values = type_colors) +
#   geom_point(aes(colour = factor(TYPE)), size = 1.5) +
#   coord_cartesian(clip = 'off') +
#   theme(plot.margin = margin(0.4, 0.8, 0.4, 0.4, "cm")) +
#   labs(fill='Type') +
#   labs(color='Type') +
#   guides(color = guide_legend(nrow = 2, byrow = TRUE))
# dev.off()
# 
# 
# # Beta-Value Clustering - Type - S_Shore
# t <- bval_sub_S_Shore_t %>% drop_na(TYPE)
# df <- bval_sub_S_Shore_t_df[intersect(row.names(t), row.names(bval_sub_S_Shore_t_df)), ]
# pca_res <- prcomp(df, scale. = TRUE)
# 
# pdf("bval-pca-plot-type-S_Shore.pdf")
# autoplot(pca_res,
#          data = t,
#          colour = 'TYPE',
#          frame = TRUE
# ) +
#   theme_bw() +
#   theme(title = element_text(size=16, face='bold'),
#         axis.title.x = element_text(size = 16),
#         axis.text.x = element_text(size = 16),
#         axis.title.y = element_text(size = 16),
#         axis.text.y = element_text(size = 16),
#         legend.text = element_text(size=14),
#         legend.position="bottom",
#         legend.title = element_text(size=14)) + 
#   ggtitle(paste("Beta value - S_Shore", prettyNum(nrow(bval_sub_S_Shore), big.mark=",", scientific=FALSE), "/", prettyNum(probe_num, big.mark=",", scientific=FALSE), "CpGs")) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_color_manual(values = type_colors) +
#   scale_fill_manual(values = type_colors) +
#   geom_point(aes(colour = factor(TYPE)), size = 1.5) +
#   coord_cartesian(clip = 'off') +
#   theme(plot.margin = margin(0.4, 0.8, 0.4, 0.4, "cm")) +
#   labs(fill='Type') +
#   labs(color='Type') +
#   guides(color = guide_legend(nrow = 2, byrow = TRUE))
# dev.off()
# 
# # M-Value Clustering - Type - S_Shore
# t <- mval_sub_S_Shore_t %>% drop_na(TYPE)
# df <- mval_sub_S_Shore_t_df[intersect(row.names(t), row.names(mval_sub_S_Shore_t_df)), ]
# pca_res <- prcomp(df, scale. = TRUE)
# 
# pdf("mval-pca-plot-type-S_Shore.pdf")
# autoplot(pca_res,
#          data = t,
#          colour = 'TYPE',
#          frame = TRUE
# ) +
#   theme_bw() +
#   theme(title = element_text(size=16, face='bold'),
#         axis.title.x = element_text(size = 16),
#         axis.text.x = element_text(size = 16),
#         axis.title.y = element_text(size = 16),
#         axis.text.y = element_text(size = 16),
#         legend.text = element_text(size=14),
#         legend.position="bottom",
#         legend.title = element_text(size=14)) + 
#   ggtitle(paste("M value - S_Shore", prettyNum(nrow(mval_sub_S_Shore), big.mark=",", scientific=FALSE), "/", prettyNum(probe_num, big.mark=",", scientific=FALSE), "CpGs")) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_color_manual(values = type_colors) +
#   scale_fill_manual(values = type_colors) +
#   geom_point(aes(colour = factor(TYPE)), size = 1.5) +
#   coord_cartesian(clip = 'off') +
#   theme(plot.margin = margin(0.4, 0.8, 0.4, 0.4, "cm")) +
#   labs(fill='Type') +
#   labs(color='Type') +
#   guides(color = guide_legend(nrow = 2, byrow = TRUE))
# dev.off()
# 
# 
# # Beta-Value Clustering - Type - N_Shelf
# t <- bval_sub_N_Shelf_t %>% drop_na(TYPE)
# df <- bval_sub_N_Shelf_t_df[intersect(row.names(t), row.names(bval_sub_N_Shelf_t_df)), ]
# pca_res <- prcomp(df, scale. = TRUE)
# 
# pdf("bval-pca-plot-type-N_Shelf.pdf")
# autoplot(pca_res,
#          data = t,
#          colour = 'TYPE',
#          frame = TRUE
# ) +
#   theme_bw() +
#   theme(title = element_text(size=16, face='bold'),
#         axis.title.x = element_text(size = 16),
#         axis.text.x = element_text(size = 16),
#         axis.title.y = element_text(size = 16),
#         axis.text.y = element_text(size = 16),
#         legend.text = element_text(size=14),
#         legend.position="bottom",
#         legend.title = element_text(size=14)) + 
#   ggtitle(paste("Beta value - N_Shelf", prettyNum(nrow(bval_sub_N_Shelf), big.mark=",", scientific=FALSE), "/", prettyNum(probe_num, big.mark=",", scientific=FALSE), "CpGs")) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_color_manual(values = type_colors) +
#   scale_fill_manual(values = type_colors) +
#   geom_point(aes(colour = factor(TYPE)), size = 1.5) +
#   coord_cartesian(clip = 'off') +
#   theme(plot.margin = margin(0.4, 0.8, 0.4, 0.4, "cm")) +
#   labs(fill='Type') +
#   labs(color='Type') +
#   guides(color = guide_legend(nrow = 2, byrow = TRUE))
# dev.off()
# 
# # M-Value Clustering - Type - N_Shelf
# t <- mval_sub_N_Shelf_t %>% drop_na(TYPE)
# df <- mval_sub_N_Shelf_t_df[intersect(row.names(t), row.names(mval_sub_N_Shelf_t_df)), ]
# pca_res <- prcomp(df, scale. = TRUE)
# 
# pdf("mval-pca-plot-type-N_Shelf.pdf")
# autoplot(pca_res,
#          data = t,
#          colour = 'TYPE',
#          frame = TRUE
# ) +
#   theme_bw() +
#   theme(title = element_text(size=16, face='bold'),
#         axis.title.x = element_text(size = 16),
#         axis.text.x = element_text(size = 16),
#         axis.title.y = element_text(size = 16),
#         axis.text.y = element_text(size = 16),
#         legend.text = element_text(size=14),
#         legend.position="bottom",
#         legend.title = element_text(size=14)) + 
#   ggtitle(paste("M value - N_Shelf", prettyNum(nrow(mval_sub_N_Shelf), big.mark=",", scientific=FALSE), "/", prettyNum(probe_num, big.mark=",", scientific=FALSE), "CpGs")) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_color_manual(values = type_colors) +
#   scale_fill_manual(values = type_colors) +
#   geom_point(aes(colour = factor(TYPE)), size = 1.5) +
#   coord_cartesian(clip = 'off') +
#   theme(plot.margin = margin(0.4, 0.8, 0.4, 0.4, "cm")) +
#   labs(fill='Type') +
#   labs(color='Type') +
#   guides(color = guide_legend(nrow = 2, byrow = TRUE))
# dev.off()
# 
# 
# # Beta-Value Clustering - Type - N_Shore
# t <- bval_sub_N_Shore_t %>% drop_na(TYPE)
# df <- bval_sub_N_Shore_t_df[intersect(row.names(t), row.names(bval_sub_N_Shore_t_df)), ]
# pca_res <- prcomp(df, scale. = TRUE)
# 
# pdf("bval-pca-plot-type-N_Shore.pdf")
# autoplot(pca_res,
#          data = t,
#          colour = 'TYPE',
#          frame = TRUE
# ) +
#   theme_bw() +
#   theme(title = element_text(size=16, face='bold'),
#         axis.title.x = element_text(size = 16),
#         axis.text.x = element_text(size = 16),
#         axis.title.y = element_text(size = 16),
#         axis.text.y = element_text(size = 16),
#         legend.text = element_text(size=14),
#         legend.position="bottom",
#         legend.title = element_text(size=14)) + 
#   ggtitle(paste("Beta value - N_Shore", prettyNum(nrow(bval_sub_N_Shore), big.mark=",", scientific=FALSE), "/", prettyNum(probe_num, big.mark=",", scientific=FALSE), "CpGs")) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_color_manual(values = type_colors) +
#   scale_fill_manual(values = type_colors) +
#   geom_point(aes(colour = factor(TYPE)), size = 1.5) +
#   coord_cartesian(clip = 'off') +
#   theme(plot.margin = margin(0.4, 0.8, 0.4, 0.4, "cm")) +
#   labs(fill='Type') +
#   labs(color='Type') +
#   guides(color = guide_legend(nrow = 2, byrow = TRUE))
# dev.off()
# 
# # M-Value Clustering - Type - N_Shore
# t <- mval_sub_N_Shore_t %>% drop_na(TYPE)
# df <- mval_sub_N_Shore_t_df[intersect(row.names(t), row.names(mval_sub_N_Shore_t_df)), ]
# pca_res <- prcomp(df, scale. = TRUE)
# 
# pdf("mval-pca-plot-type-N_Shore.pdf")
# autoplot(pca_res,
#          data = t,
#          colour = 'TYPE',
#          frame = TRUE
# ) +
#   theme_bw() +
#   theme(title = element_text(size=16, face='bold'),
#         axis.title.x = element_text(size = 16),
#         axis.text.x = element_text(size = 16),
#         axis.title.y = element_text(size = 16),
#         axis.text.y = element_text(size = 16),
#         legend.text = element_text(size=14),
#         legend.position="bottom",
#         legend.title = element_text(size=14)) + 
#   ggtitle(paste("M value - N_Shore", prettyNum(nrow(mval_sub_N_Shore), big.mark=",", scientific=FALSE), "/", prettyNum(probe_num, big.mark=",", scientific=FALSE), "CpGs")) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_color_manual(values = type_colors) +
#   scale_fill_manual(values = type_colors) +
#   geom_point(aes(colour = factor(TYPE)), size = 1.5) +
#   coord_cartesian(clip = 'off') +
#   theme(plot.margin = margin(0.4, 0.8, 0.4, 0.4, "cm")) +
#   labs(fill='Type') +
#   labs(color='Type') +
#   guides(color = guide_legend(nrow = 2, byrow = TRUE))
# dev.off()
