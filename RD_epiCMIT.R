# R Packages
# ==========
options(stringsAsFactors = F, error = NULL)
library(GenomicRanges)
library(data.table)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
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

# Change directory
setwd("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/06_epiCMIT/450k/modified_InfiniumPurify/Enmix_quantile_with_QCinfo_Harman_All/02_After_Clustering/")

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

# Read in methylation clustering information
epitypes <- read.table(file = "/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/02_Clustering/450k/modified_InfiniumPurify/Enmix_quantile_with_QCinfo_Harman_All/FL_epitypes.tsv",
                       sep = "\t", header = TRUE, na.strings = c("","NA"))

epitypes <- epitypes %>% filter(SAMPLE_ID %in% colnames(bval))

Sample_Annotations <- merge(x = Sample_Annotations,
                            y = epitypes,
                            by.x = "SAMPLE_ID",
                            by.y = "SAMPLE_ID",
                            all.x = TRUE)

dim(Sample_Annotations)
table(Sample_Annotations$TYPE)
table(Sample_Annotations$epitypes)

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

Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'Normal_GCB'] <- 'gcBC'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'bm_PC'] <- 'bm_PC'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'gcBC'] <- 'gcBC'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'memBC'] <- 'memBC'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'naiBC'] <- 'naiBC'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 't_PC'] <- 't_PC'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 't_naiBC'] <- 't_naiBC'

# Replace NAs in "epitypes" with corresponding blueprint b-cell type from "TYPE"
Sample_Annotations$epitypes[is.na(Sample_Annotations$epitypes)] <- Sample_Annotations$TYPE[is.na(Sample_Annotations$epitypes)]

table(Sample_Annotations$Batch)
table(Sample_Annotations$TYPE)
table(Sample_Annotations$epitypes)

# Rename C1, C2, C3, and UC in the epitypes column by appending their corresponding TYPE values
Sample_Annotations <- Sample_Annotations %>%
  mutate(epitypes = case_when(
    epitypes %in% c("C1", "C2", "C3", "UC") ~ paste0(epitypes, "_", TYPE),  # Rename with TYPE
    TRUE ~ epitypes  # Keep other values unchanged
  ))

table(Sample_Annotations$Batch)
table(Sample_Annotations$TYPE)
table(Sample_Annotations$epitypes)

# Filter epitypes
Sample_Annotations <- Sample_Annotations %>% filter(epitypes == "C1_FL"
                                                    | epitypes == "C2_FL"
                                                    | epitypes == "C1_t_DLBCL"
                                                    | epitypes == "C2_t_DLBCL"
                                                    | epitypes == "C1_DLBCL"
                                                    | epitypes == "C2_DLBCL"
                                                    ) %>% data.frame()

# Sample_Annotations <- Sample_Annotations %>% filter(epitypes == "C1_FL" | epitypes == "C2_FL") %>% data.frame()

table(Sample_Annotations$Batch)
table(Sample_Annotations$TYPE)
table(Sample_Annotations$epitypes)

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

# Load epiCMIT Data
load("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/epiCMIT/Estimate.epiCMIT.RData")
set.seed(1234)

DNAm.epiCMIT <- DNAm.to.epiCMIT(DNAm = data.frame(bval_sub),
                                DNAm.genome.assembly = "hg19",
                                map.DNAm.to = "Illumina.450K.epiCMIT",
                                min.epiCMIT.CpGs = 800)

# 1014 from 1348 epiCMIT-CpGs are present in your matrix - Now - Final
# 1230 from 1348 epiCMIT-CpGs are present in your matrix! # Old
# 1187 from 1348 epiCMIT-CpGs are present in your matrix! # bval
# 1186 from 1348 epiCMIT-CpGs are present in your matrix! # bval_naive_diff
# 1 from 1348 epiCMIT-CpGs are present in your matrix! # bval_naive_like

epiCMIT.Illumina <- epiCMIT(DNAm.epiCMIT = DNAm.epiCMIT,
                            return.epiCMIT.annot = FALSE,
                            export.results = TRUE,
                            export.results.dir = "",
                            export.results.name = "/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/06_epiCMIT/450k/modified_InfiniumPurify/Enmix_quantile_with_QCinfo_Harman_All/02_After_Clustering/")

head(epiCMIT.Illumina$epiCMIT.scores)
epiCMIT.Illumina$epiCMIT.run.info

# epiCMIT - epitypes
#epitypes_colors <- c(C1="#6dcc04", C2="#fc1417", C3="#7618dc", na.value="#c9c9c9")
#epitypes_colors <- c(C1="#fc1417", C2="#6dcc04", C3="#7618dc", na.value="#c9c9c9")
#epitypes_colors <- c(C1="#6dcc04", C2="#fc1417", C3="#7618dc", bm_PC="#7f8992", gcBC="#0303a8", memBC="#b5b80b", naiBC="#2480b9", t_naiBC="#6e988c", t_PC="#383838", na.value="#c9c9c9")

# Rename epitypes
table(Sample_Annotations$epitypes)

# Sample_Annotations$epitypes[Sample_Annotations$epitypes == 'C3_Normal'] <- 'Normal'
Sample_Annotations$epitypes[Sample_Annotations$epitypes == 'C1_FL'] <- 'iFL'
Sample_Annotations$epitypes[Sample_Annotations$epitypes == 'C2_FL'] <- 'aFL'

table(Sample_Annotations$epitypes)

# Define custom colors for improved visibility
epitypes_colors <- c(
  "iFL" = "#CFE7B9",
  "aFL" = "#59A14F",
  "C1_t_DLBCL" = "#FAE27F",
  "C2_t_DLBCL" = "#F28E2B",
  "C1_DLBCL" = "#F6C2C5",
  "C2_DLBCL" = "#E15759")

# # Define custom colors for improved visibility
# epitypes_colors <- c(
#   "iFL" = "#59A14F",
#   "aFL" = "#E15759"
# )

epiCMIT.Illumina$epiCMIT.scores$Samples <- gsub("^X", "", epiCMIT.Illumina$epiCMIT.scores$Samples)

long_df <- epiCMIT.Illumina$epiCMIT.scores %>%
  left_join(Sample_Annotations, by = c("Samples" = "SAMPLE_ID")) %>%
  filter(!is.na(epitypes)) %>%
  # mutate(epitypes = factor(epitypes, levels = c("t_naiBC", "naiBC", "gcBC", "t_PC", "C3", "C1", "C2", "memBC", "bm_PC"))) %>%
  # mutate(epitypes = factor(epitypes, levels = c("C3_Normal", "iFL", "aFL", "C1_t_DLBCL", "C2_t_DLBCL", "C1_DLBCL", "C2_DLBCL", "UC_FL", "UC_DLBCL", "naiBC",  "t_naiBC", "gcBC", "t_PC", "memBC", "bm_PC"))) %>% 
  mutate(epitypes = factor(epitypes, levels = c("iFL", "aFL", "C1_t_DLBCL", "C2_t_DLBCL", "C1_DLBCL", "C2_DLBCL"))) %>% 
  mutate(epiCMIT_ID = "") %>%
  select(Samples, epitypes, epiCMIT_ID, epiCMIT)

# View the first few rows of the reshaped data
head(long_df)

# Extract boxplot statistics for Tumor_Purity by TYPE
type_stats <- tapply(long_df$epiCMIT, long_df$epitypes, boxplot.stats)

# View statistics for each TYPE
type_stats

# •	Minimum (1st value): The lowest value in the data.
# •	1st Quartile (Q1, 2nd value): The 25th percentile (25% of the data points are below this value).
# •	Median (Q2, 3rd value): The 50th percentile (middle value of the data).
# •	3rd Quartile (Q3, 4th value): The 75th percentile (75% of the data points are below this value).
# •	Maximum (5th value): The highest value in the data.

# $iFL
# $iFL$stats
# [1] 0.6954173 0.7612931 0.7857660 0.8063727 0.8614118
# 
# $iFL$n
# [1] 222
# 
# $iFL$conf
# [1] 0.7809856 0.7905464
# 
# $iFL$out
# [1] 0.6743138 0.6842540 0.6879323 0.6221738 0.5899723 0.6744949
# 
# 
# $aFL
# $aFL$stats
# [1] 0.7797235 0.8144768 0.8312173 0.8467583 0.8732342
# 
# $aFL$n
# [1] 83
# 
# $aFL$conf
# [1] 0.8256188 0.8368158
# 
# $aFL$out
# numeric(0)
# 
# 
# $C1_t_DLBCL
# $C1_t_DLBCL$stats
# [1] 0.7410761 0.7650637 0.7738221 0.7891412 0.8240480
# 
# $C1_t_DLBCL$n
# [1] 5
# 
# $C1_t_DLBCL$conf
# [1] 0.7568090 0.7908352
# 
# $C1_t_DLBCL$out
# numeric(0)
# 
# 
# $C2_t_DLBCL
# $C2_t_DLBCL$stats
# [1] 0.7908360 0.8028508 0.8158959 0.8484203 0.8746534
# 
# $C2_t_DLBCL$n
# [1] 7
# 
# $C2_t_DLBCL$conf
# [1] 0.7886826 0.8431092
# 
# $C2_t_DLBCL$out
# numeric(0)
# 
# 
# $C1_DLBCL
# $C1_DLBCL$stats
# [1] 0.7027732 0.7588054 0.7804391 0.8047145 0.8388917
# 
# $C1_DLBCL$n
# [1] 21
# 
# $C1_DLBCL$conf
# [1] 0.7646103 0.7962678
# 
# $C1_DLBCL$out
# numeric(0)
# 
# 
# $C2_DLBCL
# $C2_DLBCL$stats
# [1] 0.7551361 0.8121730 0.8397696 0.8589350 0.9088802
# 
# $C2_DLBCL$n
# [1] 118
# 
# $C2_DLBCL$conf
# [1] 0.8329680 0.8465712
# 
# $C2_DLBCL$out
# [1] 0.71553

# Define custom order for sample types and Relation to Island
# custom_order <- c("C3_Normal", "iFL", "aFL", "C1_t_DLBCL", "C2_t_DLBCL", "C1_DLBCL", "C2_DLBCL", "UC_FL", "UC_DLBCL", "naiBC",  "t_naiBC", "gcBC", "t_PC", "memBC", "bm_PC")
# custom_order <- c("C3_Normal", "iFL", "aFL", "C1_t_DLBCL", "C2_t_DLBCL", "C1_DLBCL", "C2_DLBCL", "UC_FL", "UC_DLBCL")
# long_df$epitypes <- factor(long_df$epitypes, levels = custom_order)

# long_df$Relation_to_Island <- factor(long_df$Relation_to_Island, levels = c("All_Regions", "CpG_Islands", "OpenSea", "S_Shelf", "N_Shelf", "S_Shore", "N_Shore"))

# # Set up PDF output
# pdf("Proliferative_History.pdf", width = 8, height = 10)  # Adjusted width and height for better readability
# 
# # Generate grouped boxplot using custom colors
# ggplot(long_df, aes(x = epiCMIT_ID, y = epiCMIT, fill = epitypes)) + 
#   geom_boxplot(position = position_dodge(width = 0.9), size = 0.7) + 
#   facet_grid(~epiCMIT_ID, scales = "free_x", space = "free_x") +  # Facet to separate sample types
#   # facet_wrap(~Relation_to_Island, scales = "free_x", ncol = 4, nrow = 2) + 
#   scale_fill_manual(values = epitypes_colors, name = "Epitypes") + 
#   scale_y_continuous(breaks = seq(0, 1, by = 0.1)) + 
#   ggtitle("Proliferative History") + 
#   labs(x = "", y = "epiCMIT Score") + 
#   theme(
#     plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
#     legend.position = "bottom",
#     legend.text = element_text(size = 14),
#     legend.title = element_text(size = 20, face = "bold"),
#     legend.key.size = unit(0.8, "cm"),
#     legend.margin = margin(l = 40),
#     axis.ticks.x = element_blank(),
#     axis.text.x = element_text(size = 16, face = "plain"),
#     axis.text.y = element_text(size = 16, face = "plain"),
#     axis.title.x = element_text(size = 20, face = "bold"),
#     axis.title.y = element_text(size = 20, face = "bold"),
#     #plot.margin = margin(t = 20, r = 20, b = 0, l = 20),
#     panel.spacing.x = unit(1, "lines"),
#     strip.text = element_blank(),
#     strip.background = element_blank(),
#     panel.spacing = unit(2, "cm"),
#     panel.background = element_blank(),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     panel.grid.major.y = element_line(color = "gray85", linewidth = 0.5),
#     panel.grid.minor.y = element_line(color = "gray95", linewidth = 0.3)
#     # legend.position = c(0.9, 0.2),  # Adjust legend position
#     # legend.justification = c(0.8, 0),  # Align legend properly
#     # legend.direction = "horizontal"  # Make legend horizontal
#   ) + 
#   guides(fill = guide_legend(byrow = TRUE))
# 
# # Close the PDF output
# dev.off()

# Comparisons
my_comparisons <- list(c("iFL", "aFL"),
                       c("C1_t_DLBCL", "C2_t_DLBCL"),
                       c("C1_DLBCL", "C2_DLBCL"))

pdf("Proliferative_History.pdf", width = 6, height = 8)  # Adjusted width and height for better readability
ggplot(long_df, aes(x = epitypes, y = epiCMIT, fill = epitypes)) + geom_boxplot(outlier.shape = 16, outlier.size = 1, outlier.color = "black", width = 0.8, color = "black", size = 0.3) +
  # ggboxplot(df, x = "TYPE", y = "mean", color = "TYPE", add = "jitter", size = 0.3, add.params = list(size = 0.15)) +
  scale_color_manual(values = epitypes_colors) +
  scale_fill_manual(values = epitypes_colors) +
  theme_bw() +
  theme(
    title = element_text(size = 20, face = 'bold'),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # axis.title.x = element_text(size = 12),
    # axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    # axis.title.y = element_blank(),
    # axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(),
    axis.title.y = element_text(size = 20, margin = margin(r = 10)),
    axis.text.y = element_text(size = 20),
    # legend.position = "none",
    axis.line = element_line(),
    #axis.line.x = element_blank(),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20, face = "bold"),
    legend.key.size = unit(1, "cm"),
    legend.box.spacing = unit(0, "cm"),
    legend.margin = margin(10, 0, 0, 0),
    panel.spacing = unit(2, "cm"),
    panel.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray85", linewidth = 0.3),
    panel.grid.minor.y = element_line(color = "gray95", linewidth = 0.3),
    panel.border = element_blank()
  ) +
  guides(fill = guide_legend(nrow = 3, ncol = 2, byrow = TRUE)) +  # Adjust legend into 3 rows and 2 columns
  ylab("epiCMIT Score") +
  labs(fill = "Epitype") +
  ggtitle("Proliferative History") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  # Adjust the size of p-values and statistical information
  stat_compare_means(size = 6, comparisons = my_comparisons,
                     label.y = c(max(long_df$epiCMIT) - 0.03,
                                 max(long_df$epiCMIT) - 0.03,
                                 max(long_df$epiCMIT) + 0.00)) +  # Adjusted y position
  stat_compare_means(size = 6, label.x = 2, label.y = max(long_df$epiCMIT) + 0.06) +  # Adjusted label position
  
  coord_cartesian(clip = 'off') +
  theme(plot.margin = margin(0.4, 1, 0.4, 0.4, "cm")) +
  scale_x_discrete(name = "Type") +
  # scale_y_continuous(limits = c(min(df$mean) - 0.002, max(df$mean) + 0.06), breaks = seq(0, 1, by = 0.1))  # Adjusted y limits
  scale_y_continuous(limits = c(0.58, 1), breaks = seq(0, 1, by = 0.1))  # Adjusted y limits
# scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1))  # Adjusted y limits

dev.off()


# 
# pdf("epiCMIT_epitypes.pdf", width = 5, height = 6)
# epiCMIT.Illumina$epiCMIT.scores %>%
#   left_join(Sample_Annotations, by = c("Samples" = "SAMPLE_ID")) %>%
#   filter(!is.na(epitypes)) %>%
#   # mutate(epitypes = factor(epitypes, levels = c("t_naiBC", "naiBC", "gcBC", "t_PC", "C3", "C1", "C2", "memBC", "bm_PC"))) %>%
#   mutate(epitypes = factor(epitypes, levels = c("C3_Normal", "iFL", "aFL", "C1_t_DLBCL", "C2_t_DLBCL", "C1_DLBCL", "C2_DLBCL", "UC_FL", "UC_DLBCL", "naiBC",  "t_naiBC", "gcBC", "t_PC", "memBC", "bm_PC"))) %>%
#   ggpubr::ggboxplot(x = "epitypes", y = "epiCMIT", add = "jitter", color = "epitypes", size = 0.3, add.params = list(size = 0.5)) +
#   scale_color_manual(values = epitypes_colors) +
#   scale_fill_manual(values = epitypes_colors) +
#   theme_bw() +
#   theme(title = element_text(size=18, face='bold'),
#         axis.title.x = element_text(size = 18),
#         axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust=1),
#         axis.title.y = element_text(size = 18),
#         axis.text.y = element_text(size = 16),
#         legend.position = "none") +
#   ylab("epiCMIT Score") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   ggtitle("Proliferative History") +
#   stat_compare_means(size = 5, label.x = 4, label.y = max(epiCMIT.Illumina$epiCMIT.scores$epiCMIT)+0.03) +
#   coord_cartesian(clip = 'off') +
#   theme(plot.margin = margin(0.4, 0.4, 0.4, 0.4, "cm")) +
#   scale_y_continuous(limits=c(min(epiCMIT.Illumina$epiCMIT.scores$epiCMIT)-0.01,
#                               max(epiCMIT.Illumina$epiCMIT.scores$epiCMIT)+0.05)) +
#   scale_x_discrete(name ="Epitypes")
# dev.off()
