# Load R packages
library(dplyr)
library(stringr)
library(pheatmap)
library(limma)
library(lumi)

# set working directory
setwd("/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/03_MCL_Method/")

# Read in merged methylation data
bval <- readRDS("/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/02_InfiniumPurify/450k/Enmix_quantile_with_QCinfo_Harman_FlowSorted/cluster_run/modified_InfiniumPurify/beta_combat_purified.rds")
dim(bval)

# Read in sample annotations
Sample_Annotations <- read.table(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/ANN/20231031_sample_annotations.tsv", sep = "\t", header = TRUE, na.strings = c("","NA")) %>%
  subset(TYPE!="TFL") %>%
  subset(TIME_POINT!="T2") %>%
  filter(SAMPLE_ID %in% colnames(bval))

# Read in sample sheet with batch information
targets <- read.csv(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/SampleSheet/SampleSheet.csv", skip = 1, header = TRUE, na.strings = c("","NA"))
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
Benign_LN_sample_ann <- read.table(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/Benign_LN/SampleSheet/Benign_LN_sample_annotations.tsv",
                                   sep = "\t",
                                   header = TRUE,
                                   na.strings=c("","NA")) %>%
  dplyr::select(SAMPLE_ID, TYPE)

Benign_LN_sample_ann <- Benign_LN_sample_ann %>% filter(SAMPLE_ID %in% colnames(bval))
Benign_LN_sample_ann$Batch <- "Benign_LN"

# Read in the NIH sample annotation
NIH_sample_ann <- read.table(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/NIH/GSE237299/SampleSheet/GSE237299_sample_annotations.txt",
                             sep = "\t",
                             header = TRUE,
                             na.strings=c("","NA")) %>%
  dplyr::select(SAMPLE_ID, TYPE)

NIH_sample_ann$SAMPLE_ID <- gsub('-', '_', NIH_sample_ann$SAMPLE_ID)
NIH_sample_ann <- NIH_sample_ann %>% filter(SAMPLE_ID %in% colnames(bval))
NIH_sample_ann$Batch <- "GSE237299"

# Read in the EGA sample annotation
EGA_sample_ann <- read.table(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/EGA/SampleSheet/Sample_Annotations_EGAD00010001974.txt",
                             sep = "\t",
                             header = TRUE,
                             na.strings=c("","NA")) %>%
  dplyr::select(SAMPLE_ID, TYPE)

EGA_sample_ann$SAMPLE_ID <- gsub('-', '_', EGA_sample_ann$SAMPLE_ID)
EGA_sample_ann$SAMPLE_ID <- gsub('\\.', '_', EGA_sample_ann$SAMPLE_ID)
EGA_sample_ann$Batch <- "EGAD00010001974"
EGA_sample_ann <- EGA_sample_ann %>% filter(SAMPLE_ID %in% colnames(bval))

# Read in blueprint sample annotation
Blueprint_sample_ann <- read.csv(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/blueprint/Sample_Annotation_EGAS00001001196.txt",
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

# Read in clinical data
clinical_data <- read.csv(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/CLN/20231109_clinical_data.csv", 
                          header = TRUE, na.strings = c("","NA")) %>%
  dplyr::select(LY_FL_ID, ANN_ARBOR_STAGE, SEX, T_14_18, PRIM_TX_CAT, CODE_PFS, CODE_TRANSF, CODE_OS, PFS, TTT, OS)

# Read in the FlowSorted sample annotation
FlowSorted_sample_ann <- read.table(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/FlowSorted/SampleSheet/FlowSorted_Annotation.txt",
                                    sep = "\t",
                                    header = TRUE,
                                    na.strings=c("","NA")) %>%
  dplyr::select(SAMPLE_ID, TYPE, SEX)

FlowSorted_sample_ann$Batch <- 'FlowSorted'
FlowSorted_sample_ann <- FlowSorted_sample_ann %>% filter(SAMPLE_ID %in% colnames(bval))

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
Sample_Annotations <- bind_rows(Sample_Annotations, Benign_LN_sample_ann, NIH_sample_ann, EGA_sample_ann, FlowSorted_sample_ann, Blueprint_sample_ann)
dim(Sample_Annotations)

table(Sample_Annotations$TYPE)
table(Sample_Annotations$Batch)

row.names(Sample_Annotations) <- Sample_Annotations$SAMPLE_ID

# subset of merged methylation data
bval <- bval[, Sample_Annotations$SAMPLE_ID]
dim(bval)

gc()

# Write TYPE and Batch Information
df1 <- Sample_Annotations[c("SAMPLE_ID", "TYPE", "Batch")]
write.table(df1, file="TYPE_INFO.txt", sep="\t", quote = FALSE, row.names = FALSE)

# Read purity table
EpiDISH <- read.table(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/01_EpiDISH/450k/Enmix_quantile_with_QCinfo_Harman_FlowSorted/modified_InfiniumPurify/est.tsv",
                      sep = "\t",
                      header = TRUE,
                      na.strings=c("","NA"))

names(EpiDISH)[names(EpiDISH) == 'X'] <- 'SAMPLE_ID'
EpiDISH <- EpiDISH %>% filter(SAMPLE_ID %in% colnames(bval))

Sample_Annotations_Sub <- Sample_Annotations[, c("SAMPLE_ID", "TYPE", "Batch")]

EpiDISH <- merge(x = EpiDISH,
                 y = Sample_Annotations_Sub,
                 by.x = "SAMPLE_ID",
                 by.y = "SAMPLE_ID",
                 all.x = TRUE)

# read non-bcell beta
#beta_non_Bcells <- readRDS("../beta_combat.rds")
#beta_non_Bcells <- beta_non_Bcells[, c(68:116)]
#dim(beta_non_Bcells)

# Merge both methylation data
#df_beta <- merge(x = beta, y = beta_non_Bcells, by = "row.names", all.x = TRUE)
#dim(df_beta)
#bval <- df_beta[,-1]
#rownames(bval) <- df_beta[,1]
#dim(bval)

# Define tumor matrix with >= 60% B-cell purity
EpiDISH_tumor <- EpiDISH %>% filter(TYPE == "FL" | TYPE == "DLBCL" | TYPE == "DLBCL_EGA" | TYPE == "DLBCL_EN" | TYPE == "FL_A" | TYPE == "FL_B") %>% data.frame()
dim(EpiDISH_tumor) # n = 416

EpiDISH_tumor_60 <- EpiDISH_tumor %>% filter(B >= 60)
dim(EpiDISH_tumor_60) # n = 361

tumor_samples <- EpiDISH_tumor_60 %>% data.frame() %>% .$SAMPLE_ID

Y <- bval[, tumor_samples]
B <- Y
dim(B)

# Construct B2O2 term
Mono_samples <- Sample_Annotations %>% filter(TYPE == "Mono") %>% data.frame() %>% .$SAMPLE_ID
beta_Mono <- bval[, Mono_samples]
beta_Mono_mean <- as.data.frame(rowMeans(beta_Mono))
colnames(beta_Mono_mean) <- "Mono"

Neu_samples <- Sample_Annotations %>% filter(TYPE == "Neu") %>% data.frame() %>% .$SAMPLE_ID
beta_Neu <- bval[, Neu_samples]
beta_Neu_mean <- as.data.frame(rowMeans(beta_Neu))
colnames(beta_Neu_mean) <- "Neutro"

NK_samples <- Sample_Annotations %>% filter(TYPE == "NK") %>% data.frame() %>% .$SAMPLE_ID
beta_NK <- bval[, NK_samples]
beta_NK_mean <- as.data.frame(rowMeans(beta_NK))
colnames(beta_NK_mean) <- "NK"

#CD8T_samples <- Sample_Annotations %>% filter(TYPE == "CD8mem" | TYPE == "CD8nv") %>% data.frame() %>% .$SAMPLE_ID
CD8T_samples <- Sample_Annotations %>% filter(TYPE == "CD8nv") %>% data.frame() %>% .$SAMPLE_ID
beta_CD8T <- bval[, CD8T_samples]
beta_CD8T_mean <- as.data.frame(rowMeans(beta_CD8T))
colnames(beta_CD8T_mean) <- "CD8T"

#CD4T_samples <- Sample_Annotations %>% filter(TYPE == "CD4mem" | TYPE == "CD4nv") %>% data.frame() %>% .$SAMPLE_ID
CD4T_samples <- Sample_Annotations %>% filter(TYPE == "CD4nv") %>% data.frame() %>% .$SAMPLE_ID
beta_CD4T <- bval[, CD4T_samples]
beta_CD4T_mean <- as.data.frame(rowMeans(beta_CD4T))
colnames(beta_CD4T_mean) <- "CD4T"

Eosino_samples <- Sample_Annotations %>% filter(TYPE == "Eos") %>% data.frame() %>% .$SAMPLE_ID
beta_Eosino <- bval[, Eosino_samples]
beta_Eosino_mean <- as.data.frame(rowMeans(beta_Eosino))
colnames(beta_Eosino_mean) <- "Eosino"

df <- cbind(beta_Mono_mean,
            beta_Neu_mean,
            beta_NK_mean,
            beta_CD8T_mean,
            beta_CD4T_mean,
            beta_Eosino_mean)

B2 <- as.matrix(df)
dim(B2)

df1 <- EpiDISH_tumor_60[c("SAMPLE_ID", "Mono", "Neutro", "NK", "CD8T", "CD4T", "Eosino")]
df2 <- df1[,-1]
rownames(df2) <- df1[,1]

df2 <- (df2/100)

O2 <- as.matrix(df2)
dim(O2)

B2O2 <- B2 %*% t(O2)
dim(B2O2)

# Proportion of pure cell types in mixture (cancer) samples
df3 <- EpiDISH_tumor_60[c("SAMPLE_ID", "B", "TYPE")]
df4 <- df3[,-1]
rownames(df4) <- df3[,1]
df4 <- df4[c("B")]

df4 <- (df4/100)

O1 <- as.vector(df4$B)
class(O1)

# B1 term (in-silico purified betas from Y)
# Note: Getting many negative beta values in B1
# because B2O2 is larger than B, mainly due to high B2
B1 <- t(t(B-B2O2) / O1)
bval_final <- B1
dim(bval_final)

# Merge all Normal samples
normal_samples <- Sample_Annotations %>% filter(TYPE == "RLN" | TYPE == "Benign_LN" | TYPE == "LN_NN_A" | TYPE == "LN_NN_B" | 
                                                TYPE == "naiBC" | TYPE == "t_naiBC" | TYPE == "gcBC" | TYPE == "memBC" | TYPE == "bm_PC" | TYPE == "t_PC" | 
                                                TYPE == "Bas" | TYPE == "Bmem" | TYPE == "Bnv" | TYPE == "CD4mem" | TYPE == "CD4nv" | TYPE == "CD8mem" | TYPE == "CD8nv" | TYPE == "Eos" | TYPE == "MIX" | TYPE == "Mono" | TYPE == "Neu" | TYPE == "NK" | TYPE == "Treg") %>% data.frame() %>% .$SAMPLE_ID

normal <- bval[, normal_samples]
dim(normal)

# merge purified methylation data with normal
bval_subset <- merge(x = bval_final,
                     y = normal,
                     by = "row.names",
                     all.x = TRUE,
                     all.y = TRUE)

dim(bval_subset)

bval_final <- bval_subset[,-1]
rownames(bval_final) <- bval_subset[,1]
dim(bval_final)
gc()

# Convert methylation Beta-value to M-value
# Note: There were 50 or more warnings (use warnings() to see the first 50)
mval_final <- beta2m(bval_final)
dim(mval_final)

# Write purified methylation data
saveRDS(bval_final, "beta_combat_purified_MCL.rds")
saveRDS(mval_final, "mval_combat_purified_MCL.rds")
