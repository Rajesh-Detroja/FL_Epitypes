# Loading R Packages ----
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
library(patchwork)

# Change directory
setwd("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/04_Summary/450k/modified_InfiniumPurify/Enmix_quantile_with_QCinfo_Harman_All/02_After_Clustering/")

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

# Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'Normal_GCB'] <- 'gcBC'
# Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'bm_PC'] <- 'bm_PC'
# Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'gcBC'] <- 'gcBC'
# Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'memBC'] <- 'memBC'
# Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'naiBC'] <- 'naiBC'
# Sample_Annotations$TYPE[Sample_Annotations$TYPE == 't_PC'] <- 't_PC'
# Sample_Annotations$TYPE[Sample_Annotations$TYPE == 't_naiBC'] <- 't_naiBC'

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
Sample_Annotations <- Sample_Annotations %>% filter(epitypes == "C3_Normal" | epitypes == "C1_FL" | epitypes == "C2_FL") %>% data.frame()

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

# Methylation Summary
# ===================
#type_colors <- c(Normal="#7618dc", FL="#6dcc04", DLBCL="#fc1417", t_DLBCL="#F28E2B", bm_PC="#7f8992", gcBC="#0303a8", memBC="#b5b80b", naiBC="#2480b9", t_naiBC="#6e988c", t_PC="#383838", na.value="#c9c9c9")
#epitypes_colors <- c(Benign_LN="#5086F3", bm_PC="#7f8992", DLBCL="#A01641", DLBCL_EGA="#fc1417", DLBCL_EN="#f26ded", FL="#6dcc04", FL_A="#02bdaa", FL_B="#ebc602", gcBC="#0303a8", LN_NN_A="#FF7F00", LN_NN_B="#B15928", memBC="#b5b80b", naiBC="#2480b9", RLN="#7618dc", t_naiBC="#6e988c", t_PC="#383838", na.value="#c9c9c9")
#epitypes_colors <- c(C1="#6dcc04", C2="#fc1417", C3="#7618dc", bm_PC="#7f8992", gcBC="#0303a8", memBC="#b5b80b", naiBC="#2480b9", t_naiBC="#6e988c", t_PC="#383838", na.value="#c9c9c9")

# # Define custom colors for improved visibility
# epitypes_colors <- c(
#   "C3_Normal" = "#5086F3",
#   "C1_FL" = "#CFE7B9", 
#   "C2_FL" = "#59A14F",
#   "C1_t_DLBCL" = "#FAE27F", 
#   "C2_t_DLBCL" = "#F28E2B",
#   "C1_DLBCL" = "#F6C2C5", 
#   "C2_DLBCL" = "#E15759",
#   "UC_FL" = "#DEE4EA", 
#   "UC_DLBCL" = "#7f8992"
#   # "naiBC" = "#8C6D31", 
#   # "t_naiBC" = "#D1A1A6",
#   # "gcBC" = "#E4B321",
#   # "t_PC" = "#76B7B2",
#   # "memBC" = "#FF9DA7",
#   # "bm_PC" = "#6A4D8C"
# )

# Rename epitypes
table(Sample_Annotations$epitypes)

Sample_Annotations$epitypes[Sample_Annotations$epitypes == 'C3_Normal'] <- 'Normal'
Sample_Annotations$epitypes[Sample_Annotations$epitypes == 'C1_FL'] <- 'iFL'
Sample_Annotations$epitypes[Sample_Annotations$epitypes == 'C2_FL'] <- 'aFL'

table(Sample_Annotations$epitypes)

# Define custom colors for improved visibility
epitypes_colors <- c(
  "Normal" = "#5086F3",
  "iFL" = "#59A14F",
  "aFL" = "#E15759"
)

# Methylation Summary
####
#type_colors <- c(Normal="#7618dc", FL="#6dcc04", DLBCL="#fc1417", t_DLBCL="#F28E2B", bm_PC="#7f8992", gcBC="#0303a8", memBC="#b5b80b", naiBC="#2480b9", t_naiBC="#6e988c", t_PC="#383838", na.value="#c9c9c9")

# Define custom colors for improved visibility
# Define Custom Colors ----
# type_colors <- c(
# "Normal" = "#5086F3",
# "FL" = "#59A14F",
# "t_DLBCL" = "#F28E2B",
# "DLBCL" = "#E15759"
# "naiBC" = "#8C6D31", 
# "t_naiBC" = "#D1A1A6",
# "gcBC" = "#E4B321",
# "t_PC" = "#76B7B2",
# "memBC" = "#FF9DA7",
# "bm_PC" = "#6A4D8C"
# )

# Comparisons
my_comparisons <- list(c("iFL", "aFL"))

# Mean methylation types - All probes
# pdf("Mean-Methylation-All-CpGs.pdf", height = 6, width = 3)
df <- bval_sub %>%
  data.frame(check.names = FALSE) %>%
  mutate(probe = row.names(.)) %>%
  melt() %>%
  left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE")], by = c("variable" = "SAMPLE_ID")) %>%
  group_by(variable) %>%
  summarize(mean = mean(value)) %>%
  left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE", "epitypes")], by = c("variable" = "SAMPLE_ID")) %>%
  mutate(epitypes = factor(epitypes, levels = c("Normal", "iFL", "aFL"))) %>% 
  filter(epitypes != "NA")

# Filter data for only "FL" type
df_FL <- df %>% filter(TYPE == "FL")
head(df_FL)

# Create the plot
# P1 ----
P1 <- ggplot(df, aes(x = epitypes, y = mean, fill = epitypes)) + geom_boxplot(outlier.shape = 16, outlier.size = 1, outlier.color = "black", width = 0.8, color = "black", size = 0.3) +
  # ggboxplot(df, x = "TYPE", y = "mean", color = "TYPE", add = "jitter", size = 0.3, add.params = list(size = 0.15)) +
  scale_color_manual(values = epitypes_colors) +
  scale_fill_manual(values = epitypes_colors) +
  theme_bw() +
  theme(
    title = element_text(size = 14, face = 'bold'),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # axis.title.x = element_text(size = 12),
    # axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.title.y = element_text(size = 18, margin = margin(r = 10), face = 'bold'),
    axis.text.y = element_text(size = 16),
    axis.line = element_line(),
    #axis.line.x = element_blank(),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    legend.position = "none",
    panel.spacing = unit(2, "cm"),
    panel.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray85", linewidth = 0.3),
    panel.grid.minor.y = element_line(color = "gray95", linewidth = 0.3),
    panel.border = element_blank()
  ) +
  ylab("Mean Methylation Level (beta-value)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(sprintf("All\n(%s CpGs)", prettyNum(probe_num, big.mark=",", scientific=FALSE))) +
  
  # Adjust the size of p-values and statistical information
  stat_compare_means(size = 4.5, comparisons = my_comparisons,
                     label.y = c(max(df_FL$mean) + 0.01, max(df_FL$mean) + 0.05)) +  # Adjusted y position
  stat_compare_means(size = 4.5, label.x = 1.2, label.y = max(df$mean) + 0.09) +  # Adjusted label position
  
  coord_cartesian(clip = 'off') +
  theme(plot.margin = margin(0.4, 1, 0.4, 0.4, "cm")) +
  scale_x_discrete(name = "Type") +
  # scale_y_continuous(limits = c(min(df$mean) - 0.002, max(df$mean) + 0.06), breaks = seq(0, 1, by = 0.1))  # Adjusted y limits
  scale_y_continuous(limits = c(0.2, 0.8), breaks = seq(0, 1, by = 0.1))  # Adjusted y limits
# scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1))  # Adjusted y limits

# dev.off()


# Mean methylation types - Island
bval_sub_CpG_islands <- bval_sub[intersect(CpG_islands, row.names(bval_sub)),]

# pdf("Mean-Methylation-Island.pdf", height = 6, width = 3)
df <- bval_sub[intersect(CpG_islands, row.names(bval_sub)),] %>% 
  data.frame(check.names = FALSE) %>%
  mutate(probe = row.names(.)) %>%
  melt() %>%
  left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE")], by = c("variable" = "SAMPLE_ID")) %>%
  group_by(variable) %>%
  summarize(mean = mean(value)) %>%
  left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE", "epitypes")], by = c("variable" = "SAMPLE_ID")) %>%
  mutate(epitypes = factor(epitypes, levels = c("Normal", "iFL", "aFL"))) %>% 
  filter(epitypes != "NA") 

# Filter data for only "FL" type
df_FL <- df %>% filter(TYPE == "FL")
head(df_FL)

# Create the plot
# P2 ----
P2 <- ggplot(df, aes(x = epitypes, y = mean, fill = epitypes)) + geom_boxplot(outlier.shape = 16, outlier.size = 1, outlier.color = "black", width = 0.8, color = "black", size = 0.3) +
  # ggboxplot(df, x = "TYPE", y = "mean", color = "TYPE", add = "jitter", size = 0.3, add.params = list(size = 0.15)) +
  scale_color_manual(values = epitypes_colors) +
  scale_fill_manual(values = epitypes_colors) +
  theme_bw() +
  theme(
    title = element_text(size = 14, face = 'bold'),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # axis.title.x = element_text(size = 12),
    # axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    # axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    # axis.text.y = element_text(size = 14),
    axis.line = element_line(),
    axis.line.y = element_blank(),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    legend.position = "none",
    panel.spacing = unit(2, "cm"),
    panel.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray85", linewidth = 0.3),
    panel.grid.minor.y = element_line(color = "gray95", linewidth = 0.3),
    panel.border = element_blank()
  ) +
  ylab("Mean Methylation Level (Beta Value)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(sprintf("CpG Islands\n(%s CpGs)", prettyNum(nrow(bval_sub_CpG_islands), big.mark=",", scientific=FALSE))) +
  
  # Adjust the size of p-values and statistical information
  stat_compare_means(size = 4.5, comparisons = my_comparisons,
                     label.y = c(max(df_FL$mean) + 0.009, max(df_FL$mean) + 0.05)) +  # Adjusted y position
  stat_compare_means(size = 4.5, label.x = 1, label.y = max(df$mean) + 0.1) +  # Adjusted label position
  
  coord_cartesian(clip = 'off') +
  theme(plot.margin = margin(0.4, 1, 0.4, 0.4, "cm")) +
  scale_x_discrete(name = "Type") +
  # scale_y_continuous(limits = c(min(df$mean) - 0.002, max(df$mean) + 0.06), breaks = seq(0, 1, by = 0.1))  # Adjusted y limits
  scale_y_continuous(limits = c(0.2, 0.8), breaks = seq(0, 1, by = 0.1))  # Adjusted y limits
# scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1))  # Adjusted y limits

# dev.off()


# Mean methylation types - OpenSea
bval_sub_OpenSea <- bval_sub[intersect(OpenSea, row.names(bval_sub)),]

# pdf("Mean-Methylation-OpenSea.pdf", height = 6, width = 3)
df <- bval_sub[intersect(OpenSea, row.names(bval_sub)),] %>%
  data.frame(check.names = FALSE) %>%
  mutate(probe = row.names(.)) %>%
  melt() %>%
  left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE")], by = c("variable" = "SAMPLE_ID")) %>%
  group_by(variable) %>%
  summarize(mean = mean(value)) %>%
  left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE", "epitypes")], by = c("variable" = "SAMPLE_ID")) %>%
  mutate(epitypes = factor(epitypes, levels = c("Normal", "iFL", "aFL"))) %>% 
  filter(epitypes != "NA")

# Filter data for only "FL" type
df_FL <- df %>% filter(TYPE == "FL")
head(df_FL)

# Create the plot
# P3 ----
P3 <- ggplot(df, aes(x = epitypes, y = mean, fill = epitypes)) + geom_boxplot(outlier.shape = 16, outlier.size = 1, outlier.color = "black", width = 0.8, color = "black", size = 0.3) +
  # ggboxplot(df, x = "TYPE", y = "mean", color = "TYPE", add = "jitter", size = 0.3, add.params = list(size = 0.15)) +
  scale_color_manual(values = epitypes_colors) +
  scale_fill_manual(values = epitypes_colors) +
  theme_bw() +
  theme(
    title = element_text(size = 14, face = 'bold'),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # axis.title.x = element_text(size = 12),
    # axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    # axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    # axis.text.y = element_text(size = 14),
    axis.line = element_line(),
    axis.line.y = element_blank(),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    legend.position = "none",
    panel.spacing = unit(2, "cm"),
    panel.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray85", linewidth = 0.3),
    panel.grid.minor.y = element_line(color = "gray95", linewidth = 0.3),
    panel.border = element_blank()
  ) +
  ylab("Mean Methylation Level (Beta Value)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(sprintf("OpenSea\n(%s CpGs)", prettyNum(nrow(bval_sub_OpenSea), big.mark=",", scientific=FALSE))) +
  
  # Adjust the size of p-values and statistical information
  stat_compare_means(size = 4.5, comparisons = my_comparisons,
                     label.y = c(max(df_FL$mean) + 0.001, max(df_FL$mean) + 0.04)) +  # Adjusted y position
  stat_compare_means(size = 4.5, label.x = 1, label.y = max(df$mean) + 0.02) +  # Adjusted label position
  
  coord_cartesian(clip = 'off') +
  theme(plot.margin = margin(0.4, 1, 0.4, 0.4, "cm")) +
  scale_x_discrete(name = "Type") +
  # scale_y_continuous(limits = c(min(df$mean) - 0.002, max(df$mean) + 0.06), breaks = seq(0, 1, by = 0.1))  # Adjusted y limits
  scale_y_continuous(limits = c(0.2, 0.8), breaks = seq(0, 1, by = 0.1))  # Adjusted y limits
# scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1))  # Adjusted y limits

# dev.off()

# Mean methylation types - Shelf
Shelf <- c(S_Shelf, N_Shelf)
bval_sub_Shelf <- bval_sub[intersect(Shelf, row.names(bval_sub)),]

# pdf("Mean-Methylation-Shelf.pdf", height = 6, width = 3)
df <- bval_sub[intersect(Shelf, row.names(bval_sub)),] %>% 
  data.frame(check.names = FALSE) %>%
  mutate(probe = row.names(.)) %>%
  melt() %>%
  left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE")], by = c("variable" = "SAMPLE_ID")) %>%
  group_by(variable) %>%
  summarize(mean = mean(value)) %>%
  left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE", "epitypes")], by = c("variable" = "SAMPLE_ID")) %>%
  mutate(epitypes = factor(epitypes, levels = c("Normal", "iFL", "aFL"))) %>% 
  filter(epitypes != "NA")

# Filter data for only "FL" type
df_FL <- df %>% filter(TYPE == "FL")
head(df_FL)

# Create the plot
# P4 ----
P4 <- ggplot(df, aes(x = epitypes, y = mean, fill = epitypes)) + geom_boxplot(outlier.shape = 16, outlier.size = 1, outlier.color = "black", width = 0.8, color = "black", size = 0.3) +
  # ggboxplot(df, x = "TYPE", y = "mean", color = "TYPE", add = "jitter", size = 0.3, add.params = list(size = 0.15)) +
  scale_color_manual(values = epitypes_colors) +
  scale_fill_manual(values = epitypes_colors) +
  theme_bw() +
  theme(
    title = element_text(size = 14, face = 'bold'),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # axis.title.x = element_text(size = 12),
    # axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    # axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    # axis.text.y = element_text(size = 14),
    axis.line = element_line(),
    axis.line.y = element_blank(),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    legend.position = "none",
    panel.spacing = unit(2, "cm"),
    panel.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray85", linewidth = 0.3),
    panel.grid.minor.y = element_line(color = "gray95", linewidth = 0.3),
    panel.border = element_blank()
  ) +
  ylab("Mean Methylation Level (Beta Value)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(sprintf("Shelf\n(%s CpGs)", prettyNum(nrow(bval_sub_Shelf), big.mark=",", scientific=FALSE))) +
  
  # Adjust the size of p-values and statistical information
  stat_compare_means(size = 4.5, comparisons = my_comparisons,
                     label.y = c(max(df_FL$mean) + 0.001, max(df_FL$mean) + 0.05)) +  # Adjusted y position
  stat_compare_means(size = 4.5, label.x = 1, label.y = max(df$mean) + 0.01) +  # Adjusted label position
  
  coord_cartesian(clip = 'off') +
  theme(plot.margin = margin(0.4, 1, 0.4, 0.4, "cm")) +
  scale_x_discrete(name = "Type") +
  # scale_y_continuous(limits = c(min(df$mean) - 0.002, max(df$mean) + 0.06), breaks = seq(0, 1, by = 0.1))  # Adjusted y limits
  scale_y_continuous(limits = c(0.2, 0.8), breaks = seq(0, 1, by = 0.1))  # Adjusted y limits
# scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1))  # Adjusted y limits

# dev.off()

# Mean methylation types - Shore
Shore <- c(S_Shore, N_Shore)
bval_sub_Shore <- bval_sub[intersect(Shore, row.names(bval_sub)),]

# pdf("Mean-Methylation-Shore.pdf", height = 6, width = 3)
df <- bval_sub[intersect(Shore, row.names(bval_sub)),] %>%
  data.frame(check.names = FALSE) %>%
  mutate(probe = row.names(.)) %>%
  melt() %>%
  left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE")], by = c("variable" = "SAMPLE_ID")) %>%
  group_by(variable) %>%
  summarize(mean = mean(value)) %>%
  left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE", "epitypes")], by = c("variable" = "SAMPLE_ID")) %>%
  mutate(epitypes = factor(epitypes, levels = c("Normal", "iFL", "aFL"))) %>% 
  filter(epitypes != "NA")

# Filter data for only "FL" type
df_FL <- df %>% filter(TYPE == "FL")
head(df_FL)

# Create the plot
# P5 ----
P5 <- ggplot(df, aes(x = epitypes, y = mean, fill = epitypes)) + geom_boxplot(outlier.shape = 16, outlier.size = 1, outlier.color = "black", width = 0.8, color = "black", size = 0.3) +
  # ggboxplot(df, x = "TYPE", y = "mean", color = "TYPE", add = "jitter", size = 0.3, add.params = list(size = 0.15)) +
  scale_color_manual(values = epitypes_colors) +
  scale_fill_manual(values = epitypes_colors) +
  theme_bw() +
  theme(
    title = element_text(size = 14, face = 'bold'),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # axis.title.x = element_text(size = 12),
    # axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    # axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    # axis.text.y = element_text(size = 14),
    # legend.position = "none",
    axis.line = element_line(),
    axis.line.y = element_blank(),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
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
  ylab("Mean Methylation Level (Beta Value)") +
  labs(fill = "Epitype") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(sprintf("Shore\n(%s CpGs)", prettyNum(nrow(bval_sub_Shore), big.mark=",", scientific=FALSE))) +
  
  # Adjust the size of p-values and statistical information
  stat_compare_means(size = 4.5, comparisons = my_comparisons,
                     label.y = c(max(df_FL$mean) + 0.01, max(df_FL$mean) + 0.05)) +  # Adjusted y position
  stat_compare_means(size = 4.5, label.x = 1, label.y = max(df$mean) + 0.12) +  # Adjusted label position
  
  coord_cartesian(clip = 'off') +
  theme(plot.margin = margin(0.4, 1, 0.4, 0.4, "cm")) +
  scale_x_discrete(name = "Type") +
  # scale_y_continuous(limits = c(min(df$mean) - 0.002, max(df$mean) + 0.06), breaks = seq(0, 1, by = 0.1))  # Adjusted y limits
  scale_y_continuous(limits = c(0.2, 0.8), breaks = seq(0, 1, by = 0.1))  # Adjusted y limits
# scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1))  # Adjusted y limits

# dev.off()

# Combine the plots using patchwork
# Combine Plots ----
combined_plot <- P1 + P2 + P5 + P4 + P3 + plot_layout(ncol = 5, widths = c(1, 1))

pdf("Merged_Figure.pdf", width = 14, height = 6)
print(combined_plot)
dev.off()

# Save the combined plot
# ggsave("Merged_Figure.pdf", combined_plot, height = 6, width = 16)

# # Merge Box Plot
# # ==============
# # Define CpGs of interest
# # CpG_islands
# # OpenSea
# # S_Shelf
# # S_Shore
# # N_Shelf
# # N_Shore
# 
# # All_Regions
# df_All_Regions <- bval_sub %>%
#   data.frame(check.names = FALSE) %>%
#   mutate(probe = row.names(.)) %>%
#   melt() %>%
#   left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE")], by = c("variable" = "SAMPLE_ID")) %>%
#   group_by(variable) %>%
#   summarize(mean = mean(value)) %>%
#   left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE")], by = c("variable" = "SAMPLE_ID")) %>%
#   filter(TYPE != "NA") %>%
#   mutate(
#     TYPE = factor(TYPE, levels = c("Normal", "FL", "t_DLBCL", "DLBCL", "naiBC", "t_naiBC", "gcBC", "t_PC", "memBC", "bm_PC")),
#     Relation_to_Island = "All_Regions"
#   ) %>% select(variable, TYPE, Relation_to_Island, mean)
# 
# # CpG_islands
# df_CpG_islands <- bval_sub[intersect(CpG_islands, row.names(bval_sub)),] %>%
#   data.frame(check.names = FALSE) %>%
#   mutate(probe = row.names(.)) %>%
#   melt() %>%
#   left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE")], by = c("variable" = "SAMPLE_ID")) %>%
#   group_by(variable) %>%
#   summarize(mean = mean(value)) %>%
#   left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE")], by = c("variable" = "SAMPLE_ID")) %>%
#   filter(TYPE != "NA") %>%
#   mutate(
#     TYPE = factor(TYPE, levels = c("Normal", "FL", "t_DLBCL", "DLBCL", "naiBC", "t_naiBC", "gcBC", "t_PC", "memBC", "bm_PC")),
#     Relation_to_Island = "CpG_Islands"
#   ) %>% select(variable, TYPE, Relation_to_Island, mean)
# 
# # OpenSea
# df_OpenSea <- bval_sub[intersect(OpenSea, row.names(bval_sub)),] %>%
#   data.frame(check.names = FALSE) %>%
#   mutate(probe = row.names(.)) %>%
#   melt() %>%
#   left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE")], by = c("variable" = "SAMPLE_ID")) %>%
#   group_by(variable) %>%
#   summarize(mean = mean(value)) %>%
#   left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE")], by = c("variable" = "SAMPLE_ID")) %>%
#   filter(TYPE != "NA") %>%
#   mutate(
#     TYPE = factor(TYPE, levels = c("Normal", "FL", "t_DLBCL", "DLBCL", "naiBC", "t_naiBC", "gcBC", "t_PC", "memBC", "bm_PC")),
#     Relation_to_Island = "OpenSea"
#   ) %>% select(variable, TYPE, Relation_to_Island, mean)
# 
# # Shelf
# Shelf <- c(S_Shelf, N_Shelf)
# 
# df_Shelf <- bval_sub[intersect(Shelf, row.names(bval_sub)),] %>%
#   data.frame(check.names = FALSE) %>%
#   mutate(probe = row.names(.)) %>%
#   melt() %>%
#   left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE")], by = c("variable" = "SAMPLE_ID")) %>%
#   group_by(variable) %>%
#   summarize(mean = mean(value)) %>%
#   left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE")], by = c("variable" = "SAMPLE_ID")) %>%
#   filter(TYPE != "NA") %>%
#   mutate(
#     TYPE = factor(TYPE, levels = c("Normal", "FL", "t_DLBCL", "DLBCL", "naiBC", "t_naiBC", "gcBC", "t_PC", "memBC", "bm_PC")),
#     Relation_to_Island = "Shelf"
#   ) %>% select(variable, TYPE, Relation_to_Island, mean)
# 
# # Shore
# Shore <- c(S_Shore, N_Shore)
# 
# df_Shore <- bval_sub[intersect(Shore, row.names(bval_sub)),] %>%
#   data.frame(check.names = FALSE) %>%
#   mutate(probe = row.names(.)) %>%
#   melt() %>%
#   left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE")], by = c("variable" = "SAMPLE_ID")) %>%
#   group_by(variable) %>%
#   summarize(mean = mean(value)) %>%
#   left_join(Sample_Annotations[,c("SAMPLE_ID", "TYPE")], by = c("variable" = "SAMPLE_ID")) %>%
#   filter(TYPE != "NA") %>%
#   mutate(
#     TYPE = factor(TYPE, levels = c("Normal", "FL", "t_DLBCL", "DLBCL", "naiBC", "t_naiBC", "gcBC", "t_PC", "memBC", "bm_PC")),
#     Relation_to_Island = "Shore"
#   ) %>% select(variable, TYPE, Relation_to_Island, mean)
# 
# # Bind data frames
# long_df <- rbind(df_All_Regions, df_CpG_islands, df_Shore, df_Shelf, df_OpenSea)
# dim(long_df)
# 
# # Define custom colors for improved visibility
# custom_colors <- c(
#   "Normal" = "#5086F3",
#   "FL" = "#59A14F",
#   "t_DLBCL" = "#F28E2B",
#   "DLBCL" = "#E15759",
#   "naiBC" = "#8C6D31",
#   "t_naiBC" = "#D1A1A6",
#   "gcBC" = "#E4B321",
#   "t_PC" = "#76B7B2",
#   "memBC" = "#FF9DA7",
#   "bm_PC" = "#6A4D8C"
# )
# 
# # Define custom order for sample types and Relation to Island
# custom_order <- c("Normal", "FL", "t_DLBCL", "DLBCL", "naiBC", "t_naiBC", "gcBC", "t_PC", "memBC", "bm_PC")
# long_df$TYPE <- factor(long_df$TYPE, levels = custom_order)
# 
# long_df$Relation_to_Island <- factor(long_df$Relation_to_Island, levels = c("All_Regions", "CpG_Islands", "Shore", "Shelf", "OpenSea"))
# 
# # Set up PDF output
# pdf("Mean-Methylation-All-CpGs.pdf", width = 35, height = 12)  # Adjusted width and height for better readability
# 
# # Generate grouped boxplot using custom colors
# ggplot(long_df, aes(x = Relation_to_Island, y = mean, fill = TYPE)) +
#   geom_boxplot(position = position_dodge(width = 0.9), size = 0.7) +
#   facet_grid(~Relation_to_Island, scales = "free_x", space = "free_x") +  # Facet to separate sample types
#   # facet_wrap(~Relation_to_Island, scales = "free_x", ncol = 4, nrow = 2) +
#   scale_fill_manual(values = custom_colors, name = "Sample Types") +
#   scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
#   ggtitle("Methylation Levels Across Sample Types and Genomic Regions") +
#   labs(x = "Relation to Island", y = "Mean Methylation (Beta value)") +
#   scale_x_discrete(labels = c(
#     "All_Regions" = paste("All_Regions\n(", format(length(unique(rownames(bval_sub))), big.mark = ","), " CpGs)", sep = ""),
#     "CpG_Islands" = paste("CpG_Islands\n(", format(length(unique(rownames(bval_sub[intersect(CpG_islands, row.names(bval_sub)),]))), big.mark = ","), " CpGs)", sep = ""),
#     "Shore" = paste("Shore\n(", format(length(unique(rownames(bval_sub[intersect(Shore, row.names(bval_sub)),]))), big.mark = ","), " CpGs)", sep = ""),
#     "Shelf" = paste("Shelf\n(", format(length(unique(rownames(bval_sub[intersect(Shelf, row.names(bval_sub)),]))), big.mark = ","), " CpGs)", sep = ""),
#     "OpenSea" = paste("OpenSea\n(", format(length(unique(rownames(bval_sub[intersect(OpenSea, row.names(bval_sub)),]))), big.mark = ","), " CpGs)", sep = "")
#   )) +
#   theme(
#     plot.title = element_text(size = 40, face = "bold", hjust = 0.5, margin = margin(b = 42)),
#     legend.text = element_text(size = 24),
#     legend.title = element_text(size = 34, face = "bold", margin = margin(b = 30)),
#     legend.key.size = unit(1.4, "cm"),
#     legend.margin = margin(l = 40),
#     axis.ticks.x = element_blank(),
#     axis.text.x = element_text(size = 28, face = "plain"),
#     axis.text.y = element_text(size = 28, face = "plain"),
#     axis.title.x = element_text(size = 34, face = "bold", margin = margin(t = 40)),
#     axis.title.y = element_text(size = 34, face = "bold", margin = margin(r = 40)),
#     plot.margin = margin(t = 40, r = 40, b = 40, l = 40),
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
