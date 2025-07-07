# R Packages
# ==========
library(dplyr)
library(InfiniumPurify)
library(stringr)
library(limma)
library(lumi)
source("/cluster/home/t110989uhn/apps/r_packages/modified_InfiniumPurify/modified_InfiniumPurify.R")

# Change directory
setwd("/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/01_InfiniumPurify/450k/modified_InfiniumPurify/Enmix_quantile_with_QCinfo_Harman_All/")

# Read in merged methylation data
bval <- readRDS("/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/450k/Enmix_quantile_with_QCinfo_Harman_All/rds/beta_combat.rds")
dim(bval)

# Read in sample annotations
Sample_Annotations <- read.table(file = "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/db/ANN/20231031_sample_annotations.tsv", sep = "\t", header = TRUE, na.strings = c("","NA")) %>%
	subset(TYPE!="TFL") %>%
	subset(TIME_POINT!="T2") %>%
	filter(SAMPLE_ID %in% colnames(bval))

# Read in sample sheet with batch information
targets <- read.csv(file = "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/db/SampleSheet/SampleSheet.csv", skip = 1, header = TRUE, na.strings = c("","NA"))
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
Benign_LN_sample_ann <- read.table(file = "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/Benign_LN/SampleSheet/Benign_LN_sample_annotations.tsv",
                                   sep = "\t",
                                   header = TRUE,
                                   na.strings=c("","NA")) %>%
  dplyr::select(SAMPLE_ID, TYPE)

Benign_LN_sample_ann <- Benign_LN_sample_ann %>% filter(SAMPLE_ID %in% colnames(bval))
Benign_LN_sample_ann$Batch <- "Benign_LN"

# Read in the NIH sample annotation
NIH_sample_ann <- read.table(file = "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/NIH/GSE237299/SampleSheet/GSE237299_sample_annotations.txt",
                             sep = "\t",
                             header = TRUE,
                             na.strings=c("","NA")) %>%
  dplyr::select(SAMPLE_ID, TYPE)

NIH_sample_ann$SAMPLE_ID <- gsub('-', '_', NIH_sample_ann$SAMPLE_ID)
NIH_sample_ann <- NIH_sample_ann %>% filter(SAMPLE_ID %in% colnames(bval))
NIH_sample_ann$Batch <- "GSE237299"

# Read in the EGA sample annotation
EGA_sample_ann <- read.table(file = "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/EGA/SampleSheet/Sample_Annotations_EGAD00010001974.txt",
                             sep = "\t",
                             header = TRUE,
                             na.strings=c("","NA")) %>%
  dplyr::select(SAMPLE_ID, TYPE)

EGA_sample_ann$SAMPLE_ID <- gsub('-', '_', EGA_sample_ann$SAMPLE_ID)
EGA_sample_ann$SAMPLE_ID <- gsub('\\.', '_', EGA_sample_ann$SAMPLE_ID)
EGA_sample_ann$Batch <- "EGAD00010001974"
EGA_sample_ann <- EGA_sample_ann %>% filter(SAMPLE_ID %in% colnames(bval))

# Read in blueprint sample annotation
Blueprint_sample_ann <- read.csv(file = "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/db/blueprint/Sample_Annotation_EGAS00001001196.txt",
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
GEO_sample_ann <- read.table(file = "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/GEO/GEO_Annotation.tsv",
                             sep = "\t",
                             header = TRUE,
                             na.strings=c("","NA")) %>%
                             dplyr::select(SAMPLE_ID, TYPE)
                                                          
GEO_sample_ann$Batch <- 'GSE255869'
GEO_sample_ann <- GEO_sample_ann %>% filter(SAMPLE_ID %in% colnames(bval))

# Read in clinical data
clinical_data <- read.csv(file = "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/db/CLN/20240116_clinical_data.csv", 
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
dim(bval)

gc()

# Write TYPE and Batch Information
df1 <- Sample_Annotations[c("SAMPLE_ID", "TYPE", "Batch")]
write.table(df1, file="TYPE_INFO.txt", sep="\t", quote = FALSE, row.names = FALSE)

# subset methylation data into normal and tumor
normal_samples <- Sample_Annotations %>% filter(TYPE == "RLN" | TYPE == "Benign_LN" | TYPE == "LN_NN_A" | TYPE == "LN_NN_B") %>% data.frame() %>% .$SAMPLE_ID
tumor_samples <- Sample_Annotations %>% filter(TYPE == "FL" | TYPE == "FL_A" | TYPE == "FL_B" | TYPE == "DLBCL" | TYPE == "DLBCL_EGA" | TYPE == "DLBCL_EN" | TYPE == "DLBCL_GC" | TYPE == "DLBCL_nonGC" | TYPE == "t_DLBCL") %>% data.frame() %>% .$SAMPLE_ID
blueprint_samples <- Sample_Annotations %>% filter(TYPE == "naiBC" | TYPE == "t_naiBC" | TYPE == "gcBC" | TYPE == "memBC" | TYPE == "bm_PC" | TYPE == "t_PC" | TYPE == "Normal_GCB") %>% data.frame() %>% .$SAMPLE_ID

normal <- bval[, normal_samples]
dim(normal)

tumor <- bval[, tumor_samples]
dim(tumor)

blueprint <- bval[, blueprint_samples]
dim(blueprint)

# estimate tumor purity
purity <- getPurity(tumor.data = tumor,
                    normal.data = normal)

gc()

# write purity table
purity_table <- as.data.frame(purity)
purity_table <- data.frame(SAMPLE_ID = row.names(purity_table), purity_table)
rownames(purity_table) <- NULL
purity_table$purity <- round(purity_table$purity, 2)
write.table(purity_table, file="purity_table.txt", sep="\t", quote = FALSE, row.names = FALSE)

# correct tumor methylome by tumor purity
#bval_purified = InfiniumPurify(tumor.data = tumor,
#                               normal.data = normal,
#                               purity = purity)

# correct tumor methylome by tumor purity
bval_purified = modified.InfiniumPurify(tumor.data = tumor,
                                        normal.data = normal,
                                        purity = purity)

dim(bval_purified)

gc()

# merge purified methylation data
bval_subset <- merge(x = bval_purified,
                     y = normal,
                     by = "row.names",
                     all.x = TRUE,
                     all.y = TRUE)

dim(bval_subset)

bval_subset_final <- bval_subset[,-1]
rownames(bval_subset_final) <- bval_subset[,1]

dim(bval_subset_final)

bval_sub <- merge(x = bval_subset_final,
                  y = blueprint,
                  by = "row.names",
                  all.x = TRUE,
                  all.y = TRUE)

dim(bval_sub)

bval_final <- bval_sub[,-1]
rownames(bval_final) <- bval_sub[,1]

dim(bval_final)

gc()

# Convert methylation Beta-value to M-value
mval_final <- beta2m(bval_final)
dim(mval_final)

# Write purified methylation data
saveRDS(bval_final, "beta_combat_purified.rds")
saveRDS(mval_final, "mval_combat_purified.rds")
