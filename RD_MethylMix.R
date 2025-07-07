# R Packages
# ==========
library(MethylMix)
library(dplyr)
library(stringr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(annotables)
library(WGCNA)
library(parallel)
library(doParallel)
library(ggplot2)
library(cowplot)
#library(enrichR)
library(forcats)
library(DGEobj.utils)
library(edgeR)

# Change directory
setwd("/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/10_MethylMix/84k/C2_vs_C1/")

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
bval <- readRDS("/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/01_InfiniumPurify/450k/modified_InfiniumPurify/Enmix_quantile_with_QCinfo_Harman_All/beta_combat_purified.rds")
dim(bval)

# Read in methylation SNPs information
# Total CpGs: 89,678
load("/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/RNASeq/MethylMix/MethylMix_2.34.0/MethylMix/data/SNPprobes.rda")
length(SNPprobes)

# Remove SNP probes by excluding them from the bval
bval <- bval[setdiff(rownames(bval), SNPprobes), ]
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

# Read in clinical data
clinical_data <- read.csv(file = "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/db/CLN/20231109_clinical_data.csv", 
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
Sample_Annotations <- bind_rows(Sample_Annotations, Benign_LN_sample_ann, NIH_sample_ann, EGA_sample_ann, Blueprint_sample_ann)
dim(Sample_Annotations)
table(Sample_Annotations$TYPE)
table(Sample_Annotations$Batch)

# Read in methylation clustering information
epitypes <- read.table(file = "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/02_Clustering/450k/modified_InfiniumPurify/Enmix_quantile_with_QCinfo_Harman_All/FL_epitypes.tsv",
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
dim(bval)

gc()

# Create list of groups
naiBC_Annotations <- Sample_Annotations %>% filter(TYPE == "naiBC") %>% data.frame()
gcBC_Annotations <- Sample_Annotations %>% filter(TYPE == "gcBC") %>% data.frame()
memBC_Annotations <- Sample_Annotations %>% filter(TYPE == "memBC") %>% data.frame()

# Subset methylation data group
naiBC_bval <- bval[, naiBC_Annotations$SAMPLE_ID]
gcBC_bval <- bval[, gcBC_Annotations$SAMPLE_ID]
memBC_bval <- bval[, memBC_Annotations$SAMPLE_ID]

dim(naiBC_bval)
dim(gcBC_bval)
dim(memBC_bval)

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
table(keep)

# Filter Out Similar CpG Probes
bval <- bval[keep, ]
dim(bval)

gc()

# exclude blueprint samples
Sample_Annotations <- Sample_Annotations %>% filter(TYPE == "FL" | TYPE == "FL_A" | TYPE == "FL_B" | TYPE == "DLBCL" | TYPE == "DLBCL_EGA" | TYPE == "DLBCL_EN" | TYPE == "Benign_LN" | TYPE == "LN_NN_A" | TYPE == "LN_NN_B" | TYPE == "RLN")

dim(Sample_Annotations)
table(Sample_Annotations$TYPE)
table(Sample_Annotations$epitypes)

# Calculate Samples by Type per Cluster
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'FL_A'] <- 'FL'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'FL_B'] <- 'FL'

Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'DLBCL_EGA'] <- 'DLBCL'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'DLBCL_EN'] <- 'DLBCL'

Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'Benign_LN'] <- 'Normal'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'LN_NN_A'] <- 'Normal'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'LN_NN_B'] <- 'Normal'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'RLN'] <- 'Normal'

Sample_Annotations %>% dplyr::count(epitypes, TYPE)

# subset of merged methylation data
bval <- bval[, Sample_Annotations$SAMPLE_ID]
dim(bval)

gc()

# Write TYPE and Batch Information
df1 <- Sample_Annotations[c("SAMPLE_ID", "TYPE","Batch")]
write.table(df1, file="TYPE_INFO_ALL.txt", sep="\t", quote = FALSE, row.names = FALSE)

gc()

# subset matrix based on most variable probes
set.seed(1234)
probe_num = nrow(bval)
#probe_num = "10000"

# subset random probes
# bval_sub <- bval[sample(nrow(bval), probe_num), ]
# dim(bval_sub)
# bval_sub[1:5, 1:5]

# All probes
bval_sub <- bval

# Create list of groups
# FL vs. Normal
# C1 vs. Normal
# C2 vs. Normal

#DLBCL <- Sample_Annotations %>% filter(TYPE == "DLBCL") %>% data.frame() %>% .$SAMPLE_ID # n=78
FL <- Sample_Annotations %>% filter(TYPE == "FL") %>% data.frame() %>% .$SAMPLE_ID # n=338
C1 <- Sample_Annotations %>% filter(TYPE == "FL" & epitypes == "C1") %>% data.frame() %>% .$SAMPLE_ID # n=222
C2 <- Sample_Annotations %>% filter(TYPE == "FL" & epitypes == "C2") %>% data.frame() %>% .$SAMPLE_ID # n=83
Normal <- Sample_Annotations %>% filter(TYPE == "Normal") %>% data.frame() %>% .$SAMPLE_ID # n=20

# Subset methylation data group
bval_FL <- bval_sub[, FL]
dim(bval_FL)
bval_FL[1:5, 1:5]

bval_C1 <- bval_sub[, C1]
dim(bval_C1)
bval_C1[1:5, 1:5]

bval_C2 <- bval_sub[, C2]
dim(bval_C2)
bval_C2[1:5, 1:5]

bval_Normal <- bval_sub[Normal]
dim(bval_Normal)
bval_Normal[1:5, 1:5]

# Load RNA-seq data
# =================
# FL RNA-Seq, n=282

mRNA <- read.table("/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/db/RNA_Seq/2022-10-25_df_counts_adj.txt", sep = " ", header = TRUE)
passing.RNAseq.QC <- read.table("/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/db/RNA_Seq/metadata_passed_290.txt",
                                sep = "\t",
                                header = TRUE) %>% subset(phenotype == "FL")
                                
length(unique(passing.RNAseq.QC$sample_id))

mRNA <- mRNA[,passing.RNAseq.QC$id]
dim(mRNA)
mRNA[1:5, 1:5]

# Convert Gene IDs to Gene Symbol
gene_info <- read.csv("/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/db/ensembl/annotation/gene_info.csv", sep = ",", header = TRUE)
gene_coords <- gene_info[,-1]
rownames(gene_coords) <- gene_info[,1]

# Merge Gene Info with Read Count Matrix
df1 <- merge(x = mRNA,
             y = gene_coords,
             by = 'row.names',
             all.y = TRUE)
             
colnames(df1)[1] ="GENE_ID"
df1 <- df1 %>% relocate(SYMBOL, GENE_TYPE, MEAN, MEDIAN, LONGEST_ISOFORM, MERGED, .after = GENE_ID)
df1 <- df1[!is.na(df1$SYMBOL),]
df1 <- df1[!is.na(df1$GENE_TYPE),]
df1 <- na.omit(df1)

mRNA <- df1[,-1]
rownames(mRNA) <- df1[,1]
remove(df1)

# Calculate TPM
# https://rdrr.io/cran/DGEobj.utils/man/convertCounts.html

# Remove rows consistently have zero or very low counts
keep <- filterByExpr(as.matrix(mRNA[7:288]))
table(keep)
mRNA_RC <- mRNA[keep,]
dim(mRNA_RC)

mRNA_TPM <- convertCounts(as.matrix(mRNA_RC[7:288]),
                             #unit = "CPM",
                             unit = "TPM",
                             geneLength = mRNA_RC$MERGED,
                             log = FALSE,
                             normalize  = "none"
)

mRNA_TPM <- round(mRNA_TPM, digits = 2)
mRNA_TPM <- na.omit(mRNA_TPM)
mRNA_TPM[1:5, 1:5]

# Normalize TPM values to z-scores
mRNA_TPM <- t(scale(t(mRNA_TPM)))
dim(mRNA_TPM)
mRNA_TPM[1:5, 1:5]

# Merge Gene Info with TPM
df1 <- merge(x = mRNA_TPM,
             y = gene_coords,
             by = 'row.names',
             all.x = TRUE)

colnames(df1)[1] ="GENE_ID"
df1 <- df1 %>% relocate(SYMBOL, GENE_TYPE, MEAN, MEDIAN, LONGEST_ISOFORM, MERGED, .after = GENE_ID)
df1 <- df1[!is.na(df1$SYMBOL),]
df1 <- df1[!is.na(df1$GENE_TYPE),]
df1 <- na.omit(df1)

mRNA_TPM <- df1[,-1]
rownames(mRNA_TPM) <- df1[,1]
remove(df1)
mRNA_TPM[1:5, 1:10]

# Filter Matrix for Protein Coding Genes
mRNA_TPM <- mRNA_TPM %>% filter(GENE_TYPE == "protein_coding")
dim(mRNA_TPM)

# Identify duplicate genes and remove redundancy
n_occur <- data.frame(table(mRNA_TPM$SYMBOL))
n_occur[n_occur$Freq > 1,]
n_occur_df <- mRNA_TPM[mRNA_TPM$SYMBOL %in% n_occur$Var1[n_occur$Freq > 1],]
n_occur_list <- c("ENSG00000168255", "ENSG00000285053") # removed shortest isoforms
mRNA_TPM <- mRNA_TPM[!(row.names(mRNA_TPM) %in% n_occur_list),]

mRNA_TPM_Final <- mRNA_TPM[,!names(mRNA_TPM) %in% c("GENE_TYPE", "MEAN", "MEDIAN", "LONGEST_ISOFORM", "MERGED")]
colnames(mRNA_TPM_Final)[colnames(mRNA_TPM_Final) == "SYMBOL"] = "Gene"

dim(mRNA_TPM)
dim(mRNA_TPM_Final)

#write.table(mRNA_TPM_Final, file = "mRNA_TPM_Final.tsv", quote = FALSE, sep = "\t", na = "NA", row.names = FALSE)

rownames(mRNA_TPM_Final) <- mRNA_TPM_Final$Gene
mRNA_TPM_Final <- as.matrix(mRNA_TPM_Final[2:283])
dim(mRNA_TPM_Final)
mRNA_TPM_Final[1:5, 1:5]

# subset mRNA dataset
FL_RNASeq <- intersect(FL, colnames(mRNA_TPM_Final)) # n = 73 (21.59%)
mRNA_TPM_Final_FL <- mRNA_TPM_Final[, FL_RNASeq]
dim(mRNA_TPM_Final_FL)
mRNA_TPM_Final_FL[1:5, 1:5]

C1_RNASeq <- intersect(C1, colnames(mRNA_TPM_Final)) # n = 39 (11.53%) (C1: 17.56%)
mRNA_TPM_Final_C1 <- mRNA_TPM_Final[, C1_RNASeq]
dim(mRNA_TPM_Final_C1)
mRNA_TPM_Final_C1[1:5, 1:5]

C2_RNASeq <- intersect(C2, colnames(mRNA_TPM_Final)) # n = 26 (7.69%) (C2: 31.32%)
mRNA_TPM_Final_C2 <- mRNA_TPM_Final[, C2_RNASeq]
dim(mRNA_TPM_Final_C2)
mRNA_TPM_Final_C2[1:5, 1:5]

# subset methylation datasets for RNA-Seq samples
bval_FL <- bval_sub[, FL_RNASeq] # n = 73 (21.59%)
dim(bval_FL)
bval_FL[1:5, 1:5]

bval_C1 <- bval_sub[, C1_RNASeq] # n = 39 (11.53%) (C1: 17.56%)
dim(bval_C1)
bval_C1[1:5, 1:5]

bval_C2 <- bval_sub[, C2_RNASeq] # n = 26 (7.69%) (C2: 31.32%)
dim(bval_C2)
bval_C2[1:5, 1:5]

dim(bval_Normal) # n = 20 (100%)
bval_Normal[1:5, 1:5]

# Clustering probes to genes for methylation data
# -----------------------------------------------
cl <- makeCluster(30)
registerDoParallel(cl)
#res <- ClusterProbes(as.matrix(bval_FL), as.matrix(bval_Normal))
#res <- ClusterProbes(as.matrix(bval_C1), as.matrix(bval_Normal))
#res <- ClusterProbes(as.matrix(bval_C2), as.matrix(bval_Normal))
res <- ClusterProbes(as.matrix(bval_C2), as.matrix(bval_C1))
stopCluster(cl)

# Identify column names that end with "T1" and remove the last character from them
METcancer = res[[1]]
colnames(METcancer)[grep("T1$", colnames(METcancer))] <- substr(colnames(METcancer)[grep("T1$", colnames(METcancer))], 1, nchar(colnames(METcancer)[grep("T1$", colnames(METcancer))]) - 1)

METnormal = res[[2]]
colnames(METnormal)[grep("T1$", colnames(METnormal))] <- substr(colnames(METnormal)[grep("T1$", colnames(METnormal))], 1, nchar(colnames(METnormal)[grep("T1$", colnames(METnormal))]) - 1)

GEcancer = as.matrix(mRNA_TPM_Final_C2)
colnames(GEcancer)[grep("T1$", colnames(GEcancer))] <- substr(colnames(GEcancer)[grep("T1$", colnames(GEcancer))], 1, nchar(colnames(GEcancer)[grep("T1$", colnames(GEcancer))]) - 1)

saveRDS(METcancer, "METcancer.rds")
saveRDS(METnormal, "METnormal.rds")
saveRDS(GEcancer, "GEcancer.rds")

# Running MethylMix
cl <- makeCluster(30)
registerDoParallel(cl)
MethylMixResults <- MethylMix(METcancer, GEcancer, METnormal)
stopCluster(cl)

saveRDS(MethylMixResults, "MethylMixResults.rds")
