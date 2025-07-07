# R Packages
# ==========
library(minfi)
library(minfiData)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(limma)
library(maxprobes)
library(lumi)
library(dplyr)
library(stringr)
library(Harman)
library(ENmix)
library(wateRmelon)

# Read UHN IDAT Files
# ===================
# set the working directory
UHN_baseDir <- "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/db/00_FL_Samples/"
list.files(UHN_baseDir)

# Read a sample sheet
UHN_targets <- read.metharray.sheet(UHN_baseDir)

# Remove 2 machine control samples
UHN_targets <- UHN_targets %>% filter(Sample_Name != 'HCT116_DKO_methylated')

# Remove 2 samples for clinical reason
UHN_targets <- UHN_targets %>% filter(Sample_Name != 'LY_FL_311_T1' & Sample_Name != 'LY_FL_159_T1')

# Rename replicate of sample "LY_FL_159_T1"
UHN_targets$Sample_Name <- str_replace(UHN_targets$Sample_Name, "LY_FL_159_T1_rep", "LY_FL_159_T1")

# Remove 1 TFL sample
UHN_targets <- UHN_targets %>% filter(Sample_Name != 'LY_FL_127_T2')

# Remove 2 suspicious repeated samples
UHN_targets <- UHN_targets %>% filter(Sample_Name != 'LY_FL_1156_T1' & Sample_Name != 'LY_FL_571_T1')

# Read in the raw data from the IDAT files
# RGset: RGChannelSet, 1,051,539 features
UHN_RGset <- read.metharray.exp(targets = UHN_targets, force=TRUE, extended=TRUE)
dim(UHN_RGset)

# Read Benign_LN IDAT Files
# =========================
# set the working directory
Benign_LN_baseDir <- "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/Benign_LN/IDAT/"
list.files(Benign_LN_baseDir)

# Read a sample sheet
Benign_LN_targets <- read.metharray.sheet(Benign_LN_baseDir)

# Read in the raw data from the IDAT files
# RGset: RGChannelSet, 1,051,943 features
Benign_LN_RGset <- read.metharray.exp(targets = Benign_LN_targets, force=TRUE, extended=TRUE)
dim(Benign_LN_RGset)

# Read NIH IDAT Files
# ===================
# set the working directory
NIH_baseDir <- "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/NIH/GSE237299/IDAT/"
list.files(NIH_baseDir)

# Read a sample sheet
NIH_targets <- read.metharray.sheet(NIH_baseDir)
NIH_targets$Sample_Name <- gsub('-', '_', NIH_targets$Sample_Name)

# Read in the raw data from the IDAT files
# RGset: RGChannelSet, 1,051,943 features
NIH_RGset <- read.metharray.exp(targets = NIH_targets, force=TRUE, extended=TRUE)
dim(NIH_RGset)

# Read EGA 450k IDAT Files
# ========================
# set the working directory
EGA_450k_baseDir <- "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/EGA/IDAT/EGAD00010001974/450k/"
list.files(EGA_450k_baseDir)

# Read a sample sheet
EGA_450k_targets <- read.metharray.sheet(EGA_450k_baseDir)
EGA_450k_targets$Sample_Name <- gsub('-', '_', EGA_450k_targets$Sample_Name)
EGA_450k_targets$Sample_Name <- gsub('\\.', '_', EGA_450k_targets$Sample_Name)

# Read in the raw data from the IDAT files
# RGset: RGChannelSet, 622,399 features
EGA_450k_RGset <- read.metharray.exp(targets = EGA_450k_targets, force=TRUE, extended=TRUE)
dim(EGA_450k_RGset)

# Read EGA 850k IDAT Files
# ========================
# set the working directory
EGA_850k_baseDir <- "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/EGA/IDAT/EGAD00010001974/850k/"
list.files(EGA_850k_baseDir)

# Read a sample sheet
EGA_850k_targets <- read.metharray.sheet(EGA_850k_baseDir)
EGA_850k_targets$Sample_Name <- gsub('-', '_', EGA_850k_targets$Sample_Name)
EGA_850k_targets$Sample_Name <- gsub('\\.', '_', EGA_850k_targets$Sample_Name)

# Read in the raw data from the IDAT files
# RGset: RGChannelSet, 1,051,815 features
EGA_850k_RGset <- read.metharray.exp(targets = EGA_850k_targets, force=TRUE, extended=TRUE)
dim(EGA_850k_RGset)

# Read Blueprint IDAT Files
# =========================
# set the working directory
Blueprint_baseDir <- "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/db/00_Blueprint_Samples/"
list.files(Blueprint_baseDir)

# Read a sample sheet
Blueprint_targets <- read.metharray.sheet(Blueprint_baseDir)
Blueprint_targets$Sample_Name <- gsub('-', '_', Blueprint_targets$Sample_Name)
Blueprint_targets$Batch <- 'Blueprint'

# Read in the raw data from the IDAT files
# RGset: RGChannelSet, 622,399 features
Blueprint_RGset <- read.metharray.exp(targets = Blueprint_targets, force=TRUE, extended=TRUE)
dim(Blueprint_RGset)

# Read GEO IDAT Files
# ===================
# set the working directory
GEO_baseDir <- "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/GEO/IDAT/"
list.files(GEO_baseDir)

# Read a sample sheet
GEO_targets <- read.metharray.sheet(GEO_baseDir)

# Read in the raw data from the IDAT files
# RGset: RGChannelSet, 1,051,943 features
GEO_RGset <- read.metharray.exp(targets = GEO_targets, force=TRUE, extended=TRUE)
dim(GEO_RGset)

# Combine UHN, Benign_LN, NIH RGset, EGA_RGset and Blueprint_RGset
# RGset: RGChannelSet, 574,981 features
RGset_R1 <- combineArrays(UHN_RGset, Benign_LN_RGset, outType = "IlluminaHumanMethylationEPIC")
dim(RGset_R1)

RGset_R2 <- combineArrays(RGset_R1, NIH_RGset, outType = "IlluminaHumanMethylationEPIC")
dim(RGset_R2)

RGset_R3 <- combineArrays(RGset_R2, EGA_850k_RGset, outType = "IlluminaHumanMethylationEPIC")
dim(RGset_R3)

RGset_R4 <- combineArrays(RGset_R3, GEO_RGset, outType = "IlluminaHumanMethylationEPIC")
dim(RGset_R4)

RGset_R5 <- combineArrays(RGset_R4, EGA_450k_RGset, outType = "IlluminaHumanMethylation450k")
dim(RGset_R5)

RGset <- combineArrays(RGset_R5, Blueprint_RGset, outType = "IlluminaHumanMethylation450k")
dim(RGset)

remove(RGset_R1, RGset_R2, RGset_R3, RGset_R4, RGset_R5, UHN_RGset, Benign_LN_RGset, NIH_RGset, EGA_850k_RGset, GEO_RGset, EGA_450k_RGset, Blueprint_RGset)
gc()

# Combine sample sheets
targets <- bind_rows(UHN_targets, Benign_LN_targets, NIH_targets, EGA_850k_targets, GEO_targets, EGA_450k_targets, Blueprint_targets)
dim(targets)

# Give the samples descriptive names
targets$ID <- gsub(" ", "", paste(targets$Sample_Name,".",targets$Array))
sampleNames(RGset) <- targets$ID
RGset
dim(RGset)
colnames(RGset)

# See which annotation package being used by
annotation(RGset)

# Read in the UHN sample annotation
UHN_sample_ann <- read.table(file = "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/db/ANN/20231031_sample_annotations.tsv",
                             sep = "\t",
                             header = TRUE,
                             na.strings=c("","NA"))

# Read in the Benign_LN sample annotation
Benign_LN_sample_ann <- read.table(file = "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/Benign_LN/SampleSheet/Benign_LN_sample_annotations.tsv",
                                   sep = "\t",
                                   header = TRUE,
                                   na.strings=c("","NA")) %>%
                                   dplyr::select(SAMPLE_ID, TYPE)

# Read in the NIH sample annotation
NIH_sample_ann <- read.table(file = "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/NIH/GSE237299/SampleSheet/GSE237299_sample_annotations.txt",
                             sep = "\t",
                             header = TRUE,
                             na.strings=c("","NA")) %>% 
                             dplyr::select(SAMPLE_ID, TYPE)

NIH_sample_ann$SAMPLE_ID <- gsub('-', '_', NIH_sample_ann$SAMPLE_ID)

# Read in the EGA sample annotation
EGA_sample_ann <- read.table(file = "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/EGA/SampleSheet/Sample_Annotations_EGAD00010001974.txt",
                             sep = "\t",
                             header = TRUE,
                             na.strings=c("","NA")) %>%
                             dplyr::select(SAMPLE_ID, TYPE)

EGA_sample_ann$SAMPLE_ID <- gsub('-', '_', EGA_sample_ann$SAMPLE_ID)
EGA_sample_ann$SAMPLE_ID <- gsub('\\.', '_', EGA_sample_ann$SAMPLE_ID)

# Read in the Blueprint sample annotation
Blueprint_sample_ann <- read.table(file = "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/db/blueprint/Sample_Annotation_EGAS00001001196.txt",
                                   sep = "\t",
                                   header = TRUE,
                                   na.strings=c("","NA")) %>%
                                   dplyr::select(Sample_Name, Source_Name, CELL_TYPE, TISSUE_TYPE, DONOR_SEX)

colnames(Blueprint_sample_ann) <- c("SAMPLE_ID", "LY_FL_ID", "TYPE", "SITE_BIOPSY", "SEX")
Blueprint_sample_ann$SAMPLE_ID <- gsub('-', '_', Blueprint_sample_ann$SAMPLE_ID)
Blueprint_sample_ann$TYPE <- gsub('-', '_', Blueprint_sample_ann$TYPE)
Blueprint_sample_ann$SEX <- gsub('Male', 'M', Blueprint_sample_ann$SEX)
Blueprint_sample_ann$SEX <- gsub('Female', 'F', Blueprint_sample_ann$SEX)

# Read in the GEO sample annotation
GEO_sample_ann <- read.table(file = "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/GEO/GEO_Annotation.tsv",
                                   sep = "\t",
                                   header = TRUE,
                                   na.strings=c("","NA")) %>%
                                   dplyr::select(SAMPLE_ID, TYPE)
                                   
# Combine sample annotations
sample_ann <- bind_rows(UHN_sample_ann, Benign_LN_sample_ann, NIH_sample_ann, EGA_sample_ann, GEO_sample_ann, Blueprint_sample_ann)

# Left outer join sample sheet and sample annotation
targets <- merge(x = targets,
                 y = sample_ann,
                 by.x = "Sample_Name",
                 by.y = "SAMPLE_ID",
                 all.x = TRUE)
                 
dim(targets)

gc()

## Set Working Current Directory
## =============================
setwd("/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/450k/Enmix_quantile_with_QCinfo_Harman_All/")

# Filter Samples By p-value
# =========================
# Calculate the detection p-values
detP <- detectionP(RGset)
head(detP)

# Examine mean detection p-values across all samples to identify any failed samples
pval <- (colMeans(detP)) %>% as.data.frame()
pval <- cbind(rownames(pval), data.frame(pval, row.names=NULL))
colnames(pval) <- c("Sample_Name", "mean_pval")
write.table(pval, file = "pval.tsv", sep="\t", quote = FALSE, row.names = FALSE)

# Examine mean detection p-values across all samples to identify any failed samples
pdf("pval.pdf", height = 40)
par(mar=c(6,6,2,2)) # Increase margin size
barplot(pval$mean_pval,
        names=pval$Sample_Name,
        col='#454545',
        xlab='Mean detection p-values',
        horiz = TRUE,
        las=1,
        cex.names=0.4,
        xlim = c(0, 0.6),
        border = NA
)
abline(v = 0.01, col = "red")
dev.off()

# remove poor quality samples
keep <- colMeans(detP) < 0.01
table(keep)

RGset <- RGset[,keep]
dim(RGset)

# remove poor quality samples from detection p-value table
detP <- detP[,keep]
dim(detP)

gc()

# Re-Examine mean detection p-values across all samples to identify any failed samples
pval <- (colMeans(detP)) %>% as.data.frame()
pval <- cbind(rownames(pval), data.frame(pval, row.names=NULL))
colnames(pval) <- c("Sample_Name", "mean_pval")
write.table(pval, file = "pval_updated.tsv", sep="\t", quote = FALSE, row.names = FALSE)

# Re-Examine mean detection p-values across all samples to identify any failed samples
pdf("pval_updated.pdf", height = 40)
par(mar=c(6,6,2,2)) # Increase margin size
barplot(pval$mean_pval,
        names=pval$Sample_Name,
        col='#454545',
        xlab='Mean detection p-values',
        horiz = TRUE,
        las=1,
        cex.names=0.4,
        xlim = c(0, 0.6),
        border = NA
)
abline(v = 0.01, col = "red")
dev.off()

# Filter known repeated samples by a # of significant probes
# ==========================================================
# Examine significant number of probes in each samples
sig_probes <- (colSums(detP <= 0.01)) %>% as.data.frame()
sig_probes <- cbind(rownames(sig_probes), data.frame(sig_probes, row.names=NULL))
colnames(sig_probes) <- c("Sample_Name", "significant_probes")
write.table(sig_probes, file = "significant_probes.tsv", sep="\t", quote = FALSE, row.names = FALSE)

# Examine significant number of probes in each samples # Total probes: 452,453
pdf("significant_probes.pdf", height = 40)
par(mar=c(6,6,2,2)) # Increase margin size
barplot(sig_probes$significant_probes,
        names=sig_probes$Sample_Name,
        col='#454545',
        xlab = "Significant number of probes",
        horiz = TRUE,
        las=1,
        cex.names=0.4,
        xlim = c(0,485512),
        border = NA
)
abline(v = 407208, col = "red")
dev.off()

# Identify known repeated sample with minimum # of significant probes
sig_probes_t <- data.frame(sig_probes[,-1], row.names=sig_probes[,1])
colnames(sig_probes_t) <- c("significant_probes")

# Compare and find known repeated sample names with minimum # of significant probes
sample_1 <- c("LY_FL_529_T1.R01C01",
              "LY_FL_527_T1.R01C01",
              "LY_FL_525_T1.R03C01",
              "LY_FL_524_T1.R02C01",
              "LY_FL_523_T1.R01C01",
              "LY_FL_488_T1.R06C01",
              "LY_FL_479_T1.R02C01",
              "LY_FL_158_T1.R01C01",
              "LY_FL_498_T1.R06C01",
              "LY_FL_535_T1.R06C01")

sample_2 <- c("LY_FL_529_T1.R05C01",
              "LY_FL_527_T1.R04C01",
              "LY_FL_525_T1.R02C01",
              "LY_FL_524_T1.R03C01",
              "LY_FL_523_T1.R04C01",
              "LY_FL_488_T1.R04C01",
              "LY_FL_479_T1.R07C01",
              "LY_FL_158_T1.R08C01",
              "LY_FL_498_T1.R05C01",
              "LY_FL_535_T1.R08C01")

rep_samples <- data.frame(sample_1, sample_2)

sample_min <- c()

for (i in 1:10) {
  df <- subset(sig_probes_t, rownames(sig_probes_t) %in% c(rep_samples$sample_1[i], rep_samples$sample_2[i]))
  sample_min <- append(rownames(df)[which(df == min(df), arr.ind = TRUE)[ , 1]], sample_min)
}

sample_min
#[1] "LY_FL_535_T1.R08C01" "LY_FL_498_T1.R06C01" "LY_FL_158_T1.R01C01" "LY_FL_479_T1.R02C01" "LY_FL_488_T1.R06C01"
#[6] "LY_FL_523_T1.R01C01" "LY_FL_524_T1.R02C01" "LY_FL_525_T1.R03C01" "LY_FL_527_T1.R01C01" "LY_FL_529_T1.R01C01"

# remove repeated 10 samples from mean detection p-values
pval <- pval[!(pval$Sample_Name %in% sample_min),]

# filter sample sheet with final remaining samples
targets <- targets[(targets$ID %in% pval$Sample_Name),]
dim(targets)

# subset "UHN_targets" sample sheet
UHN_targets_final <- targets %>% filter(grepl('Genome-Quebec|PM-Genomics-Centre|PM-OICR-TGL', Batch))
dim(UHN_targets_final)

# Read in clinical data
clinical_data <- read.csv(file = "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/db/CLN/20240116_clinical_data.csv",
                          header = TRUE, na.strings = c("","NA")) %>% 
                          dplyr::select(LY_FL_ID, ANN_ARBOR_STAGE, SEX, T_14_18)

# Left outer join sample sheet and clinical data
UHN_targets_final <- merge(x = UHN_targets_final,
                           y = clinical_data,
                           by.x = "LY_FL_ID",
                           by.y = "LY_FL_ID",
                           all.x = TRUE)

# subset "Benign_LN_targets" sample sheet
Benign_LN_targets_final <- targets %>% filter(grepl('Benign_LN', Batch))
dim(Benign_LN_targets_final)

# subset "NIH" sample sheet
NIH_targets_final <- targets %>% filter(grepl('GSE237299', Batch))
dim(NIH_targets_final)

# Read in NIH clinical data
NIH_clinical_data <- read.csv(file = "/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/NIH/GSE237299/SampleSheet/GSE237299_sample_annotations.txt", 
                              sep = "\t", header = TRUE, na.strings = c("","NA")) %>%
                              dplyr::select(SAMPLE_ID, SEX)

NIH_clinical_data$SAMPLE_ID <- gsub('-', '_', NIH_clinical_data$SAMPLE_ID)

# Left outer join NIH sample sheet and clinical data
NIH_targets_final <- NIH_targets_final[1:(length(NIH_targets_final)-1)]

NIH_targets_final <- merge(x = NIH_targets_final,
                           y = NIH_clinical_data,
                           by.x = "Sample_Name",
                           by.y = "SAMPLE_ID",
                           all.x = TRUE)

# subset "EGA" sample sheet
EGA_targets_final <- targets %>% filter(grepl('EGAD00010001974', Batch))
dim(EGA_targets_final)

EGA_850k_targets_final <- EGA_targets_final %>% filter(grepl('850k', Basename))
EGA_450k_targets_final <- EGA_targets_final %>% filter(grepl('450k', Basename))

# subset "Blueprint_targets" sample sheet
Blueprint_targets_final <- targets %>% filter(grepl('Blueprint', Batch))
dim(Blueprint_targets_final)

# subset "GEO_targets" sample sheet
GEO_targets_final <- targets %>% filter(grepl('GSE255869', Batch))
dim(GEO_targets_final)

# Re-generate Final RGset, detP, pval, sig_probes with final unique samples
# =========================================================================
# Remove previous
remove(RGset, detP)
gc()

# Read in the raw data from the IDAT files
# RGset: RGChannelSet, 1,051,539 features
UHN_RGset <- read.metharray.exp(targets = UHN_targets_final, force=TRUE, extended=TRUE)
dim(UHN_RGset)

# Read in the raw data from the IDAT files
# RGset: RGChannelSet, 1,051,943 features
Benign_LN_RGset <- read.metharray.exp(targets = Benign_LN_targets_final, force=TRUE, extended=TRUE)
dim(Benign_LN_RGset)

# Read in the raw data from the IDAT files
# RGset: RGChannelSet, 1,051,943 features
NIH_RGset <- read.metharray.exp(targets = NIH_targets_final, force=TRUE, extended=TRUE)
dim(NIH_RGset)

# Read in the raw data from the IDAT files
# RGset: RGChannelSet, 1,051,815 features
EGA_850k_RGset <- read.metharray.exp(targets = EGA_850k_targets_final, force=TRUE, extended=TRUE)
dim(EGA_850k_RGset)

# Read in the raw data from the IDAT files
# RGset: RGChannelSet, 622,399 features
EGA_450k_RGset <- read.metharray.exp(targets = EGA_450k_targets_final, force=TRUE, extended=TRUE)
dim(EGA_450k_RGset)

# Read in the raw data from the IDAT files
# RGset: RGChannelSet, 622,399 features
Blueprint_RGset <- read.metharray.exp(targets = Blueprint_targets_final, force=TRUE, extended=TRUE)
dim(Blueprint_RGset)

# Read in the raw data from the IDAT files
# RGset: RGChannelSet, 1,051,943 features
GEO_RGset <- read.metharray.exp(targets = GEO_targets_final, force=TRUE, extended=TRUE)
dim(GEO_RGset)

# Combine UHN, Benign_LN, NIH RGset, EGA_RGset and Blueprint_RGset
# RGset: RGChannelSet, 574,981 features
RGset_R1 <- combineArrays(UHN_RGset, Benign_LN_RGset, outType = "IlluminaHumanMethylationEPIC")
dim(RGset_R1)

RGset_R2 <- combineArrays(RGset_R1, NIH_RGset, outType = "IlluminaHumanMethylationEPIC")
dim(RGset_R2)

RGset_R3 <- combineArrays(RGset_R2, EGA_850k_RGset, outType = "IlluminaHumanMethylationEPIC")
dim(RGset_R3)

RGset_R4 <- combineArrays(RGset_R3, GEO_RGset, outType = "IlluminaHumanMethylationEPIC")
dim(RGset_R4)

RGset_R5 <- combineArrays(RGset_R4, EGA_450k_RGset, outType = "IlluminaHumanMethylation450k")
dim(RGset_R5)

RGset <- combineArrays(RGset_R5, Blueprint_RGset, outType = "IlluminaHumanMethylation450k")
dim(RGset)

remove(RGset_R1, RGset_R2, RGset_R3, RGset_R4, RGset_R5, UHN_RGset, Benign_LN_RGset, NIH_RGset, EGA_850k_RGset, GEO_RGset, EGA_450k_RGset, Blueprint_RGset)
gc()

# Combine sample sheets
targets <- bind_rows(UHN_targets_final, Benign_LN_targets_final, NIH_targets_final, EGA_850k_targets_final, GEO_targets_final, EGA_450k_targets_final, Blueprint_targets_final)
dim(targets)

# Give the samples original and unique descriptive names
sampleNames(RGset) <- targets$Sample_Name
RGset

# Calculate the detection p-values
detP <- detectionP(RGset)
head(detP)

## Set Working Current Directory
## =============================
setwd("/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/450k/Enmix_quantile_with_QCinfo_Harman_All/")

# Examine mean detection p-values across all samples to identify any failed samples
pval <- (colMeans(detP)) %>% as.data.frame()
pval <- cbind(rownames(pval), data.frame(pval, row.names=NULL))
colnames(pval) <- c("Sample_Name", "mean_pval")
write.table(pval, file = "pval_final.tsv", sep="\t", quote = FALSE, row.names = FALSE)

# Examine mean detection p-values across all samples to identify any failed samples
pdf("pval_final.pdf", height = 40)
par(mar=c(6,6,2,2)) # Increase margin size
pal <- c("#5086F3", "#7f8992", "#D02060", "#fc1417", "#f26ded", "#D87087", "#9A3B8A", "#6dcc04", "#02bdaa", "#ebc602", "#0303a8", "#FF7F00", "#B15928", "#b5b80b", "#2480b9", "#2E2EAA", "#7618dc", "#6E0A0A", "#6e988c", "#383838")
barplot(pval$mean_pval,
        names=pval$Sample_Name,
        col=pal[factor(targets$TYPE)],
        xlab='Mean detection p-values',
        horiz = TRUE,
        las=1,
        cex.names=0.4,
        xlim = c(0, 0.6),
        border = NA
)
abline(v = 0.01, col = "red")
dev.off()

# Examine significant number of probes in each samples
sig_probes <- (colSums(detP <= 0.01)) %>% as.data.frame()
sig_probes <- cbind(rownames(sig_probes), data.frame(sig_probes, row.names=NULL))
colnames(sig_probes) <- c("Sample_Name", "significant_probes")
write.table(sig_probes, file = "significant_probes_final.tsv", sep="\t", quote = FALSE, row.names = FALSE)

# Examine significant number of probes in each samples
pdf("significant_probes_final.pdf", height = 40)
par(mar=c(6,6,2,2)) # Increase margin size
pal <- c("#5086F3", "#7f8992", "#D02060", "#fc1417", "#f26ded", "#D87087", "#9A3B8A", "#6dcc04", "#02bdaa", "#ebc602", "#0303a8", "#FF7F00", "#B15928", "#b5b80b", "#2480b9", "#2E2EAA", "#7618dc", "#6E0A0A", "#6e988c", "#383838")
barplot(sig_probes$significant_probes,
        names=sig_probes$Sample_Name,
        col=pal[factor(targets$TYPE)],
        xlab = "Significant number of probes",
        horiz = TRUE,
        las=1,
        cex.names=0.4,
        xlim = c(0,485512),
        border = NA
)
abline(v = 407208, col = "red")
dev.off()

# Also take a look, If any samples got significant probes <= 407208
filter(sig_probes, significant_probes <= 407208)

# Writing data sets
# Save RGset
saveRDS(RGset, file = "RGset.rds")

## QC information
## ==============
qc <- QCinfo(RGset)
qc$badsample
qc$badCpG
qc$outlier_sample
gc()

# Pre-Processing And Filtering
# ============================
# create a MethylSet object from the raw data for plotting
# Mset: MethylSet
# Mset <- preprocessQuantile(RGset)
Mset_enmix <- preprocessENmix(RGset, nCores = 6, QCinfo = qc)
Mset_enmix_quantile <- norm.quantile(Mset_enmix, method = "quantile1")
Mset <- mapToGenome(Mset_enmix_quantile)
Mset
dim(Mset)

gc()

# ensure probes are in the same order in the Mset and detP objects
detP <- detP[match(featureNames(Mset),rownames(detP)),]
dim(detP)

# Writing data sets
# Save detP
saveRDS(detP, file = "detP.rds")

# keep only probes that have passed in 100% of the samples
per_samples <- round(((length(colnames(Mset)) * 100) / 100))
keep <- rowSums(detP < 0.01) >= per_samples
table(keep)

MsetFlt <- Mset[keep,]
message ("After keeping only probes that have passed in 100% of the samples")
dim(MsetFlt)

# get the 450k annotation data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)

# if your data includes males and females, remove probes on the sex chromosomes
keep <- !(featureNames(MsetFlt) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
table(keep)
MsetFlt <- MsetFlt[keep,]
message ("After removing probes on the sex chromosomes")
dim(MsetFlt)

# remove probes with SNPs at CpG site
MsetFlt <- dropLociWithSnps(MsetFlt)
message ("After removing probes with SNPs at CpG site")
dim(MsetFlt)

# exclude cross reactive probes
MsetFlt <- maxprobes::dropXreactiveLoci(MsetFlt)
message ("After excluding cross-reactive probes")
dim(MsetFlt)

# Combine intensity values
Mset <- MethylSet(Meth = getMeth(MsetFlt), Unmeth = getUnmeth(MsetFlt), annotation = annotation(MsetFlt), preprocessMethod = preprocessMethod(MsetFlt))
Mset

# Writing data sets
# Save Mset
saveRDS(Mset, file = "Mset.rds")

# probe-type bias adjustment
beta <- rcp(Mset, qcscore = qc)
mval <- beta2m(beta)
dim(beta)
dim(mval)

# Extracting Beta and M-values
#beta <- getBeta(MsetFlt)
#mval <- getM(MsetFlt)

remove(qc, Mset_enmix, Mset_enmix_quantile, Mset, MsetFlt, detP)
gc()

# Read in methylation SNPs information
# Total CpGs: 89,678
load("/cluster/home/t110989uhn/kridelgroup/rajesh/00_Methylation_Analysis/00_Minfi/RNASeq/MethylMix/MethylMix_2.34.0/MethylMix/data/SNPprobes.rda")
length(SNPprobes)

# Remove SNP probes by excluding them from the beta and mval matrices
beta <- beta[setdiff(rownames(beta), SNPprobes), ]
mval <- mval[setdiff(rownames(mval), SNPprobes), ]

dim(beta)
dim(mval)

# Writing data sets
# Beta and M values 
saveRDS(beta, "beta.rds")
saveRDS(mval, "mval.rds")

# ensure samples are in the same order in the beta and mval
targets <- targets[match(colnames(beta), targets$Sample_Name),]

# Writing data sets
# Save sample annotations
write.table(targets, file = "targets.tsv", sep="\t", quote = FALSE, row.names = FALSE)
gc()

# Batch effects correction with Harman
table(targets$TYPE)
table(targets$Batch)

# Copy TYPE to TYPES
targets$TYPES <- targets$TYPE
table(targets$TYPES)

# Update TYPE based on COO while keeping other rows unchanged
targets$TYPES[targets$COO == "GCB"] <- "DLBCL_GCB"
targets$TYPES[targets$COO == "ABC"] <- "DLBCL_ABC"

table(targets$TYPES)
table(targets$Batch)

# Rename Type
targets$TYPES[targets$TYPES == 'FL'] <- 'FL'
targets$TYPES[targets$TYPES == 'FL_A'] <- 'FL_A'
targets$TYPES[targets$TYPES == 'FL_B'] <- 'FL_B'

# targets$TYPES[targets$TYPES == 'DLBCL'] <- 'DLBCL'
targets$TYPES[targets$TYPES == 'DLBCL_ABC'] <- 'DLBCL_ABC'
targets$TYPES[targets$TYPES == 'DLBCL_GCB'] <- 'DLBCL_GCB'
targets$TYPES[targets$TYPES == 'DLBCL_EGA'] <- 'DLBCL_EGA'
targets$TYPES[targets$TYPES == 'DLBCL_EN'] <- 'DLBCL_EN'
targets$TYPES[targets$TYPES == 'DLBCL_GC'] <- 'DLBCL_GC'
targets$TYPES[targets$TYPES == 'DLBCL_nonGC'] <- 'DLBCL_nonGC'
targets$TYPES[targets$TYPES == 't_DLBCL'] <- 't_DLBCL'

targets$TYPES[targets$TYPES == 'Benign_LN'] <- 'Benign_LN'
targets$TYPES[targets$TYPES == 'LN_NN_A'] <- 'LN_NN_A'
targets$TYPES[targets$TYPES == 'LN_NN_B'] <- 'LN_NN_B'
targets$TYPES[targets$TYPES == 'RLN'] <- 'RLN'

targets$TYPES[targets$TYPES == 'Normal_GCB'] <- 'Normal_GCB'
targets$TYPES[targets$TYPES == 'bm_PC'] <- 'bm_PC'
targets$TYPES[targets$TYPES == 'gcBC'] <- 'gcBC'
targets$TYPES[targets$TYPES == 'memBC'] <- 'memBC'
targets$TYPES[targets$TYPES == 'naiBC'] <- 'naiBC'
targets$TYPES[targets$TYPES == 't_PC'] <- 't_PC'
targets$TYPES[targets$TYPES == 't_naiBC'] <- 't_naiBC'

table(targets$TYPES)
table(targets$Batch)

methHarman <- harman(mval, expt = targets$TYPES, batch = targets$Batch)
summary(methHarman)

combat_mval <- reconstructData(methHarman)
combat_beta <- m2beta(combat_mval)

dim(combat_mval)
dim(combat_beta)

remove(methHarman)
gc()

# Writing out M and beta matrices
saveRDS(combat_mval, "mval_combat.rds")
saveRDS(combat_beta, "beta_combat.rds")

# Data exploration 1 - Before and After normalization
# ===================================================
type_colors <- c("#5086F3", "#7f8992", "#D02060", "#fc1417", "#f26ded", "#D87087", "#9A3B8A", "#6dcc04", "#02bdaa", "#ebc602", "#0303a8", "#FF7F00", "#B15928", "#b5b80b", "#2480b9", "#2E2EAA", "#7618dc", "#6E0A0A", "#6e988c", "#383838")
batch_colors <- c("#383838", "#0303a8", "#02bdaa", "#A01641", "#6dcc04", "#FF7F00", "#7818e1", "#b5b80b")

# factor:Type, Batch - visualize what the data looks like before and after normalization
pdf("Density-Raw-Beta-RGset-Type.pdf")
pal <- type_colors
densityPlot(RGset,sampGroups=targets$TYPE, main="Beta-Value - Raw Data - Type", legend=FALSE, pal=pal)
legend("top", legend = levels(factor(targets$TYPE)),
       text.col=pal)
dev.off()

pdf("Density-Raw-Beta-RGset-Batch.pdf")
pal <- batch_colors
densityPlot(RGset, sampGroups=targets$Batch, main="Beta-Value - Raw Data - Batch", legend=FALSE, pal=pal)
legend("top", legend = levels(factor(targets$Batch)),
       text.col=pal)
dev.off()

pdf("Density-Norm-Beta-MsetFlt-Type.pdf")
pal <- type_colors
densityPlot(beta, sampGroups=targets$TYPE, main="Beta-Value - Normalized & Filtered Data - Type", legend = FALSE, pal=pal)
legend("top", legend = levels(factor(targets$TYPE)), 
       text.col=pal)
dev.off()

pdf("Density-Norm-Beta-MsetFlt-Batch.pdf")
pal <- batch_colors
densityPlot(beta, sampGroups=targets$Batch, main="Beta-Value - Normalized and Filtered - Batch", legend = FALSE, pal=pal)
legend("top", legend = levels(factor(targets$Batch)),
       text.col=pal)
dev.off()

# MDS plots to look at largest sources of variation
pdf("MDS-Norm-Beta-MsetFlt-Type.pdf")
pal <- type_colors
plotMDS(beta,
        top=1000, gene.selection="common",
        col=pal[factor(targets$TYPE)],
        pch=19,
        main="Beta-Value - Normalized & Filtered Data - Type"
)
legend("topleft",
       legend=levels(factor(targets$TYPE)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()

# MDS plots to look at largest sources of variation
pdf("MDS-Norm-Beta-MsetFlt-Batch.pdf")
pal <- batch_colors
plotMDS(beta,
        top=1000, gene.selection="common",
        col=pal[factor(targets$Batch)],
        pch=19,
        main="Beta-Value - Normalized & Filtered Data - Batch"
)
legend("topleft",
       legend=levels(factor(targets$Batch)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()

# MDS plots to look at largest sources of variation
pdf("MDS-Norm-M-MsetFlt-Type.pdf")
pal <- type_colors
plotMDS(mval,
        top=1000, gene.selection="common",
        col=pal[factor(targets$TYPE)],
        pch=19,
        main="M-Value - Normalized & Filtered Data - Type"
)
legend("topleft",
       legend=levels(factor(targets$TYPE)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()

# MDS plots to look at largest sources of variation
pdf("MDS-Norm-M-MsetFlt-Batch.pdf")
pal <- batch_colors
plotMDS(mval,
        top=1000, gene.selection="common",
        col=pal[factor(targets$Batch)],
        pch=19,
        main="M-Value - Normalized & Filtered Data - Batch"
)
legend("topleft",
       legend=levels(factor(targets$Batch)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()

gc()

# Data exploration 2 - after batch effects correction 
# ===================================================
# factor: Type, Batch - Visualize what the data looks like after batch effects correction
pdf("Density-Norm-Combat-Beta-Type.pdf")
pal <- type_colors
densityPlot(combat_beta, sampGroups=targets$TYPE, main="Beta-Value - Normalized, Filtered, Batch Corrected Data - Type", legend = FALSE, pal=pal)
legend("top", legend = levels(factor(targets$TYPE)), 
       text.col=pal)
dev.off()

pdf("Density-Norm-Combat-Beta-Batch.pdf")
pal <- batch_colors
densityPlot(combat_beta, sampGroups=targets$Batch, main="Beta-Value - Normalized, Filtered, Batch Corrected Data - Batch", legend = FALSE, pal=pal)
legend("top", legend = levels(factor(targets$Batch)),
       text.col=pal)
dev.off()

# MDS plots to look at largest sources of variation
pdf("MDS-Norm-Combat-Beta-Type.pdf")
pal <- type_colors
plotMDS(combat_beta,
        top=1000, gene.selection="common",
        col=pal[factor(targets$TYPE)],
        pch=19,
        main="Beta-Value - Normalized, Filtered, Batch Corrected Data - Type"
)
legend("topleft",
       legend=levels(factor(targets$TYPE)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()

# MDS plots to look at largest sources of variation
pdf("MDS-Norm-Combat-Beta-Batch.pdf")
pal <- batch_colors
plotMDS(combat_beta,
        top=1000, gene.selection="common",
        col=pal[factor(targets$Batch)],
        pch=19,
        main="Beta-Value - Normalized, Filtered, Batch Corrected Data - Batch"
)
legend("topleft",
       legend=levels(factor(targets$Batch)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()

# MDS plots to look at largest sources of variation
pdf("MDS-Norm-Combat-Beta-Stage.pdf")
pal <- c("#5086F3", "#0303a8", "#9A3B8A", "#6E0A0A")
plotMDS(combat_beta,
        top=1000, gene.selection="common",
        col=pal[factor(targets$ANN_ARBOR_STAGE)],
        pch=19,
        main="Beta-Value - Normalized, Filtered, Batch Corrected Data - Stage"
)
legend("topleft",
       legend=levels(factor(targets$ANN_ARBOR_STAGE)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()

# MDS plots to look at largest sources of variation
pdf("MDS-Norm-Combat-M-Type.pdf")
pal <- type_colors
plotMDS(combat_mval,
        top=1000, gene.selection="common",
        col=pal[factor(targets$TYPE)],
        pch=19,
        main="M-Value - Normalized, Filtered, Batch Corrected Data - Type"
)
legend("topleft",
       legend=levels(factor(targets$TYPE)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()

# MDS plots to look at largest sources of variation
pdf("MDS-Norm-Combat-M-Batch.pdf")
pal <- batch_colors
plotMDS(combat_mval,
        top=1000, gene.selection="common",
        col=pal[factor(targets$Batch)],
        pch=19,
        main="M-Value - Normalized, Filtered, Batch Corrected Data - Batch"
)
legend("topleft",
       legend=levels(factor(targets$Batch)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()

# MDS plots to look at largest sources of variation
pdf("MDS-Norm-Combat-M-Stage.pdf")
pal <- c("#5086F3", "#0303a8", "#9A3B8A", "#6E0A0A")
plotMDS(combat_mval,
        top=1000, gene.selection="common",
        col=pal[factor(targets$ANN_ARBOR_STAGE)],
        pch=19,
        main="M-Value - Normalized, Filtered, Batch Corrected Data - Stage"
)
legend("topleft",
       legend=levels(factor(targets$ANN_ARBOR_STAGE)),
       text.col=pal,
       bg="white",
       cex=0.7)
dev.off()

gc()

# Detailed Information About Current R Session
sessionInfo()
