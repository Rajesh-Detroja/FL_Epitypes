# R Packages
# ==========
library(dplyr)
library(limma)
library(missMethyl)
library(DMRcate)
library(ggplot2)
library(RColorBrewer)
library(edgeR)
library(ggrepel)
library(stringr)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(IlluminaHumanMethylationEPICv2manifest)
library(DMRcatedata)

# Change directory
setwd("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/07_DMP_DMR/450k/hg38/FL_vs_DLBCL/")

# Read EPICv2_hg38 annotation
EPICv2_hg38 <- read.table(file = "/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/ANN_EPIC_v2/EPICv2_hg38.tsv", sep = "\t", header = TRUE, row.names = 1, na.strings = c("","NA"))
dim(EPICv2_hg38)

# Read in methylation data
bval <- readRDS("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/01_InfiniumPurify/450k/modified_InfiniumPurify/Enmix_quantile_with_QCinfo_Harman_All/beta_combat_purified.rds")
mval <- readRDS("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/01_InfiniumPurify/450k/modified_InfiniumPurify/Enmix_quantile_with_QCinfo_Harman_All/mval_combat_purified.rds")

dim(bval)
dim(mval)

# Filter in methylation data
bval <- bval[intersect(EPICv2_hg38$EPICv1_Loci, rownames(bval)) ,]
mval <- mval[intersect(EPICv2_hg38$EPICv1_Loci, rownames(mval)) ,]

dim(bval)
dim(mval)

EPICv2_hg38 <- EPICv2_hg38 %>% filter(EPICv1_Loci %in% rownames(bval))
dim(EPICv2_hg38)

gc()

# Rename CpG probes from EPICv1 to EPICv2
bval[1:5, 1:5]
mval[1:5, 1:5]

# Create a named vector with new row names from EPICv2_hg38
new_row_names <- EPICv2_hg38$Name
names(new_row_names) <- EPICv2_hg38$EPICv1_Loci
length(new_row_names)

# Match and replace row names in methylation data
matched_names <- new_row_names[row.names(bval)]
row.names(bval) <- matched_names[row.names(bval)]

matched_names <- new_row_names[row.names(mval)]
row.names(mval) <- matched_names[row.names(mval)]

# Check the class of methylation data
bval[1:5, 1:5]
mval[1:5, 1:5]

class(bval)
class(mval)

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

# Filter samples
Sample_Annotations <- Sample_Annotations %>% filter(TYPE == "FL" | TYPE == "FL_A" | TYPE == "FL_B" | TYPE == "DLBCL" | TYPE == "DLBCL_EGA" | TYPE == "DLBCL_EN" | TYPE == "DLBCL_GC" | TYPE == "DLBCL_nonGC" | TYPE == "t_DLBCL" | TYPE == "Benign_LN" | TYPE == "LN_NN_A" | TYPE == "LN_NN_B" | TYPE == "RLN") %>% data.frame()

Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'FL_A'] <- 'FL'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'FL_B'] <- 'FL'

Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'DLBCL_EGA'] <- 'DLBCL'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'DLBCL_EN'] <- 'DLBCL'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'DLBCL_GC'] <- 'DLBCL'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'DLBCL_nonGC'] <- 'DLBCL'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 't_DLBCL'] <- 'DLBCL'

Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'Benign_LN'] <- 'Normal'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'LN_NN_A'] <- 'Normal'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'LN_NN_B'] <- 'Normal'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'RLN'] <- 'Normal'

table(Sample_Annotations$TYPE)
table(Sample_Annotations$Batch)

row.names(Sample_Annotations) <- Sample_Annotations$SAMPLE_ID

# subset of merged methylation data
bval <- bval[, Sample_Annotations$SAMPLE_ID]
mval <- mval[, Sample_Annotations$SAMPLE_ID]
dim(bval)
dim(mval)

# Write TYPE and Batch Information
df1 <- Sample_Annotations[c("SAMPLE_ID", "TYPE","Batch")]
write.table(df1, file="TYPE_INFO_ALL.txt", sep="\t", quote = FALSE, row.names = FALSE)

gc()

# Create list of groups
FL <- Sample_Annotations %>% filter(TYPE == "FL") %>% data.frame() %>% .$SAMPLE_ID # n = 338
length(unique(FL))

DLBCL <- Sample_Annotations %>% filter(TYPE == "DLBCL") %>% data.frame() %>% .$SAMPLE_ID # n = 157
length(unique(DLBCL))

Normal <- Sample_Annotations %>% filter(TYPE == "Normal") %>% data.frame() %>% .$SAMPLE_ID # n = 20
length(unique(Normal))

# subset annotations for group
ref <- "FL" # CHANGE
epi_types <- "DLBCL" # CHANGE
group_samples <- c(get(ref), get(epi_types))
Group_Annotations <- Sample_Annotations %>% filter(SAMPLE_ID %in% group_samples)
table(Group_Annotations$TYPE)

# Set Variables
dmp_tsv <- gsub(" ", "", paste("DMP_",ref,"_",epi_types,".tsv"))
dmp_pdf <- gsub(" ", "", paste("DMP_",ref,"_",epi_types,".pdf"))
dmr_tsv_005 <- gsub(" ", "", paste("DMR_",ref,"_",epi_types,"_FDR_005",".tsv"))
dmr_tsv_001 <- gsub(" ", "", paste("DMR_",ref,"_",epi_types,"_FDR_001",".tsv"))
dmp_vol_005 <- gsub(" ", "", paste("DMP_VOL_005_",ref,"_",epi_types,".pdf"))
dmp_vol_001 <- gsub(" ", "", paste("DMP_VOL_001_",ref,"_",epi_types,".pdf"))
diff_type <- gsub("  ", " ", paste(ref," vs." , epi_types))

# Subset methylation data group
bval_group <- bval[, Group_Annotations$SAMPLE_ID]
mval_group <- mval[, Group_Annotations$SAMPLE_ID]
dim(bval_group)
dim(mval_group)

# Making sure order of samples in `Group_Annotations` and `matrices`
table(colnames(bval_group) %in% Group_Annotations$SAMPLE_ID)
table(colnames(mval_group) %in% Group_Annotations$SAMPLE_ID)

all(row.names(Group_Annotations$SAMPLE_ID) == colnames(bval_group))
all(row.names(Group_Annotations$SAMPLE_ID) == colnames(mval_group))

## Differentially methylated probes (DMPs) analysis
## ------------------------------------------------
# check levels
Group_Annotations$TYPE <- factor(Group_Annotations$TYPE)
Group_Annotations$TYPE
# Levels: DLBCL FL (Reference level is DLBCL)

Group_Annotations$TYPE <- relevel(Group_Annotations$TYPE, ref = ref)
Group_Annotations$TYPE
# Levels: FL DLBCL (Reference level is FL)

# Model Design
design <- model.matrix(~ Group_Annotations$TYPE)

# fit the linear model
fit <- lmFit(mval_group, design)
fit2 <- eBayes(fit)

# clear RAM
gc()

# look at the numbers of DMPs at FDR < 0.05
summary(decideTests(fit2, p.value = 0.05, lfc = 0))

# look at the numbers of DMPs at FDR < 0.01
summary(decideTests(fit2, p.value = 0.01, lfc = 1))

# Generate table of DMPs with annotation
EPICv2_hg38_sub <- EPICv2_hg38[match(rownames(mval_group), EPICv2_hg38$Name), c(2:4,42,19,20,30,37)]
dim(EPICv2_hg38_sub)

# extract significant DMPs
deff.meth <- topTable(fit2,
                      coef = ncol(design),
                      p.value = 0.05,
                      sort.by = "p",
                      number = nrow(mval_group),
                      genelist = EPICv2_hg38_sub)

# add `HYPER` and `HYPO` methylation status for FDR < 0.05
deff.meth$HYPER_HYPO_005 <- "NO_CHANGE"

deff.meth$HYPER_HYPO_005[deff.meth$logFC > 0 & deff.meth$adj.P.Val < 0.05] <- "HYPER"
deff.meth$HYPER_HYPO_005[deff.meth$logFC < 0 & deff.meth$adj.P.Val < 0.05] <- "HYPO"
table(deff.meth$HYPER_HYPO_005)

# add `HYPER` and `HYPO` methylation status for FDR < 0.01
deff.meth$HYPER_HYPO_001 <- "NO_CHANGE"

deff.meth$HYPER_HYPO_001[deff.meth$logFC > 1 & deff.meth$adj.P.Val < 0.01] <- "HYPER"
deff.meth$HYPER_HYPO_001[deff.meth$logFC < -1 & deff.meth$adj.P.Val < 0.01] <- "HYPO"
table(deff.meth$HYPER_HYPO_001)

deff.meth <- deff.meth %>% relocate(HYPER_HYPO_005, .before = chr)
deff.meth <- deff.meth %>% relocate(HYPER_HYPO_001, .before = chr)

write.table(deff.meth, file = dmp_tsv, quote = FALSE, sep = '\t', col.names = NA)

## Visualization
## -------------
# plot top 10 most significant DMPs
pdf(dmp_pdf)
par(mfrow = c(2,5))
sapply(rownames(deff.meth)[1:10], function(cpg){
  plotCpg(bval_group, cpg = cpg, pheno = Group_Annotations$TYPE, ylab = "Beta values")
})
dev.off()

gc()

# generate a volcano plot
# -----------------------
# FDR: 0.05, logFC: 0
dat <- dplyr::select(deff.meth, c(logFC, adj.P.Val))
dat$mstatus <- "NO_CHANGE"

dat$mstatus[dat$logFC > 0 & dat$adj.P.Val < 0.05] <- "HYPER"
dat$mstatus[dat$logFC < 0 & dat$adj.P.Val < 0.05] <- "HYPO"
table(dat$mstatus)

HYPER_COUNTS <- sum(dat$mstatus == "HYPER")
HYPO_COUNTS <- sum(dat$mstatus == "HYPO")

dat$mlabel <- NA  # Initialize mlabel with NA
# Optionally, specify which points to label (e.g., top significant ones)
# dat$mlabel[dat$adj.P.Val < 1e-10] <- rownames(dat)[dat$adj.P.Val < 1e-10]

pdf(dmp_vol_005)
cols <- c("HYPO" = "blue", "HYPER" = "red")
ggplot(data = dat, aes(x = logFC, y = -log10(adj.P.Val), col = mstatus, label = mlabel)) +
  geom_point() +
  theme_gray() +
  geom_text_repel(show.legend = FALSE, na.rm = TRUE) +  # na.rm = TRUE to remove NAs
  scale_colour_manual(values = cols) +
  ggtitle(paste(diff_type, "| HYPER: ", HYPER_COUNTS, "| HYPO: ", HYPO_COUNTS, "| FDR: 0.05 | logFC: 0")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +
  xlab("Log (Fold Change)") +
  ylab("-log10 (adj.P.Val)") +
  labs(color = "CpGs") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  geom_vline(xintercept = 0, colour = "black", linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.05), colour="black", linetype = "dashed")
dev.off()


# FDR: 0.01, logFC: 1
dat <- dplyr::select(deff.meth, c(logFC, adj.P.Val))
dat$mstatus <- "NO_CHANGE"

dat$mstatus[dat$logFC > 1 & dat$adj.P.Val < 0.01] <- "HYPER"
dat$mstatus[dat$logFC < -1 & dat$adj.P.Val < 0.01] <- "HYPO"

dat$mlabel <- NA  # Initialize mlabel with NA
# Optionally, specify which points to label (e.g., top significant ones)
# dat$mlabel[dat$adj.P.Val < 1e-10] <- rownames(dat)[dat$adj.P.Val < 1e-10]

HYPER_COUNTS <- sum(dat$mstatus == "HYPER")
HYPO_COUNTS <- sum(dat$mstatus == "HYPO")
NC_COUNTS <- sum(dat$mstatus == "NO_CHANGE")

table(dat$mstatus)

pdf(dmp_vol_001)
cols <- c("HYPO" = "blue", "HYPER" = "red", "NO_CHANGE" = "grey")
ggplot(data = dat, aes(x = logFC, y = -log10(adj.P.Val), col = mstatus, label = mlabel)) +
  geom_point() +
  theme_gray() +
  geom_text_repel(show.legend = FALSE, na.rm = TRUE) +  # na.rm = TRUE to remove NAs
  scale_colour_manual(values = cols) +
  ggtitle(paste(diff_type, "| HYPER: ", HYPER_COUNTS, " | HYPO: ", HYPO_COUNTS, " | NO_CHANGE: ", NC_COUNTS, " | FDR: 0.01 | logFC: 1")) +
  theme(plot.title = element_text(hjust = 0.5, size = 9)) +
  xlab("Log (Fold Change)") +
  ylab("-log10 (adj.P.Val)") +
  labs(color = "CpGs") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  geom_vline(xintercept = 1, colour = "black", linetype = "dashed") + 
  geom_vline(xintercept = -1, colour = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), colour = "black", linetype = "dashed")
dev.off()


# clear RAM
gc()

## Differentially methylated regions (DMRs) analysis
## =================================================

## DMR Analysis - FDR = 0.05
## =========================
# Only include DMPs
hyper_hypo_rownames <- rownames(deff.meth[deff.meth$HYPER_HYPO_005 %in% c("HYPER", "HYPO"), ])

mval_group_005 <- mval_group[intersect(hyper_hypo_rownames, rownames(mval_group)) ,]
dim(mval_group_005)

# If not a matrix, convert to matrix
if (!inherits(mval_group, "matrix")) {
  mval_group_005 <- as.matrix(mval_group_005)
}

class(mval_group_005)

# Remove remaining datasets
rm(EPICv2_hg38, EPICv2_hg38_sub, bval, mval, bval_group, dat, fit, fit2)

gc()

# annotation settings
myAnnotation <- cpg.annotate(object = mval_group_005,
                             datatype = "array",
                             what = "M",
                             analysis.type = "differential",
                             design = design,
                             contrasts = FALSE,
                             coef = ncol(design),
                             arraytype = "EPICv2",
                             fdr = 0.05,
                             annotation=c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
                             )
str(myAnnotation)

# clear RAM
gc()

# DMR analysis
DMRs <- dmrcate(myAnnotation, lambda = 1000, C = 2)
results.ranges <- extractRanges(DMRs)
results.ranges

# visualization
dmr.table <- data.frame(results.ranges)

# Assuming `meanDiff` indicates the direction of methylation change
dmr.table$HYPER_HYPO_005 <- "NO_CHANGE"

dmr.table$HYPER_HYPO_005[dmr.table$meandiff > 0] <- "HYPER"
dmr.table$HYPER_HYPO_005[dmr.table$meandiff < 0] <- "HYPO"
table(dmr.table$HYPER_HYPO_005)

dmr.table <- dmr.table %>% relocate(HYPER_HYPO_005, .after = no.cpgs)

write.table(dmr.table, file = dmr_tsv_005, quote = FALSE, sep = '\t', row.names = FALSE)

# Remove remaining datasets
rm(myAnnotation, DMRs, results.ranges, dmr.table)

# clear RAM
gc()

## DMR Analysis - FDR = 0.01
## =========================
# Only include DMPs
hyper_hypo_rownames <- rownames(deff.meth[deff.meth$HYPER_HYPO_001 %in% c("HYPER", "HYPO"), ])

mval_group_001 <- mval_group[intersect(hyper_hypo_rownames, rownames(mval_group)) ,]
dim(mval_group_001)

# If not a matrix, convert to matrix
if (!inherits(mval_group_001, "matrix")) {
  mval_group_001 <- as.matrix(mval_group_001)
}

class(mval_group_001)

# annotation settings
myAnnotation <- cpg.annotate(object = mval_group_001,
                             datatype = "array",
                             what = "M",
                             analysis.type = "differential",
                             design = design,
                             contrasts = FALSE,
                             coef = ncol(design),
                             arraytype = "EPICv2",
                             fdr = 0.01,
                             annotation=c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
)
str(myAnnotation)

# clear RAM
gc()

# DMR analysis
DMRs <- dmrcate(myAnnotation, lambda = 1000, C = 2)
results.ranges <- extractRanges(DMRs)
results.ranges

# visualization
dmr.table <- data.frame(results.ranges)

# Assuming `meanDiff` indicates the direction of methylation change
dmr.table$HYPER_HYPO_001 <- "NO_CHANGE"

dmr.table$HYPER_HYPO_001[dmr.table$meandiff > 0] <- "HYPER"
dmr.table$HYPER_HYPO_001[dmr.table$meandiff < 0] <- "HYPO"
table(dmr.table$HYPER_HYPO_001)

dmr.table <- dmr.table %>% relocate(HYPER_HYPO_001, .after = no.cpgs)

write.table(dmr.table, file = dmr_tsv_001, quote = FALSE, sep = '\t', row.names = FALSE)

# Remove remaining datasets
rm(myAnnotation, DMRs, results.ranges, dmr.table)

# clear RAM
gc()
