# R Packages
# ==========
library(dplyr)
library(stringr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(ggplot2)
library(DGEobj.utils)
library(edgeR)
library(ggpubr)
library(DESeq2)
library(EnhancedVolcano)
library(apeglm)

# Change directory
setwd("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/18_DGE_Analysis/edgeR/")

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

# exclude blueprint samples
Sample_Annotations <- Sample_Annotations %>% filter(TYPE == "FL" | TYPE == "FL_A" | TYPE == "FL_B" | TYPE == "t_DLBCL" | TYPE == "DLBCL" | TYPE == "DLBCL_EGA" | TYPE == "DLBCL_EN" | TYPE == "DLBCL_GC" | TYPE == "DLBCL_nonGC" | TYPE == "Benign_LN" | TYPE == "LN_NN_A" | TYPE == "LN_NN_B" | TYPE == "RLN")

dim(Sample_Annotations)
table(Sample_Annotations$TYPE)
table(Sample_Annotations$epitypes)

# Calculate Samples by Type per Cluster
Sample_Annotations$TYPE[Sample_Annotations$TYPE == "FL"] <- "FL"
Sample_Annotations$TYPE[Sample_Annotations$TYPE == "FL_A"] <- "FL"
Sample_Annotations$TYPE[Sample_Annotations$TYPE == "FL_B"] <- "FL"

Sample_Annotations$TYPE[Sample_Annotations$TYPE == "DLBCL"] <- "DLBCL"
# Sample_Annotations$TYPE[Sample_Annotations$TYPE == "t_DLBCL"] <- "DLBCL"
Sample_Annotations$TYPE[Sample_Annotations$TYPE == "DLBCL_EGA"] <- "DLBCL"
Sample_Annotations$TYPE[Sample_Annotations$TYPE == "DLBCL_EN"] <- "DLBCL"
Sample_Annotations$TYPE[Sample_Annotations$TYPE == "DLBCL_GC"] <- "DLBCL"
Sample_Annotations$TYPE[Sample_Annotations$TYPE == "DLBCL_nonGC"] <- "DLBCL"

Sample_Annotations$TYPE[Sample_Annotations$TYPE == "Benign_LN"] <- "Normal"
Sample_Annotations$TYPE[Sample_Annotations$TYPE == "LN_NN_A"] <- "Normal"
Sample_Annotations$TYPE[Sample_Annotations$TYPE == "LN_NN_B"] <- "Normal"
Sample_Annotations$TYPE[Sample_Annotations$TYPE == "RLN"] <- "Normal"

Sample_Annotations %>% dplyr::count(epitypes, TYPE)

# Rename epitypes details
Sample_Annotations <- Sample_Annotations %>%
  filter((epitypes == "C1" & TYPE == "FL") |
           (epitypes == "C2" & TYPE == "FL") |
           (epitypes == "C3" & TYPE == "Normal")) %>%
  mutate(epitypes = case_when(
    epitypes == "C1" & TYPE == "FL"     ~ "iFL",
    epitypes == "C2" & TYPE == "FL"     ~ "aFL",
    epitypes == "C3" & TYPE == "Normal" ~ "Normal"
  ))

Sample_Annotations %>% dplyr::count(epitypes, TYPE)

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
# iFL vs. aFL
FL <- Sample_Annotations %>% filter(TYPE == "FL") %>% data.frame() %>% .$SAMPLE_ID # n=305
iFL <- Sample_Annotations %>% filter(TYPE == "FL" & epitypes == "iFL") %>% data.frame() %>% .$SAMPLE_ID # n=222
aFL <- Sample_Annotations %>% filter(TYPE == "FL" & epitypes == "aFL") %>% data.frame() %>% .$SAMPLE_ID # n=83
Normal <- Sample_Annotations %>% filter(TYPE == "Normal") %>% data.frame() %>% .$SAMPLE_ID # n=20

# Load RNA-seq data
# =================
# FL RNA-Seq, n=282

mRNA <- read.table("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/RNA_Seq/2022-10-25_df_counts_adj.txt", sep = " ", header = TRUE)
passing.RNAseq.QC <- read.table("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/RNA_Seq/metadata_passed_290.txt",
                                sep = "\t",
                                header = TRUE) %>% subset(phenotype == "FL")
                                
length(unique(passing.RNAseq.QC$sample_id))

mRNA <- mRNA[,passing.RNAseq.QC$id]
dim(mRNA)
mRNA[1:5, 1:5]

# Remove rows consistently have zero or very low counts
keep <- filterByExpr(mRNA)
table(keep)

mRNA <- mRNA[keep,]
mRNA[1:5, 1:5]
dim(mRNA)

# Convert Gene IDs to Gene Symbol
gene_info <- read.csv("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/ensembl/annotation/gene_info.csv", sep = ",", header = TRUE)
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

# Filter Matrix for Protein Coding Genes
mRNA <- mRNA %>% filter(GENE_TYPE == "protein_coding")
dim(mRNA)

# Identify duplicate genes and remove redundancy
# n_occur <- data.frame(table(mRNA$SYMBOL))
# n_occur[n_occur$Freq > 1,]
# n_occur_df <- mRNA[mRNA$SYMBOL %in% n_occur$Var1[n_occur$Freq > 1], ]
# n_occur_list <- c("ENSG00000168255", "ENSG00000285053") # removed shortest isoforms
# mRNA <- mRNA[!(row.names(mRNA) %in% n_occur_list),]

# Filter mRNA to exclude rows with duplicated SYMBOL
n_occur <- data.frame(table(mRNA$SYMBOL))
dup_genes <- n_occur[n_occur$Freq > 1, "Var1"]
mRNA <- mRNA[!(mRNA$SYMBOL %in% dup_genes), ]
dim(mRNA)
mRNA[1:10, 1:10]

# Set SYMBOL as row names
rownames(mRNA) <- mRNA$SYMBOL

# Remove columns
mRNA <- mRNA[,!names(mRNA) %in% c("SYMBOL", "GENE_TYPE", "MEAN", "MEDIAN", "LONGEST_ISOFORM", "MERGED")]
dim(mRNA)

# Prepare phenotype sample data
coldata <- Sample_Annotations %>% select(epitypes)
# coldata <- coldata %>% mutate(type = "paired-end")
names(coldata)[names(coldata) == "epitypes"] <- "condition"
coldata <- coldata[rownames(coldata) %in% colnames(mRNA), , drop = FALSE]
coldata$condition <- factor(coldata$condition, levels = c("iFL", "aFL"))
table(coldata$condition)
head(coldata)
dim(coldata)

# Filter read counts file
cts <- mRNA[, colnames(mRNA) %in% rownames(coldata)]
cts[1:5, 1:5]
dim(cts)

# check sample order
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))

# Create a sample information for the count data
sample_info <- coldata$condition
class(sample_info) # factor

# Create DGEList data class for count and sample information
dge <- DGEList(counts = cts, group = sample_info)
class(dge)
dge

# Filter out the genes with low counts
keep <- filterByExpr(y = dge)
table(keep)
dge <- dge[keep, ,keep.lib.sizes = FALSE]
dge

# Normalization and effective library sizes
dge <- calcNormFactors(object = dge)
dge

# Model fitting and estimating dispersions
dge <- estimateDisp(y = dge)
dge

# Testing for differential gene expression
et <- exactTest(object = dge)
et

# Get a summary DGE table - log fold change at least 1 and adjusted p value < 0.05
# summary(decideTests(object = et, lfc = 1))

# Extract the table with adjusted p values (FDR)
top_degs <- topTags(object = et, n = "Inf")
top_degs
dim(top_degs)

# Export differential gene expression analysis table to CSV file
write.table(as.data.frame(top_degs), file = "DEG_aFL_vs_iFL_All.tsv", quote = FALSE, sep = "\t")

# Get summary of differential gene expression with adjusted p value cut-off (FDR) at 0.05
top_degs05 <- subset(top_degs$table, top_degs$table$FDR < 0.05)
dim(top_degs05)
write.table(top_degs05, file = "DEG_aFL_vs_iFL_top_degs05.tsv", quote = FALSE, sep = "\t")

# up-regulated genes
top_degs05up <- subset(top_degs05, logFC >= 0.1)
dim(top_degs05up)
write.table(top_degs05up, file = "DEG_aFL_vs_iFL_top_degs05up.tsv", quote = FALSE, sep = "\t")

# down-regulated genes
top_degs05down <- subset(top_degs05, logFC <= -0.1)
dim(top_degs05down)
write.table(top_degs05down, file = "DEG_aFL_vs_iFL_top_degs05down.tsv", quote = FALSE, sep = "\t")

# Visualize the shrinkage estimation of LFCs with MA plot and compare it without shrinkage of LFCs
pdf("MAplot.pdf", height = 5, width = 10)
plotMD(object = et)
dev.off()

# Volcano Plot
# png("volcano-plot.png", units="px", width=4500, height=4000, res=600)
maxY <- max(-log10(top_degs$table$FDR), na.rm = TRUE)

pdf("volcano-plot.pdf", height = 8, width = 8)
EnhancedVolcano(top_degs$table,
                lab = rownames(top_degs$table),
                x = "logFC",
                y = "FDR",
                pCutoff = 0.05,
                FCcutoff = 0.1,
                selectLab = c("NA"),
                ylim = c(0, maxY + 0.1))
dev.off()
