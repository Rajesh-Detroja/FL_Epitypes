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
library(conumee2)
library(fixr)
library(readr)

# Change directory
setwd("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/09_CNV/450k/01_C1_FL/")

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

# subset of merged methylation data
bval <- bval[, Sample_Annotations$SAMPLE_ID]
mval <- mval[, Sample_Annotations$SAMPLE_ID]
dim(bval)
dim(mval)

# Write TYPE and Batch Information
df1 <- Sample_Annotations[c("SAMPLE_ID", "TYPE","Batch")]
write.table(df1, file="TYPE_INFO.txt", sep="\t", quote = FALSE, row.names = FALSE)

gc()

# Load the RDS file
Mset <- readRDS("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/450k/Enmix_quantile_with_QCinfo_Harman_All/rds/Mset.rds")
Mset
dim(Mset)

# Subset Mset 
common_cpgs <- intersect(featureNames(Mset), rownames(bval))
Mset <- Mset[common_cpgs, ]

common_samples <- intersect(colnames(Mset), colnames(bval))
Mset <- Mset[, common_samples]

Mset
dim(Mset)

## CNV Analysis
## ============
# Create annotation object
data(exclude_regions)
data(detail_regions)
head(detail_regions, n = 2)

anno <- CNV.create_anno(array_type = "450k",
                        exclude_regions = exclude_regions,
                        detail_regions = detail_regions,
                        chrXY = FALSE)

anno

# Filter overlapped probes
# https://support.bioconductor.org/p/97909/
Mset <- mapToGenome(Mset)
Mset

anno@probes <- subsetByOverlaps(anno@probes, granges(Mset))
anno

# Prepare Mset contains only the probes found in the annotation
methylset_probes <- rownames(Mset)
anno_probes <- names(anno@probes)

common_probes <- intersect(methylset_probes, anno_probes)
length(common_probes)

MsetFlt <- Mset[common_probes, ]
MsetFlt

Mset <- MethylSet(Meth = getMeth(MsetFlt), Unmeth = getUnmeth(MsetFlt), annotation = annotation(MsetFlt), preprocessMethod = preprocessMethod(MsetFlt))
Mset

# Combine intensity values
load.data <- CNV.load(Mset)
load.data

# Generate tumor and normal Intensity
# normal_samples <- Sample_Annotations %>% filter(TYPE == "RLN" | TYPE == "Benign_LN" | TYPE == "LN_NN_A" | TYPE == "LN_NN_B") %>% data.frame() %>% .$Sample_Name
# tumor_samples <- Sample_Annotations %>% filter(TYPE == "FL" | TYPE == "DLBCL" | TYPE == "DLBCL_EGA" | TYPE == "DLBCL_EN" | TYPE == "FL_A" | TYPE == "FL_B") %>% data.frame() %>% .$Sample_Name
# tumor_samples <- Sample_Annotations %>% filter(SAMPLE_ID == "LY_DLC_001") %>% data.frame() %>% .$SAMPLE_ID

normal_samples <- Sample_Annotations %>% filter(epitypes == "C3_Normal") %>% data.frame() %>% .$SAMPLE_ID
tumor_samples <- Sample_Annotations %>% filter(epitypes == "C1_FL") %>% data.frame() %>% .$SAMPLE_ID

normal <- load.data[normal_samples, ]
normal

tumor <- load.data[tumor_samples, ]
tumor

x <- CNV.fit(tumor, normal, anno)
x <- CNV.bin(x)
x <- CNV.detail(x)

# Remove both negative values and NA values from weights
total_col <- (length(tumor_samples))
total_col

for (i in 1:total_col) {
  cat("Removing negative and NAs values from sample: ", tumor_samples[i], "\n")
  
  # Calculate weights and remove bins with NA or negative values for CNV.segment analysis
  weights_x <- (1/x@bin$variance[[i]][names(x@anno@bins)])
  cat("`weights_x` before removing: ", length(weights_x), "\n")
  
  # Remove NA and negative values
  #weights_x_names <- names(weights_x)
  #print(table(is.na(weights_x)))  # Check for NA values
  
  cat("`weights_x` with NA or negative values: ", sum(is.na(weights_x) | weights_x < 0), "\n")
  cat("`weights_x` without NA and negative values: ", sum(!is.na(weights_x) & weights_x > 0), "\n")
  
  # Identify valid weights (non-NA and non-negative)
  valid_weights_x <- weights_x[!is.na(weights_x) & weights_x > 0]
  valid_weights_x_names <- names(valid_weights_x)
  
  # Check that invalid weights (NA or negative) are removed
  #print(table(is.na(valid_weights_x)))  # No NA should remain
  #print(table(valid_weights_x > 0))     # All values should be positive
  
  cat("`valid_weights_x` with NA or negative values: ", sum(is.na(valid_weights_x) | valid_weights_x < 0), "\n")
  cat("`valid_weights_x` without NA and negative values: ", sum(!is.na(valid_weights_x) & valid_weights_x > 0), "\n")
  
  # Update x@anno@bins to only keep bins with valid weights
  x@anno@bins <- x@anno@bins[names(x@anno@bins) %in% valid_weights_x_names]
  cat("`x@anno@bins` after removing NA or negative values: ", length(x@anno@bins), "\n")
  
  # Update x@bin$variance and x@bin$ratio with valid weights
  x@bin$variance[[i]] <- x@bin$variance[[i]][names(x@bin$variance[[i]]) %in% valid_weights_x_names]
  cat("`x@bin$variance` after removing NA or negative values: ", length(x@bin$variance[[i]]), "\n")
  
  x@bin$ratio[[i]] <- x@bin$ratio[[i]][names(x@bin$ratio[[i]]) %in% valid_weights_x_names]
  cat("`x@bin$ratio` after removing NA or negative values: ", length(x@bin$ratio[[i]]), "\n\n")
}

x <- CNV.segment(x)
show(x)
#saveRDS(x, "x.rds")
write_rds(x, "x.rds")

xf <- CNV.focal(x)
show(xf)
# saveRDS(xf, "xf.rds")
write_rds(xf, "xf.rds")

gc()

# Output generation
sample_name <- names(xf)
length(sample_name)

# CNV.genomeplot - All
# CNV.genomeplot(xf[1])
# CNV.genomeplot(xf, chr = 'chr6')

for (i in seq_along(sample_name)){
  pdf_filename <- paste0(sample_name[i], ".pdf")
  print(pdf_filename)
  
  pdf(pdf_filename, height = 14, width = 24)
  CNV.genomeplot(xf[i])
  dev.off()
}

# CNV.detailplot - All
#CNV.detailplot(xf, name = 'PTEN', output = "pdf", directory = dir)
#CNV.detailplot_wrap(xf)
CNV.detailplot_wrap(xf, output = "pdf")
# Error:
#DL45_D1173
#Error in xy.coords(x, y) : 'x' and 'y' lengths differ

# CNV.summaryplot
# CNV.summaryplot(xf, threshold = 0.15)
pdf("CNV.summaryplot.pdf", height = 10, width = 20)
CNV.summaryplot(xf)
dev.off()

# CNV.heatmap
# CNV.heatmap(xf)
pdf("CNV.heatmap.pdf", height = 10, width = 20)
CNV.heatmap(xf)
dev.off()

# output text files
segments <- CNV.write(xf, what = 'segments')
detail <- CNV.write(xf, what = 'detail')
bins <- CNV.write(xf, what = 'bins')
probes <- CNV.write(xf, what = 'probes')


# Downstream analysis
# -------------------
x <- readRDS("x.rds")
xf <- readRDS("xf.rds")

# output text files
segments <- CNV.write(xf, what = 'segments')
detail <- CNV.write(xf, what = 'detail')
bins <- CNV.write(xf, what = 'bins')
probes <- CNV.write(xf, what = 'probes')

# Write segments
segments_df <- rbind.fill(segments)
colnames(segments_df)[colnames(segments_df) == "ID"] <- "Sample_ID"
write.table(segments_df, file = "C1_FL_Segments.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Initialize result table
cnv_counts <- data.frame(SampleID = names(segments), 
                         Amplifications = integer(length(names(segments))), 
                         Deletions = integer(length(names(segments))),
                         stringsAsFactors = FALSE)

# Loop through each sample in `segments`
for (sample in names(segments)) {
  cnv_counts[cnv_counts$SampleID == sample, "Amplifications"] <- sum(segments[[sample]]$seg.mean > 0, na.rm = TRUE)
  cnv_counts[cnv_counts$SampleID == sample, "Deletions"] <- sum(segments[[sample]]$seg.mean < 0, na.rm = TRUE)
}

# Print final table
print(cnv_counts)
write.table(cnv_counts, file = "cnv_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# # Count positive values
# positive_count <- sum(segments$HP_AR70$seg.mean > 0)
# 
# # Count negative values
# negative_count <- sum(segments$HP_AR70$seg.mean < 0)
# 
# # Print results
# cat("Positive values:", positive_count, "\n")
# cat("Negative values:", negative_count, "\n")

# CNV.summaryplot
# CNV.summaryplot(xf, threshold = 0.15)
pdf("CNV_summaryplot_C1_FL.pdf", height = 5, width = 10)
CNV.summaryplot(xf, main = "C1_FL", res = 500)
dev.off()
