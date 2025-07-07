# R Packages
# ==========
library(EpiDISH)
library(dplyr)
library(stringr)
library(pheatmap)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(ggplot2)
library(reshape2)
library(tidyr)

# Change directory
setwd("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/03_EpiDISH/450k/02_after_purity/")

# Read in methylation data
bval <- readRDS("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/01_InfiniumPurify/450k/modified_InfiniumPurify/Enmix_quantile_with_QCinfo_Harman_All/beta_combat_purified.rds")
dim(bval)

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
  dplyr::select(LY_FL_ID, ANN_ARBOR_STAGE, SEX, T_14_18, PRIM_TX_CAT, LDH_ELEVATED, HEMOGLOBIN_LESS_120, FLIPI_BINARY, AGE_AT_DIAGNOSIS, AGE_CAT, MORE_4_NODAL_SITES, CODE_PFS, CODE_TRANSF, CODE_OS, PFS, TTT, OS)

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

# Remove Blueprint samples before clustering
# Sample_Annotations <- Sample_Annotations %>% filter(TYPE == "FL" | TYPE == "FL_A" | TYPE == "FL_B" | TYPE == "DLBCL" | TYPE == "DLBCL_EGA" | TYPE == "DLBCL_EN" | TYPE == "DLBCL_GC" | TYPE == "DLBCL_nonGC" | TYPE == "t_DLBCL" | TYPE == "Benign_LN" | TYPE == "LN_NN_A" | TYPE == "LN_NN_B" | TYPE == "RLN") %>% data.frame()

# Calculate Samples by Type per Cluster
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'FL_A'] <- 'FL'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'FL_B'] <- 'FL'

# Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'DLBCL_ABC'] <- 'DLBCL'
# Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'DLBCL_GCB'] <- 'DLBCL'
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

dim(Sample_Annotations)
table(Sample_Annotations$TYPE)

# subset of merged methylation data
bval <- bval[, Sample_Annotations$SAMPLE_ID]
dim(bval)

# Write TYPE and Batch Information
df1 <- Sample_Annotations[c("SAMPLE_ID", "TYPE", "Batch")]
write.table(df1, file="TYPE_INFO.txt", sep="\t", quote = FALSE, row.names = FALSE)

gc()

# EpiDISH
# =======
data(centDHSbloodDMC.m)
dim(centDHSbloodDMC.m) # 333   7

out.l <- epidish(beta.m = bval,
                 ref.m = centDHSbloodDMC.m,
                 method = "RPC")

# the output list. estF is the matrix of estimated cell-type fractions
est <- out.l$estF

pdf("EpiDISH.pdf")
boxplot(est)
dev.off()

write.table(round(est*100,1), file='est.tsv', quote=FALSE, sep='\t', col.names = NA)

# ref is the reference centroid matrix used
dim(out.l$ref)

# dataREF is the subset of the input data matrix over the probes defined in the reference matrix
dim(out.l$dataREF)

# HEpiDISH
# ========
# Ref-1
data(centEpiFibIC.m)
dim(centEpiFibIC.m) # 716   3

# Ref-2
data(centBloodSub.m)
dim(centBloodSub.m) # 188   7

frac.m <- hepidish(beta.m = bval,
                   ref1.m = centEpiFibIC.m,
                   ref2.m = centBloodSub.m,
                   h.CT.idx = 3,
                   method = 'RPC')
frac.m

pdf("HEpiDISH.pdf")
boxplot(frac.m)
dev.off()

write.table(round(frac.m*100,1), file='frac.m.tsv', quote=FALSE, sep='\t', col.names = NA)

# Generate a plot
# ===============
est_df <- as.data.frame(est)
dim(est_df)

est_df$SAMPLE_ID <- rownames(est_df)

Sample_Annotations_Merged <- merge(x = Sample_Annotations,
                                   y = est_df,
                                   by.x = "SAMPLE_ID",
                                   by.y = "SAMPLE_ID",
                                   all.x = TRUE)

Sample_Annotations_Merged <- Sample_Annotations_Merged %>% dplyr::select(SAMPLE_ID, TYPE, B, NK, CD4T, CD8T, Mono, Neutro, Eosino)

# Reshape data: Convert wide format to long format
long_df <- Sample_Annotations_Merged %>%
  pivot_longer(cols = -c(SAMPLE_ID, TYPE),  # Keep SAMPLE_ID & TYPE
               names_to = "Cell_Type", 
               values_to = "Proportion")

# Define custom colors for improved visibility
custom_colors <- c(
  "Normal" = "#5086F3",
  "FL" = "#59A14F",
  "t_DLBCL" = "#F28E2B",
  "DLBCL" = "#E15759",
  "naiBC" = "#8C6D31", 
  "t_naiBC" = "#D1A1A6",
  "gcBC" = "#EDC949",
  "t_PC" = "#76B7B2",
  "memBC" = "#FF9DA7",
  "bm_PC" = "#6A4D8C"
)

# Define custom order for sample types and cell types
custom_order <- c("Normal", "FL", "t_DLBCL", "DLBCL", "naiBC", "t_naiBC", "gcBC", "t_PC", "memBC", "bm_PC")
long_df$TYPE <- factor(long_df$TYPE, levels = custom_order)

long_df$Cell_Type <- factor(long_df$Cell_Type, levels = c("B", "CD4T", "CD8T", "Mono", "NK", "Eosino", "Neutro"))

# Set up PDF output
pdf("CellType_Boxplot_After.pdf", width = 35, height = 12)  # Adjusted width and height for better readability

# Generate grouped boxplot using custom colors
ggplot(long_df, aes(x = Cell_Type, y = Proportion, fill = TYPE)) + 
  geom_boxplot(position = position_dodge(width = 0.9), size = 0.7) + 
  facet_grid(~Cell_Type, scales = "free_x", space = "free_x") +  # Facet to separate sample types
  scale_fill_manual(values = custom_colors, name = "Sample Type") + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) + 
  ggtitle("Cell Type Proportions by Sample Type") + 
  labs(x = "Cell Type", y = "Proportion") + 
  theme(
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5, margin = margin(b = 42)),
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 34, face = "bold", margin = margin(b = 30)),
    legend.key.size = unit(1.8, "cm"),
    legend.margin = margin(l = 40),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 30, face = "bold"),
    axis.text.y = element_text(size = 24, face = "plain"),
    axis.title.x = element_text(size = 34, face = "bold", margin = margin(t = 40)),
    axis.title.y = element_text(size = 34, face = "bold", margin = margin(r = 40)),
    plot.margin = margin(t = 40, r = 40, b = 40, l = 40),
    panel.spacing.x = unit(1, "lines"),
    strip.text = element_blank(),
    strip.background = element_blank(),
    panel.spacing = unit(2, "cm"),
    panel.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray85", linewidth = 0.5),
    panel.grid.minor.y = element_line(color = "gray95", linewidth = 0.3)
  ) + 
  guides(fill = guide_legend(byrow = TRUE))

# Close the PDF output
dev.off()

# Filter only FL and DLBCL within the B-cell category
filtered_data <- long_df %>% 
  filter(Cell_Type == "B" & TYPE %in% c("FL", "DLBCL"))
table(filtered_data)

# Perform Wilcoxon test
wilcox_test_result <- wilcox.test(Proportion ~ TYPE, data = filtered_data)

# Print p-value
wilcox_test_result$p.value #p-value: 1.791517e-16
