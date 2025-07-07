# R Packages
# ==========
library(dplyr)
library(ggplot2)
library(reshape)
library(ggpubr)
library(factoextra)
library(readxl)
library(FactoMineR)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(stringr)
library(tidyr)
library(RColorBrewer)
library(cluster)
library(pheatmap)
library(GGally)

# Change directory
setwd("/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/01_Clustering/850k/02_hclust/")

# get the 850k annotation data
ann850k <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
EPIC_Annotations <- data.frame(ann850k)
remove(ann850k)

# Define CpGs of interest
CpG_islands <- EPIC_Annotations %>% filter(Relation_to_Island == "Island") %>% data.frame() %>% .$Name
OpenSea <- EPIC_Annotations %>% filter(Relation_to_Island == "OpenSea") %>% data.frame() %>% .$Name

# Read in methylation data
bval <- readRDS("/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/850k/preprocessQuantile/rds/beta_combat.rds")
mval <- readRDS("/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/850k/preprocessQuantile/rds/mval_combat.rds")
dim(bval)
dim(mval)

# Read in sample annotations
Sample_Annotations <- read.table(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/ANN/20230518_sample_annotations.tsv",
                                 sep = "\t", header = TRUE, na.strings = c("","NA")) %>%
  subset(TYPE!="TFL") %>%
  subset(TIME_POINT!="T2") %>%
  filter(SAMPLE_ID %in% colnames(bval))

# Read in ClusterAIC information
flexmix_clust <- read.table(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/ClusterAIC/2023-05-06_14_GMM_Cluster_Labels_flexmix_clusters.txt",
                            sep = "\t", header = TRUE, na.strings = c("","NA")) %>%
  dplyr::select(SAMPLE_ID, ClusterAIC)
flexmix_clust <- flexmix_clust %>% filter(SAMPLE_ID %in% colnames(bval))

# Read in sample sheet with batch information
targets <- read.csv(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/SampleSheet/SampleSheet.csv", 
                    skip = 1, header = TRUE, na.strings = c("","NA"))

targets$ID <- gsub(" ", "", paste(targets$Sample_Name,".",targets$Sentrix_Position))

# Remove repeated samples to obtain correct batch information
rep_rm_samples <- c("LY_FL_158_T1.R01C01", "LY_FL_479_T1.R02C01", "LY_FL_488_T1.R06C01",
                    "LY_FL_498_T1.R06C01", "LY_FL_523_T1.R01C01", "LY_FL_524_T1.R02C01",
                    "LY_FL_525_T1.R03C01", "LY_FL_527_T1.R04C01", "LY_FL_529_T1.R01C01",
                    "LY_FL_1156_T1.R06C01", "LY_FL_571_T1.R01C01", "HCT116_DKO_methylated.R08C01",
                    "HCT116_DKO_methylated.R04C01", "LY_FL_311_T1.R04C01", "LY_FL_159_T1.R08C01",
                    "LY_FL_535_T1.R08C01", "LY_FL_536_T1.R01C01")

targets <- targets[ ! targets$ID %in% rep_rm_samples, ]
targets$Sample_Name <- str_replace(targets$Sample_Name, "LY_FL_159_T1_rep", "LY_FL_159_T1")
targets <- targets %>% dplyr::select(Sample_Name, Batch)
colnames(targets)[1]  <- "SAMPLE_ID"
targets <- targets %>% filter(SAMPLE_ID %in% colnames(bval))

# Read in clinical data
clinical_data <- read.csv(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/CLN/20230526_clinical_data.csv", 
                          header = TRUE, na.strings = c("","NA")) %>%
  dplyr::select(LY_FL_ID, ANN_ARBOR_STAGE)

# Merge data frames
df <- merge(x = targets,
            y = flexmix_clust,
            by.x = "SAMPLE_ID",
            by.y = "SAMPLE_ID",
            all.x = TRUE,
            all.y = TRUE)

df1 <- merge(x = Sample_Annotations,
             y = df,
             by.x = "SAMPLE_ID",
             by.y = "SAMPLE_ID",
             all.x = TRUE)

df2 <- merge(x = df1,
             y = clinical_data,
             by.x = "LY_FL_ID",
             by.y = "LY_FL_ID",
             all.x = TRUE)

remove(Sample_Annotations)
Sample_Annotations <- df2
remove(df, df1, df2)

# Create a additional column with Clinical Stage, RLN, DLBCL
Sample_Annotations$STAGE_TYPE <- Sample_Annotations$STAGE
Sample_Annotations <- Sample_Annotations %>% mutate(STAGE_TYPE = coalesce(STAGE_TYPE, TYPE))
Sample_Annotations$STAGE_TYPE[Sample_Annotations$STAGE_TYPE == "FL"] <- NA

# Create a additional column with Ann Arbor Stage, RLN, DLBCL
Sample_Annotations$ANN_ARBOR_STAGE_TYPE <- Sample_Annotations$ANN_ARBOR_STAGE
Sample_Annotations <- Sample_Annotations %>% mutate(ANN_ARBOR_STAGE_TYPE = coalesce(as.character(ANN_ARBOR_STAGE_TYPE), TYPE))
Sample_Annotations$ANN_ARBOR_STAGE_TYPE[Sample_Annotations$ANN_ARBOR_STAGE_TYPE == "FL"] <- NA

# Create a additional column with ClusterAIC, RLN, DLBCL
Sample_Annotations$ClusterAIC_TYPE <- Sample_Annotations$ClusterAIC
Sample_Annotations <- Sample_Annotations %>% mutate(ClusterAIC_TYPE = coalesce(ClusterAIC_TYPE, TYPE))
Sample_Annotations$ClusterAIC_TYPE[Sample_Annotations$ClusterAIC_TYPE == "FL"] <- NA

row.names(Sample_Annotations) <- Sample_Annotations$SAMPLE_ID

# Subset of merged methylation data
bval <- bval[, Sample_Annotations$SAMPLE_ID]
mval <- mval[, Sample_Annotations$SAMPLE_ID]
dim(bval)
dim(mval)

bval <- bval[,1:328] # subset FL samples
mval <- mval[,1:328] # subset FL samples
dim(bval)
dim(mval)

# subset matrix based on most variable probes
set.seed(1234)
probe_num = "20000"

# using `standard deviation`
bval_sd <- apply(bval, 1, sd)
bval_sub <- bval[match(names(sort(bval_sd, decreasing = TRUE)[1:probe_num]), rownames(bval)), ]
mval_sub <- mval[intersect(rownames(bval_sub), rownames(mval)) ,]

df <- t(bval_sub)
dim(df)


## Hierarchical Clustering

# Find the Linkage Method to Use
# ------------------------------
# define linkage methods
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")
m

# function to compute agglomerative coefficient
ac <- function(x) {
  agnes(df, method = x)$ac
}

# calculate agglomerative coefficient for each clustering linkage method
sapply(m, ac)

# perform hierarchical clustering using Ward's minimum variance
clust <- agnes(df, method = "ward")

# produce dendrogram
pdf("FL_epitypes.pdf")
pltree(clust, cex = 0.6, hang = -1, main = "Dendrogram")
dev.off()


# Determine the Optimal Number of Clusters
# ----------------------------------------
# Find the Optimal Number of Clusters - Elbow Method
pdf("fviz_nbclust_elbow_method.pdf")
fviz_nbclust(df, hcut, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)
dev.off()

# Find the Optimal Number of Clusters - silhouette
pdf("fviz_nbclust_silhouette_score.pdf")
fviz_nbclust(df, hcut, method = "silhouette")
dev.off()


# Apply Cluster Labels to Original Dataset
# ----------------------------------------
# compute distance matrix
d <- dist(df, method = "euclidean")

# perform hierarchical clustering using Ward's method
final_clust <- hclust(d, method = "ward.D2")

# cut the dendrogram into 2 clusters
groups <- cutree(final_clust, k = 4)

# find number of observations in each cluster
table(groups)

kclust_epitypes <- data.frame(epitypes = groups)

kclust_epitypes$epitypes[kclust_epitypes$epitypes == '1'] <- 'C1'
kclust_epitypes$epitypes[kclust_epitypes$epitypes == '2'] <- 'C2'
kclust_epitypes$epitypes[kclust_epitypes$epitypes == '3'] <- 'C3'
kclust_epitypes$epitypes[kclust_epitypes$epitypes == '4'] <- 'C4'

table(kclust_epitypes$epitypes)

kclust_epitypes$SAMPLE_ID <- rownames(kclust_epitypes)
kclust_epitypes <- kclust_epitypes %>% relocate(SAMPLE_ID, .before = epitypes)

write.table(kclust_epitypes, file = "FL_epitypes.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

# Visualize silhouette information
set.seed(101)

df_dm <- daisy(df)

pdf("silhouette_info.pdf")
plot(silhouette(groups, df_dm), col = c("#f58231","#4363d8"), border=NA)
dev.off()
