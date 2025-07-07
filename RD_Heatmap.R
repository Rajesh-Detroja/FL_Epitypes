# R Packages
# ==========
library(dplyr)
library(pheatmap)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(stringr)
library(tidyr)
library(RColorBrewer)

# Change directory
setwd("/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/01_Clustering/850k/02_hclust/Heatmap/")

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

# Read in methylation clustering information
epitypes <- read.table(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/01_Clustering/850k/04_ConsensusClusterPlus/FL_epitypes.tsv",
                       sep = "\t", header = TRUE, na.strings = c("","NA"))

epitypes <- epitypes %>% filter(SAMPLE_ID %in% colnames(bval))

# Read in EZH2_STATUS
EZH2_STATUS <- read.csv(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/EZH2/2023-03-21_CAPSEQ_EZH2_variant_status_mat.tsv",
                        sep = "\t", header = TRUE, na.strings=c("","NA"))

# Read in IGH_STATUS
IGH_STATUS <- read.csv(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/IGH/2023-06-06_clonotypes_igh_fl_vdjgrouped_dominant_edited.txt",
                       sep = "\t", header = TRUE, na.strings=c("","NA"))

colnames(IGH_STATUS)[1] ="SAMPLE_ID"
colnames(IGH_STATUS)[2] ="IGH_STATUS"

IGH_STATUS <- IGH_STATUS %>% dplyr::select(SAMPLE_ID, IGH_STATUS)

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

df3 <- merge(x = df2,
             y = epitypes,
             by.x = "SAMPLE_ID",
             by.y = "SAMPLE_ID",
             all.x = TRUE)

df4 <- merge(x = df3,
             y = EZH2_STATUS,
             by.x = "SAMPLE_ID",
             by.y = "SAMPLE_ID",
             all.x = TRUE)

df5 <- merge(x = df4,
             y = IGH_STATUS,
             by.x = "SAMPLE_ID",
             by.y = "SAMPLE_ID",
             all.x = TRUE)

remove(Sample_Annotations)
Sample_Annotations <- df5
remove(df, df1, df2, df3, df4, df5)

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

# Create a additional column with epitypes, RLN, DLBCL
Sample_Annotations$epitypes_TYPE <- Sample_Annotations$epitypes
Sample_Annotations <- Sample_Annotations %>% mutate(epitypes_TYPE = coalesce(epitypes_TYPE, TYPE))
Sample_Annotations$epitypes_TYPE[Sample_Annotations$epitypes_TYPE == "FL"] <- NA

row.names(Sample_Annotations) <- Sample_Annotations$SAMPLE_ID

Sample_Annotations <- Sample_Annotations %>% filter(TYPE == "FL")

# Subset of merged methylation data
bval <- bval[, Sample_Annotations$SAMPLE_ID]
mval <- mval[, Sample_Annotations$SAMPLE_ID]
dim(bval)
dim(mval)

# Subset of methylation data for FL samples
#bval_FL <- bval[, c(11:338)]
#mval_FL <- mval[, c(11:338)]
#dim(bval_FL)
#dim(mval_FL)

# subset matrix based on random probes
#set.seed(12345)

#bval_sub_df <- bval[sample(nrow(bval), 20000), ]
#bval_sub_df <- sample_n(as.data.frame(bval), 20000)
#dim(bval_sub_df)

# subset matrix based on most variable probes
#set.seed(1234)
probe_num = "20000"

# using `standard deviation`
#bval_sd <- apply(bval_FL, 1, sd)
#bval_sub <- bval[match(names(sort(bval_sd, decreasing = TRUE)[1:probe_num]), rownames(bval)), ]
#mval_sub <- mval[intersect(rownames(bval_sub), rownames(mval)) ,]

# read FL 20k probes
sd_probes_20k <- read.table(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/01_Clustering/850k/04_ConsensusClusterPlus/sd_probes_20k.tsv", 
                            sep = "\t", header = TRUE, na.strings = c("","NA"))

bval_sub <- bval[intersect(sd_probes_20k$PROBE_ID, rownames(bval)),]
mval_sub <- mval[intersect(sd_probes_20k$PROBE_ID, rownames(mval)),]
dim(bval_sub)
dim(mval_sub)

# bval - All CpG
bval_sub_t <- t(bval_sub)
bval_sub_t <- merge(bval_sub_t, Sample_Annotations, by=0)
names(bval_sub_t)[1] <- ""
bval_sub_t <- data.frame(bval_sub_t[,-1], row.names=bval_sub_t[,1])
bval_sub_t_df <- bval_sub_t[1:probe_num]

# mval - All CpG
mval_sub_t <- t(mval_sub)
mval_sub_t <- merge(mval_sub_t, Sample_Annotations, by=0)
names(mval_sub_t)[1] <- ""
mval_sub_t <- data.frame(mval_sub_t[,-1], row.names=mval_sub_t[,1])
mval_sub_t_df <- mval_sub_t[1:probe_num]

# bval - CpG Islands
bval_sub_CpG_islands <- bval_sub[intersect(CpG_islands, row.names(bval_sub)), ]
bval_sub_CpG_islands_t <- t(bval_sub_CpG_islands)
bval_sub_CpG_islands_t <- merge(bval_sub_CpG_islands_t, Sample_Annotations, by=0)
names(bval_sub_CpG_islands_t)[1] <- ""
bval_sub_CpG_islands_t <- data.frame(bval_sub_CpG_islands_t[,-1], row.names=bval_sub_CpG_islands_t[,1])
bval_sub_CpG_islands_t_df <- bval_sub_CpG_islands_t[1:ncol(t(bval_sub_CpG_islands))]

# mval - CpG Islands
mval_sub_CpG_islands <- mval_sub[intersect(CpG_islands, row.names(mval_sub)), ]
mval_sub_CpG_islands_t <- t(mval_sub_CpG_islands)
mval_sub_CpG_islands_t <- merge(mval_sub_CpG_islands_t, Sample_Annotations, by=0)
names(mval_sub_CpG_islands_t)[1] <- ""
mval_sub_CpG_islands_t <- data.frame(mval_sub_CpG_islands_t[,-1], row.names=mval_sub_CpG_islands_t[,1])
mval_sub_CpG_islands_t_df <- mval_sub_CpG_islands_t[1:ncol(t(mval_sub_CpG_islands))]

# bval - OpenSea
bval_sub_OpenSea <- bval_sub[intersect(OpenSea, row.names(bval_sub)), ]
bval_sub_OpenSea_t <- t(bval_sub_OpenSea)
bval_sub_OpenSea_t <- merge(bval_sub_OpenSea_t, Sample_Annotations, by=0)
names(bval_sub_OpenSea_t)[1] <- ""
bval_sub_OpenSea_t <- data.frame(bval_sub_OpenSea_t[,-1], row.names=bval_sub_OpenSea_t[,1])
bval_sub_OpenSea_t_df <- bval_sub_OpenSea_t[1:ncol(t(bval_sub_OpenSea))]

# mval - OpenSea
mval_sub_OpenSea <- mval_sub[intersect(OpenSea, row.names(mval_sub)), ]
mval_sub_OpenSea_t <- t(mval_sub_OpenSea)
mval_sub_OpenSea_t <- merge(mval_sub_OpenSea_t, Sample_Annotations, by=0)
names(mval_sub_OpenSea_t)[1] <- ""
mval_sub_OpenSea_t <- data.frame(mval_sub_OpenSea_t[,-1], row.names=mval_sub_OpenSea_t[,1])
mval_sub_OpenSea_t_df <- mval_sub_OpenSea_t[1:ncol(t(mval_sub_OpenSea))]

## epitypes
## ========
# Beta-Value Clustering - All CpGs
t <- bval_sub_t %>% drop_na(epitypes_TYPE)
df <- bval_sub_t_df[intersect(row.names(t), row.names(bval_sub_t_df)), ]

ann = data.frame("Epitypes" = t$epitypes_TYPE,
                 "Stage" = t$STAGE,
                 "Ann_Arbor_Stage" = as.character(t$ANN_ARBOR_STAGE),
                 "Cluster" = t$ClusterAIC,
                 "EZH2_Status" = t$EZH2_STATUS
                 )

rownames(ann) = rownames(df)

ann_colors = list(
  Epitypes = c("RLN"="#911eb4", "C1"="#F8766D", "C2"="#7CAE00", "C3"="#00BFC4", "C4"="#C77CFF", "DLBCL"="#c7284f"),
  Stage = c("LIMITED"="#028796", "ADVANCED"="#7d2233"),
  Ann_Arbor_Stage = c("1"="#3dd1d1", "2"="#028796", "3"="#e0abab", "4"="#7d2233"),
  Cluster = c(CS="#F8766D", TT="#A3A500", GM="#00BF7D", Q="#00B0F6", AR="#E76BF3"),
  EZH2_Status = c(EZH2_Mut_Pos= "#7d2233", EZH2_Mut_Neg = "#3dd1d1")
)

png(filename = "bval-heatmap-Type-All-CpGs.png")
pheatmap(t(df),
         scale="none",
         clustering_method = "ward.D2",
         main="Beta Value - All CpGs",
         annotation_col = ann,
         annotation_colors = ann_colors,
         # color = rev(brewer.pal(n = 10, name = 'RdBu')),
         color=colorRampPalette(c("navy", "white", "red"))(10),
         cutree_rows = 4,
         cutree_cols = 4,
         fontsize = 8,
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()


#hcl = hclust(dist(t(df)))
#plot(hcl, labels = FALSE, hang= -1)
#rect.hclust(hcl, h = 12, border = "red")

#clu.k3 = cutree(hcl, k = 4) # cut tree so that there are 3 clusters
#clu.h12 = cutree(hcl, h = 12) # cut tree/dendrogram from height 12
#table(clu.k3) # number of samples for each cluster
#table(clu.h12) # number of samples for each cluster

# M-Value Clustering - Clinical Stage - All CpGs
t <- mval_sub_t %>% drop_na(TYPE)
df <- mval_sub_t_df[intersect(row.names(t), row.names(mval_sub_t_df)), ]

ann = data.frame("Type"=t$TYPE,
                 # "Batch"=t$Batch,
                 "Stage"=t$STAGE,
                 "Ann_Arbor_Stage"=as.character(t$ANN_ARBOR_STAGE),
                 "Cluster"=t$ClusterAIC,
                 "EZH2_Status"=t$EZH2_STATUS
)

rownames(ann) = rownames(df)

ann_colors = list(
  Type = c("RLN"="#8b25fa", "FL"="#6dcc04", "DLBCL"="#fc1417", "naiBC"="#B15928", "t_naiBC"="#fc51a8", "gcBC"="#0303a8", "memBC"="#1B9E77", "t_PC"="#383838", "bm_PC"="#FF7F00"),
  # Batch = c("Blueprint"="#0303a8", "PM-OICR-TGL"="#fc51a8", "PM-Genomics-Centre"="#6dcc04", "Genome-Quebec"="#FF7F00"),
  Stage = c("LIMITED"="#028796", "ADVANCED"="#7d2233"),
  Ann_Arbor_Stage = c("1"="#3dd1d1", "2"="#028796", "3"="#e0abab", "4"="#7d2233"),
  Cluster = c(CS="#F8766D", TT="#A3A500", GM="#00BF7D", Q="#00B0F6", AR="#E76BF3"),
  EZH2_Status = c(EZH2_Mut_Pos= "#7d2233", EZH2_Mut_Neg = "#3dd1d1")
)

pdf("mval-heatmap-Type-All-CpGs.pdf")
pheatmap(t(df),
         scale="none",
         clustering_method = "ward.D2",
         main="M Value - All CpGs",
         annotation_col = ann,
         annotation_colors = ann_colors,
         # color = rev(brewer.pal(n = 10, name = 'RdBu')),
         color=colorRampPalette(c("navy", "white", "red"))(10),
         cutree_rows = 2,
         cutree_cols = 3,
         fontsize = 8,
         show_rownames = FALSE,
         show_colnames = FALSE)
graphics.off()


# Beta-Value Clustering - Clinical Stage - CpG Islands
t <- bval_sub_CpG_islands_t %>% drop_na(TYPE)
df <- bval_sub_CpG_islands_t_df[intersect(row.names(t), row.names(bval_sub_CpG_islands_t_df)), ]

ann = data.frame("Type"=t$TYPE,
                 # "Batch"=t$Batch,
                 "Stage"=t$STAGE,
                 "Ann_Arbor_Stage"=as.character(t$ANN_ARBOR_STAGE),
                 "Cluster"=t$ClusterAIC,
                 "EZH2_Status"=t$EZH2_STATUS
)

rownames(ann) = rownames(df)

ann_colors = list(
  Type = c("RLN"="#8b25fa", "FL"="#6dcc04", "DLBCL"="#fc1417", "naiBC"="#B15928", "t_naiBC"="#fc51a8", "gcBC"="#0303a8", "memBC"="#1B9E77", "t_PC"="#383838", "bm_PC"="#FF7F00"),
  # Batch = c("Blueprint"="#0303a8", "PM-OICR-TGL"="#fc51a8", "PM-Genomics-Centre"="#6dcc04", "Genome-Quebec"="#FF7F00"),
  Stage = c("LIMITED"="#028796", "ADVANCED"="#7d2233"),
  Ann_Arbor_Stage = c("1"="#3dd1d1", "2"="#028796", "3"="#e0abab", "4"="#7d2233"),
  Cluster = c(CS="#F8766D", TT="#A3A500", GM="#00BF7D", Q="#00B0F6", AR="#E76BF3"),
  EZH2_Status = c(EZH2_Mut_Pos= "#7d2233", EZH2_Mut_Neg = "#3dd1d1")
)

pdf("bval-heatmap-Type-CpG-Islands.pdf")
pheatmap(t(df),
         scale="none",
         clustering_method = "ward.D2",
         main="Beta Value - CpG Islands",
         annotation_col = ann,
         annotation_colors = ann_colors,
         # color = rev(brewer.pal(n = 10, name = 'RdBu')),
         color=colorRampPalette(c("navy", "white", "red"))(10),
         cutree_rows = 2,
         cutree_cols = 3,
         fontsize = 8,
         show_rownames = FALSE,
         show_colnames = FALSE)
graphics.off()


# M-Value Clustering - Clinical Stage - CpG Islands
t <- mval_sub_CpG_islands_t %>% drop_na(TYPE)
df <- mval_sub_CpG_islands_t_df[intersect(row.names(t), row.names(mval_sub_CpG_islands_t_df)), ]

ann = data.frame("Type"=t$TYPE,
                 # "Batch"=t$Batch,
                 "Stage"=t$STAGE,
                 "Ann_Arbor_Stage"=as.character(t$ANN_ARBOR_STAGE),
                 "Cluster"=t$ClusterAIC,
                 "EZH2_Status"=t$EZH2_STATUS
)

rownames(ann) = rownames(df)

ann_colors = list(
  Type = c("RLN"="#8b25fa", "FL"="#6dcc04", "DLBCL"="#fc1417", "naiBC"="#B15928", "t_naiBC"="#fc51a8", "gcBC"="#0303a8", "memBC"="#1B9E77", "t_PC"="#383838", "bm_PC"="#FF7F00"),
  # Batch = c("Blueprint"="#0303a8", "PM-OICR-TGL"="#fc51a8", "PM-Genomics-Centre"="#6dcc04", "Genome-Quebec"="#FF7F00"),
  Stage = c("LIMITED"="#028796", "ADVANCED"="#7d2233"),
  Ann_Arbor_Stage = c("1"="#3dd1d1", "2"="#028796", "3"="#e0abab", "4"="#7d2233"),
  Cluster = c(CS="#F8766D", TT="#A3A500", GM="#00BF7D", Q="#00B0F6", AR="#E76BF3"),
  EZH2_Status = c(EZH2_Mut_Pos= "#7d2233", EZH2_Mut_Neg = "#3dd1d1")
)

pdf("mval-heatmap-Type-CpG-Islands.pdf")
pheatmap(t(df),
         scale="none",
         clustering_method = "ward.D2",
         main="M Value - CpG Islands",
         annotation_col = ann,
         annotation_colors = ann_colors,
         # color = rev(brewer.pal(n = 10, name = 'RdBu')),
         color=colorRampPalette(c("navy", "white", "red"))(10),
         cutree_rows = 2,
         cutree_cols = 3,
         fontsize = 8,
         show_rownames = FALSE,
         show_colnames = FALSE)
graphics.off()


# Beta-Value Clustering - Clinical Stage - OpenSea
t <- bval_sub_OpenSea_t %>% drop_na(TYPE)
df <- bval_sub_OpenSea_t_df[intersect(row.names(t), row.names(bval_sub_OpenSea_t_df)), ]

ann = data.frame("Type"=t$TYPE,
                 # "Batch"=t$Batch,
                 "Stage"=t$STAGE,
                 "Ann_Arbor_Stage"=as.character(t$ANN_ARBOR_STAGE),
                 "Cluster"=t$ClusterAIC,
                 "EZH2_Status"=t$EZH2_STATUS
)

rownames(ann) = rownames(df)

ann_colors = list(
  Type = c("RLN"="#8b25fa", "FL"="#6dcc04", "DLBCL"="#fc1417", "naiBC"="#B15928", "t_naiBC"="#fc51a8", "gcBC"="#0303a8", "memBC"="#1B9E77", "t_PC"="#383838", "bm_PC"="#FF7F00"),
  # Batch = c("Blueprint"="#0303a8", "PM-OICR-TGL"="#fc51a8", "PM-Genomics-Centre"="#6dcc04", "Genome-Quebec"="#FF7F00"),
  Stage = c("LIMITED"="#028796", "ADVANCED"="#7d2233"),
  Ann_Arbor_Stage = c("1"="#3dd1d1", "2"="#028796", "3"="#e0abab", "4"="#7d2233"),
  Cluster = c(CS="#F8766D", TT="#A3A500", GM="#00BF7D", Q="#00B0F6", AR="#E76BF3"),
  EZH2_Status = c(EZH2_Mut_Pos= "#7d2233", EZH2_Mut_Neg = "#3dd1d1")
)

pdf("bval-heatmap-Type-OpenSea.pdf")
pheatmap(t(df),
         scale="none",
         clustering_method = "ward.D2",
         main="Beta Value - OpenSea",
         annotation_col = ann,
         annotation_colors = ann_colors,
         # color = rev(brewer.pal(n = 10, name = 'RdBu')),
         color=colorRampPalette(c("navy", "white", "red"))(10),
         cutree_rows = 2,
         cutree_cols = 3,
         fontsize = 8,
         show_rownames = FALSE,
         show_colnames = FALSE)
graphics.off()


# M-Value Clustering - Clinical Stage - OpenSea
t <- mval_sub_OpenSea_t %>% drop_na(TYPE)
df <- mval_sub_OpenSea_t_df[intersect(row.names(t), row.names(mval_sub_OpenSea_t_df)), ]

ann = data.frame("Type"=t$TYPE,
                 # "Batch"=t$Batch,
                 "Stage"=t$STAGE,
                 "Ann_Arbor_Stage"=as.character(t$ANN_ARBOR_STAGE),
                 "Cluster"=t$ClusterAIC,
                 "EZH2_Status"=t$EZH2_STATUS
)

rownames(ann) = rownames(df)

ann_colors = list(
  Type = c("RLN"="#8b25fa", "FL"="#6dcc04", "DLBCL"="#fc1417", "naiBC"="#B15928", "t_naiBC"="#fc51a8", "gcBC"="#0303a8", "memBC"="#1B9E77", "t_PC"="#383838", "bm_PC"="#FF7F00"),
  # Batch = c("Blueprint"="#0303a8", "PM-OICR-TGL"="#fc51a8", "PM-Genomics-Centre"="#6dcc04", "Genome-Quebec"="#FF7F00"),
  Stage = c("LIMITED"="#028796", "ADVANCED"="#7d2233"),
  Ann_Arbor_Stage = c("1"="#3dd1d1", "2"="#028796", "3"="#e0abab", "4"="#7d2233"),
  Cluster = c(CS="#F8766D", TT="#A3A500", GM="#00BF7D", Q="#00B0F6", AR="#E76BF3"),
  EZH2_Status = c(EZH2_Mut_Pos= "#7d2233", EZH2_Mut_Neg = "#3dd1d1")
)

pdf("mval-heatmap-Type-OpenSea.pdf")
pheatmap(t(df),
         scale="none",
         clustering_method = "ward.D2",
         main="M Value - OpenSea",
         annotation_col = ann,
         annotation_colors = ann_colors,
         # color = rev(brewer.pal(n = 10, name = 'RdBu')),
         color=colorRampPalette(c("navy", "white", "red"))(10),
         cutree_rows = 2,
         cutree_cols = 3,
         fontsize = 8,
         show_rownames = FALSE,
         show_colnames = FALSE)
graphics.off()
