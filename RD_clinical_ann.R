# R Packages
# ==========
library(dplyr)
library(ggplot2)
library(reshape)
library(ggpubr)
library(stringr)
library(tidyr)
library(RColorBrewer)

# Change directory
setwd("/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/04_Clinical_Analysis/450k/04_ConsensusClusterPlus/Enmix_quantile_with_QCinfo_Harman/modified_InfiniumPurify/")

# Read in methylation data
bval <- readRDS("/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/02_InfiniumPurify/450k/Enmix_quantile_with_QCinfo_Harman/cluster_run/20_normal/modified_InfiniumPurify/beta_combat_purified.rds")
dim(bval)

# Read in sample annotations
Sample_Annotations <- read.table(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/ANN/20231031_sample_annotations.tsv", sep = "\t", header = TRUE, na.strings = c("","NA")) %>%
                      subset(TYPE!="TFL") %>%
                      subset(TIME_POINT!="T2") %>%
                      filter(SAMPLE_ID %in% colnames(bval))

# Read in ClusterAIC information
flexmix_clust <- read.table(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/ClusterAIC/2023-05-06_14_GMM_Cluster_Labels_flexmix_clusters.txt",
                            sep = "\t", header = TRUE, na.strings = c("","NA")) %>%
                            dplyr::select(SAMPLE_ID, ClusterAIC)

flexmix_clust <- flexmix_clust %>% filter(SAMPLE_ID %in% colnames(bval))

# Read in sample sheet with batch information
targets <- read.csv(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/SampleSheet/SampleSheet.csv", skip = 1, header = TRUE, na.strings = c("","NA"))
targets$ID <- gsub(" ", "", paste(targets$Sample_Name,".",targets$Sentrix_Position))

# Remove repeated samples to obtain correct batch information
rep_rm_samples <- c("LY_FL_158_T1.R01C01", "LY_FL_479_T1.R02C01", "LY_FL_488_T1.R06C01",
                    "LY_FL_498_T1.R06C01", "LY_FL_523_T1.R01C01", "LY_FL_524_T1.R02C01",
                    "LY_FL_525_T1.R03C01", "LY_FL_527_T1.R04C01", "LY_FL_529_T1.R01C01",
                    "LY_FL_1156_T1.R06C01", "LY_FL_571_T1.R01C01", "HCT116_DKO_methylated.R08C01",
                    "HCT116_DKO_methylated.R04C01", "LY_FL_311_T1.R04C01", "LY_FL_159_T1.R08C01",
                    "LY_FL_535_T1.R08C01", "LY_FL_536_T1.R01C01", "LY_FL_127_T2.R01C01") # check details on "LY_FL_527_T1.R01C01"

targets <- targets[!targets$ID %in% rep_rm_samples,]
targets$Sample_Name <- str_replace(targets$Sample_Name, "LY_FL_159_T1_rep", "LY_FL_159_T1")
targets <- targets %>% dplyr::select(Sample_Name, Batch)
colnames(targets)[1]  <- "SAMPLE_ID"
targets <- targets %>% filter(SAMPLE_ID %in% colnames(bval))

# Read in the Benign_LN sample annotation
Benign_LN_sample_ann <- read.table(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/Benign_LN/SampleSheet/Benign_LN_sample_annotations.tsv",
                                   sep = "\t",
                                   header = TRUE,
                                   na.strings=c("","NA")) %>%
                                   dplyr::select(SAMPLE_ID, TYPE)

Benign_LN_sample_ann <- Benign_LN_sample_ann %>% filter(SAMPLE_ID %in% colnames(bval))
Benign_LN_sample_ann$Batch <- "Benign_LN"

# Read in the NIH sample annotation
NIH_sample_ann <- read.table(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/NIH/GSE237299/SampleSheet/GSE237299_sample_annotations.txt",
                             sep = "\t",
                             header = TRUE,
                             na.strings=c("","NA")) %>%
                             dplyr::select(SAMPLE_ID, TYPE)

NIH_sample_ann$SAMPLE_ID <- gsub('-', '_', NIH_sample_ann$SAMPLE_ID)
NIH_sample_ann <- NIH_sample_ann %>% filter(SAMPLE_ID %in% colnames(bval))
NIH_sample_ann$Batch <- "GSE237299"

# Read in the EGA sample annotation
EGA_sample_ann <- read.table(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/EGA/SampleSheet/Sample_Annotations_EGAD00010001974.txt", sep = "\t",
                             header = TRUE,
                             na.strings=c("","NA")) %>%
                             dplyr::select(SAMPLE_ID, TYPE)

EGA_sample_ann$SAMPLE_ID <- gsub('-', '_', EGA_sample_ann$SAMPLE_ID)
EGA_sample_ann$SAMPLE_ID <- gsub('\\.', '_', EGA_sample_ann$SAMPLE_ID)

EGA_sample_ann <- EGA_sample_ann %>% filter(SAMPLE_ID %in% colnames(bval))
EGA_sample_ann$Batch <- "EGAD00010001974"

# Read in blueprint sample annotation
Blueprint_sample_ann <- read.csv(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/blueprint/Sample_Annotation_EGAS00001001196.txt",
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
Blueprint_sample_ann$SEX[Blueprint_sample_ann$SEX == 'Mixed'] <- NA
Blueprint_sample_ann <- Blueprint_sample_ann %>% filter(SAMPLE_ID %in% colnames(bval))

# Read in clinical data
clinical_data <- read.csv(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/CLN/20231109_clinical_data.csv", 
                          header = TRUE, na.strings = c("","NA")) %>%
                          dplyr::select(LY_FL_ID, ANN_ARBOR_STAGE, SEX, T_14_18, PRIM_TX_CAT, CODE_PFS, CODE_TRANSF, CODE_OS, PFS, TTT, OS)

# Read in methylation clustering information
epitypes <- read.table(file = "/home/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/03_Clustering/450k/04_ConsensusClusterPlus/Enmix_quantile_with_QCinfo_Harman/modified_InfiniumPurify/00_clustering_runs/RD_ConsensusClusterPlus_01_Rep_500/FL_epitypes.tsv",
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
             y = EZH2_STATUS,
             by.x = "SAMPLE_ID",
             by.y = "SAMPLE_ID",
             all.x = TRUE)

df4 <- merge(x = df3,
             y = IGH_STATUS,
             by.x = "SAMPLE_ID",
             by.y = "SAMPLE_ID",
             all.x = TRUE)

df5 <- bind_rows(df4, Benign_LN_sample_ann, NIH_sample_ann, EGA_sample_ann, Blueprint_sample_ann)

df6 <- merge(x = df5,
             y = epitypes,
             by.x = "SAMPLE_ID",
             by.y = "SAMPLE_ID",
             all.x = TRUE)

remove(Sample_Annotations)
Sample_Annotations <- df6
remove(df, df1, df2, df3, df4, df5, df6)

# Create a additional column with Ann Arbor Stage, RLN, DLBCL
#Sample_Annotations$ANN_ARBOR_STAGE_TYPE <- Sample_Annotations$ANN_ARBOR_STAGE
#Sample_Annotations <- Sample_Annotations %>% mutate(ANN_ARBOR_STAGE_TYPE = coalesce(as.character(ANN_ARBOR_STAGE_TYPE), TYPE))
#Sample_Annotations$ANN_ARBOR_STAGE_TYPE[Sample_Annotations$ANN_ARBOR_STAGE_TYPE == "FL"] <- NA

# Create a additional column with ClusterAIC, RLN, DLBCL
#Sample_Annotations$ClusterAIC_TYPE <- Sample_Annotations$ClusterAIC
#Sample_Annotations <- Sample_Annotations %>% mutate(ClusterAIC_TYPE = coalesce(ClusterAIC_TYPE, TYPE))
#Sample_Annotations$ClusterAIC_TYPE[Sample_Annotations$ClusterAIC_TYPE == "FL"] <- NA

# Create a additional column with epitypes, RLN, DLBCL
#Sample_Annotations$epitypes_TYPE <- Sample_Annotations$epitypes
#Sample_Annotations <- Sample_Annotations %>% mutate(epitypes_TYPE = coalesce(epitypes_TYPE, TYPE))
#Sample_Annotations$epitypes_TYPE[Sample_Annotations$epitypes_TYPE == "FL"] <- NA

row.names(Sample_Annotations) <- Sample_Annotations$SAMPLE_ID

#Sample_Annotations <- Sample_Annotations %>% filter(TYPE == "FL") %>% data.frame()
#Sample_Annotations <- Sample_Annotations %>% filter(TYPE == "FL" | TYPE == "FL_A" | TYPE == "FL_B") %>% data.frame()
Sample_Annotations <- Sample_Annotations %>% filter(TYPE == "FL" | TYPE == "FL_A" | TYPE == "FL_B" | TYPE == "DLBCL" | TYPE == "DLBCL_EGA" | TYPE == "DLBCL_EN" | TYPE == "Benign_LN" | TYPE == "LN_NN_A" | TYPE == "LN_NN_B" | TYPE == "RLN") %>% data.frame()

## Rename
#Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'FL_A'] <- 'FL'
#Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'FL_B'] <- 'FL'
#Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'DLBCL_EGA'] <- 'DLBCL'
#Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'DLBCL_EN'] <- 'DLBCL'
#Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'Benign_LN'] <- 'Normal_LN'
#Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'LN_NN_A'] <- 'Normal_LN'
#Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'LN_NN_B'] <- 'Normal_LN'
#Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'RLN'] <- 'Normal_LN'

# Subset of merged methylation data
bval <- bval[, Sample_Annotations$SAMPLE_ID]
dim(bval)

# Annotation - SEX
# https://www.geeksforgeeks.org/grouped-stacked-and-percent-stacked-barplot-in-ggplot2/
pdf("SEX_TEST.pdf")
d <- Sample_Annotations %>% dplyr::count(epitypes, SEX) %>% mutate(pct=n/sum(n))
ggplot(d, aes(fill = SEX, 
              y = n,
              x = epitypes))+ 
  geom_bar(position = position_fill(), stat="identity") +
  geom_text(aes(label=paste0(sprintf("%1.1f", pct*100),"%, ", n)),
            position=position_fill(vjust=0.5), size=5) +
  ggtitle("SEX")+ 
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Annotation - SEX
pdf("SEX.pdf")
ggplot(Sample_Annotations %>% dplyr::count(epitypes, SEX) %>%
         mutate(pct=n/sum(n)),
       aes(x = as.character(epitypes), n, fill=SEX)) +
  geom_bar(position = position_fill(), stat="identity") +
  #geom_text(aes(label=paste0(sprintf("%1.1f", pct*100),"%, ", n)),
  geom_text(aes(label=paste0(n)),
            position=position_fill(vjust=0.5), size=5) +
  scale_fill_brewer(palette="Dark2", na.value="grey") +
  ggtitle("Sex") +
  theme(title = element_text(size=16, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.title = element_text(size=16)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("M"="#8199f0", "F"="#ffa8b9"), na.value="#c9c9c9") +
  scale_fill_manual(values = c("M"="#8199f0", "F"="#ffa8b9"), na.value="#c9c9c9") +
  coord_cartesian(clip = 'off') +
  xlab("FL epitypes") +
  ylab("FL Patient Counts") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  guides(fill=guide_legend(title="Sex"))
dev.off()

# Annotation - BATCH
pdf("BATCH.pdf")
ggplot(Sample_Annotations %>% dplyr::count(epitypes, Batch) %>%
         mutate(pct=n/sum(n)),
       aes(x = as.character(epitypes), n, fill=Batch)) +
  geom_bar(position = position_fill(), stat="identity") +
  #geom_text(aes(label=paste0(sprintf("%1.1f", pct*100),"%, ", n)),
  geom_text(aes(label=paste0(n)),
            position=position_fill(vjust=0.5), size=5) +
  scale_fill_brewer(palette="Dark2", na.value="grey") +
  ggtitle("Batch") +
  theme(title = element_text(size=16, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.title = element_text(size=16)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("Benign_LN"="#7f8992", "Blueprint"="#4E4EC2", "EGAD00010001974"="#02bdaa", "Genome-Quebec"="#6dcc04", "GSE237299"="#FF7F00", "PM-Genomics-Centre"="#7818e1", "PM-OICR-TGL"="#b5b80b")) +
  scale_fill_manual(values = c("Benign_LN"="#7f8992", "Blueprint"="#4E4EC2", "EGAD00010001974"="#02bdaa", "Genome-Quebec"="#6dcc04", "GSE237299"="#FF7F00", "PM-Genomics-Centre"="#7818e1", "PM-OICR-TGL"="#b5b80b")) +
  coord_cartesian(clip = 'off') +
  xlab("FL epitypes") +
  ylab("FL Patient Counts") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  guides(fill=guide_legend(title="Batch", nrow = 4))
dev.off()


# Annotation - TYPE
pdf("TYPE.pdf")
ggplot(Sample_Annotations %>% dplyr::count(epitypes, TYPE) %>%
         mutate(pct=n/sum(n)),
       aes(x = as.character(epitypes), n, fill=TYPE)) +
  #geom_bar(stat="identity") +
  #geom_bar(position = position_fill(), stat="identity") +
  #geom_text(aes(label=paste0(n)),
            #geom_text(aes(label=paste0(sprintf("%1.1f", pct*100),"%, ", n)),
            #position=position_stack(vjust=0.5), size=5) +
  geom_bar(position = position_fill(), stat="identity") +
  #geom_text(aes(label=paste0(sprintf("%1.1f", pct*100),"%, ", n)),
  geom_text(aes(label=paste0(n)),
            position=position_fill(vjust=0.5), size=3.5) +
  scale_fill_brewer(palette="Dark2", na.value="grey") +
  ggtitle("Total Samples") +
  theme(title = element_text(size=16, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text=element_text(size=14),
        legend.position="bottom",
        legend.title = element_text(size=16)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c(Benign_LN="#5086F3", bm_PC="#7f8992", DLBCL="#A01641", DLBCL_EGA="#fc1417", DLBCL_EN="#f26ded", FL="#6dcc04", FL_A="#02bdaa", FL_B="#ebc602", gcBC="#4E4EC2", LN_NN_A="#FF7F00", LN_NN_B="#B15928", memBC="#b5b80b", naiBC="#2480b9", RLN="#7618dc", t_naiBC="#6e988c", t_PC="#383838", na.value="#c9c9c9")) +
  scale_fill_manual(values = c(Benign_LN="#5086F3", bm_PC="#7f8992", DLBCL="#A01641", DLBCL_EGA="#fc1417", DLBCL_EN="#f26ded", FL="#6dcc04", FL_A="#02bdaa", FL_B="#ebc602", gcBC="#4E4EC2", LN_NN_A="#FF7F00", LN_NN_B="#B15928", memBC="#b5b80b", naiBC="#2480b9", RLN="#7618dc", t_naiBC="#6e988c", t_PC="#383838", na.value="#c9c9c9")) +
  #scale_color_manual(values = c(Benign_LN="#5086F3", bm_PC="#7f8992", DLBCL="#fc1417", FL="#6dcc04", gcBC="#0303a8", LN_NN_A="#FF7F00", LN_NN_B="#B15928", memBC="#b5b80b", naiBC="#2480b9", RLN="#7618dc", t_naiBC="#6e988c", t_PC="#383838", na.value="#c9c9c9")) +
  #scale_fill_manual(values = c(Benign_LN="#5086F3", bm_PC="#7f8992", DLBCL="#fc1417", FL="#6dcc04", gcBC="#0303a8", LN_NN_A="#FF7F00", LN_NN_B="#B15928", memBC="#b5b80b", naiBC="#2480b9", RLN="#7618dc", t_naiBC="#6e988c", t_PC="#383838", na.value="#c9c9c9")) +
  coord_cartesian(clip = 'off') +
  xlab("epitypes") +
  ylab("Sample Counts") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  guides(fill=guide_legend(title="Type", nrow = 3))
dev.off()


# Annotation - ANN_ARBOR_STAGE
pdf("ANN_ARBOR_STAGE.pdf")
ggplot(Sample_Annotations %>% dplyr::count(epitypes, ANN_ARBOR_STAGE) %>%
         mutate(pct=n/sum(n)),
       aes(x = as.character(epitypes), n, fill = as.character(ANN_ARBOR_STAGE))) +
  geom_bar(position = position_fill(), stat="identity") +
  #geom_text(aes(label=paste0(sprintf("%1.1f", pct*100),"%, ", n)),
  geom_text(aes(label=paste0(n)),
            position=position_fill(vjust=0.5), size=4) +
  scale_fill_brewer(palette="Dark2", na.value="grey") +
  ggtitle("Ann Arbor Stage") +
  theme(title = element_text(size=16, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.title = element_text(size=16)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("1"="#8199f0", "2"="#5370d4", "3"="#ffa8b9", "4"="#db7488"), na.value="#c9c9c9") +
  scale_fill_manual(values = c("1"="#8199f0", "2"="#5370d4", "3"="#ffa8b9", "4"="#db7488"), na.value="#c9c9c9") +
  coord_cartesian(clip = 'off') +
  xlab("FL epitypes") +
  ylab("Relative proportion") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  guides(fill=guide_legend(title="Ann Arbor Stage"))
dev.off()

# Annotation - ClusterAIC
pdf("ClusterAIC.pdf")
ggplot(Sample_Annotations %>% dplyr::count(epitypes, ClusterAIC) %>%
         mutate(pct=n/sum(n)),
       aes(x = as.character(epitypes), n, fill=ClusterAIC)) +
  geom_bar(position = position_fill(), stat="identity") +
  #geom_text(aes(label=paste0(sprintf("%1.1f", pct*100),"%, ", n)),
  geom_text(aes(label=paste0(n)),
            position=position_fill(vjust=0.5), size=4) +
  scale_fill_brewer(palette="Dark2", na.value="grey") +
  ggtitle("ClusterAIC") +
  theme(title = element_text(size=16, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.title = element_text(size=16)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c(CS="#ffa9a3", TT="#e7e8b0", GM="#addece", Q="#a2cfe0", AR="#e7c0eb"), na.value="#c9c9c9") +
  scale_fill_manual(values = c(CS="#ffa9a3", TT="#e7e8b0", GM="#addece", Q="#a2cfe0", AR="#e7c0eb"), na.value="#c9c9c9") +
  coord_cartesian(clip = 'off') +
  xlab("FL epitypes") +
  ylab("FL Patient Counts") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  guides(fill=guide_legend(title="ClusterAIC"))
dev.off()


# Annotation - EZH2
pdf("EZH2.pdf")
ggplot(Sample_Annotations %>% dplyr::count(epitypes, EZH2_STATUS) %>%
         mutate(pct=n/sum(n)),
       aes(x = as.character(epitypes), n, fill=EZH2_STATUS)) +
  geom_bar(position = position_fill(), stat="identity") +
  #geom_text(aes(label=paste0(sprintf("%1.1f", pct*100),"%, ", n)),
  geom_text(aes(label=paste0(n)),
            position=position_fill(vjust=0.5), size=4) +
  scale_fill_brewer(palette="Dark2", na.value="grey") +
  ggtitle("EZH2 Mutation") +
  theme(title = element_text(size=16, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.title = element_text(size=16)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("EZH2_Mut_Pos"="#ffa8b9", "EZH2_Mut_Neg"="#8199f0"), na.value="#c9c9c9") +
  scale_fill_manual(values = c("EZH2_Mut_Pos"="#ffa8b9", "EZH2_Mut_Neg"="#8199f0"), na.value="#c9c9c9") +
  coord_cartesian(clip = 'off') +
  xlab("FL epitypes") +
  ylab("FL Patient Counts") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  guides(fill=guide_legend(title="EZH2 Status"))
dev.off()


# Annotation - IGH
pdf("IGH.pdf")
ggplot(Sample_Annotations %>% dplyr::count(epitypes, IGH_STATUS) %>%
         mutate(pct=n/sum(n)),
       aes(x = as.character(epitypes), n, fill=IGH_STATUS)) +
  geom_bar(position = position_fill(), stat="identity") +
  #geom_text(aes(label=paste0(sprintf("%1.1f", pct*100),"%, ", n)),
  geom_text(aes(label=paste0(n)),
            position=position_fill(vjust=0.5), size=4) +
  scale_fill_brewer(palette="Dark2", na.value="grey") +
  ggtitle("IGH") +
  theme(title = element_text(size=16, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.title = element_text(size=16)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c(IGHA="#e7c0eb", IGHG="#addece", IGHM="#ffa9a3"), na.value="#c9c9c9") +
  scale_fill_manual(values = c(IGHA="#e7c0eb", IGHG="#addece", IGHM="#ffa9a3"), na.value="#c9c9c9") +
  coord_cartesian(clip = 'off') +
  xlab("FL epitypes") +
  ylab("FL Patient Counts") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  guides(fill=guide_legend(title="IGH Status"))
dev.off()

# Annotation - T_14_18
pdf("T_14_18.pdf")
ggplot(Sample_Annotations %>% dplyr::count(epitypes, T_14_18) %>%
         mutate(pct=n/sum(n)),
       aes(x = as.character(epitypes), n, fill= as.character(T_14_18))) +
  geom_bar(position = position_fill(), stat="identity") +
  #geom_text(aes(label=paste0(sprintf("%1.1f", pct*100),"%, ", n)),
  geom_text(aes(label=paste0(n)),
            position=position_fill(vjust=0.5), size=4) +
  scale_fill_brewer(palette="Dark2", na.value="grey") +
  ggtitle("T_14_18") +
  theme(title = element_text(size=16, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.title = element_text(size=16)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("1"="#ffa8b9", "0"="#8199f0"), na.value="#c9c9c9") +
  scale_fill_manual(values = c("1"="#ffa8b9", "0"="#8199f0"), na.value="#c9c9c9") +
  coord_cartesian(clip = 'off') +
  xlab("FL epitypes") +
  ylab("FL Patient Counts") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  guides(fill=guide_legend(title="T_14_18 Status"))
dev.off()

# Annotation - INSTITUTION
pdf("INSTITUTION.pdf")
ggplot(Sample_Annotations %>% dplyr::count(epitypes, INSTITUTION) %>%
         mutate(pct=n/sum(n)),
       aes(x = as.character(epitypes), n, fill=INSTITUTION)) +
  geom_bar(position = position_fill(), stat="identity") +
  #geom_text(aes(label=paste0(sprintf("%1.1f", pct*100),"%, ", n)),
  geom_text(aes(label=paste0(n)),
            position=position_fill(vjust=0.5), size=4) +
  scale_fill_brewer(palette="Dark2", na.value="grey") +
  ggtitle("INSTITUTION") +
  theme(title = element_text(size=16, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.title = element_text(size=16)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c(AARHUS="#ffa9a3", BCCA="#e7e8b0", E2408="#addece", JGH="#a2cfe0", KINGSTON="#e7c0eb", OSLO="#8199f0", UHN="#cdb6f2", na.value="#c9c9c9")) +
  scale_fill_manual(values = c(AARHUS="#ffa9a3", BCCA="#e7e8b0", E2408="#addece", JGH="#a2cfe0", KINGSTON="#e7c0eb", OSLO="#8199f0", UHN="#cdb6f2", na.value="#c9c9c9")) +
  coord_cartesian(clip = 'off') +
  xlab("FL epitypes") +
  ylab("FL Patient Counts") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  guides(fill = guide_legend(title="INSTITUTION", nrow = 3, byrow = TRUE))
dev.off()


# Annotation - TIME_POINT
pdf("TIME_POINT.pdf")
ggplot(Sample_Annotations %>% dplyr::count(epitypes, TIME_POINT) %>%
         mutate(pct=n/sum(n)),
       aes(x = as.character(epitypes), n, fill= as.character(TIME_POINT))) +
  geom_bar(position = position_fill(), stat="identity") +
  #geom_text(aes(label=paste0(sprintf("%1.1f", pct*100),"%, ", n)),
  geom_text(aes(label=paste0(n)),
            position=position_fill(vjust=0.5), size=4) +
  scale_fill_brewer(palette="Dark2", na.value="grey") +
  ggtitle("TIME_POINT") +
  theme(title = element_text(size=16, face='bold'),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.text=element_text(size=16),
        legend.position="bottom",
        legend.title = element_text(size=16)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("T2"="#ffa8b9", "T1"="#8199f0"), na.value="#c9c9c9") +
  scale_fill_manual(values = c("T2"="#ffa8b9", "T1"="#8199f0"), na.value="#c9c9c9") +
  coord_cartesian(clip = 'off') +
  xlab("FL epitypes") +
  ylab("FL Patient Counts") +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  guides(fill=guide_legend(title="TIME_POINT"))
dev.off()
