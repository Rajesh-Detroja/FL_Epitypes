# R Packages
# ==========
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(stringr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(GenomicRanges)
library(purrr)

# Set working directory
setwd("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/08_GO_Pathway/450k/chromatin_state/")

# Load DMPs
DMPs <- read.delim("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/07_DMP_DMR/450k/hg19/04_naiBC_vs_C1/DMP_naiBC_C1_FL.tsv",
                   header = TRUE)[, 1:12]
DMPs <- data.frame(DMPs[,-1], row.names = DMPs[,1])

# DMPs <- read.delim("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/07_DMP_DMR/450k/hg19/07_naiBC_vs_C2/DMP_naiBC_C2_FL.tsv", 
#                    header = TRUE)[, 1:12]
# DMPs <- data.frame(DMPs[,-1], row.names = DMPs[,1])

# Filter DMPs based on HYPER_HYPO_001 column
DMPs_hyper <- DMPs %>% filter(HYPER_HYPO_001 == "HYPER")
DMPs_hypo <- DMPs %>% filter(HYPER_HYPO_001 == "HYPO")

# Define epitypes (modify as necessary)
epitypes <- "C1"

# Create GRanges for hyper and hypo DMPs
DMPs_hyper_GR <- makeGRangesFromDataFrame(DMPs_hyper %>% select(chr = chr, start = pos, end = pos))
DMPs_hypo_GR <- makeGRangesFromDataFrame(DMPs_hypo %>% select(chr = chr, start = pos, end = pos))

# Load chromatin states data
chrom_states <- read.table("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/chromatin_state/wgEncodeBroadHmmGm12878HMM_hg19.bed",
                           sep = "\t", na.strings = c("", "NA"))

# Extract unique chromatin state names
chrom_states_names <- unique(chrom_states$V4)

# Function to process chromatin state overlaps
count_overlap <- function(DMPs_GR, chrom_states_df) {
  count_overlaps <- map(chrom_states_names, function(state) {
    chrom_state <- chrom_states_df %>% filter(V4 == state) %>% select(chr = V1, start = V2, end = V3)
    chrom_state_GR <- makeGRangesFromDataFrame(chrom_state)
    overlap_count <- countOverlaps(DMPs_GR, chrom_state_GR)
    overlap_table <- table(overlap_count)
    if (length(overlap_table) < 2) overlap_table[2] <- 0
    return(c(state = state, overlap_table))
  })
  
  # Combine results into a dataframe
  overlap_data <- do.call(rbind, count_overlaps)
  overlap_df <- as.data.frame(overlap_data)
  return(overlap_df)
}

# Process hypo and hyper chromatin states
chrom_states_list_hypo <- count_overlap(DMPs_hypo_GR, chrom_states)
chrom_states_list_hyper <- count_overlap(DMPs_hyper_GR, chrom_states)

# Combine results and calculate proportions
chrom_states_res <- bind_rows(
  mutate(chrom_states_list_hypo, direction = "hypo"),
  mutate(chrom_states_list_hyper, direction = "hyper")
)

# Calculate percentage of overlap
chrom_states_res$overlap_prop <- (as.numeric(chrom_states_res$`1`) / 
                                    (as.numeric(chrom_states_res$`0`) + as.numeric(chrom_states_res$`1`))) * 100

# Add labels for not_overlap and overlap to the result for writing to a file
chrom_states_res <- chrom_states_res %>%
  mutate(not_overlap = as.numeric(chrom_states_res$`0`),
         overlap = as.numeric(chrom_states_res$`1`))

# Write results to file with additional labels
chrom_states_tsv <- paste0("chrom_states_res_", epitypes, ".tsv")
write.table(chrom_states_res, file = chrom_states_tsv, quote = FALSE, sep = '\t', row.names = FALSE)

# Plot chromatin states
chrom_states_pdf <- paste0("RD_chromatin_states_", epitypes, ".pdf")
pdf(chrom_states_pdf)
chrom_states_res %>%
  mutate(direction = factor(direction, levels = c("hypo", "hyper"))) %>%
  mutate(chrom_state_class = factor(state, levels = c("15_Repetitive/CNV",
                                                      "14_Repetitive/CNV",
                                                      "13_Heterochrom/lo",
                                                      "12_Repressed",
                                                      "11_Weak_Txn",
                                                      "10_Txn_Elongation",
                                                      "9_Txn_Transition",
                                                      "8_Insulator",
                                                      "7_Weak_Enhancer", 
                                                      "6_Weak_Enhancer",
                                                      "5_Strong_Enhancer",
                                                      "4_Strong_Enhancer",
                                                      "3_Poised_Promoter",
                                                      "2_Weak_Promoter",
                                                      "1_Active_Promoter"))) %>%
  ggplot(aes(direction, chrom_state_class, fill = overlap_prop)) +
  ggtitle(paste("Methylation Alterations - Chromatin States - ", epitypes)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_raster() +
  geom_text(aes(label = round(overlap_prop, 1)), color = "white") +
  scale_fill_gradient(low = "#4393c3", high = "#b2182b") +
  theme_bw() +
  theme(title = element_text(size = 12, face = 'bold'),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))
dev.off()
