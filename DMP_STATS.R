# R Packages
# ==========
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/07_DMP_DMR/450k/hg19/DMP_STATS/")

# Read in 01_C3_vs_C1 (35,828)
temp <- read.delim("../01_C3_vs_C1/DMP_C3_Normal_C1_FL.tsv", header = TRUE)[, 1:12]
C3_vs_C1 <- temp[,-1]
rownames(C3_vs_C1) <- temp[,1]
dim(C3_vs_C1)

# Read in 02_C3_vs_C2 (39,152)
temp <- read.delim("../02_C3_vs_C2/DMP_C3_Normal_C2_FL.tsv", header = TRUE)[, 1:12]
C3_vs_C2 <- temp[,-1]
rownames(C3_vs_C2) <- temp[,1]
dim(C3_vs_C2)

# Read in 03_C1_vs_C2 (29,364)
temp <- read.delim("../03_C1_vs_C2/DMP_C1_FL_C2_FL.tsv", header = TRUE)[, 1:12]
C1_vs_C2 <- temp[,-1]
rownames(C1_vs_C2) <- temp[,1]
dim(C1_vs_C2)

# Read in 04_naiBC_vs_C1 (30,584)
temp <- read.delim("../04_naiBC_vs_C1/DMP_naiBC_C1_FL.tsv", header = TRUE)[, 1:12]
naiBC_vs_C1 <- temp[,-1]
rownames(naiBC_vs_C1) <- temp[,1]
dim(naiBC_vs_C1)

# Read in 05_gcBC_vs_C1 (30,556)
temp <- read.delim("../05_gcBC_vs_C1/DMP_gcBC_C1_FL.tsv", header = TRUE)[, 1:12]
gcBC_vs_C1 <- temp[,-1]
rownames(gcBC_vs_C1) <- temp[,1]
dim(gcBC_vs_C1)

# Read in 06_memBC_vs_C1 (13,809)
temp <- read.delim("../06_memBC_vs_C1/DMP_memBC_C1_FL.tsv", header = TRUE)[, 1:12]
memBC_vs_C1 <- temp[,-1]
rownames(memBC_vs_C1) <- temp[,1]
dim(memBC_vs_C1)

# Read in 07_naiBC_vs_C2 (33,847)
temp <- read.delim("../07_naiBC_vs_C2/DMP_naiBC_C2_FL.tsv", header = TRUE)[, 1:12]
naiBC_vs_C2 <- temp[,-1]
rownames(naiBC_vs_C2) <- temp[,1]
dim(naiBC_vs_C2)

# Read in 08_gcBC_vs_C2 (34,122)
temp <- read.delim("../08_gcBC_vs_C2/DMP_gcBC_C2_FL.tsv", header = TRUE)[, 1:12]
gcBC_vs_C2 <- temp[,-1]
rownames(gcBC_vs_C2) <- temp[,1]
dim(gcBC_vs_C2)

# Read in 09_memBC_vs_C2 (24,901)
temp <- read.delim("../09_memBC_vs_C2/DMP_memBC_C2_FL.tsv", header = TRUE)[, 1:12]
memBC_vs_C2 <- temp[,-1]
rownames(memBC_vs_C2) <- temp[,1]
dim(memBC_vs_C2)

remove(temp)

# List of all data frames
df_list <- list(
  "C1_vs_C2" = C1_vs_C2,
  "C3_vs_C1" = C3_vs_C1,
  "C3_vs_C2" = C3_vs_C2,
  "naiBC_vs_C1" = naiBC_vs_C1,
  "naiBC_vs_C2" = naiBC_vs_C2,
  "gcBC_vs_C1" = gcBC_vs_C1,
  "gcBC_vs_C2" = gcBC_vs_C2,
  "memBC_vs_C1" = memBC_vs_C1,
  "memBC_vs_C2" = memBC_vs_C2
)

# Process each data frame and store results
combined_table <- bind_rows(
  lapply(names(df_list), function(df_name) {
    df_list[[df_name]] %>%
      filter(HYPER_HYPO_001 %in% c("HYPER", "HYPO")) %>%  # Keep only HYPER and HYPO
      mutate(
        # Combine N_Shelf and S_Shelf into Shelf
        Relation_to_Island = recode(Relation_to_Island, 
                                    "N_Shelf" = "Shelf", 
                                    "S_Shelf" = "Shelf", 
                                    "N_Shore" = "Shore", 
                                    "S_Shore" = "Shore"),
        Comparison = factor(df_name, 
                            levels = c("C1_vs_C2", "C3_vs_C1", "C3_vs_C2", 
                                       "naiBC_vs_C1", "naiBC_vs_C2", 
                                       "gcBC_vs_C1", "gcBC_vs_C2", 
                                       "memBC_vs_C1", "memBC_vs_C2"))
      ) %>%
      group_by(Relation_to_Island, HYPER_HYPO_001, Comparison) %>%
      summarise(Count = n(), .groups = "drop") %>%
      pivot_wider(names_from = HYPER_HYPO_001, values_from = Count, values_fill = list(Count = 0))  # Pivot table
  })
)

# Rename Comparison Text
combined_table <- combined_table %>%
  mutate(Comparison = case_when(
    Comparison == "C1_vs_C2" ~ "aFL vs. iFL", #"C2_FL vs. C1_FL", # "C1_FL vs. C2_FL",
    Comparison == "C3_vs_C1" ~ "iFL vs. C3_", # "C1_FL vs. C3_", # "C3_ vs. C1_FL",
    Comparison == "C3_vs_C2" ~ "aFL vs. C3_", # "C2_FL vs. C3_", # "C3_ vs. C2_FL",
    Comparison == "naiBC_vs_C1" ~ "iFL vs. naiBC", # "C1_FL vs. naiBC", # "naiBC vs. C1_FL",
    Comparison == "naiBC_vs_C2" ~ "aFL vs. naiBC", # "C2_FL vs. naiBC", # "naiBC vs. C2_FL",
    Comparison == "gcBC_vs_C1" ~  "iFL vs. gcBC", # "C1_FL vs. gcBC", # "gcBC vs. C1_FL",
    Comparison == "gcBC_vs_C2" ~ "aFL vs. gcBC", # "C2_FL vs. gcBC", # "gcBC vs. C2_FL",
    Comparison == "memBC_vs_C1" ~ "iFL vs. memBC", # "C1_FL vs. memBC", # "memBC vs. C1_FL",
    Comparison == "memBC_vs_C2" ~ "aFL vs. memBC", # "C2_FL vs. memBC", # "memBC vs. C2_FL",
    TRUE ~ as.character(Comparison)
  ))

combined_table <- combined_table %>%
  mutate(Comparison = factor(Comparison, levels = c(
    # "C1_FL vs. C2_FL", "C3_ vs. C1_FL", "C3_ vs. C2_FL", 
    # "naiBC vs. C1_FL", "naiBC vs. C2_FL", 
    # "gcBC vs. C1_FL", "gcBC vs. C2_FL", 
    # "memBC vs. C1_FL", "memBC vs. C2_FL"
    
    # "C2_FL vs. C1_FL", "C1_FL vs. C3_", "C2_FL vs. C3_",
    # "C1_FL vs. naiBC", "C2_FL vs. naiBC",
    # "C1_FL vs. gcBC", "C2_FL vs. gcBC",
    # "C1_FL vs. memBC", "C2_FL vs. memBC"
    
    "aFL vs. iFL", "iFL vs. C3_", "aFL vs. C3_",
    "iFL vs. naiBC", "aFL vs. naiBC",
    "iFL vs. gcBC", "aFL vs. gcBC",
    "iFL vs. memBC", "aFL vs. memBC"
  ))) %>%
  arrange(Comparison)

# Reorder the Relation_to_Island factor levels
combined_table$Relation_to_Island <- factor(combined_table$Relation_to_Island, 
                                            levels = c("Island", "Shore", "Shelf", "OpenSea"))

# View combined table
print(combined_table)

# Summarize totals for each comparison
summary_totals <- combined_table %>%
  group_by(Comparison) %>%
  summarise(
    Total_HYPER = sum(HYPER),
    Total_HYPO = sum(HYPO)
  )

print(summary_totals)

# Save plot to PDF
pdf("Relation_to_Island.pdf", width = 9, height = 6)

# Reshape data to separate HYPER and HYPO into long format
combined_table_long <- combined_table %>%
  pivot_longer(cols = c("HYPER", "HYPO"), names_to = "HYPER_HYPO_001", values_to = "Count")

# Plot with facets for HYPER and HYPO, showing counts
ggplot(combined_table_long, aes(x = Comparison, y = Count, fill = Relation_to_Island)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +  # Thinner black outline
  labs(
    x = "C1: C1_FL, C2: C2_FL",
    y = "Differentially Methylated CpGs",
    fill = "Relation to Island"
  ) +
  theme_minimal() +
  theme(
    axis.line = element_line(),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    # axis.title.x = element_text(margin = margin(t = 20)),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Angle and adjust x-axis text alignment
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    strip.text = element_text(size = 16, face = "bold"),  # Make facet labels (HYPER, HYPO) bold and larger
    legend.position = "bottom",  # Position the legend at the bottom
    legend.box = "horizontal",  # Ensure the legend is horizontal
    legend.box.spacing = unit(0.5, "cm")  # Adjust spacing between legend items
  ) +
  scale_fill_brewer(palette = "Set2") +  # Use the Set2 palette for vibrant colors
  scale_x_discrete(labels = function(x) gsub("C3_", "Normal", x)) +  # Replace "C3_" with "Normal_" on x-axis
  facet_wrap(~ HYPER_HYPO_001, scales = "free_y") +  # Facet by HYPER_HYPO_001 (HYPER vs HYPO)
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))  # Ensure legend is in a single row for the 'fill' aesthetic

dev.off()

# UCSC_RefGene_Group
# ==================
# # List of all data frames
# df_list <- list(
#   "C3_vs_C1" = C3_vs_C1,
#   "C3_vs_C2" = C3_vs_C2,
#   "C1_vs_C2" = C1_vs_C2,
#   "naiBC_vs_C1" = naiBC_vs_C1,
#   "gcBC_vs_C1" = gcBC_vs_C1,
#   "memBC_vs_C1" = memBC_vs_C1,
#   "naiBC_vs_C2" = naiBC_vs_C2,
#   "gcBC_vs_C2" = gcBC_vs_C2,
#   "memBC_vs_C2" = memBC_vs_C2
# )
# 
# # Function to clean UCSC_RefGene_Group column
# tidy_UCSC_RefGene_Group <- function(group) {
#   sapply(strsplit(group, ";"), `[`, 1)  # Take the first annotation
# }
# 
# # Process each data frame and store results
# combined_table <- bind_rows(
#   lapply(names(df_list), function(df_name) {
#     df_list[[df_name]] %>%
#       filter(HYPER_HYPO_001 %in% c("HYPER", "HYPO")) %>%  # Keep only HYPER and HYPO
#       mutate(UCSC_RefGene_Group = tidy_UCSC_RefGene_Group(UCSC_RefGene_Group)) %>%  # Clean column
#       drop_na(UCSC_RefGene_Group) %>%  # Remove NA entries
#       group_by(UCSC_RefGene_Group, HYPER_HYPO_001) %>%
#       summarise(Count = n(), .groups = "drop") %>%
#       pivot_wider(names_from = HYPER_HYPO_001, values_from = Count, values_fill = list(Count = 0)) %>%
#       mutate(Comparison = df_name)  # Add column to identify comparison
#   })
# )
# 
# # View combined table
# print(combined_table)
# 
# pdf("combined_plot.pdf", width = 12, height = 5)
# # Reshape data to separate HYPER and HYPO into long format
# combined_table_long <- combined_table %>%
#   pivot_longer(cols = c("HYPER", "HYPO"), names_to = "HYPER_HYPO_001", values_to = "Count")
# 
# # Generate the plot
# num_categories <- length(unique(combined_table_long$UCSC_RefGene_Group))
# color_palette <- if (num_categories > 12) {
#   scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(num_categories))
# } else {
#   scale_fill_brewer(palette = "Set3")
# }
# 
# ggplot(combined_table_long, aes(x = Comparison, y = Count, fill = UCSC_RefGene_Group)) +
#   geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +  # Thinner black outline
#   labs(
#     x = "Comparison",
#     y = "Differentially Methylated CpGs",
#     fill = "UCSC RefGene Group"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.title = element_text(size = 14, face = "bold"),
#     axis.text = element_text(size = 12),
#     axis.text.x = element_text(angle = 45, hjust = 1),  # Angle and adjust x-axis text alignment
#     legend.title = element_text(size = 14, face = "bold"),
#     legend.text = element_text(size = 12),
#     plot.title = element_text(size = 16, face = "bold"),
#     plot.margin = margin(1, 1, 1, 1, "cm")
#   ) +
#   ggtitle("HYPER and HYPO by UCSC RefGene Group (DMP)") +
#   color_palette +  # Dynamically chosen color palette
#   facet_wrap(~ HYPER_HYPO_001, scales = "free_y")  # Facet by HYPER_HYPO_001 (HYPER vs HYPO)
# 
# dev.off()

# List of all data frames
df_list <- list(
  "C3_vs_C1" = C3_vs_C1,
  "C3_vs_C2" = C3_vs_C2,
  "C1_vs_C2" = C1_vs_C2,
  "naiBC_vs_C1" = naiBC_vs_C1,
  "gcBC_vs_C1" = gcBC_vs_C1,
  "memBC_vs_C1" = memBC_vs_C1,
  "naiBC_vs_C2" = naiBC_vs_C2,
  "gcBC_vs_C2" = gcBC_vs_C2,
  "memBC_vs_C2" = memBC_vs_C2
)

# Function to clean UCSC_RefGene_Group column
tidy_UCSC_RefGene_Group <- function(group) {
  annotations <- strsplit(group, ";")[[1]]
  priority_regions <- c("1stExon", "3'UTR", "5'UTR", "TSS1500", "TSS200")
  prioritized_annotation <- intersect(annotations, priority_regions)
  if (length(prioritized_annotation) > 0) return(prioritized_annotation[1])
  return(annotations[1])
}

# Combine and process data frames
combined_table <- bind_rows(
  lapply(names(df_list), function(df_name) {
    df_list[[df_name]] %>%
      filter(HYPER_HYPO_001 %in% c("HYPER", "HYPO")) %>%
      mutate(UCSC_RefGene_Group = sapply(UCSC_RefGene_Group, tidy_UCSC_RefGene_Group)) %>%
      drop_na(UCSC_RefGene_Group) %>%
      group_by(UCSC_RefGene_Group, HYPER_HYPO_001) %>%
      summarise(Count = n(), .groups = "drop") %>%
      pivot_wider(names_from = HYPER_HYPO_001, values_from = Count, values_fill = list(Count = 0)) %>%
      mutate(Comparison = df_name)
  })
)

# Reorder the UCSC_RefGene_Group factor
combined_table$UCSC_RefGene_Group <- factor(combined_table$UCSC_RefGene_Group,
                                            levels = c("Body", setdiff(unique(combined_table$UCSC_RefGene_Group), "Body")))

# Reorder the Comparison factor
combined_table$Comparison <- factor(combined_table$Comparison, 
                                    levels = c("C1_vs_C2", "C3_vs_C1", "C3_vs_C2", 
                                               "naiBC_vs_C1", "naiBC_vs_C2", 
                                               "gcBC_vs_C1", "gcBC_vs_C2", 
                                               "memBC_vs_C1", "memBC_vs_C2"))
# Rename Comparison Text
combined_table <- combined_table %>%
  mutate(Comparison = case_when(
    Comparison == "C1_vs_C2" ~ "aFL vs. iFL", #"C2_FL vs. C1_FL", # "C1_FL vs. C2_FL",
    Comparison == "C3_vs_C1" ~ "iFL vs. C3_", # "C1_FL vs. C3_", # "C3_ vs. C1_FL",
    Comparison == "C3_vs_C2" ~ "aFL vs. C3_", # "C2_FL vs. C3_", # "C3_ vs. C2_FL",
    Comparison == "naiBC_vs_C1" ~ "iFL vs. naiBC", # "C1_FL vs. naiBC", # "naiBC vs. C1_FL",
    Comparison == "naiBC_vs_C2" ~ "aFL vs. naiBC", # "C2_FL vs. naiBC", # "naiBC vs. C2_FL",
    Comparison == "gcBC_vs_C1" ~  "iFL vs. gcBC", # "C1_FL vs. gcBC", # "gcBC vs. C1_FL",
    Comparison == "gcBC_vs_C2" ~ "aFL vs. gcBC", # "C2_FL vs. gcBC", # "gcBC vs. C2_FL",
    Comparison == "memBC_vs_C1" ~ "iFL vs. memBC", # "C1_FL vs. memBC", # "memBC vs. C1_FL",
    Comparison == "memBC_vs_C2" ~ "aFL vs. memBC", # "C2_FL vs. memBC", # "memBC vs. C2_FL",
    TRUE ~ as.character(Comparison)
  ))

combined_table <- combined_table %>%
  mutate(Comparison = factor(Comparison, levels = c(
    # "C1_FL vs. C2_FL", "C3_ vs. C1_FL", "C3_ vs. C2_FL", 
    # "naiBC vs. C1_FL", "naiBC vs. C2_FL", 
    # "gcBC vs. C1_FL", "gcBC vs. C2_FL", 
    # "memBC vs. C1_FL", "memBC vs. C2_FL"
    
    # "C2_FL vs. C1_FL", "C1_FL vs. C3_", "C2_FL vs. C3_",
    # "C1_FL vs. naiBC", "C2_FL vs. naiBC",
    # "C1_FL vs. gcBC", "C2_FL vs. gcBC",
    # "C1_FL vs. memBC", "C2_FL vs. memBC"
    
    "aFL vs. iFL", "iFL vs. C3_", "aFL vs. C3_",
    "iFL vs. naiBC", "aFL vs. naiBC",
    "iFL vs. gcBC", "aFL vs. gcBC",
    "iFL vs. memBC", "aFL vs. memBC"
  ))) %>%
  arrange(Comparison)

# Set the colors for the plot
color_palette <- scale_fill_manual(
  values = c(
    "Body" = "white", "1stExon" = "#80A8F5", "3'UTR" = "#81C78E", 
    "5'UTR" = "#9B7BB6", "TSS1500" = "#F17B77", "TSS200" = "#F2C85A", 
    "Other" = "#A8D2D0"
  )
)

# Summarize totals for each comparison
summary_totals <- combined_table %>%
  group_by(Comparison) %>%
  summarise(
    Total_HYPER = sum(HYPER),
    Total_HYPO = sum(HYPO)
  )

print(summary_totals)

# Plot and save to a PDF
pdf("UCSC_RefGene_Group.pdf", width = 9, height = 6)

# Reshape data and generate the plot
combined_table_long <- combined_table %>%
  pivot_longer(cols = c("HYPER", "HYPO"), names_to = "HYPER_HYPO_001", values_to = "Count")

ggplot(combined_table_long, aes(x = Comparison, y = Count, fill = UCSC_RefGene_Group)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
  labs(x = "C1: C1_FL, C2: C2_FL", y = "Differentially Methylated CpGs", fill = "UCSC RefGene Group") +
  theme_minimal() +
  theme(
    axis.line = element_line(),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    # axis.title.x = element_text(margin = margin(t = 20)),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.spacing = unit(0.5, "cm"),
    strip.text = element_text(size = 16, face = "bold")
  ) +
  color_palette +
  scale_x_discrete(labels = function(x) gsub("C3_", "Normal", x)) +
  facet_wrap(~ HYPER_HYPO_001, scales = "free_y") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

dev.off()
