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
library(ggplot2)
library(tidyr)
library(ggpubr)  

# Change directory
setwd("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/09_CNV/450k/CNV_Mean_Plot/")

# C1
# ==
# CNV.summaryplot
setMethod("CNV.summaryplot", signature(object = "CNV.analysis"), function(object,
                                                                          set_par = TRUE, main = NULL, output = "local", directory = getwd(), width = 12, height = 6, res = 720, threshold = 0.1,...) {
  
  if (set_par) {
    mfrow_original <- par()$mfrow
    mar_original <- par()$mar
    oma_original <- par()$oma
  }
  
  if(ncol(object@fit$ratio) <= 1) {
    stop("Please use multiple query samples to create a summaryplot")
  }
  
  # Output setup (PDF or PNG) -- removed if not needed for final output
  if(output == "pdf"){
    p_names <- paste(directory,"/","genome_summaryplot",".pdf",sep="")
    pdf(p_names, width = width, height = height)
    par(mfrow = c(1, 1), mar = c(4, 4, 4, 4), oma = c(0, 0, 0, 0))
  }
  
  if(output == "png"){
    p_names <- paste(directory,"/", "genome_summaryplot",".png",sep="")
    png(p_names, units = "in", width = width, height = height, res = res)
    par(mfrow = c(1, 1), mar = c(4, 4, 4, 4), oma = c(0, 0, 0, 0))
  }
  
  # Generate thresholded data for CNV analysis
  y <- CNV.write(object, what = "threshold", threshold = threshold)
  
  message("creating summaryplot")
  
  segments.i <- GRanges(seqnames=y$Chromosome,ranges=IRanges(y$Start_Position, y$End_Position))
  segments_chromosomes <- GRanges(seqnames = object@anno@genome$chr, ranges = IRanges(start = 1, end = object@anno@genome$size))
  segments <- c(segments.i, segments_chromosomes)
  d_segments <- as.data.frame(GenomicRanges::disjoin(segments))
  
  overview <- as.data.frame(matrix(nrow = 0, ncol = 4))
  for (i in 1:nrow(d_segments)) {
    x <- d_segments[i,]
    involved_segements <- y[y$Chromosome == x$seqnames & y$Start_Position <= x$start & y$End_Position >= x$end,]
    balanced <- sum(involved_segements$Alteration == "balanced")
    gain <- sum(involved_segements$Alteration == "gain")
    loss <- sum(involved_segements$Alteration == "loss")
    overview <- rbind(overview, c(as.character(x$seqnames), balanced, gain, loss))
  }
  
  colnames(overview) <- c("disjoined_segment", "count_balanced", "count_gains", "count_losses")
  overview$count_balanced <- as.numeric(overview$count_balanced)
  overview$count_gains <- as.numeric(overview$count_gains)
  overview$count_losses <- as.numeric(overview$count_losses)
  
  d_segments$gains <- overview$count_gains / length(unique(y$Sample)) * 100
  d_segments$losses <- overview$count_losses / length(unique(y$Sample)) * 100
  d_segments$balanced <- overview$count_balanced / length(unique(y$Sample)) * 100
  
  segments_pl <- d_segments[rep(1:nrow(d_segments), each = 2), ]
  odd_indexes <- seq(1, nrow(segments_pl), 2)
  even_indexes <- seq(2, nrow(segments_pl), 2)
  segments_pl$xpos <- NA
  segments_pl$xpos[odd_indexes] <- segments_pl$start[odd_indexes]
  segments_pl$xpos[even_indexes] <- segments_pl$end[even_indexes]
  
  segments_pl$seqnames <- factor(segments_pl$seqnames, levels = object@anno@genome$chr)
  segments_pl <- segments_pl[order(segments_pl$seqnames, segments_pl$start),]
  
  # Write the final segments_pl to CSV (Retained as required)
  write.csv(segments_pl, file = paste(directory, "/C1_segments_pl.csv", sep = ""), row.names = FALSE)
  message("Final segments_pl written to: ", paste(directory, "/C1_segments_pl.csv", sep = ""))
  
  # Remove the plotting and final processing code after writing segments_pl to CSV
  if (set_par)
    par(mfrow = mfrow_original, mar = mar_original, oma = oma_original)
})

# x <- readRDS("../01_C1_FL/x.rds")
xf <- readRDS("../01_C1_FL/xf.rds")
CNV.summaryplot(xf)


# C2
# ==
# CNV.summaryplot
setMethod("CNV.summaryplot", signature(object = "CNV.analysis"), function(object,
                                                                          set_par = TRUE, main = NULL, output = "local", directory = getwd(), width = 12, height = 6, res = 720, threshold = 0.1,...) {
  
  if (set_par) {
    mfrow_original <- par()$mfrow
    mar_original <- par()$mar
    oma_original <- par()$oma
  }
  
  if(ncol(object@fit$ratio) <= 1) {
    stop("Please use multiple query samples to create a summaryplot")
  }
  
  # Output setup (PDF or PNG) -- removed if not needed for final output
  if(output == "pdf"){
    p_names <- paste(directory,"/","genome_summaryplot",".pdf",sep="")
    pdf(p_names, width = width, height = height)
    par(mfrow = c(1, 1), mar = c(4, 4, 4, 4), oma = c(0, 0, 0, 0))
  }
  
  if(output == "png"){
    p_names <- paste(directory,"/", "genome_summaryplot",".png",sep="")
    png(p_names, units = "in", width = width, height = height, res = res)
    par(mfrow = c(1, 1), mar = c(4, 4, 4, 4), oma = c(0, 0, 0, 0))
  }
  
  # Generate thresholded data for CNV analysis
  y <- CNV.write(object, what = "threshold", threshold = threshold)
  
  message("creating summaryplot")
  
  segments.i <- GRanges(seqnames=y$Chromosome,ranges=IRanges(y$Start_Position, y$End_Position))
  segments_chromosomes <- GRanges(seqnames = object@anno@genome$chr, ranges = IRanges(start = 1, end = object@anno@genome$size))
  segments <- c(segments.i, segments_chromosomes)
  d_segments <- as.data.frame(GenomicRanges::disjoin(segments))
  
  overview <- as.data.frame(matrix(nrow = 0, ncol = 4))
  for (i in 1:nrow(d_segments)) {
    x <- d_segments[i,]
    involved_segements <- y[y$Chromosome == x$seqnames & y$Start_Position <= x$start & y$End_Position >= x$end,]
    balanced <- sum(involved_segements$Alteration == "balanced")
    gain <- sum(involved_segements$Alteration == "gain")
    loss <- sum(involved_segements$Alteration == "loss")
    overview <- rbind(overview, c(as.character(x$seqnames), balanced, gain, loss))
  }
  
  colnames(overview) <- c("disjoined_segment", "count_balanced", "count_gains", "count_losses")
  overview$count_balanced <- as.numeric(overview$count_balanced)
  overview$count_gains <- as.numeric(overview$count_gains)
  overview$count_losses <- as.numeric(overview$count_losses)
  
  d_segments$gains <- overview$count_gains / length(unique(y$Sample)) * 100
  d_segments$losses <- overview$count_losses / length(unique(y$Sample)) * 100
  d_segments$balanced <- overview$count_balanced / length(unique(y$Sample)) * 100
  
  segments_pl <- d_segments[rep(1:nrow(d_segments), each = 2), ]
  odd_indexes <- seq(1, nrow(segments_pl), 2)
  even_indexes <- seq(2, nrow(segments_pl), 2)
  segments_pl$xpos <- NA
  segments_pl$xpos[odd_indexes] <- segments_pl$start[odd_indexes]
  segments_pl$xpos[even_indexes] <- segments_pl$end[even_indexes]
  
  segments_pl$seqnames <- factor(segments_pl$seqnames, levels = object@anno@genome$chr)
  segments_pl <- segments_pl[order(segments_pl$seqnames, segments_pl$start),]
  
  # Write the final segments_pl to CSV (Retained as required)
  write.csv(segments_pl, file = paste(directory, "/C2_segments_pl.csv", sep = ""), row.names = FALSE)
  message("Final segments_pl written to: ", paste(directory, "/C2_segments_pl.csv", sep = ""))
  
  # Remove the plotting and final processing code after writing segments_pl to CSV
  if (set_par)
    par(mfrow = mfrow_original, mar = mar_original, oma = oma_original)
})

# x <- readRDS("../02_C2_FL/x.rds")
xf <- readRDS("../02_C2_FL/xf.rds")
CNV.summaryplot(xf)

# For statistical comparisons
# Step 1: Read the data
C1_data <- read.csv("C1_segments_pl.csv")  # Replace with the actual file path
C2_data <- read.csv("C2_segments_pl.csv")  # Replace with the actual file path

# Step 2: Add a new column indicating the dataset (C1 or C2)
C1_data$Dataset <- "C1_FL"
C2_data$Dataset <- "C2_FL"

# Step 3: Compute total (sum of gains and losses) and append it to the dataset
C1_data$Total <- C1_data$gains + C1_data$losses
C2_data$Total <- C2_data$gains + C2_data$losses

# Step 4: Combine both datasets
combined_data <- bind_rows(C1_data, C2_data)

# Step 5: Reshape the data to a long format for easy plotting
combined_data_long <- combined_data %>%
  gather(key = "Type", value = "Value", Total, gains, losses)  # Reordered: Total → Gains → Losses

# Step 6: Reorder the "Type" column to ensure "Total" comes first, followed by "Gains" and "Losses"
combined_data_long$Type <- factor(combined_data_long$Type, levels = c("Total", "gains", "losses"))

# Step 7: Compute the maximum mean value of "Total" for y-axis limit
max_total_mean <- combined_data_long %>%
  filter(Type == "Total") %>%
  group_by(Dataset) %>%
  summarise(mean_value = mean(Value, na.rm = TRUE)) %>%
  summarise(max_mean = max(mean_value)) %>%
  pull(max_mean)

y_limit <- max_total_mean * 1.5  # Adding 20% padding

# Rename epitypes
table(combined_data_long$Dataset)

combined_data_long$Dataset[combined_data_long$Dataset == 'C1_FL'] <- 'iFL'
combined_data_long$Dataset[combined_data_long$Dataset == 'C2_FL'] <- 'aFL'
combined_data_long$Dataset <- factor(combined_data_long$Dataset, levels = c("iFL", "aFL"))

table(combined_data_long$Dataset)

# Step 8: Create the mean bar plot with p-values
pdf("mean_plot_with_total_reordered.pdf", width = 4, height = 5.2)

ggplot(combined_data_long, aes(x = Type, y = Value, fill = Dataset)) +
  stat_summary(geom = "bar", fun = "mean", position = position_dodge(width = 0.7), width = 0.6) +  # Bar plot for mean
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(width = 0.7), width = 0.3) +  # Error bars
  labs(x = "CNV Type", 
       y = "Average percentage of samples \n exhibiting genome-wide CNVs",
       fill = "Epitypes") +  # Legend title
  scale_fill_manual(values = c("iFL" = "#59A14F", 
                               "aFL" = "#E15759")) + # Modern, more vibrant colors ("#81C78E", "#F17B77")
  scale_x_discrete(labels = c("Total" = "Total", "gains" = "Gains", "losses" = "Losses")) +  # Capitalize labels
  theme_light() +  # Light and modern theme
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_text(size = 12, face = "bold", color = "#333333"),  # Bold legend title
    legend.text = element_text(size = 12, color = "#333333"),  # Synchronize legend text size
    panel.grid.major = element_line(color = "#E0E0E0"),  # Soft grid lines
    panel.grid.minor = element_blank(),  # No minor grid lines
    axis.text = element_text(size = 12, color = "#333333"),  # Dark axis labels
    axis.title = element_text(size = 12, face = "bold", color = "#333333"),  # Bold axis titles
    plot.title = element_text(size = 18, face = "bold", color = "#333333", hjust = 0.5),  # Centered title
    axis.title.x = element_text(margin = margin(t = 10)),  # Margin between x-axis and x-label
    axis.title.y = element_text(margin = margin(r = 10), lineheight = 1.2),  # Margin between y-axis and y-label
    axis.line = element_line(),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    panel.border = element_blank(),  # Remove top and right borders
    strip.text = element_text(size = 16)  # Increase size of strip text
  ) +
  
  # Adjust y-axis limits based on "Total"
  coord_cartesian(ylim = c(0, y_limit)) +
  
  # Add p-values at 90% of max "Total" mean, replacing "=" with ":"
  stat_compare_means(aes(group = Dataset), method = "t.test", label = "p.signif", 
                     label.y = y_limit * 0.9, size = 6, label.x = 1.5)  # Increased p-value size and label adjustments

dev.off()
