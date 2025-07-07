# load R packages:
# ===============
library(pacman)
pacman::p_load(tidyr,
               tidyverse,
               data.table,
               ComplexHeatmap, readxl, patchwork, gridExtra, cowplot, survival, survminer, ggpubr, ggplotify)

# set directory
# =============
setwd("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/12_SNP/Scripts/01_Davidson/")

# set the colors for the variant types in the oncoprint
col = c('Missense'='#1ca900', 
        'Frameshift'='#203bff',
        'Nonsense'='#ff9700',
        'INDEL'='#c758ff', 
        'Splice'='#ff0204')

# format the bars for each variant type to allow for overlapping variant types
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  Missense = function(x, y, w, h) grid.rect(x, y, w*0.9, h*1, 
                                            gp = gpar(fill=col['Missense'])),
  Frameshift = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.75,
                                              gp = gpar(fill=col['Frameshift'])),
  Nonsense = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.6,
                                            gp = gpar(fill=col['Nonsense'])),
  INDEL = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4,
                                         gp = gpar(fill=col['INDEL'])),
  Splice = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.2,
                                          gp = gpar(fill=col['Splice'])))

# Load epitypes
# =============
epitypes_df <- read.delim('/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/02_Clustering/450k/modified_InfiniumPurify/Enmix_quantile_with_QCinfo_Harman_All/FL_epitypes.tsv')

# preview epitypes_df
head(epitypes_df)

# verify the number of samples in each epitype
table(epitypes_df$epitypes)

# Load FLOMICS mutation data
# ==========================
# get list of samples in FLOMICS paper. Not all samples had detected mutations
# FLOMICS_samples <- read_excel("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/12_SNP/Scripts/01_Davidson/data/41408_2024_1111_MOESM2_ESM-2.xlsx", sheet = 7) %>% 
#   pull(SAMPLE_ID)

FLOMICS_samples <- read.delim("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/ClusterAIC/2023-05-06_14_GMM_Cluster_Labels_flexmix_clusters.txt") %>% 
  pull(SAMPLE_ID)

# subset epitypes_df for samples with FLOMICS mutation data, remove UC samples
epitypes_FL_FLOMICS_df <- epitypes_df %>% 
  filter(SAMPLE_ID %in% FLOMICS_samples) %>% 
  filter(!epitypes == 'UC')

# get the list of C1 samples with FLOMICS mutation data
FLOMICS_C1_samples <- epitypes_FL_FLOMICS_df %>% 
  filter(epitypes %in% c('C1')) %>% 
  pull(SAMPLE_ID)

# get the list of C2 samples with FLOMICS mutation data
FLOMICS_C2_samples <- epitypes_FL_FLOMICS_df %>% 
  filter(epitypes %in% c('C2')) %>% 
  pull(SAMPLE_ID)

table(epitypes_df$epitypes, dnn = c('Epitypes (all samples)'))

table(epitypes_FL_FLOMICS_df$epitypes, dnn = c('Epitypes for samples with FLOMICS data'))

# Load the mutation data
mutations_long_df <- read_excel("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/12_SNP/Scripts/01_Davidson/data/41408_2024_1111_MOESM2_ESM-2.xlsx", sheet = 5) %>% 
  filter(SAMPLE_ID %in% epitypes_FL_FLOMICS_df$SAMPLE_ID) %>% 
  mutate(Variant_Classification_2 = case_when(Variant_Classification %in% c('Missense_Mutation') ~ 'Missense',
                                              Variant_Classification %in% c('Frame_Shift_Del', 'Frame_Shift_Ins') ~ 'Frameshift', 
                                              Variant_Classification %in% c('In_Frame_Del', 'In_Frame_Ins') ~ 'INDEL', 
                                              Variant_Classification %in% c('Splice_Site') ~ 'Splice', 
                                              Variant_Classification %in% c('Nonsense_Mutation') ~ 'Nonsense'))

# preview mutations_long_df
head(mutations_long_df)

# verify my variant classification (Variant_Classification_2) is accurate
table(mutations_long_df$Variant_Classification_2, mutations_long_df$Variant_Classification)

# check the number of FLOMICS epitype samples with/without any detected mutation
summary(epitypes_FL_FLOMICS_df$SAMPLE_ID %in% mutations_long_df$SAMPLE_ID)

# get the list of 7 samples without any FLOMICS mutation - need to add to final mutations df
samples_no_mut <- epitypes_FL_FLOMICS_df$SAMPLE_ID[!epitypes_FL_FLOMICS_df$SAMPLE_ID %in% mutations_long_df$SAMPLE_ID]
samples_no_mut

# set up the final mutations df - keep all variant types for each gene per sample for the plot
mutations_wide <- mutations_long_df %>% 
  distinct(SAMPLE_ID, Hugo_Symbol, Variant_Classification_2) %>% 
  pivot_wider(names_from = Hugo_Symbol, values_from = Variant_Classification_2, values_fn = ~ paste(.x, collapse = ';'))

mutations_wide[is.na(mutations_wide)] <- ''
mutations_wide <- data.frame(mutations_wide[,-1], row.names=mutations_wide$SAMPLE_ID)

mutations_transposed <- transpose(mutations_wide)
colnames(mutations_transposed) <- rownames(mutations_wide)
rownames(mutations_transposed) <- colnames(mutations_wide)

# add the samples without detected FLOMICS mutation to mutations_transposed
for (i in samples_no_mut) {
  mutations_transposed[[i]] <- ''
}

# ensure a match between the number of samples in mutations df and the expected number
ncol(mutations_transposed) == epitypes_FL_FLOMICS_df$SAMPLE_ID %>% length()

# Set up metadata
# ===============
# get patient metadata
metadata_df <- read.csv("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/CLN/20240116_clinical_data.csv") %>% 
  filter(SAMPLE_ID %in% epitypes_FL_FLOMICS_df$SAMPLE_ID) %>% 
  mutate(epitype = case_when(SAMPLE_ID %in% FLOMICS_C1_samples ~ 'iFL', 
                             SAMPLE_ID %in% FLOMICS_C2_samples ~ 'aFL'))

# Row order of metadata needs to match with the column order in the mutations dataframe for oncoplot
order_samples_mutation_df <- colnames(mutations_transposed)

metadata_df <- metadata_df %>%
  arrange(match(SAMPLE_ID, order_samples_mutation_df))

# check that all epitype samples with FLOMICS mutation data have metadata
summary(epitypes_FL_FLOMICS_df$SAMPLE_ID %in% metadata_df$SAMPLE_ID)

# Reorder epitypes names
table(metadata_df$epitype)
metadata_df$epitype <- factor(metadata_df$epitype, levels = c("iFL", "aFL"))
table(metadata_df$epitype)

# oncoPrint
# =========
# order samples in the oncoplot by epitype and ClusterAIC
FLOMICS_C1_sample_order <- metadata_df %>% 
  filter(epitype == 'iFL') %>% 
  mutate(ClusterAIC = factor(ClusterAIC, levels = c('AR','CS','GM','Q','TT'))) %>% 
  arrange(ClusterAIC) %>% 
  pull(SAMPLE_ID)

FLOMICS_C2_sample_order <- metadata_df %>% 
  filter(epitype == 'aFL') %>% 
  mutate(ClusterAIC = factor(ClusterAIC, levels = c('AR','CS','GM','Q','TT'))) %>% 
  arrange(ClusterAIC) %>% 
  pull(SAMPLE_ID)

order_samples <- c(FLOMICS_C1_sample_order, FLOMICS_C2_sample_order)


ha = HeatmapAnnotation(Epitype = metadata_df[,'epitype'],
                       ClusterAIC = metadata_df[,'ClusterAIC'],
                       col=list(Epitype = c('iFL'='#59A14F','aFL'='#E15759'), 
                                ClusterAIC = c('CS'='#FFA9A3','TT'='#E7E8B0','GM'='#ADDECE','Q'='#A2CFE0','AR'='#E7C0EB')),
                       annotation_height = unit(c(1), "cm"),
                       annotation_legend_param=list(Epitype=list(title='Epitype')),
                       annotation_name_gp = gpar(fontsize=7))

p <- oncoPrint(mutations_transposed,
               show_column_names = TRUE,
               alter_fun = alter_fun,
               col = col,
               column_order = order_samples,
               alter_fun_is_vectorized = FALSE,
               column_title = NULL,
               # heatmap_legend_param = heatmap_legend_param,
               pct_side = 'right',
               row_names_side = 'left',
               right_annotation = rowAnnotation(rbar=anno_oncoprint_barplot(show_fraction=FALSE)),
               bottom_annotation = ha,
               row_title_gp = gpar(fontsize=7),
               column_names_gp = gpar(fontsize=2),
               row_names_gp = gpar(fontsize=12, fontface='italic')) %>% 
  draw() %>% 
  grid.grabExpr() %>% 
  as.ggplot()

pdf("oncoPrint.pdf", width = 15, height = 10)
p
dev.off()

# Frequency of mutated genes across epitypes
# ==========================================
# get list of recurrently mutated genes (n >= 10)
recurrent_mutations <- mutations_long_df %>% 
  distinct(SAMPLE_ID, Hugo_Symbol) %>% 
  group_by(Hugo_Symbol) %>%
  # count() %>%
  summarise(n = n()) %>%
  filter(n >= 10) %>% 
  arrange(desc(n)) %>%
  pull(Hugo_Symbol)

# comparison of mutated gene rate across groups using Fisher's exact test
mutations_comparison_df <- mutations_wide
mutations_comparison_df[mutations_comparison_df == ''] <- 'WT'
# mutations_comparison_df[mutations_comparison_df == ''] <- 'WT'
mutations_comparison_df[mutations_comparison_df != 'WT'] <- 'MUT'

# modify the dataframe for tidyverse format
mutations_comparison_df$SAMPLE_ID <- rownames(mutations_wide)
rownames(mutations_comparison_df) <- NULL

# add on the samples without any mutations
mutations_comparison_df <- mutations_comparison_df %>% 
  bind_rows(tibble(SAMPLE_ID = samples_no_mut))
mutations_comparison_df[is.na(mutations_comparison_df)] <- 'WT'

mutations_comparison_df <- mutations_comparison_df %>% 
  mutate(epitype = case_when(SAMPLE_ID %in% FLOMICS_C1_samples ~ 'C1', 
                             SAMPLE_ID %in% FLOMICS_C2_samples ~ 'C2'))

# preview mutations_comparison_df
head(mutations_comparison_df)

# check the tail of mutations_comparison_df should be all WT (no mut samples)
tail(mutations_comparison_df)

# confirm the number of samples in mutations_comparison_df
mutations_comparison_df %>% nrow() == epitypes_FL_FLOMICS_df$SAMPLE_ID %>% length()

# create vectors to store the gene and p-value for each iteration of Fisher's test
# restrict analysis to recurrent mutations
gene <- NULL
Fisher_test_p_value <- NULL
for (i in recurrent_mutations){
  i.gene <- i
  i.Fisher_test_p_value <- fisher.test(mutations_comparison_df$epitype, mutations_comparison_df[[i]])$p.value
  if(is.null(gene)){gene <- i.gene} else {gene <- c(gene, i.gene)}
  if(is.null(Fisher_test_p_value)){Fisher_test_p_value <- i.Fisher_test_p_value} else {Fisher_test_p_value <- c(Fisher_test_p_value, i.Fisher_test_p_value)}
}

# Create Fisher test dataframe and perform bonferroni correction (multiply by number of comparisons)
comp_mutated_gene_df <- data.frame(gene = gene,
                                   Fisher_test_p_value = Fisher_test_p_value) %>% 
  mutate(FDR_bonferroni = Fisher_test_p_value*length(recurrent_mutations))

# Preview the Fisher test comparison dataframe
head(comp_mutated_gene_df)

# Get the genes with significant difference across epitypes (adj. P-value < 0.05)
significant_gene_comp <- comp_mutated_gene_df %>% 
  filter(FDR_bonferroni < 0.05) %>% 
  pull(gene)

# Get the max value in plot for those significant genes (for offsetting the asterisks)
max_values <- mutations_long_df %>% 
  distinct(SAMPLE_ID, Hugo_Symbol) %>% 
  mutate(epitype = case_when(SAMPLE_ID %in% FLOMICS_C1_samples ~ 'C1', 
                             SAMPLE_ID %in% FLOMICS_C2_samples ~ 'C2')) %>% 
  group_by(epitype, Hugo_Symbol) %>% 
  # count() %>% 
  summarise(n = n()) %>% 
  mutate(prop = case_when(epitype %in% 'C1' ~ n/length(FLOMICS_C1_samples), 
                          epitype %in% 'C2' ~ n/length(FLOMICS_C2_samples))) %>% 
  ungroup() %>% 
  filter(Hugo_Symbol %in% significant_gene_comp) %>% 
  arrange(desc(prop)) %>% 
  distinct(Hugo_Symbol, .keep_all = TRUE) %>% 
  pull(prop)

# Set up asterisks significant dataframe
sig_data <- data.frame(x_pos = match(significant_gene_comp, recurrent_mutations),
                       y_pos = max_values)

p_barplot <- mutations_long_df %>% 
  distinct(SAMPLE_ID, Hugo_Symbol) %>% 
  filter(Hugo_Symbol %in% recurrent_mutations) %>%
  mutate(epitype = case_when(SAMPLE_ID %in% FLOMICS_C1_samples ~ 'iFL', 
                             SAMPLE_ID %in% FLOMICS_C2_samples ~ 'aFL')) %>% 
  group_by(epitype, Hugo_Symbol) %>% 
  # count() %>% 
  summarise(n = n()) %>%
  mutate(prop = case_when(epitype %in% 'iFL' ~ n/length(FLOMICS_C1_samples), 
                          epitype %in% 'aFL' ~ n/length(FLOMICS_C2_samples))) %>%
  mutate(Hugo_Symbol = factor(Hugo_Symbol, levels = recurrent_mutations)) %>% 
  mutate(epitype = factor(epitype, levels = c('iFL','aFL'))) %>% 
  ggplot(aes(x = Hugo_Symbol, y = prop, fill = epitype)) + 
  geom_bar(position = 'dodge', stat = 'identity') + 
  scale_fill_manual(values = c('iFL'='#59A14F','aFL'='#E15759')) + 
  labs(y = 'Proportion of samples with mutated gene',
       x = '', 
       fill = 'Epitypes') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "italic"),
        axis.title.y = element_text(face = "bold"),
        panel.border = element_blank(), 
        axis.line = element_line(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.title = element_text(face = "bold"), # R
        legend.position = "bottom", # R
        legend.spacing.y = unit(0.1, "cm"), # R
        legend.margin = margin(t = -15, r = 0, b = 0, l = 0), # R
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2)) +  # R
  # geom_segment(data = sig_data, aes(x = x_pos - 0.2, xend = x_pos + 0.2, y = y_pos + 0.02, yend = y_pos + 0.02),
  #              size = 0.5, inherit.aes = FALSE) + 
  geom_segment(data = sig_data, aes(x = x_pos - 0.2, xend = x_pos + 0.2, y = y_pos + 0.02, yend = y_pos + 0.02),
               linewidth = 0.5, inherit.aes = FALSE) + 
  geom_text(data = sig_data, aes(x = x_pos, y = y_pos + 0.03, label = "*"), 
            size = 6, inherit.aes = FALSE)

pdf("Barplot.pdf", width = 10, height = 5)
p_barplot
dev.off()

# Compare number of mutated genes across epitype groups
# =====================================================
# create a dataframe with the samples with no mutations, to be rbinded
FLOMICS_no_mutations_df <- data.frame(SAMPLE_ID = samples_no_mut, n = rep(0, length(samples_no_mut)))

# preview FLOMICS_no_mutations_df
head(FLOMICS_no_mutations_df)

p_boxplot_mutated_genes <- mutations_long_df %>% 
  distinct(SAMPLE_ID, Hugo_Symbol) %>% 
  group_by(SAMPLE_ID) %>% 
  # count() %>% 
  summarise(n = n()) %>%
  rbind(FLOMICS_no_mutations_df) %>% 
  mutate(epitype = case_when(SAMPLE_ID %in% FLOMICS_C1_samples ~ 'C1', 
                             SAMPLE_ID %in% FLOMICS_C2_samples ~ 'C2')) %>% 
  mutate(epitype = factor(epitype, levels = c('C1','C2'))) %>% 
  ggplot(aes(x = epitype, y = n, fill = epitype)) + 
  geom_boxplot(outliers = TRUE) + 
  stat_compare_means(comparisons = list(c('C1','C2')), method = 'wilcox', label = 'p.format', label.y = 14.5, size = 3) + 
  scale_fill_manual(values = c('C1'='#59A14F','C2'='#E15759')) + 
  theme_bw() + 
  theme(legend.position = 'none', 
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 8, face = "bold"),
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  labs(y = 'Number of mutated genes per sample') + 
  ylim(c(0,16))

pdf("Mutated_Genes.pdf", width = 3, height = 5)
p_boxplot_mutated_genes
dev.off()


# Compare number of mutations across epitype groups
# =================================================
# create a dataframe with the samples with no mutations, to be rbinded
FLOMICS_no_mutations_df <- data.frame(SAMPLE_ID = samples_no_mut, n = rep(0, length(samples_no_mut)))

# preview FLOMICS_no_mutations_df
head(FLOMICS_no_mutations_df)

p_boxplot_mutations <- mutations_long_df %>% 
  group_by(SAMPLE_ID) %>% 
  # count() %>% 
  summarise(n = n()) %>%
  rbind(FLOMICS_no_mutations_df) %>% 
  mutate(epitype = case_when(SAMPLE_ID %in% FLOMICS_C1_samples ~ 'C1', 
                             SAMPLE_ID %in% FLOMICS_C2_samples ~ 'C2')) %>% 
  mutate(epitype = factor(epitype, levels = c('C1','C2'))) %>% 
  ggplot(aes(x = epitype, y = n, fill = epitype)) + 
  geom_boxplot(outliers = TRUE) + 
  stat_compare_means(comparisons = list(c('C1','C2')), method = 'wilcox', label = 'p.format', label.y = 20.5, size = 3) + 
  scale_fill_manual(values = c('C1'='#59A14F','C2'='#E15759')) + 
  theme_bw() + 
  theme(legend.position = 'none', 
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 8, face = "bold"),
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  labs(y = 'Number of mutations per sample') + 
  ylim(c(0,23))

pdf("Total_Mutations.pdf", width = 3, height = 5)
p_boxplot_mutations
dev.off()

# Compare number of somatic hypermutations across clusters
# ========================================================
pacman::p_load("dplyr", "ggplot2", "survival", "ggpubr", "stringr", "tidyr",
               "gridExtra", "GenomicRanges", "BSgenome.Hsapiens.UCSC.hg19", "magrittr",
               "data.table")

hg19_genome <- BSgenome.Hsapiens.UCSC.hg19
seqnames(hg19_genome) <- gsub("chr", "", seqnames(hg19_genome))

rgyw_motif <- DNAString("RGYW")
wrcy_motif <- DNAString("WRCY")

snvs_rgyw_gr <- mutations_long_df %>%
  dplyr::select(chr = Chromosome, pos = Start_Position,
                ref = Reference_Allele) %>%
  filter(ref == "G") %>%
  distinct(chr, pos, ref) %$%
  GRanges(chr, IRanges(start = pos - 1, end = pos + 2),
          id = str_c(chr, ":", pos))

snvs_wrcy_gr <- mutations_long_df %>%
  dplyr::select(chr = Chromosome, pos = Start_Position,
                ref = Reference_Allele) %>%
  filter(ref == "C") %>%
  distinct(chr, pos, ref) %$%
  GRanges(chr, IRanges(start = pos - 2, end = pos + 1),
          id = str_c(chr, ":", pos))

snvs_rgyw_seq <- getSeq(hg19_genome, snvs_rgyw_gr)
names(snvs_rgyw_seq) <- mcols(snvs_rgyw_gr)[, "id"]

snvs_wrcy_seq <- getSeq(hg19_genome, snvs_wrcy_gr)
names(snvs_wrcy_seq) <- mcols(snvs_wrcy_gr)[, "id"]

snvs_rgyw_seq_df <- data.table(posID = names(snvs_rgyw_seq),
                               RGYW.Motif = as.character(snvs_rgyw_seq))
snvs_wrcy_seq_df <- data.table(posID = names(snvs_wrcy_seq),
                               WRCY.Motif = as.character(snvs_wrcy_seq))

# Check if any of the SNVs overlap with a SHM motifs
snvs_rgyw_vmatch <- vmatchPattern(rgyw_motif, snvs_rgyw_seq, fixed = FALSE)
snvs_wrcy_vmatch <- vmatchPattern(wrcy_motif, snvs_wrcy_seq, fixed = FALSE)

mutations_maf_motif <- mutations_long_df %>%
  mutate(posID = str_c(Chromosome, ":", Start_Position)) %>%
  mutate(motif = ifelse(posID %in% names(unlist(snvs_rgyw_vmatch)), "RGYW",
                        ifelse(posID %in% names(unlist(snvs_wrcy_vmatch)),
                               "WRCY", "None"))) %>%
  left_join(snvs_rgyw_seq_df) %>%
  left_join(snvs_wrcy_seq_df)

mutations_maf_motif <- mutations_maf_motif %>%
  arrange(SAMPLE_ID, Chromosome, Start_Position)

# column motif will be used to identify SHM. confirm motif annotation is correct.
#  R   G   Y    W
# A/G  G  C/T  A/T
table(mutations_maf_motif$RGYW.Motif, mutations_maf_motif$motif)

#  W    R   C   Y
# A/T  A/G  C  C/T
table(mutations_maf_motif$WRCY.Motif, mutations_maf_motif$motif)

number_SHM_df <- mutations_maf_motif %>% 
  filter(!motif == 'None') %>% 
  group_by(SAMPLE_ID) %>% 
  # count() %>% 
  summarise(n = n()) %>%
  ungroup() %>% 
  as.data.frame()

colnames(number_SHM_df)[2] <- 'n_SHM'

# Create dataframe for samples without SHM
samples_no_SHM <- epitypes_FL_FLOMICS_df$SAMPLE_ID[!epitypes_FL_FLOMICS_df$SAMPLE_ID %in% number_SHM_df$SAMPLE_ID]
samples_no_SHM_df <- data.frame(SAMPLE_ID = samples_no_SHM,
                                n_SHM = rep(0, length(samples_no_SHM)))

# Combine dataframes
number_SHM_df <- number_SHM_df %>% 
  rbind(samples_no_SHM_df) %>% 
  mutate(epitype = case_when(SAMPLE_ID %in% FLOMICS_C1_samples ~ 'C1', 
                             SAMPLE_ID %in% FLOMICS_C2_samples ~ 'C2'))

# Confirm that the correct number of samples are included 
number_SHM_df$SAMPLE_ID %>% length() == epitypes_FL_FLOMICS_df$SAMPLE_ID %>% length()


p_boxplot_SHM <- number_SHM_df %>% 
  ggplot(aes(x = epitype, y = n_SHM, fill = epitype)) + 
  geom_boxplot(outliers = FALSE) +
  geom_jitter(height = 0) + 
  stat_compare_means(comparisons = list(c('C1','C2')), method = 'wilcox', label = 'p.format', label.y = 9, size = 3) + 
  scale_fill_manual(values = c('C1'='#59A14F','C2'='#E15759')) + 
  theme_bw() + 
  theme(legend.position = 'none', 
        # axis.text.x = element_blank(),
        axis.title.y = element_text(size = 8),
        # axis.title.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  labs(y = 'Number of SHM mutations / sample') + 
  ylim(c(0,10))

pdf("Boxplot_SHM.pdf", width = 3, height = 5)
p_boxplot_SHM
dev.off()

p_dotplot_SHM <- number_SHM_df %>% 
  ggplot(aes(x = epitype, y = n_SHM, fill = epitype)) + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.2) + 
  stat_compare_means(comparisons = list(c('C1','C2')), method = 'wilcox', label = 'p.format', label.y = 9, size = 3) + 
  scale_fill_manual(values = c('C1'='#59A14F','C2'='#E15759')) + 
  theme_bw() + 
  theme(legend.position = 'none', 
        # axis.text.x = element_blank(),
        axis.title.y = element_text(size = 8, face = "bold"),
        # axis.title.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  labs(y = 'Number of SHM mutations per sample') + 
  ylim(c(0,10))

pdf("Dotplot_SHM.pdf", width = 3, height = 5)
p_dotplot_SHM
dev.off()

# Compare the proportion of samples with at least 1 SHM across epitypes
number_SHM_df <- number_SHM_df %>% 
  mutate(has_SHM = case_when(n_SHM > 0 ~ 'Y', 
                             n_SHM == 0 ~ 'N'))

table(number_SHM_df$epitype, number_SHM_df$has_SHM)

fisher.test(number_SHM_df$epitype, number_SHM_df$has_SHM)

# Combine barplot with boxplots
# =============================
# Define layout using area(): area(row_start, col_start, row_end, col_end)
layout <- c(area(1, 1, 3.5, 15), 
            area(1, 10, 2, 12), 
            area(1, 13, 2, 15))

# Arrange the plots with overlay
pdf("Final_Plot.pdf", width = 10, height = 5)
p_barplot + p_boxplot_mutations + p_boxplot_mutated_genes + plot_layout(design = layout)
dev.off()

# CREBBP KAT missense mutations
# =============================
# ENSG00000005339:ENST00000262367:exon26:c.C4336T:p.R1446C, isoform 1
# ENSG00000005339:ENST00000382070:exon25:c.C4222T:p.R1408C, isoform 2

# R1446 isoform 1 (longer) = R1408 isoform 2 (shorter), chr16:3788617-3788618
# Based on Garcia Ramirez HAT is from 1342 to 1649 in isoform 1
# V1342 isoform 1 (longer) = V1304 isoform 2 (shorter), chr16:3790509
# I1649 isoform 1 (longer) = I1611 isoform 2 (shorter), chr16:3781418

# For isoform 2, HAT domain is from 1304 to 1611
# Here: mutations are called based on isoform 2 (R1408 is present in the FLOMICS data and a common variant)

# https://www.pnas.org/doi/10.1073/pnas.1501199112
# (R1408 in isoform b, R1446 in isoform a; SI Appendix, Fig. S9)
KAT_domain_AAs_range <- seq(1304, 1611)

# Load the mutation data
CREBBP_KAT_domain_mutations <- read_excel("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/12_SNP/Scripts/01_Davidson/data/41408_2024_1111_MOESM2_ESM-2.xlsx", sheet = 5) %>% 
  filter(SAMPLE_ID %in% epitypes_FL_FLOMICS_df$SAMPLE_ID) %>% 
  filter(Variant_Classification == 'Missense_Mutation') %>% 
  filter(Hugo_Symbol == 'CREBBP') %>% 
  filter(grepl(paste(KAT_domain_AAs_range, collapse = '|'), Protein_Change))

# Check for the most common CREBBP KAT mutation (R1408) - isoform 2
CREBBP_KAT_domain_mutations$Protein_Change %>% table() %>% sort()

# Check that KAT domain missense mutations in the correct range
KAT_domain_AAs <- CREBBP_KAT_domain_mutations %>% 
  mutate(AA = str_sub(Protein_Change, 2)) %>% 
  mutate(AA = str_sub(AA, 1, -2))
KAT_domain_AAs$AA %>% str_sort(numeric = TRUE)

# Check that non-KAT domain missense mutations are excluded from KAT domain range
non_KAT_domain_AAs <- read_excel("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/12_SNP/Scripts/01_Davidson/data/41408_2024_1111_MOESM2_ESM-2.xlsx", sheet = 5) %>% 
  filter(SAMPLE_ID %in% epitypes_FL_FLOMICS_df$SAMPLE_ID) %>% 
  filter(Variant_Classification == 'Missense_Mutation') %>% 
  filter(Hugo_Symbol == 'CREBBP') %>% 
  filter(!grepl(paste(KAT_domain_AAs_range, collapse = '|'), Protein_Change)) %>% 
  mutate(AA = str_sub(Protein_Change, 2)) %>% 
  mutate(AA = str_sub(AA, 1, -2))
non_KAT_domain_AAs$AA %>% str_sort(numeric = TRUE)

epitypes_FL_FLOMICS_CREBBP_KAT_df <- epitypes_FL_FLOMICS_df %>% 
  mutate(CREBBP_KAT_mutation = case_when(SAMPLE_ID %in% CREBBP_KAT_domain_mutations$SAMPLE_ID ~ 'Y', 
                                         TRUE ~ 'N'))

table(epitypes_FL_FLOMICS_CREBBP_KAT_df$epitypes, epitypes_FL_FLOMICS_CREBBP_KAT_df$CREBBP_KAT_mutation)

fisher.test(epitypes_FL_FLOMICS_CREBBP_KAT_df$epitypes, epitypes_FL_FLOMICS_CREBBP_KAT_df$CREBBP_KAT_mutation)

# Package versions
# ================
sessionInfo()
