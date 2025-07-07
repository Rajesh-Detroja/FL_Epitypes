# R Packages
# ==========
library(MethylMix)
library(dplyr)
library(stringr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(annotables)
library(WGCNA)
library(parallel)
library(doParallel)
library(ggplot2)
library(cowplot)
library(enrichR)
library(forcats)
library(DGEobj.utils)
library(edgeR)

# Change directory
setwd("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/10_MethylMix/84110_CpGs/cluster_run/C2_vs_C1/RD_MethylMix_Analysis/")

# MethylMix Analysis
# ==================
# MethylationDrivers: Genes identified as transcriptionally predictive and differentially methylated by MethylMix.
# NrComponents: The number of methylation states found for each driver gene.
# MixtureStates: A list with the DM-values for each driver gene.
# MethylationStates: Matrix with DM-values for all driver genes (rows) and all samples (columns).
# Classifications: Matrix with integers indicating to which mixture component each cancer sample was assigned to, for each gene.
# Models: Beta mixture model parameters for each driver gene.
# Differential Methylation values (DM-values) are defined as the difference between the methylation mean in one mixture component of cancer samples and the methylation mean in the normal samples, for a given gene.

METcancer <- readRDS("../METcancer.rds")
class(METcancer)

METnormal <- readRDS("../METnormal.rds")
class(METnormal)

GEcancer <- readRDS("../GEcancer.rds")
class(GEcancer)

MethylMixResults <- readRDS("../MethylMixResults.rds")
class(MethylMixResults)

# The output from the MethylMix
MethylMixResults$MethylationDrivers

results <- data.frame(MethylationDrivers = MethylMixResults$MethylationDrivers,
                      Symbol = gsub("(.*)---.*", "\\1", MethylMixResults$MethylationDrivers),
                      Nb_components = MethylMixResults$NrComponents,
                      Meth_aFL = rowSums(METcancer[MethylMixResults$MethylationDrivers,])/ncol(METcancer),
                      Meth_iFL = rowSums(METnormal[MethylMixResults$MethylationDrivers,])/ncol(METnormal))

results <- results %>% mutate(Methylation_status_aFL = ifelse(Meth_aFL > Meth_iFL, "HYPER", "HYPO"))
write.table(results, file="MethylMix_Results.tsv", sep="\t", quote = FALSE, row.names = FALSE)

hypermethylated_genes <- unique(results %>% filter(Methylation_status_aFL == "HYPER") %>% .$Symbol)
hypomethylated_genes <- unique(results %>% filter(Methylation_status_aFL == "HYPO") %>% .$Symbol)
length(hypermethylated_genes) # n = 28
length(hypomethylated_genes) # n = 21

write.table(hypermethylated_genes, file="HYPER_Genes.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(hypomethylated_genes, file="HYPO_Genes.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Perform pathway analysis
dbs <- c("BioPlanet_2019",
         "GO_Biological_Process_2021", "GO_Molecular_Function_2021",
         "WikiPathway_2021_Human", "KEGG_2021_Human", "Reactome_2016")

enrichr_results_hyper <- enrichr(hypermethylated_genes, databases = dbs)
enrichr_results_hyper <- do.call("rbind", enrichr_results_hyper)
colnames(enrichr_results_hyper) <- c("Term",	"Overlap",	"P-value",	"Adjusted P-value",	"Old P-value",	"Old Adjusted P-value",	"Odds Ratio",	"Combined Score", "Genes")
enrichr_results_hyper <- enrichr_results_hyper %>%
  mutate(nb.genes.overlap = stringr::str_extract(enrichr_results_hyper$Overlap, "[^/]+")) %>%
  mutate(nb.genes.overlap = as.numeric(nb.genes.overlap)) %>%
  filter(nb.genes.overlap > 2) %>%
  select(-nb.genes.overlap)
write.table(enrichr_results_hyper, file="HYPER_EnrichR.tsv", sep="\t", quote = FALSE, row.names = FALSE)

enrichr_results_hypo <- enrichr(hypomethylated_genes, databases = dbs)
enrichr_results_hypo <- do.call("rbind", enrichr_results_hypo)
colnames(enrichr_results_hypo) <- c("Term",	"Overlap",	"P-value",	"Adjusted P-value",	"Old P-value",	"Old Adjusted P-value",	"Odds Ratio",	"Combined Score", "Genes")
enrichr_results_hypo <- enrichr_results_hypo %>%
  mutate(nb.genes.overlap = stringr::str_extract(enrichr_results_hypo$Overlap, "[^/]+")) %>%
  mutate(nb.genes.overlap = as.numeric(nb.genes.overlap)) %>%
  filter(nb.genes.overlap > 2) %>%
  select(-nb.genes.overlap)
write.table(enrichr_results_hypo, file="HYPO_EnrichR.tsv", sep="\t", quote = FALSE, row.names = FALSE)

# pdf("HYPER_EnrichR.pdf", width = 15, height = 5)
# enrichr_results_hyper %>%
p1 <- enrichr_results_hyper %>%
  select(Term, adj.P.value = "Adjusted P-value") %>%
  mutate(adj.P.value = -log10(as.numeric(adj.P.value))) %>%
  group_by(Term) %>%
  top_n(1, adj.P.value) %>%
  ungroup() %>%
  mutate(Term = fct_reorder(Term, adj.P.value)) %>%
  top_n(20, adj.P.value) %>%
  ggplot(aes(x = Term, y = adj.P.value, colour = adj.P.value, size = adj.P.value)) + 
  geom_point(stat = "identity") +
  coord_flip() +
  ylab("-log10(adjusted P value)") +
  theme(axis.title.x = element_blank()) +
  ggtitle("Hyper Methylated Genes") +
  ylim(0,20) +
  theme_bw() 
# dev.off()

# pdf("HYPO_EnrichR.pdf", width = 15, height = 5)
#enrichr_results_hypo %>%
p2 <- enrichr_results_hypo %>%
  select(Term, adj.P.value = "Adjusted P-value") %>%
  mutate(adj.P.value = -log10(as.numeric(adj.P.value))) %>%
  group_by(Term) %>%
  top_n(1, adj.P.value) %>%
  ungroup() %>%
  mutate(Term = fct_reorder(Term, adj.P.value)) %>%
  top_n(20, adj.P.value) %>%
  ggplot(aes(x = Term, y = adj.P.value, colour = adj.P.value, size = adj.P.value)) + 
  geom_point(stat = "identity") +
  coord_flip() +
  ylab("-log10(adjusted P value)") +
  theme(axis.title.x = element_blank()) + 
  ggtitle("Hypo Methylated Genes") +
  ylim(0,20) +
  theme_bw()
# dev.off()

p <- plot_grid(p1, p2, labels = c("A", "B"), ncol = 1, nrow = 2)
ggsave("EnrichR_Plot.pdf", width = 12, height = 10)

# Gene analysis
# suppressWarnings(
# plots <- MethylMix_PlotModel("LMO2---Cluster2", MethylMixResults, METcancer, GEcancer, METnormal)
# )
# 
# p1 <- plots$MixtureModelPlot
# p2 <- plots$CorrelationPlot
#  
# p <- plot_grid(p1, p2, labels = c("A", "B"), ncol = 2)
# ggsave("MethylMix_LMO2_Cluster2.pdf", width = 12, height = 6)

# # Plot all functional and differential genes
# for (gene in MethylMixResults$MethylationDrivers) {
#   MethylMix_PlotModel(gene, MethylMixResults, METcancer, METnormal = METnormal)
# }
