# Load libraries
# ==============
pacman::p_load(tidyverse, readxl, survival, survminer, ggpubr, patchwork, scales)

# Change directory
setwd("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/11_Survival/450K/modified_InfiniumPurify/Enmix_quantile_with_QCinfo_Harman_All/")

# Load the data
clinical <- read.csv('/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/00_Minfi/db/CLN/20240116_clinical_data.csv')
FL_epitypes <- read.delim('/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/04_Summary/450k/modified_InfiniumPurify/Enmix_quantile_with_QCinfo_Harman_All/03_Sample_Overview/Epitypes_Info.txt')
# FL_epitypes <- read.delim('/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/02_Clustering/450k/modified_InfiniumPurify/Enmix_quantile_with_QCinfo_Harman_All/FL_epitypes.tsv')

FL_epitypes <- subset(FL_epitypes, TYPE == "FL" & epitypes %in% c("C1", "C2"))
FL_epitypes$epitypes <- gsub("C1", "iFL", FL_epitypes$epitypes)
FL_epitypes$epitypes <- gsub("C2", "aFL", FL_epitypes$epitypes)
FL_epitypes$epitypes <- factor(FL_epitypes$epitypes, levels = c("iFL", "aFL")) # reorder
table(FL_epitypes$epitypes)

FL_epitypes_clinical <- merge(x = FL_epitypes, y = clinical,
                              by.x = 'SAMPLE_ID', by.y = 'SAMPLE_ID', 
                              all.x = FALSE) %>% filter(epitypes %in% c('iFL', 'aFL'))
table(FL_epitypes_clinical$epitypes)

# KM Plots - All Patients - OS
fit <- survfit(Surv(OS, CODE_OS) ~ epitypes, data = FL_epitypes_clinical)
summary(fit, time = 5)

names(fit$strata) <- gsub('epitypes=', '', names(fit$strata))
p <-  ggsurvplot(fit,
                 pval = TRUE,
                 pval.size = 5,
                 risk.table = TRUE,
                 title = 'All Patients',
                 font.title = c(16, "bold"),
                 font.x = c(14, "bold"),
                 font.y = c(14, "bold"),
                 font.tickslab = c(12, "bold"),
                 font.legend = list(size = 12),
                 xlim = c(0,10),
                 break.x.by = 2,
                 palette = 'lancet',
                 risk.table.y.text = FALSE,
                 risk.table.col = 'strata',
                 risk.table.height = 0.2,
                 tables.theme = theme_cleantable(),
                 legend.title = 'Epitypes',
                 xlab = "Time (years)",
                 ylab = "Overall Survival",
                 ggtheme = theme_classic2(base_size=12), 
                 theme = theme_minimal() + 
                   theme(plot.margin = unit(c(0, 0, 0, 0), 'cm')))

# Centering the title
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5))

pdf("KM_OS_All.pdf", width = 5, height = 5)
print(p, newpage = FALSE)
dev.off()

### All Patients ###
# KM Plots - All Patients - PFS
fit <- survfit(Surv(PFS, CODE_PFS) ~ epitypes, data = FL_epitypes_clinical)
summary(fit, time = 5)

names(fit$strata) <- gsub('epitypes=', '', names(fit$strata))
p <-  ggsurvplot(fit,
                 pval = TRUE,
                 pval.size = 5,
                 risk.table = TRUE,
                 title = 'All Patients',
                 font.title = c(16, "bold"),
                 font.x = c(14, "bold"),
                 font.y = c(14, "bold"),
                 font.tickslab = c(12, "bold"),
                 font.legend = list(size = 12),
                 xlim = c(0,10),
                 break.x.by = 2,
                 palette = 'lancet',
                 risk.table.y.text = FALSE,
                 risk.table.col = 'strata',
                 risk.table.height = 0.2,
                 tables.theme = theme_cleantable(),
                 legend.title = 'Epitypes',
                 xlab = "Time (years)",
                 ylab = "Progression Free Survival",
                 ggtheme = theme_classic2(base_size=12), 
                 theme = theme_minimal() + 
                   theme(plot.margin = unit(c(0, 0, 0, 0), 'cm')))

# Centering the title
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5))

pdf("KM_PFS_All.pdf", width = 5, height = 5)
print(p, newpage = FALSE)
dev.off()


# KM Plots - All Patients - TTT
fit <- survfit(Surv(TTT, CODE_TRANSF) ~ epitypes, data = FL_epitypes_clinical)
summary(fit, time = 5)

names(fit$strata) <- gsub('epitypes=', '', names(fit$strata))
p <-  ggsurvplot(fit,
                 pval = TRUE,
                 pval.size = 5,
                 risk.table = TRUE,
                 title = 'All Patients',
                 font.title = c(16, "bold"),
                 font.x = c(14, "bold"),
                 font.y = c(14, "bold"),
                 font.tickslab = c(12, "bold"),
                 font.legend = list(size = 12),
                 xlim = c(0,10),
                 break.x.by = 2,
                 fun = 'event',
                 palette = 'lancet',
                 risk.table.y.text = FALSE,
                 risk.table.col = 'strata',
                 risk.table.height = 0.2,
                 tables.theme = theme_cleantable(),
                 legend.title = 'Epitypes',
                 xlab = "Time (years)",
                 ylab = "Time to Transformation",
                 ggtheme = theme_classic2(base_size=12), 
                 theme = theme_minimal() + 
                   theme(plot.margin = unit(c(0, 0, 0, 0), 'cm')))

# Centering the title
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5))

pdf("KM_TTT_All.pdf", width = 5, height = 5)
print(p, newpage = FALSE)
dev.off()


# POD24 status
# ============
# Define POD24
FL_epitypes_clinical <- FL_epitypes_clinical %>%
  mutate(
    POD24 = case_when(
      CODE_PFS == "1" & PFS < 2 ~ "YES",
      PFS >= 2 ~ "NO",
      TRUE ~ NA_character_
      )
    )

# Check POD24 counts per epitype
table(FL_epitypes_clinical$POD24)
table(FL_epitypes_clinical$epitypes, FL_epitypes_clinical$POD24, useNA = "ifany")

# Test for association
# Use a Chi-square test (if all expected counts are ≥5) or Fisher’s Exact Test (if some are small)

# Chi-square test
chisq.test(table(FL_epitypes_clinical$epitypes, FL_epitypes_clinical$POD24))

# Pearson's Chi-squared test with Yates' continuity correction
# data:  table(FL_epitypes_clinical$epitypes, FL_epitypes_clinical$POD24)
# X-squared = 6.3953, df = 1, p-value = 0.01144

# fisher.test(table(df$Epitype, df$POD24))
fisher.test(table(FL_epitypes_clinical$epitypes, FL_epitypes_clinical$POD24))

# Fisher's Exact Test for Count Data
# data:  table(FL_epitypes_clinical$epitypes, FL_epitypes_clinical$POD24)
# p-value = 0.009048
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.190035 5.091118
# sample estimates:
# odds ratio 
#   2.475432

# Visualize it
# Save plot as PDF
plot_data <- FL_epitypes_clinical %>%
  filter(!is.na(POD24)) %>%
  group_by(epitypes, POD24) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(epitypes) %>%
  mutate(
    Proportion = Count / sum(Count),
    Proportion_Label = paste0(Count, " (", percent(Proportion, accuracy = 0.1), ")")
  )

pdf("POD24_All.pdf", width = 6.5, height = 6)
ggplot(plot_data, aes(x = epitypes, y = Proportion, fill = POD24)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = Proportion_Label), 
            position = position_stack(vjust = 0.5), size = 4.5, color = "white") +
  ylab("Proportion") +
  ggtitle("All Patients") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    axis.line = element_line(),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(0.8, "cm"),
    legend.box.spacing = unit(0.2, "cm"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.margin = margin(0.5, 0.8, 0.5, 0.5, "cm")
  )
dev.off()

### R-Chemo ###
# FL_epitypes_clinical_R_CHEMO <- FL_epitypes_clinical %>% 
#   filter(PRIM_TX_CAT %in% c('BR', 'R-CHOP', 'R-CVP'))

# check the treatment grouping:
table(FL_epitypes_clinical$PRIM_TX_CAT, FL_epitypes_clinical$SURVIVAL_COHORT)

FL_epitypes_clinical_R_CHEMO <- FL_epitypes_clinical %>% 
  filter(SURVIVAL_COHORT == 'R_CHEMO')

# KM Plots - Chemotherapy - OS
fit <- survfit(Surv(OS, CODE_OS) ~ epitypes, data = FL_epitypes_clinical_R_CHEMO)
summary(fit, time = 5)

names(fit$strata) <- gsub('epitypes=', '', names(fit$strata))
p <-  ggsurvplot(fit,
                 pval = TRUE,
                 pval.size = 5,
                 risk.table = TRUE,
                 title = 'Chemoimmunotherapy Patients',
                 font.title = c(16, "bold"),
                 font.x = c(14, "bold"),
                 font.y = c(14, "bold"),
                 font.tickslab = c(12, "bold"),
                 font.legend = list(size = 12),
                 xlim = c(0,10),
                 break.x.by = 2,
                 palette = 'lancet',
                 risk.table.y.text = FALSE,
                 risk.table.col = 'strata',
                 risk.table.height = 0.2,
                 tables.theme = theme_cleantable(),
                 legend.title = 'Epitypes',
                 xlab = "Time (years)",
                 ylab = "Overall Survival",
                 ggtheme = theme_classic2(base_size=12), 
                 theme = theme_minimal() + 
                   theme(plot.margin = unit(c(0, 0, 0, 0), 'cm')))

# Centering the title
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5))

pdf("KM_OS_Chemo.pdf", width = 5, height = 5)
print(p, newpage = FALSE)
dev.off()


# KM Plots - All Patients - PFS
fit <- survfit(Surv(PFS, CODE_PFS) ~ epitypes, data = FL_epitypes_clinical_R_CHEMO)
summary(fit, time = 5)

names(fit$strata) <- gsub('epitypes=', '', names(fit$strata))
p <-  ggsurvplot(fit,
                 pval = TRUE,
                 pval.size = 5,
                 risk.table = TRUE,
                 title = 'Chemoimmunotherapy Patients',
                 font.title = c(16, "bold"),
                 font.x = c(14, "bold"),
                 font.y = c(14, "bold"),
                 font.tickslab = c(12, "bold"),
                 font.legend = list(size = 12),
                 xlim = c(0,10),
                 break.x.by = 2,
                 palette = 'lancet',
                 risk.table.y.text = FALSE,
                 risk.table.col = 'strata',
                 risk.table.height = 0.2,
                 tables.theme = theme_cleantable(),
                 legend.title = 'Epitypes',
                 xlab = "Time (years)",
                 ylab = "Progression Free Survival",
                 ggtheme = theme_classic2(base_size=12), 
                 theme = theme_minimal() + 
                   theme(plot.margin = unit(c(0, 0, 0, 0), 'cm')))

# Centering the title
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5))

pdf("KM_PFS_Chemo.pdf", width = 5, height = 5)
print(p, newpage = FALSE)
dev.off()


# KM Plots - All Patients - TTT
fit <- survfit(Surv(TTT, CODE_TRANSF) ~ epitypes, data = FL_epitypes_clinical_R_CHEMO)
summary(fit, time = 5)

names(fit$strata) <- gsub('epitypes=', '', names(fit$strata))
p <-  ggsurvplot(fit,
                 pval = TRUE,
                 pval.size = 5,
                 risk.table = TRUE,
                 title = 'Chemoimmunotherapy Patients',
                 font.title = c(16, "bold"),
                 font.x = c(14, "bold"),
                 font.y = c(14, "bold"),
                 font.tickslab = c(12, "bold"),
                 font.legend = list(size = 12),
                 xlim = c(0,10),
                 break.x.by = 2,
                 fun = 'event',
                 palette = 'lancet',
                 risk.table.y.text = FALSE,
                 risk.table.col = 'strata',
                 risk.table.height = 0.2,
                 tables.theme = theme_cleantable(),
                 legend.title = 'Epitypes',
                 xlab = "Time (years)",
                 ylab = "Time to Transformation",
                 ggtheme = theme_classic2(base_size=12), 
                 theme = theme_minimal() + 
                   theme(plot.margin = unit(c(0, 0, 0, 0), 'cm')))

# Centering the title
p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5))

pdf("KM_TTT_Chemo.pdf", width = 5, height = 5)
print(p, newpage = FALSE)
dev.off()


# POD24 status
# ============
# Define POD24
FL_epitypes_clinical_R_CHEMO <- FL_epitypes_clinical_R_CHEMO %>%
  mutate(
    POD24 = case_when(
      CODE_PFS == "1" & PFS < 2 ~ "YES",
      PFS >= 2 ~ "NO",
      TRUE ~ NA_character_
    )
  )

# Check POD24 counts per epitype
table(FL_epitypes_clinical_R_CHEMO$POD24)
table(FL_epitypes_clinical_R_CHEMO$epitypes, FL_epitypes_clinical_R_CHEMO$POD24, useNA = "ifany")

# Test for association
# Use a Chi-square test (if all expected counts are ≥5) or Fisher’s Exact Test (if some are small)

# Chi-square test
chisq.test(table(FL_epitypes_clinical_R_CHEMO$epitypes, FL_epitypes_clinical_R_CHEMO$POD24))

# Pearson's Chi-squared test with Yates' continuity correction
# data:  table(FL_epitypes_clinical_R_CHEMO$epitypes, FL_epitypes_clinical_R_CHEMO$POD24)
# X-squared = 7.1832, df = 1, p-value = 0.007359

# fisher.test
fisher.test(table(FL_epitypes_clinical_R_CHEMO$epitypes, FL_epitypes_clinical_R_CHEMO$POD24))
# Fisher's Exact Test for Count Data
# data:  table(FL_epitypes_clinical_R_CHEMO$epitypes, FL_epitypes_clinical_R_CHEMO$POD24)
# p-value = 0.005563
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.307896 7.549518
# sample estimates:
#   odds ratio 
# 3.125124

# Visualize it
# Save plot as PDF
plot_data <- FL_epitypes_clinical_R_CHEMO %>%
  filter(!is.na(POD24)) %>%
  group_by(epitypes, POD24) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(epitypes) %>%
  mutate(
    Proportion = Count / sum(Count),
    Proportion_Label = paste0(Count, " (", percent(Proportion, accuracy = 0.1), ")")
  )

pdf("POD24_Chemo.pdf", width = 6.5, height = 6)
ggplot(plot_data, aes(x = epitypes, y = Proportion, fill = POD24)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = Proportion_Label), 
            position = position_stack(vjust = 0.5), size = 4.5, color = "white") +
  ylab("Proportion") +
  ggtitle("Chemoimmunotherapy Patients") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    axis.line = element_line(),
    axis.line.x.top = element_blank(),
    axis.line.y.right = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(0.8, "cm"),
    legend.box.spacing = unit(0.2, "cm"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.margin = margin(0.5, 0.8, 0.5, 0.5, "cm")
  )
dev.off()

# Forest plots - MV analysis for all patients
# ===========================================

# OS Forest Plot - FL Epitypes and FLIPI
# --------------------------------------
fit <- coxph(Surv(OS, CODE_OS) ~ epitypes + FLIPI_BINARY, data = FL_epitypes_clinical)
summary(fit)

coefficients <- summary(fit)$coefficients %>% as.data.frame()
coefficients$Feature <- rownames(coefficients)
colnames(coefficients)[colnames(coefficients) == 'exp(coef)'] <- 'HR'

conf_int <- summary(fit)$conf.int %>% as.data.frame()
conf_int$Feature <- rownames(conf_int)

COX_df <- merge(coefficients, conf_int, 
                by.x = 'Feature', by.y = 'Feature') %>% 
  mutate(HR_low = `lower .95`) %>% 
  mutate(HR_high = `upper .95`) %>% 
  mutate(significant = case_when(`Pr(>|z|)` < 0.05 ~ 'Yes', 
                                 TRUE ~ 'No')) %>% 
  mutate(Feature = case_when(Feature == 'epitypesaFL' ~ 'Epitype - aFL', 
                             Feature == 'FLIPI_BINARY' ~ 'FLIPI - High')) %>% 
  rename(p.value = `Pr(>|z|)`)

# wrangle results into pre-plotting table form
res_plot <- COX_df %>%
  mutate(across(c(HR, HR_low, HR_high), ~sprintf('%.2f', .x))) %>% 
  mutate(estimate_lab = paste0(HR)) %>% 
  mutate(p.value = case_when(p.value < .001 ~ '<0.001', 
                             TRUE ~ str_pad(as.character(round(p.value, 3)), 
                                            width = 4, pad = '0', side = 'right'))) |>
  bind_rows(data.frame(Feature = 'Feature', estimate_lab = 'HR', HR_low = '',
                       HR_high = '', p.value = 'p-value')) %>% 
  mutate(Feature = factor(Feature, levels = Feature))

p_forest <- COX_df %>% 
  mutate(Feature = factor(Feature, levels = levels(res_plot$Feature))) %>% 
  ggplot(aes(y = Feature, color = significant)) + 
  geom_point(aes(x = HR), shape = 15, size = 3) + 
  geom_linerange(aes(xmin = HR_low, xmax = HR_high)) + 
  geom_errorbar(aes(xmin = HR_low, xmax = HR_high), width = 0.2) + 
  coord_cartesian(ylim = c(1, nrow(COX_df)+1), xlim = c(-2, 12)) + 
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'grey') + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(), 
        axis.line.x = element_line(), 
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(), 
        legend.position = 'none') + 
  scale_color_manual(values = c('Yes' = 'red', 'No' = 'black'))
p_left <- res_plot %>% 
  ggplot(aes(y = Feature)) + 
  geom_text(aes(x = 0, label = Feature), hjust = 0, fontface = 'bold') +
  geom_text(aes(x = 1.8, label = estimate_lab), hjust = 0, fontface = ifelse(res_plot$estimate_lab == 'HR', 'bold', 'plain')) +
  theme_void() + 
  coord_cartesian(xlim = c(0,2.9))

# right side of plot - pvalues
p_right <- res_plot %>% 
  ggplot() +
  geom_text(aes(x = 0, y = Feature, label = p.value), 
            hjust=0, 
            fontface = ifelse(res_plot$p.value == 'p-value', 'bold', 'plain')) +
  theme_void() 

layout <- c(area(t = 0, l = 0, b = 30, r = 4),      
            area(t = 0, l = 4.5, b = 30, r = 6.5), 
            area(t = 0, l = 6.5, b = 30, r = 8))

p_MV <- p_left + p_forest + p_right + plot_layout(design = layout) +
  plot_annotation(title = 'Multivariable Analysis') & 
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, face = "bold"), plot.margin = unit(c(0.1, -0.5, 0.1, -0.1), "cm"))

pdf("MV_OS_All.pdf", width = 5, height = 1.5)
print(p_MV, newpage = FALSE)
dev.off()

rm(p_MV, p_left, p_forest, p_right, layout, res_plot, COX_df, coefficients, conf_int)

# PFS Forest Plot - FL Epitypes and FLIPI
# ---------------------------------------
fit <- coxph(Surv(PFS, CODE_PFS) ~ epitypes + FLIPI_BINARY, data = FL_epitypes_clinical)
summary(fit)

coefficients <- summary(fit)$coefficients %>% as.data.frame()
coefficients$Feature <- rownames(coefficients)
colnames(coefficients)[colnames(coefficients) == 'exp(coef)'] <- 'HR'

conf_int <- summary(fit)$conf.int %>% as.data.frame()
conf_int$Feature <- rownames(conf_int)

COX_df <- merge(coefficients, conf_int, 
                by.x = 'Feature', by.y = 'Feature') %>% 
  mutate(HR_low = `lower .95`) %>% 
  mutate(HR_high = `upper .95`) %>% 
  mutate(significant = case_when(`Pr(>|z|)` < 0.05 ~ 'Yes', 
                                 TRUE ~ 'No')) %>% 
  mutate(Feature = case_when(Feature == 'epitypesaFL' ~ 'Epitype - aFL', 
                             Feature == 'FLIPI_BINARY' ~ 'FLIPI - High')) %>% 
  rename(p.value = `Pr(>|z|)`)

# wrangle results into pre-plotting table form
res_plot <- COX_df %>%
  mutate(across(c(HR, HR_low, HR_high), ~sprintf('%.2f', .x))) %>% 
  mutate(estimate_lab = paste0(HR)) %>% 
  mutate(p.value = case_when(p.value < .001 ~ '<0.001', 
                             TRUE ~ str_pad(as.character(round(p.value, 3)), 
                                            width = 4, pad = '0', side = 'right'))) |>
  bind_rows(data.frame(Feature = 'Feature', estimate_lab = 'HR', HR_low = '',
                       HR_high = '', p.value = 'p-value')) %>% 
  mutate(Feature = factor(Feature, levels = Feature))

p_forest <- COX_df %>% 
  mutate(Feature = factor(Feature, levels = levels(res_plot$Feature))) %>% 
  ggplot(aes(y = Feature, color = significant)) + 
  geom_point(aes(x = HR), shape = 15, size = 3) + 
  geom_linerange(aes(xmin = HR_low, xmax = HR_high)) + 
  geom_errorbar(aes(xmin = HR_low, xmax = HR_high), width = 0.2) + 
  coord_cartesian(ylim = c(1, nrow(COX_df)+1), xlim = c(-2, 12)) + 
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'grey') + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(), 
        axis.line.x = element_line(), 
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(), 
        legend.position = 'none') + 
  scale_color_manual(values = c('Yes' = 'red', 'No' = 'black'))
p_left <- res_plot %>% 
  ggplot(aes(y = Feature)) + 
  geom_text(aes(x = 0, label = Feature), hjust = 0, fontface = 'bold') +
  geom_text(aes(x = 1.8, label = estimate_lab), hjust = 0, fontface = ifelse(res_plot$estimate_lab == 'HR', 'bold', 'plain')) +
  theme_void() + 
  coord_cartesian(xlim = c(0,2.9))

# right side of plot - pvalues
p_right <- res_plot %>% 
  ggplot() +
  geom_text(aes(x = 0, y = Feature, label = p.value), 
            hjust=0, 
            fontface = ifelse(res_plot$p.value == 'p-value', 'bold', 'plain')) +
  theme_void() 

layout <- c(area(t = 0, l = 0, b = 30, r = 4),      
            area(t = 0, l = 4.5, b = 30, r = 6.5), 
            area(t = 0, l = 6.5, b = 30, r = 8))

p_MV <- p_left + p_forest + p_right + plot_layout(design = layout) +
  plot_annotation(title = 'Multivariable Analysis') & 
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, face = "bold"), plot.margin = unit(c(0.1, -0.5, 0.1, -0.1), "cm"))

pdf("MV_PFS_All.pdf", width = 5, height = 1.5)
print(p_MV, newpage = FALSE)
dev.off()

rm(p_MV, p_left, p_forest, p_right, layout, res_plot, COX_df, coefficients, conf_int)

# TTT Forest Plot - FL Epitypes and FLIPI
# ---------------------------------------
fit <- coxph(Surv(TTT, CODE_TRANSF) ~ epitypes + FLIPI_BINARY, data = FL_epitypes_clinical)
summary(fit)

coefficients <- summary(fit)$coefficients %>% as.data.frame()
coefficients$Feature <- rownames(coefficients)
colnames(coefficients)[colnames(coefficients) == 'exp(coef)'] <- 'HR'

conf_int <- summary(fit)$conf.int %>% as.data.frame()
conf_int$Feature <- rownames(conf_int)

COX_df <- merge(coefficients, conf_int, 
                by.x = 'Feature', by.y = 'Feature') %>% 
  mutate(HR_low = `lower .95`) %>% 
  mutate(HR_high = `upper .95`) %>% 
  mutate(significant = case_when(`Pr(>|z|)` < 0.05 ~ 'Yes', 
                                 TRUE ~ 'No')) %>% 
  mutate(Feature = case_when(Feature == 'epitypesaFL' ~ 'Epitype - aFL', 
                             Feature == 'FLIPI_BINARY' ~ 'FLIPI - High')) %>% 
  rename(p.value = `Pr(>|z|)`)

# wrangle results into pre-plotting table form
res_plot <- COX_df %>%
  mutate(across(c(HR, HR_low, HR_high), ~sprintf('%.2f', .x))) %>% 
  mutate(estimate_lab = paste0(HR)) %>% 
  mutate(p.value = case_when(p.value < .001 ~ '<0.001', 
                             TRUE ~ str_pad(as.character(round(p.value, 3)), 
                                            width = 4, pad = '0', side = 'right'))) |>
  bind_rows(data.frame(Feature = 'Feature', estimate_lab = 'HR', HR_low = '',
                       HR_high = '', p.value = 'p-value')) %>% 
  mutate(Feature = factor(Feature, levels = Feature))

p_forest <- COX_df %>% 
  mutate(Feature = factor(Feature, levels = levels(res_plot$Feature))) %>% 
  ggplot(aes(y = Feature, color = significant)) + 
  geom_point(aes(x = HR), shape = 15, size = 3) + 
  geom_linerange(aes(xmin = HR_low, xmax = HR_high)) + 
  geom_errorbar(aes(xmin = HR_low, xmax = HR_high), width = 0.2) + 
  coord_cartesian(ylim = c(1, nrow(COX_df)+1), xlim = c(-2, 12)) + 
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'grey') + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(), 
        axis.line.x = element_line(), 
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(), 
        legend.position = 'none') + 
  scale_color_manual(values = c('Yes' = 'red', 'No' = 'black'))
p_left <- res_plot %>% 
  ggplot(aes(y = Feature)) + 
  geom_text(aes(x = 0, label = Feature), hjust = 0, fontface = 'bold') +
  geom_text(aes(x = 1.8, label = estimate_lab), hjust = 0, fontface = ifelse(res_plot$estimate_lab == 'HR', 'bold', 'plain')) +
  theme_void() + 
  coord_cartesian(xlim = c(0,2.9))

# right side of plot - pvalues
p_right <- res_plot %>% 
  ggplot() +
  geom_text(aes(x = 0, y = Feature, label = p.value), 
            hjust=0, 
            fontface = ifelse(res_plot$p.value == 'p-value', 'bold', 'plain')) +
  theme_void() 

layout <- c(area(t = 0, l = 0, b = 30, r = 4),      
            area(t = 0, l = 4.5, b = 30, r = 6.5), 
            area(t = 0, l = 6.5, b = 30, r = 8))

p_MV <- p_left + p_forest + p_right + plot_layout(design = layout) +
  plot_annotation(title = 'Multivariable Analysis') & 
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, face = "bold"), plot.margin = unit(c(0.1, -0.5, 0.1, -0.1), "cm"))

pdf("MV_TTT_All.pdf", width = 5, height = 1.5)
print(p_MV, newpage = FALSE)
dev.off()

rm(p_MV, p_left, p_forest, p_right, layout, res_plot, COX_df, coefficients, conf_int)


# Forest Plots - MV Analysis for Chemotherapy Patients
# ====================================================

# OS Forest Plot - FL Epitypes and FLIPI
# --------------------------------------
fit <- coxph(Surv(OS, CODE_OS) ~ epitypes + FLIPI_BINARY, data = FL_epitypes_clinical_R_CHEMO)
summary(fit)

coefficients <- summary(fit)$coefficients %>% as.data.frame()
coefficients$Feature <- rownames(coefficients)
colnames(coefficients)[colnames(coefficients) == 'exp(coef)'] <- 'HR'

conf_int <- summary(fit)$conf.int %>% as.data.frame()
conf_int$Feature <- rownames(conf_int)

COX_df <- merge(coefficients, conf_int, 
                by.x = 'Feature', by.y = 'Feature') %>% 
  mutate(HR_low = `lower .95`) %>% 
  mutate(HR_high = `upper .95`) %>% 
  mutate(significant = case_when(`Pr(>|z|)` < 0.05 ~ 'Yes', 
                                 TRUE ~ 'No')) %>% 
  mutate(Feature = case_when(Feature == 'epitypesaFL' ~ 'Epitype - aFL', 
                             Feature == 'FLIPI_BINARY' ~ 'FLIPI - High')) %>% 
  rename(p.value = `Pr(>|z|)`)

# wrangle results into pre-plotting table form
res_plot <- COX_df %>%
  mutate(across(c(HR, HR_low, HR_high), ~sprintf('%.2f', .x))) %>% 
  mutate(estimate_lab = paste0(HR)) %>% 
  mutate(p.value = case_when(p.value < .001 ~ '<0.001', 
                             TRUE ~ str_pad(as.character(round(p.value, 3)), 
                                            width = 4, pad = '0', side = 'right'))) |>
  bind_rows(data.frame(Feature = 'Feature', estimate_lab = 'HR', HR_low = '',
                       HR_high = '', p.value = 'p-value')) %>% 
  mutate(Feature = factor(Feature, levels = Feature))

p_forest <- COX_df %>% 
  mutate(Feature = factor(Feature, levels = levels(res_plot$Feature))) %>% 
  ggplot(aes(y = Feature, color = significant)) + 
  geom_point(aes(x = HR), shape = 15, size = 3) + 
  geom_linerange(aes(xmin = HR_low, xmax = HR_high)) + 
  geom_errorbar(aes(xmin = HR_low, xmax = HR_high), width = 0.2) + 
  coord_cartesian(ylim = c(1, nrow(COX_df)+1), xlim = c(-2, 12)) + 
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'grey') + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(), 
        axis.line.x = element_line(), 
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(), 
        legend.position = 'none') + 
  scale_color_manual(values = c('Yes' = 'red', 'No' = 'black'))
p_left <- res_plot %>% 
  ggplot(aes(y = Feature)) + 
  geom_text(aes(x = 0, label = Feature), hjust = 0, fontface = 'bold') +
  geom_text(aes(x = 1.8, label = estimate_lab), hjust = 0, fontface = ifelse(res_plot$estimate_lab == 'HR', 'bold', 'plain')) +
  theme_void() + 
  coord_cartesian(xlim = c(0,2.9))

# right side of plot - pvalues
p_right <- res_plot %>% 
  ggplot() +
  geom_text(aes(x = 0, y = Feature, label = p.value), 
            hjust=0, 
            fontface = ifelse(res_plot$p.value == 'p-value', 'bold', 'plain')) +
  theme_void() 

layout <- c(area(t = 0, l = 0, b = 30, r = 4),      
            area(t = 0, l = 4.5, b = 30, r = 6.5), 
            area(t = 0, l = 6.5, b = 30, r = 8))

p_MV <- p_left + p_forest + p_right + plot_layout(design = layout) +
  plot_annotation(title = 'Multivariable Analysis') & 
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, face = "bold"), plot.margin = unit(c(0.1, -0.5, 0.1, -0.1), "cm"))

pdf("MV_OS_Chemo.pdf", width = 5, height = 1.5)
print(p_MV, newpage = FALSE)
dev.off()

rm(p_MV, p_left, p_forest, p_right, layout, res_plot, COX_df, coefficients, conf_int)


# PFS Forest Plot - FL Epitypes and FLIPI
# ---------------------------------------
fit <- coxph(Surv(PFS, CODE_PFS) ~ epitypes + FLIPI_BINARY, data = FL_epitypes_clinical_R_CHEMO)
summary(fit)

coefficients <- summary(fit)$coefficients %>% as.data.frame()
coefficients$Feature <- rownames(coefficients)
colnames(coefficients)[colnames(coefficients) == 'exp(coef)'] <- 'HR'

conf_int <- summary(fit)$conf.int %>% as.data.frame()
conf_int$Feature <- rownames(conf_int)

COX_df <- merge(coefficients, conf_int, 
                by.x = 'Feature', by.y = 'Feature') %>% 
  mutate(HR_low = `lower .95`) %>% 
  mutate(HR_high = `upper .95`) %>% 
  mutate(significant = case_when(`Pr(>|z|)` < 0.05 ~ 'Yes', 
                                 TRUE ~ 'No')) %>% 
  mutate(Feature = case_when(Feature == 'epitypesaFL' ~ 'Epitype - aFL', 
                             Feature == 'FLIPI_BINARY' ~ 'FLIPI - High')) %>% 
  rename(p.value = `Pr(>|z|)`)

# wrangle results into pre-plotting table form
res_plot <- COX_df %>%
  mutate(across(c(HR, HR_low, HR_high), ~sprintf('%.2f', .x))) %>% 
  mutate(estimate_lab = paste0(HR)) %>% 
  mutate(p.value = case_when(p.value < .001 ~ '<0.001', 
                             TRUE ~ str_pad(as.character(round(p.value, 3)), 
                                            width = 4, pad = '0', side = 'right'))) |>
  bind_rows(data.frame(Feature = 'Feature', estimate_lab = 'HR', HR_low = '',
                       HR_high = '', p.value = 'p-value')) %>% 
  mutate(Feature = factor(Feature, levels = Feature))

p_forest <- COX_df %>% 
  mutate(Feature = factor(Feature, levels = levels(res_plot$Feature))) %>% 
  ggplot(aes(y = Feature, color = significant)) + 
  geom_point(aes(x = HR), shape = 15, size = 3) + 
  geom_linerange(aes(xmin = HR_low, xmax = HR_high)) + 
  geom_errorbar(aes(xmin = HR_low, xmax = HR_high), width = 0.2) + 
  coord_cartesian(ylim = c(1, nrow(COX_df)+1), xlim = c(-2, 12)) + 
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'grey') + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(), 
        axis.line.x = element_line(), 
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(), 
        legend.position = 'none') + 
  scale_color_manual(values = c('Yes' = 'red', 'No' = 'black'))
p_left <- res_plot %>% 
  ggplot(aes(y = Feature)) + 
  geom_text(aes(x = 0, label = Feature), hjust = 0, fontface = 'bold') +
  geom_text(aes(x = 1.8, label = estimate_lab), hjust = 0, fontface = ifelse(res_plot$estimate_lab == 'HR', 'bold', 'plain')) +
  theme_void() + 
  coord_cartesian(xlim = c(0,2.9))

# right side of plot - pvalues
p_right <- res_plot %>% 
  ggplot() +
  geom_text(aes(x = 0, y = Feature, label = p.value), 
            hjust=0, 
            fontface = ifelse(res_plot$p.value == 'p-value', 'bold', 'plain')) +
  theme_void() 

layout <- c(area(t = 0, l = 0, b = 30, r = 4),      
            area(t = 0, l = 4.5, b = 30, r = 6.5), 
            area(t = 0, l = 6.5, b = 30, r = 8))

p_MV <- p_left + p_forest + p_right + plot_layout(design = layout) +
  plot_annotation(title = 'Multivariable Analysis') & 
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, face = "bold"), plot.margin = unit(c(0.1, -0.5, 0.1, -0.1), "cm"))

pdf("MV_PFS_Chemo.pdf", width = 5, height = 1.5)
print(p_MV, newpage = FALSE)
dev.off()

rm(p_MV, p_left, p_forest, p_right, layout, res_plot, COX_df, coefficients, conf_int)


# TTT Forest Plot - FL Epitypes and FLIPI
# ---------------------------------------
fit <- coxph(Surv(TTT, CODE_TRANSF) ~ epitypes + FLIPI_BINARY, data = FL_epitypes_clinical_R_CHEMO)
summary(fit)

coefficients <- summary(fit)$coefficients %>% as.data.frame()
coefficients$Feature <- rownames(coefficients)
colnames(coefficients)[colnames(coefficients) == 'exp(coef)'] <- 'HR'

conf_int <- summary(fit)$conf.int %>% as.data.frame()
conf_int$Feature <- rownames(conf_int)

COX_df <- merge(coefficients, conf_int, 
                by.x = 'Feature', by.y = 'Feature') %>% 
  mutate(HR_low = `lower .95`) %>% 
  mutate(HR_high = `upper .95`) %>% 
  mutate(significant = case_when(`Pr(>|z|)` < 0.05 ~ 'Yes', 
                                 TRUE ~ 'No')) %>% 
  mutate(Feature = case_when(Feature == 'epitypesaFL' ~ 'Epitype - aFL', 
                             Feature == 'FLIPI_BINARY' ~ 'FLIPI - High')) %>% 
  rename(p.value = `Pr(>|z|)`)

# wrangle results into pre-plotting table form
res_plot <- COX_df %>%
  mutate(across(c(HR, HR_low, HR_high), ~sprintf('%.2f', .x))) %>% 
  mutate(estimate_lab = paste0(HR)) %>% 
  mutate(p.value = case_when(p.value < .001 ~ '<0.001', 
                             TRUE ~ str_pad(as.character(round(p.value, 3)), 
                                            width = 4, pad = '0', side = 'right'))) |>
  bind_rows(data.frame(Feature = 'Feature', estimate_lab = 'HR', HR_low = '',
                       HR_high = '', p.value = 'p-value')) %>% 
  mutate(Feature = factor(Feature, levels = Feature))

p_forest <- COX_df %>% 
  mutate(Feature = factor(Feature, levels = levels(res_plot$Feature))) %>% 
  ggplot(aes(y = Feature, color = significant)) + 
  geom_point(aes(x = HR), shape = 15, size = 3) + 
  geom_linerange(aes(xmin = HR_low, xmax = HR_high)) + 
  geom_errorbar(aes(xmin = HR_low, xmax = HR_high), width = 0.2) + 
  coord_cartesian(ylim = c(1, nrow(COX_df)+1), xlim = c(-2, 12)) + 
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'grey') + 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(), 
        axis.line.x = element_line(), 
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(), 
        legend.position = 'none') + 
  scale_color_manual(values = c('Yes' = 'red', 'No' = 'black'))
p_left <- res_plot %>% 
  ggplot(aes(y = Feature)) + 
  geom_text(aes(x = 0, label = Feature), hjust = 0, fontface = 'bold') +
  geom_text(aes(x = 1.8, label = estimate_lab), hjust = 0, fontface = ifelse(res_plot$estimate_lab == 'HR', 'bold', 'plain')) +
  theme_void() + 
  coord_cartesian(xlim = c(0,2.9))

# right side of plot - pvalues
p_right <- res_plot %>% 
  ggplot() +
  geom_text(aes(x = 0, y = Feature, label = p.value), 
            hjust=0, 
            fontface = ifelse(res_plot$p.value == 'p-value', 'bold', 'plain')) +
  theme_void() 

layout <- c(area(t = 0, l = 0, b = 30, r = 4),      
            area(t = 0, l = 4.5, b = 30, r = 6.5), 
            area(t = 0, l = 6.5, b = 30, r = 8))

p_MV <- p_left + p_forest + p_right + plot_layout(design = layout) +
  plot_annotation(title = 'Multivariable Analysis') & 
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, face = "bold"), plot.margin = unit(c(0.1, -0.5, 0.1, -0.1), "cm"))

pdf("MV_TTT_Chemo.pdf", width = 5, height = 1.5)
print(p_MV, newpage = FALSE)
dev.off()

rm(p_MV, p_left, p_forest, p_right, layout, res_plot, COX_df, coefficients, conf_int)
