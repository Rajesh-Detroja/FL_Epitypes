# R Packages
# ==========
library(dplyr)
library(ggplot2)
library(reshape)
library(ggpubr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(stringr)
library(limma)
library(tidyr)
library(weights)
library(glmnet)
library(caret)
library(doParallel)
library(MLmetrics)
# library(randomForest)

# Change directory
setwd("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/03_Feature_Selection/84110_CpGs/")

# get the 450k annotation data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k <- as.data.frame(ann450k)
head(ann450k)

# Define CpGs of interest
CpG_islands <- ann450k %>% filter(Relation_to_Island == "Island") %>% data.frame() %>% .$Name
OpenSea <- ann450k %>% filter(Relation_to_Island == "OpenSea") %>% data.frame() %>% .$Name
S_Shelf <- ann450k %>% filter(Relation_to_Island == "S_Shelf") %>% data.frame() %>% .$Name
S_Shore <- ann450k %>% filter(Relation_to_Island == "S_Shore") %>% data.frame() %>% .$Name
N_Shelf <- ann450k %>% filter(Relation_to_Island == "N_Shelf") %>% data.frame() %>% .$Name
N_Shore <- ann450k %>% filter(Relation_to_Island == "N_Shore") %>% data.frame() %>% .$Name

# Read in methylation data
bval <- readRDS("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/01_InfiniumPurify/450k/modified_InfiniumPurify/Enmix_quantile_with_QCinfo_Harman_All/beta_combat_purified.rds")
mval <- readRDS("/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/01_InfiniumPurify/450k/modified_InfiniumPurify/Enmix_quantile_with_QCinfo_Harman_All/mval_combat_purified.rds")

dim(bval)
dim(mval)

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
  dplyr::select(LY_FL_ID, ANN_ARBOR_STAGE, SEX, T_14_18, PRIM_TX_CAT, CODE_PFS, CODE_TRANSF, CODE_OS, PFS, TTT, OS)

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

# Read in methylation clustering information
epitypes <- read.table(file = "/Users/rajesh/Dropbox/Documents/CAN_IMG/01_PMCRC/Projects/00_Projects/01_FLOMICS/00_Methylation_Analysis/02_Clustering/450k/modified_InfiniumPurify/Enmix_quantile_with_QCinfo_Harman_All/FL_epitypes.tsv",
                       sep = "\t", header = TRUE, na.strings = c("","NA"))

epitypes <- epitypes %>% filter(SAMPLE_ID %in% colnames(bval))

Sample_Annotations <- merge(x = Sample_Annotations,
                            y = epitypes,
                            by.x = "SAMPLE_ID",
                            by.y = "SAMPLE_ID",
                            all.x = TRUE)

dim(Sample_Annotations)
table(Sample_Annotations$TYPE)
table(Sample_Annotations$epitypes)

row.names(Sample_Annotations) <- Sample_Annotations$SAMPLE_ID

# subset of merged methylation data
bval <- bval[, Sample_Annotations$SAMPLE_ID]
mval <- mval[, Sample_Annotations$SAMPLE_ID]
dim(bval)
dim(mval)

gc()

# Create list of groups
naiBC_Annotations <- Sample_Annotations %>% filter(TYPE == "naiBC") %>% data.frame()
gcBC_Annotations <- Sample_Annotations %>% filter(TYPE == "gcBC") %>% data.frame()
memBC_Annotations <- Sample_Annotations %>% filter(TYPE == "memBC") %>% data.frame()

# Subset methylation data group
naiBC_bval <- bval[, naiBC_Annotations$SAMPLE_ID]
gcBC_bval <- bval[, gcBC_Annotations$SAMPLE_ID]
memBC_bval <- bval[, memBC_Annotations$SAMPLE_ID]

naiBC_mval <- mval[, naiBC_Annotations$SAMPLE_ID]
gcBC_mval <- mval[, gcBC_Annotations$SAMPLE_ID]
memBC_mval <- mval[, memBC_Annotations$SAMPLE_ID]

dim(naiBC_bval)
dim(naiBC_mval)

dim(gcBC_bval)
dim(gcBC_mval)

dim(memBC_bval)
dim(memBC_mval)

# Calculate Mean Methylation for Each CpG in Naive, GC, and Memory B Cells
naiBC_means <- rowMeans(naiBC_bval)
gcBC_means <- rowMeans(gcBC_bval)
memBC_means <- rowMeans(memBC_bval)

length(naiBC_means)
length(gcBC_means)
length(memBC_means)

# Identify CpG Probes with Similar Methylation Profiles
threshold <- 0.1

# Identify CpGs with similar methylation profiles to naive B cells
gcBC_naiBC <- abs(gcBC_means - naiBC_means) <= threshold
table(gcBC_naiBC)

memBC_naiBC <- abs(memBC_means - naiBC_means) <= threshold
table(memBC_naiBC)

# Combine the results to find CpGs similar to naive in both GC and memory B cells
naiBC_gcBC_memBC <- gcBC_naiBC & memBC_naiBC
table(naiBC_gcBC_memBC)

# CpG probes to keep (those that are NOT similar across all B cell types)
keep <- !naiBC_gcBC_memBC
remove <- naiBC_gcBC_memBC
table(keep)
table(remove)

# Subset Methylation Data
bval <- bval[keep, ]
mval <- mval[keep, ]

dim(bval)
dim(mval)

gc()

# exclude blueprint samples
Sample_Annotations <- Sample_Annotations %>% filter(TYPE == "FL" | TYPE == "FL_A" | TYPE == "FL_B" | TYPE == "DLBCL" | TYPE == "DLBCL_EGA" | TYPE == "DLBCL_EN" | TYPE == "Benign_LN" | TYPE == "LN_NN_A" | TYPE == "LN_NN_B" | TYPE == "RLN")

dim(Sample_Annotations)
table(Sample_Annotations$TYPE)
table(Sample_Annotations$epitypes)

# Calculate Samples by Type per Cluster
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'FL_A'] <- 'FL'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'FL_B'] <- 'FL'

Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'DLBCL_EGA'] <- 'DLBCL'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'DLBCL_EN'] <- 'DLBCL'

Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'Benign_LN'] <- 'Normal'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'LN_NN_A'] <- 'Normal'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'LN_NN_B'] <- 'Normal'
Sample_Annotations$TYPE[Sample_Annotations$TYPE == 'RLN'] <- 'Normal'

Sample_Annotations %>% dplyr::count(epitypes, TYPE)

# Rename epitypes details
Sample_Annotations <- Sample_Annotations %>%
  filter((epitypes == "C1" & TYPE == "FL") |
           (epitypes == "C2" & TYPE == "FL") |
           (epitypes == "C3" & TYPE == "Normal")) %>%
  mutate(epitypes = case_when(
    epitypes == "C1" & TYPE == "FL"     ~ "iFL",
    epitypes == "C2" & TYPE == "FL"     ~ "aFL",
    epitypes == "C3" & TYPE == "Normal" ~ "Normal"
  ))

Sample_Annotations %>% dplyr::count(epitypes, TYPE)

# subset of merged methylation data
bval <- bval[, Sample_Annotations$SAMPLE_ID]
mval <- mval[, Sample_Annotations$SAMPLE_ID]
dim(bval)
dim(mval)

# Write TYPE and Batch Information
df1 <- Sample_Annotations[c("SAMPLE_ID", "TYPE","Batch")]
write.table(df1, file="TYPE_INFO_ALL.txt", sep="\t", quote = FALSE, row.names = FALSE)

gc()

# subset matrix
set.seed(1234)
probe_num = nrow(bval)
#probe_num = "10000"

# # subset random probes
# bval_sub <- bval[sample(nrow(bval), probe_num), ]
# mval_sub <- mval[intersect(rownames(bval_sub), rownames(mval)), ]

# All probes
bval_sub <- bval
mval_sub <- mval

dim(bval_sub)
bval_sub[1:5, 1:5]

dim(mval_sub)
mval_sub[1:5, 1:5]

gc()

# Transform data sets
bval_sub <- t(bval_sub)
dim(bval_sub)
bval_sub[1:5, 1:5]

mval_sub <- t(mval_sub)
dim(mval_sub)
mval_sub[1:5, 1:5]

gc()

# Preparing input data
bval_sub <- data.matrix(bval_sub)
class(bval_sub)
dim(bval_sub)

mval_sub <- data.matrix(mval_sub)
class(mval_sub)
dim(mval_sub)

# Ensure class labels are factors
class_labels <- Sample_Annotations$epitypes
class_labels <- as.factor(class_labels)
class_labels

# Normalize data
methylation_data_normalized <- scale(mval_sub)
dim(methylation_data_normalized)

# Split data into training and testing sets
set.seed(123)  # For reproducibility
trainIndex <- createDataPartition(class_labels, p = .8, list = FALSE, times = 1)
trainData <- methylation_data_normalized[trainIndex, ]
testData  <- methylation_data_normalized[-trainIndex, ]
trainLabels <- class_labels[trainIndex]
testLabels  <- class_labels[-trainIndex]

# Check for missing values
if (any(is.na(trainData))) {
  stop("Training data contains missing values.")
}

# Set up parallel processing
cl <- makePSOCKcluster(6)
registerDoParallel(cl)

# Train control with cross-validation
train_control <- trainControl(method = "cv", 
                              # method = "repeatedcv", # Alternatively, explore repeated cross-validation
                              number = 10,
                              allowParallel = TRUE,
                              classProbs = TRUE,
                              # summaryFunction = twoClassSummary,
                              summaryFunction = multiClassSummary,
                              # sampling = "up", # Handle class imbalance by oversampling the minority class
                              verboseIter = TRUE)  # Verbose output for more info

# Train the elastic net model
set.seed(123)
elastic_net_model <- train(x = trainData, y = trainLabels,
                           method = "glmnet",
                           # method = "rf", # Alternatively, explore Random Forest
                           trControl = train_control,
                           # metric = "ROC",  # Optimize for ROC since it's binary classification (twoClassSummary)
                           tuneLength = 10,
                           #tuneLength = 20 # If required, Increase from 10 to 20
                           )

# Stop parallel processing
stopCluster(cl)
registerDoSEQ()

# Extract the optimal lambda value from the trained model
best_lambda <- elastic_net_model$bestTune$lambda

# Extract the coefficients of the best model for each class
coefficients <- coef(elastic_net_model$finalModel, s = best_lambda)

# For multiClassSummary Only
# ==========================
# Combine selected features for all classes
combined_features <- list()

for (i in 1:length(levels(trainLabels))) {
  class <- levels(trainLabels)[i]
  selected_features <- as.data.frame(as.matrix(coefficients[[i]]))
  colnames(selected_features) <- class
  selected_features <- selected_features[selected_features[, 1] != 0, , drop = FALSE]
  selected_features$feature <- rownames(selected_features)
  combined_features[[class]] <- selected_features
}

# Ensure the row names are preserved when merging
combined_features_df <- Reduce(function(x, y) merge(x, y, by = "feature", all = TRUE), combined_features)
rownames(combined_features_df) <- combined_features_df$feature
combined_features_df <- combined_features_df[, -1]
combined_features_df[is.na(combined_features_df)] <- 0  # Replace NA with 0 for clarity

# Display combined features
combined_features_df1 <- as.data.frame(combined_features_df)
head(combined_features_df1)
dim(combined_features_df1)
write.table(combined_features_df1, file="combined_features.txt", sep="\t", quote = FALSE, row.names = TRUE)

# Writing out M and beta matrices
beta_value <- bval[intersect(rownames(combined_features_df), rownames(bval)) ,]
mval_value <- mval[intersect(rownames(combined_features_df), rownames(mval)) ,]

dim(beta_value)
dim(mval_value)

beta_value[1:5, 1:5]
mval_value[1:5, 1:5]

saveRDS(beta_value, "beta_combat_purified_features.rds")
saveRDS(mval_value, "mval_combat_purified_features.rds")

# Predict on test data
predictions <- predict(elastic_net_model, newdata = testData)

# Confusion matrix to evaluate performance
confusion_matrix <- confusionMatrix(predictions, testLabels)
print(confusion_matrix)

# # Only RF Mode
# # ============
# # Get feature importance
# importance_rf <- varImp(elastic_net_model, scale = FALSE)
# 
# # Display importance
# print(importance_rf)
# 
# pdf("importance_rf.pdf")
# # For caret model
# ggplot(importance_rf, aes(x = reorder(Features, Overall), y = Overall)) +
#   geom_bar(stat = "identity") +
#   coord_flip() +
#   theme_minimal() +
#   xlab("Features") +
#   ylab("Importance") +
#   ggtitle("Feature Importance from Random Forest")
# dev.off()
# 
# # Convert importance to a data frame and add rownames as a column
# importance_rf_df <- as.data.frame(importance_rf$importance)
# importance_rf_df$feature <- rownames(importance_rf_df)
# 
# # Set threshold value
# threshold_value <- 0.1
# 
# # Filter features based on threshold
# important_features <- importance_rf_df[importance_rf_df$Overall > threshold_value, ]
# 
# # View the filtered features with row names (features) included
# print(important_features)
# dim(important_features)
# 
# # Write selected features to a file
# write.table(important_features, file = "important_features.txt", sep = "\t", quote = FALSE, row.names = TRUE)
# 
# # Writing out M and beta matrices
# beta_value <- bval[intersect(rownames(important_features), rownames(bval)) ,]
# mval_value <- mval[intersect(rownames(important_features), rownames(mval)) ,]
# 
# dim(beta_value)
# dim(mval_value)
# 
# beta_value[1:5, 1:5]
# mval_value[1:5, 1:5]
# 
# saveRDS(beta_value, "beta_combat_purified_features.rds")
# saveRDS(mval_value, "mval_combat_purified_features.rds")


# # For twoClassSummary Only
# # ========================
# # Convert sparse matrix to a standard data frame
# selected_features <- as.data.frame(as.matrix(coefficients))
# colnames(selected_features) <- "Coefficient"  # Assign a meaningful column name
# 
# # Filter for non-zero coefficients (selected features)
# selected_features <- selected_features[selected_features$Coefficient != 0, , drop = FALSE]
# 
# # Add feature names as a column
# selected_features$feature <- rownames(selected_features)
# 
# # Rearrange columns to ensure clarity
# selected_features <- selected_features[, c("feature", "Coefficient")]
# 
# # Ensure row names are preserved
# rownames(selected_features) <- selected_features$feature
# 
# # Replace NA values with 0 if present (optional safeguard)
# selected_features[is.na(selected_features)] <- 0
# 
# # Display selected features
# head(selected_features)
# dim(selected_features)
# 
# # Write selected features to a file
# write.table(selected_features, file = "selected_features.txt", sep = "\t", quote = FALSE, row.names = TRUE)
# 
# # Writing out M and beta matrices
# beta_value <- bval[intersect(rownames(selected_features), rownames(bval)) ,]
# mval_value <- mval[intersect(rownames(selected_features), rownames(mval)) ,]
# 
# dim(beta_value)
# dim(mval_value)
# 
# beta_value[1:5, 1:5]
# mval_value[1:5, 1:5]
# 
# saveRDS(beta_value, "beta_combat_purified_features.rds")
# saveRDS(mval_value, "mval_combat_purified_features.rds")
# 
# # Predict on test data
# predictions <- predict(elastic_net_model, newdata = testData)
# 
# # Confusion matrix to evaluate performance
# confusion_matrix <- confusionMatrix(predictions, testLabels)
# print(confusion_matrix)
