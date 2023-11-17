
## Set working directory
setwd("C:/Users/ERDP5G/Desktop/Manuscripts/SNFProtocol/")

## Read in data set
in_ibd = readRDS("IBD_Example3.rds")

## Interrogate data set
str(in_ibd)
dim(in_ibd)

## Distribution of variable of interest: 
  ## IBD diagnosis, 1 = ulcerative colitis, 2 = Crohn's disease
table(in_ibd$other_dx)


    ### 
    ###   SPLIT DATA SET: 14 TEST SAMPLES 
    ###   

## Use splitstackshape package to split stratified by diagnosis
library("splitstackshape")

## Set seed so split is replicable
set.seed(4321)

## Split data with 7 samples from each group (7 UC, 7 CD)
ibd_test = stratified(in_ibd,c("other_dx"), 7)

## Look at test data
str(ibd_test)

## Isolate training data (exclude test data)
ibd_train = in_ibd[!(in_ibd$MRN %in% ibd_test$MRN),]

## identify features with only 1 class
apply(ibd_train, 2, table)
  
    ####
    ### DEFINE FEATURE SPLITS BY SOURCE
    ####
hist.features <- c("Granuloma","focal_chronic_duodenitis", "focal_active_colitis", "FEG", "ileitis_mild_cecum", "pattern_involvement_worse.distally")
endosc.features <- c("classic_backwash", "ileal_inflammation", "reverse_gradient", "small._ulcers_SB", "X5_small_ulcers_colon", "less_5_ulcer_colon", "skip_lesion", "relative_patchiness")
clin.hist.features <- c("inflammed_tag", "non_bloody_diarrhea")
otherhist.features <-c("basal_plasma_cells", "activity", "gastritis", "duodenitis", "crypt_distortion", "chronic_inflammation")

    ###
    ### CREATE MANHATTAN OF INPUT DATA 
    ###

hist_assoc = sapply(hist.features, function(x){
  fisher.test(ibd_train[,x],ibd_train$other_dx)$p.val
})

endosc_assoc = sapply(endosc.features, function(x){
  if(length(unique(ibd_train[,x])) > 1){
    fisher.test(ibd_train[,x],ibd_train$other_dx)$p.val
  } else{
    NA
  }
})

clin.hist_assoc = sapply(clin.hist.features, function(x){
  if(length(unique(ibd_train[,x])) > 1){
    fisher.test(ibd_train[,x],ibd_train$other_dx)$p.val
  } else{
    NA
  }
})

otherhist_assoc = sapply(otherhist.features, function(x){
  if(length(unique(ibd_train[,x])) > 1){
    fisher.test(ibd_train[,x],ibd_train$other_dx)$p.val
  } else{
    NA
  }
})

p_vals = c(hist_assoc,endosc_assoc,clin.hist_assoc,otherhist_assoc)

neglog10p_val = sapply(p_vals, function(x){
  -log10(x)
})

p_val_weights = sapply(neglog10p_val, function(x){
  x/sum(na.omit(neglog10p_val))
})

var_group = c(rep("Histology", length(hist_assoc)),
              rep("Endoscopy", length(endosc_assoc)),
              rep("Clinical Observation", length(clin.hist_assoc)),
              rep("OtherHistology", length(otherhist_assoc)))

p_out = data.frame(names(p_vals), p_vals, var_group)
names(p_out) = c("Feature","p.value","FeatureGroup")
p_out$neglog10p = -log10(p_out$p.value)
p_out$Feature = factor(p_out$Feature, levels = p_out$Feature)

## Create nicer column names
p_out$Features = c("Granuloma", "Focal Chronic Duodenitis",
                   "Focal Active Colitis","FEG", "Ileitis Mild Cecum",
                   "Pattern Involvement Worst Distally", 
                   "Classic Backwash","Ileal Inflammation", 
                   "Reverse Gradient", "Small Ulcers SB", "Small Colonic Ulcers",
                   "<5 Colon Ulcers", "Skip Lesion",
                   "Relative Patchiness", "Inflammed Tag",
                   "Non-Bloody Diarrhea","Basal Plasma Cells", 
                   "Activity","Gastritis", "Duodenitis",
                   "Crypt Distortion","Chronic Inflammation")

p_out$Features = factor(p_out$Features, levels = p_out$Features)

library(ggplot2)

theme_set(
  theme_bw(base_size = 15)
)
ggplot(p_out, aes(x = Features, y = neglog10p, col = FeatureGroup)) + 
  geom_point(size = 3) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ylab("-log(p-value)") + geom_hline(yintercept = -log10(0.05), col = "red", lty = 2, lwd = 1)


### split data and update code below

library(metasnf)

ibd_dl_results = generate_data_list(
  list(
    data = ibd_train[,c("MRN",hist.features)],
    name = "HistologicalFeatures",
    domain = "Histology",
    type = "categorical"
  ),
  list(
    data = ibd_train[,c("MRN",endosc.features)],
    name = "EndoscopicFeatures",
    domain = "Endoscopy",
    type = "categorical"
  ),
  list(
    data = ibd_train[,c("MRN",clin.hist.features)],
    name = "ClinicalFeatures",
    domain = "ClinObservation",
    type = "categorical"
  ),
  list(
    data = ibd_train[,c("MRN",otherhist.features)],
    name = "OtherHistologicalFeatures",
    domain = "Histology",
    type = "categorical"
  ),
  uid = "MRN"
)

ibd_data_list <- ibd_dl_results

summarize_dl(data_list = ibd_data_list)

ibd_settings_matrix <- generate_settings_matrix(
  ibd_dl_results,
  nrow = 20,
  min_k = 10,
  max_k = 30,
  seed = 42
)

weights_matrix = matrix(nrow = 20, ncol = length(p_val_weights))

for(i in 1:20){
  weights_matrix[i,] = p_val_weights
}

names(weights_matrix) = names(p_val_weights)

ibd_solutions_matrix = batch_snf(ibd_dl_results,
                             ibd_settings_matrix, 
                             weights_matrix = weights_matrix)


ibd_cluster_solutions = get_cluster_solutions(ibd_solutions_matrix)


ibd_solutions_matrix_aris <- calc_om_aris(ibd_solutions_matrix)

adjusted_rand_index_heatmap(ibd_solutions_matrix_aris)

ibd_meta_cluster_order <- get_heatmap_order(ibd_solutions_matrix_aris)

settings_matrix_heatmap(ibd_settings_matrix, order = ibd_meta_cluster_order)

ibd_batch_snf_results <- batch_snf(
  ibd_data_list,
  ibd_settings_matrix,
  return_similarity_matrices = TRUE
)

ibd_solutions_matrix <- ibd_batch_snf_results$"solutions_matrix"
ibd_similarity_matrices <- ibd_batch_snf_results$"similarity_matrices"

ibd_silhouette_scores <- calculate_silhouettes(ibd_solutions_matrix, 
                                               ibd_similarity_matrices)

## Note: Package "clv" must be installed to use this function.
ibd_dunn_indices <- calculate_dunn_indices(ibd_solutions_matrix, 
                                           ibd_similarity_matrices)

ibd_db_indices <- calculate_db_indices(ibd_solutions_matrix, 
                                       ibd_similarity_matrices)

ibd_data_list_subsamples <- subsample_data_list(
  ibd_data_list,
  n_subsamples = 30, # calculate 30 subsamples
  subsample_fraction = 0.8 # for each subsample, use random 80% of patients
)

### how does this interact with weights/snf solution 
  ### -- maybe don't use this in this example? 
ibd_pairwise_aris <- subsample_pairwise_aris(
  ibd_data_list_subsamples,
  ibd_settings_matrix
)

str(ibd_pairwise_aris)

ibd_fraction_together <- fraction_clustered_together(
  ibd_data_list_subsamples,
  ibd_settings_matrix,
  ibd_solutions_matrix
)

str(ibd_fraction_together)

### UPDATE BELOW

names(ibd_train)

ibd_outcomes = c("other_dx", "other_ordinal_dx", "Porto_dx", "Porto_ordinal_dx")
ibd_confounders = c("age", "provider")

### split these up by data type
ibd_target_list <- generate_target_list(
  list(ibd_train[,c("MRN", ibd_outcomes)], "IBDOutcomes", "categorical"),
  list(ibd_train[,c("MRN", ibd_confounders)], "IBDConfounders", "categorical"),
  uid = "MRN"
)

summarize_target_list(ibd_target_list)

ibd_extended_solutions_matrix <- extend_solutions(ibd_solutions_matrix, 
                                                  ibd_target_list, cat_test = "fisher_exact")

ibd_target_pvals <- p_val_select(ibd_extended_solutions_matrix)

head(ibd_target_pvals)

pvals_heatmap(ibd_target_pvals, order = ibd_meta_cluster_order)


###
###     SNF HEATMAP
###

similarity_matrix_heatmap(
  similarity_matrix = ibd_similarity_matrices[[17]],
  cluster_solution = ibd_cluster_solutions$"17",
  scale_diag = "mean",
  log_graph = TRUE,
  data_list = ibd_target_list,
  left_hm = list(
    "Diagnosis" = "other_dx"
  ),
  heatmap_height = grid::unit(10, "cm"),
  heatmap_width = grid::unit(10, "cm")
)


###
###     LABEL PROPAGATION
###

ibd_ids_vec = c(ibd_train$MRN, ibd_test$MRN)
ibd_traintest_vec = c(rep("train",length(ibd_train$MRN)),
                      rep("test",length(ibd_test$MRN)))
ibd_assigned_splits = data.frame(ibd_ids_vec, ibd_traintest_vec)
names(ibd_assigned_splits) = c("subjectkey","split")

ibd_full_data_list = generate_data_list(
  list(
    data = in_ibd[match(ibd_ids_vec,in_ibd$MRN),
                  c("MRN",hist.features)],
    name = "HistologicalFeatures",
    domain = "Histology",
    type = "categorical"
  ),
  list(
    data = in_ibd[match(ibd_ids_vec,in_ibd$MRN),
                  c("MRN",endosc.features)],
    name = "EndoscopicFeatures",
    domain = "Endoscopy",
    type = "categorical"
  ),
  list(
    data = in_ibd[match(ibd_ids_vec,in_ibd$MRN),
                  c("MRN",clin.hist.features)],
    name = "ClinicalFeatures",
    domain = "ClinObservation",
    type = "categorical"
  ),
  list(
    data = in_ibd[match(ibd_ids_vec,in_ibd$MRN),
                  c("MRN",otherhist.features)],
    name = "OtherHistologicalFeatures",
    domain = "Histology",
    type = "categorical"
  ),
  uid = "MRN",
  assigned_splits = ibd_assigned_splits
)

ibd_top_row = ibd_extended_solutions_matrix[17, ]

# We start with the order of training subjects in the solutions matrix
train_subjects <- colnames(subs(ibd_top_row))[-1] # the -1 removes the "row_id" column

# This is a vector containing unordered training and testing subjects
train_and_test_subjects <- dl_subjects(ibd_full_data_list, remove_prefix = FALSE)

# This is a vector containing just test subjects
test_subjects <- train_and_test_subjects[!train_and_test_subjects %in% train_subjects]

# The full subject order will be the training order (matching the order in the
#  solutions matrix) + the test subjects in any order
valid_subject_order <- c(train_subjects, test_subjects)

# Reorder the data list so that patients are ordered by valid_subject_order
full_data_list_reordered <- lapply(ibd_full_data_list,
                                   function(x) {
                                     x$"data" <- x$"data"[match(valid_subject_order, x$"data"$"subjectkey"), ]
                                     return(x)
                                   }
)

# Now the label propagation works
ibd_propagated_labels <- lp_row(ibd_top_row, full_data_list_reordered)

# Extract original IDs
ibd_propagated_labels$MRN = sub(pattern = "subject_",replacement = "",x = ibd_propagated_labels$subjectkey)
str(ibd_propagated_labels)
names(ibd_propagated_labels)[3] = "SNF_group"

## Check structure of outcomes and confounders
str(in_ibd[,c("MRN", ibd_outcomes)])
str(in_ibd[,c("MRN", ibd_confounders)])

## Merging cluster label, outcomes, and confounders to test and plot these individually
merged_ibd_outcomes = merge(ibd_propagated_labels, in_ibd[,c("MRN", ibd_outcomes)],by = "MRN")
merged_ibd_outcomes$Cluster = factor(merged_ibd_outcomes$SNF_group, levels = c(1,2), labels = c("Group 1", "Group 2"))
merged_ibd_confounders = merge(ibd_propagated_labels, in_ibd[,c("MRN", ibd_confounders)], by = "MRN")

barplot(table(merged_ibd_outcomes$other_dx[merged_ibd_outcomes$group == "train"],
              merged_ibd_outcomes$Cluster[merged_ibd_outcomes$group == "train"]),
        beside = TRUE, col = c("firebrick","deeppink2"))

barplot(table(merged_ibd_outcomes$other_dx[merged_ibd_outcomes$group == "test"],
              merged_ibd_outcomes$Cluster[merged_ibd_outcomes$group == "test"]),
        beside = TRUE, col = c("firebrick","deeppink2"))


####
###   MANHATTAN PLOT OF EACH INTEGRATED FEATURE VS. THE SELECTED CLUSTER SOLUTION 
####

## Test selected 
ibd_features_extended_solutions_matrix <- extend_solutions(ibd_extended_solutions_matrix[17, ], 
                                                           ibd_data_list, cat_test = "fisher_exact")

ibd_feature_snf_pvals <- p_val_select(ibd_features_extended_solutions_matrix)

manhattan_plot(ibd_feature_snf_pvals, threshold = 0.05, bonferroni_line = FALSE)

manhattan_plot(target_pvals, threshold = 0.05)


####
###   ADDITIONAL ANALYSIS COMPARING FEATURE SIGNFICANCE BETWEEN TRAIN AND TEST
####

merged_ibd_all = merge(ibd_propagated_labels, in_ibd,by = "MRN")
merged_ibd_all$Cluster = factor(merged_ibd_outcomes$SNF_group, levels = c(1,2), labels = c("Group 1", "Group 2"))
int_features = c(hist.features,endosc.features,clin.hist.features,otherhist.features)

train_pvals = rep(NA,length(int_features))

for(i in 1:length(train_pvals)){
  train_pvals[i] = chisq.test(merged_ibd_all$SNF_group[merged_ibd_all$group == "train"],
                              merged_ibd_all[merged_ibd_all$group == "train",int_features[i]])$p.value
}

train_pvals = unlist(train_pvals)
names(train_pvals) = int_features

test_pvals = rep(NA,length(int_features))

for(i in 1:length(train_pvals)){
  test_pvals[i] = chisq.test(merged_ibd_all$SNF_group[merged_ibd_all$group == "test"],
                             merged_ibd_all[merged_ibd_all$group == "test",int_features[i]])$p.value
}

test_pvals = unlist(test_pvals)
names(test_pvals) = int_features


tr_test_pval_df = data.frame(Features = names(test_pvals),
                             Train_p = train_pvals,
                             Test_p = test_pvals)
tr_test_pval_df = tr_test_pval_df[!is.na(tr_test_pval_df$Features),]

ggplot(tr_test_pval_df, aes(x = Train_p, y = Test_p)) + 
  geom_point(size=5) + ylim(0,1)


## which test feature was most significant? 

tr_test_pval_df[as.numeric(tr_test_pval_df$Test_p) < 0.99 & !is.na(tr_test_pval_df$Test_p),]

