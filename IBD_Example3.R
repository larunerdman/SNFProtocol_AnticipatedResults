
## Set working directory
setwd("C:/Users/ERDP5G/Desktop/Manuscripts/SNFProtocol/")

### install metasnf if needed
# devtools::install_github("BRANCHlab/metasnf", force = TRUE)
library(metasnf)
library(SNFtool)
library(ggplot2)

## Read in data set
in_ibd = readRDS("IBD_Example3.rds")

## Interrogate data set
str(in_ibd)
dim(in_ibd)

## Distribution of variable of interest: 
  ## IBD diagnosis, 1 = ulcerative colitis, 2 = Crohn's disease
table(in_ibd$IBD_dx)

    ### 
    ###   SPLIT DATA SET: 14 TEST SAMPLES 
    ###   

## Use splitstackshape package to split stratified by diagnosis
library("splitstackshape")

## Set seed so split is replicable
set.seed(4321)

## Split data with 7 samples from each group (7 UC, 7 CD)
ibd_test = stratified(in_ibd,c("IBD_dx"), 7)

## Look at test data
str(ibd_test)

## Isolate training data (exclude test data)
ibd_train = in_ibd[!(in_ibd$ID %in% ibd_test$ID),]

## identify features with only 1 class
apply(ibd_train, 2, table)
  
    ####
    ### DEFINE FEATURE SPLITS BY SOURCE
    ####

hist.features <- c("Granuloma","focal_chronic_duodenitis", "focal_active_colitis", "FEG", 
                   "ileitis_mild_cecum", "pattern_involvement_worse.distally",
                   "basal_plasma_cells", "activity", "gastritis", "duodenitis", 
                   "crypt_distortion", "chronic_inflammation")
endosc.features <- c("classic_backwash", "ileal_inflammation", "reverse_gradient", "small._ulcers_SB",
                     "X5_small_ulcers_colon", "less_5_ulcer_colon", "skip_lesion", "relative_patchiness")
clin.obs.features <- c("inflammed_tag", "non_bloody_diarrhea")

    ###
    ### CREATE MANHATTAN OF INPUT DATA 
    ###

hist_assoc = sapply(hist.features, function(x){
  fisher.test(ibd_train[,x],ibd_train$IBD_dx)$p.val
})

endosc_assoc = sapply(endosc.features, function(x){
  if(length(unique(ibd_train[,x])) > 1){
    fisher.test(ibd_train[,x],ibd_train$IBD_dx)$p.val
  } else{
    NA
  }
})

clin.obs_assoc = sapply(clin.obs.features, function(x){
  if(length(unique(ibd_train[,x])) > 1){
    fisher.test(ibd_train[,x],ibd_train$IBD_dx)$p.val
  } else{
    NA
  }
})

p_vals = c(hist_assoc,endosc_assoc,clin.obs_assoc)

neglog10p_val = sapply(p_vals, function(x){
  -log10(x)
})

var_group = c(rep("Histology", length(hist_assoc)),
              rep("Endoscopy", length(endosc_assoc)),
              rep("Clinical Observation", length(clin.obs_assoc)))

p_out = data.frame(names(p_vals), p_vals, var_group)
names(p_out) = c("Feature","p.value","FeatureGroup")
p_out$neglog10p = -log10(p_out$p.value)
p_out$Feature = factor(p_out$Feature, levels = p_out$Feature)

## Create nicer column names
p_out$Features = c("Granuloma", "Focal Chronic Duodenitis",
                   "Focal Active Colitis","FEG", "Ileitis Mild Cecum",
                   "Pattern Involvement Worst Distally", 
                   "Basal Plasma Cells", "Activity",
                   "Gastritis", "Duodenitis",
                   "Crypt Distortion","Chronic Inflammation",
                   "Classic Backwash","Ileal Inflammation",
                   "Reverse Gradient", "Small Ulcers SB", "Small Colonic Ulcers",
                   "<5 Colon Ulcers", "Skip Lesion",
                   "Relative Patchiness", "Inflammed Tag",
                   "Non-Bloody Diarrhea")

p_out$Features = factor(p_out$Features, levels = p_out$Features)

library(ggplot2)

theme_set(
  theme_bw(base_size = 15)
)
ggplot(p_out, aes(x = Features, y = neglog10p, col = FeatureGroup)) + 
  geom_point(size = 3) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ylab("-log(p-value)") + geom_hline(yintercept = -log10(0.05), col = "red", lty = 2, lwd = 1) + 
  guides(color = guide_legend(title = "Data Source"))


  ###
  ###   DEFINE WEIGHTS
  ###

p_val_weights = sapply(neglog10p_val, function(x){
  x/sum(na.omit(neglog10p_val))
})

  ###
  ###   metaSNF
  ###

ibd_dl_results = generate_data_list(
  list(
    data = ibd_train[,c("ID",hist.features)],
    name = "HistologicalFeatures",
    domain = "Histology",
    type = "categorical"
  ),
  list(
    data = ibd_train[,c("ID",endosc.features)],
    name = "EndoscopicFeatures",
    domain = "Endoscopy",
    type = "categorical"
  ),
  list(
    data = ibd_train[,c("ID",clin.obs.features)],
    name = "ClinicalFeatures",
    domain = "ClinObservation",
    type = "categorical"
  ),
  uid = "ID"
)

ibd_data_list <- ibd_dl_results

summarize_dl(data_list = ibd_data_list)

ibd_settings_matrix <- generate_settings_matrix(
  ibd_dl_results,
  nrow = 20,
  min_k = 15,
  max_k = 25,
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

solutions_matrix_aris <- calc_aris(ibd_solutions_matrix)

meta_cluster_order <- get_matrix_order(solutions_matrix_aris)

# adjusted_rand_index_heatmap(
#   solutions_matrix_aris,
#   order = meta_cluster_order
# )

heatmap_colours <- colorRampPalette(
  rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu"))
)(100)

adjusted_rand_index_heatmap(
  solutions_matrix_aris,
  order = meta_cluster_order,
  col = heatmap_colours,
  show_row_names = TRUE,
  show_column_names = TRUE,
  rect_gp = grid::gpar(col = "gray", lwd = 1)
)

settings_matrix_heatmap(ibd_settings_matrix, order = meta_cluster_order)


ibd_batch_snf_results <- batch_snf(
  ibd_data_list,
  ibd_settings_matrix,
  return_similarity_matrices = TRUE
)

ibd_solutions_matrix <- ibd_batch_snf_results$"solutions_matrix"
ibd_similarity_matrices <- ibd_batch_snf_results$"similarity_matrices"

names(ibd_train)

ibd_outcomes = c("IBD_dx")
ibd_confounders_cont = c("age")
ibd_confounders_cat = c("provider")


### split these up by data type
ibd_target_list <- generate_data_list(
  list(ibd_train[,c("ID", ibd_outcomes)], "IBDOutcomes", "IBDOutcomes", "categorical"),
  list(ibd_train[,c("ID", ibd_confounders_cat)], "IBDConfoundersCat", "IBDConfoundersCat", "categorical"),
  list(ibd_train[,c("ID", ibd_confounders_cont)], "IBDConfoundersCont","IBDConfoundersCont", "numeric"),
  uid = "ID"
)

summarize_dl(ibd_target_list)

ibd_extended_solutions_matrix <- extend_solutions(ibd_solutions_matrix, 
                                                  ibd_target_list, cat_test = "fisher_exact")

ibd_target_pvals <- get_pvals(ibd_extended_solutions_matrix)

head(ibd_target_pvals)

ibd_target_pvals_sub = ibd_target_pvals[,1:4]

pval_heatmap(ibd_target_pvals_sub, order = meta_cluster_order)


  ###
  ###     SNF HEATMAP
  ###

similarity_matrix_heatmap(
  similarity_matrix = ibd_similarity_matrices[[9]],
  cluster_solution = ibd_cluster_solutions$"9",
  scale_diag = "mean",
  log_graph = TRUE,
  data_list = ibd_target_list,
  left_hm = list(
    "IBDDiagnosis" = "IBD_dx"
  ),
  heatmap_height = grid::unit(10, "cm"),
  heatmap_width = grid::unit(10, "cm")
)


  ####
  ###   MANHATTAN PLOT OF EACH INTEGRATED FEATURE, CONFOUNDER, AND TARGET VS. THE SELECTED CLUSTER SOLUTION 
  ####

###  Extract clustering from top solution
select_cluster_sol = ibd_cluster_solutions$"9"

### Find p-value associations 
hist_clust_assoc = sapply(hist.features, function(x){
  fisher.test(ibd_train[,x],select_cluster_sol)$p.val
})

endosc_clust_assoc = sapply(endosc.features, function(x){
  if(length(unique(ibd_train[,x])) > 1){
    fisher.test(ibd_train[,x],select_cluster_sol)$p.val
  } else{
    NA
  }
})

clin.obs_clust_assoc = sapply(clin.obs.features, function(x){
  if(length(unique(ibd_train[,x])) > 1){
    fisher.test(ibd_train[,x],select_cluster_sol)$p.val
  } else{
    NA
  }
})


age_clust_assoc = kruskal.test(ibd_train$age, select_cluster_sol)$p.val

provider_clust_assoc = fisher.test(ibd_train$provider, select_cluster_sol)$p.val

ibd_clust_assoc = fisher.test(ibd_train$IBD_dx, select_cluster_sol)$p.val

p_vals_clust = c(hist_clust_assoc,endosc_clust_assoc,clin.obs_clust_assoc, age_clust_assoc,
                 provider_clust_assoc, ibd_clust_assoc)
names(p_vals_clust)[23:25] = c("Age_at_dx","Provider","IBD_dx")

neglog10p_clust_val = sapply(p_vals_clust, function(x){
  -log10(x)
})

var_group = c(rep("Histology", length(hist_clust_assoc)),
              rep("Endoscopy", length(endosc_clust_assoc)),
              rep("Clinical Observation", length(clin.obs_clust_assoc)),
              rep("Confounder", 2),
              "Target")

p_clust_out = data.frame(names(p_vals_clust), p_vals_clust, var_group)
names(p_clust_out) = c("Feature","p.value","FeatureGroup")
p_clust_out$neglog10p = -log10(p_clust_out$p.value)
p_clust_out$Feature = factor(p_clust_out$Feature, levels = p_clust_out$Feature)

## Create nicer column names
p_clust_out$Features = c("Granuloma", "Focal Chronic Duodenitis",
                   "Focal Active Colitis","FEG", "Ileitis Mild Cecum",
                   "Pattern Involvement Worst Distally", 
                   "Basal Plasma Cells", "Activity",
                   "Gastritis", "Duodenitis",
                   "Crypt Distortion","Chronic Inflammation",
                   "Classic Backwash","Ileal Inflammation",
                   "Reverse Gradient", "Small Ulcers SB", "Small Colonic Ulcers",
                   "<5 Colon Ulcers", "Skip Lesion",
                   "Relative Patchiness", "Inflammed Tag",
                   "Non-Bloody Diarrhea", "Age_at_dx","Provider","IBD_dx")

p_clust_out$Features = factor(p_clust_out$Features, levels = p_clust_out$Features)

library(ggplot2)

theme_set(
  theme_bw(base_size = 15)
)
ggplot(p_clust_out, aes(x = Features, y = neglog10p, col = FeatureGroup)) + 
  geom_point(size = 3) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ylab("-log(p-value)") + geom_hline(yintercept = -log10(0.05), col = "red", lty = 2, lwd = 1) + 
  guides(color = guide_legend(title = "Data Source")) + 
  xlab("Features, Confounders, Target")


  ###
  ###     LABEL PROPAGATION
  ###

## Create full dataset with training and testing data
ibd_ids_vec = c(ibd_train$ID, ibd_test$ID)
ibd_traintest_vec = c(rep("train",length(ibd_train$ID)),
                      rep("test",length(ibd_test$ID)))
ibd_assigned_splits = data.frame(ibd_ids_vec, ibd_traintest_vec)
names(ibd_assigned_splits) = c("subjectkey","split")

## Create full data list object
ibd_full_data_list = generate_data_list(
  list(
    data = in_ibd[match(ibd_ids_vec,in_ibd$ID),
                  c("ID",hist.features)],
    name = "HistologicalFeatures",
    domain = "Histology",
    type = "categorical"
  ),
  list(
    data = in_ibd[match(ibd_ids_vec,in_ibd$ID),
                  c("ID",endosc.features)],
    name = "EndoscopicFeatures",
    domain = "Endoscopy",
    type = "categorical"
  ),
  list(
    data = in_ibd[match(ibd_ids_vec,in_ibd$ID),
                  c("ID",clin.obs.features)],
    name = "ClinicalFeatures",
    domain = "ClinObservation",
    type = "categorical"
  ),
  uid = "ID"
)


train_solutions_matrix <- batch_snf(
  ibd_data_list,
  ibd_settings_matrix
)

# Generated predicted clustering for all 
extended_solutions_matrix <- lp_solutions_matrix(
  train_solutions_matrix,
  ibd_full_data_list
)

ibd_propagated_labels = extended_solutions_matrix[,c("subjectkey","group","9")]
names(ibd_propagated_labels)[3] = "SNF_group"
ibd_propagated_labels$ID = sub(pattern = "subject_",replacement = "",x = ibd_propagated_labels$subjectkey)

  ####
  ###   CLUSTER ~ TARGET ASSOCIATIONS
  ####

## Merging cluster label, outcomes, and confounders to test and plot these individually
merged_ibd_outcomes = merge(ibd_propagated_labels, in_ibd[,c("ID", ibd_outcomes)],by = "ID")
merged_ibd_outcomes$Cluster = factor(merged_ibd_outcomes$SNF_group, levels = c(1,2), labels = c("Group 1", "Group 2"))
merged_ibd_confounders = merge(ibd_propagated_labels, in_ibd[,c("ID", ibd_outcomes)], by = "ID")

    ## Training sample graph
barplot(table(merged_ibd_outcomes$IBD_dx[merged_ibd_outcomes$group == "train"],
              merged_ibd_outcomes$Cluster[merged_ibd_outcomes$group == "train"]),
        beside = TRUE, col = c("firebrick","deeppink2"))

    ## Test/held-out sample graph
barplot(table(merged_ibd_outcomes$IBD_dx[merged_ibd_outcomes$group == "test"],
              merged_ibd_outcomes$Cluster[merged_ibd_outcomes$group == "test"]),
        beside = TRUE, col = c("firebrick","deeppink2"))



