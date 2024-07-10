### To use each of these packages, one must ensure they are installed: the installation step need only be performed one time
### We will reference the packages used for each function in the steps below
### However, the example section will assume all 3 packages are installed and called into the workspace using the library function

# install.packages("SNFtool")
# install.packages("devtools")
# devtools::install_github("pfruan/abSNF")
# devtools::install_github("BRANCHlab/metasnf", force = TRUE)

library(SNFtool)
library(abSNF)
library(metasnf)

####
### creating + importing data to demonstrate steps
###     These are not part of the SNF steps, rather they are a set up for the example data
####

### Data for SNF integration
data(Data1)

## Adding missing data 
set.seed(1234)
D1_miss = runif(3, min = 1, max = nrow(Data1))

Data1[D1_miss, 2] = NA

data(Data2)

## Adding missing data
set.seed(1122)
D2_miss = runif(3, min = 1, max = nrow(Data2))

Data2[D2_miss, 1] = NA

all_data = cbind(Data1,Data2)
features_vec1 = c("V1", "V2")
features_vec2 = c("V3", "V4")

data(label)
age = 25*(label) + rnorm(n = length(label), mean = 0, sd = 5)

  ###
  ###    STEP 1: (Optional: Divide data into training and test sets)
  ###

## identify 25% of data set samples to hold out
test_subject_indices = sample(1:nrow(all_data), size = 50,replace = FALSE)
## keep remaining samples to develop SNF model with
train_subject_indices = (1:nrow(all_data))[!((1:nrow(all_data)) %in% test_subject_indices)]

train_data = all_data[train_subject_indices,]
test_data = all_data[test_subject_indices,]

train_label = label[train_subject_indices]
test_label = label[test_subject_indices]

train_age = age[train_subject_indices]
test_age = age[test_subject_indices]

  ###
  ###    STEP 2: Account for missing data 
  ###

# Here we are removing the few missing samples and ensuring that 
  # our label and confounder data include the same individuals

complete_train_data = which(complete.cases(train_data))
complete_test_data = which(complete.cases(test_data))

train_data = train_data[complete_train_data,]
train_label = train_label[complete_train_data]
train_age = train_age[complete_train_data]

test_data = test_data[complete_test_data,]
test_label = test_label[complete_test_data]
test_age = test_age[complete_test_data]

  
  ###
  ###    STEP 3: Division of features by data source
  ###

  ## Standard data source division
train_data1 = train_data[,features_vec1]
train_data2 = train_data[,features_vec2]

test_data1 = test_data[,features_vec1]
test_data2 = test_data[,features_vec2]

  ## Data source division for metaSNF
train_data1_msnf = train_data1
train_data2_msnf = train_data2
train_label_msnf = data.frame(train_label)
names(train_label_msnf) = "Label"
train_age_msnf = data.frame(train_age)
names(train_age_msnf) = "Age"

train_data1_msnf$sample_ids = rownames(train_data1_msnf)
train_data2_msnf$sample_ids = rownames(train_data2_msnf)
train_label_msnf$sample_ids = rownames(train_data2_msnf)
train_age_msnf$sample_ids = rownames(train_data2_msnf)

data1_msnf = rbind(train_data1, test_data1)
data1_msnf$sample_ids = rownames(data1_msnf)
data2_msnf = rbind(train_data2, test_data2)
data2_msnf$sample_ids = rownames(data2_msnf)

identical(data1_msnf$sample_ids, data2_msnf$sample_ids)

  ###
  ###    VISUALIZATION: HEATMAP OF CORRELATION P-VALUES 
  ###

# prepare your data into a metaSNF object 
heatmap_data_list = metasnf::generate_data_list(
  list(
    data = train_data1_msnf,
    name = "Data1",
    domain = "Data1",
    type = "continuous"
  ),
  list(
    data = train_data2_msnf,
    name = "Data2",
    domain = "Data2",
    type = "continuous"
  ),
  list(
    data = train_label_msnf,
    name = "Label",
    domain = "Label",
    type = "ordinal"
  ),
  list(
    data = train_age_msnf,
    name = "Age",
    domain = "Age",
    type = "continuous"
  ),
  uid = "sample_ids"
)

summarize_dl(heatmap_data_list)


assoc_pval_matrix <- calc_assoc_pval_matrix(heatmap_data_list)

assoc_pval_heatmap(
  assoc_pval_matrix,
  confounders = list(
    "Age" = "Age"
  ),
  out_of_models = list(
    "Label" = "Label"
  ),
  row_km = 3,
  column_km = 3
)


  ###
  ###    STEP 4: (OPTIONAL) NORMALIZE DATA
  ###

standardized_data1 = SNFtool::standardNormalization(train_data1)
standardized_data2 = SNFtool::standardNormalization(train_data2)

  ###
  ###    STEP 5: (OPTIONAL) GENERATE WEIGHTS FEATURE-BY-FEATURE
  ###

### function to transform p-values to weights 
make_weights = function(p_vec){out_vec = -log10(p_vec)/sum(-log10(p_vec)); return(out_vec)}

### uniform weights (this is equivalent to no weighting)
unif_p = rep(0.5, length(c(features_vec1, features_vec2)))
unif_weights = make_weights(p_vec = unif_p)
names(unif_weights) = c(features_vec1, features_vec2)

### suppose weights associated with a categorical outcome
p_vals = unlist(apply(train_data, 2, function(x){kruskal.test(x, train_label)$p.value}))
p_weights = make_weights(p_vec = p_vals)
names(p_weights) = c(features_vec1, features_vec2)

  ###
  ### VISUALIZATION: Manhattan plot of p-values
  ###

library(reshape2)

pval_df = melt(p_vals)
pval_df$Feature = rownames(pval_df)
names(pval_df)[1] = "pvalue"
pval_df$neglog10p = -log10(pval_df$pvalue)

library(ggplot2)

theme_set(
  theme_classic(base_size = 20)
)
ggplot(pval_df, aes(x = Feature, y = neglog10p)) + 
  geom_point(size = 4, col = "blue") + ylab("-log10(p)")


  ###
  ###    STEP 6: GENERATING DISTANCE MATRICES
  ###

data_dist1 = SNFtool::dist2(X = standardized_data1, C = standardized_data1) 
data_dist2 = SNFtool::dist2(X = standardized_data2, C = standardized_data2)

  ###
  ###    STEP 7: CREATE AFFINITY MATRICES
  ###

##Set the hyper-parameters defaults for your SNF run : 
K=20 ##number of neighbors,must be greater than 1. usually(10~30) 
alpha=0.5 ##hyperparameter, usually (0.3~0.8) 
T=20 ###Number of Iterations, usually (10~50)

data_aff1 = SNFtool::affinityMatrix(data_dist1, K = K, sigma = alpha)
data_aff2 = SNFtool::affinityMatrix(data_dist2, K = K, sigma = alpha)

  ###
  ###    STEP 8: SIMILARITY NETWORK FUSION STEP
  ###

snf_out = SNFtool::SNF(W = list(data_aff1, data_aff2), K = K, t = T)

  ###
  ###    STEP 9: CLUSTER SIMILARITY NETWORK SOLUTION
  ###

estimated_snf_clusters = estimateNumberOfClustersGivenGraph(snf_out, NUMC = 2:6)
clusters_snf = spectralClustering(snf_out, estimated_snf_clusters[[1]])

### Clustering data source-specific similarity network
  ## Data set 1
estimated_data1_clusters = estimateNumberOfClustersGivenGraph(data_aff1, NUMC = 2:6)
clusters_data1 = spectralClustering(data_aff1, estimated_data1_clusters[[1]])

  ## Data set 2
estimated_data2_clusters = estimateNumberOfClustersGivenGraph(data_aff2, NUMC = 2:6)
clusters_data2 = spectralClustering(data_aff1, estimated_data1_clusters[[1]])

  ###
  ###    VISUALIZATION: Heatmaps 
  ###

library(RColorBrewer)

## Data 1
displayClustersWithHeatmap(data_aff1, group = clusters_data1,
                           col = brewer.pal(n = 9, name = "Blues"))
## Data 2
displayClustersWithHeatmap(data_aff2, group = clusters_data2,
                           col = brewer.pal(n = 9, name = "Purples"))

## SNF
displayClustersWithHeatmap(snf_out, group = clusters_snf,
                           col = brewer.pal(n = 9, name = "Greens"))


  ###
  ###    STEP 10: (OPTIONAL) RESAMPLE DATA AND ITERATE STEPS 7-10 TO ASSESS SNF STABILITY
  ###

# prepare your data using a metaSNF 
my_data_list = metasnf::generate_data_list(
  list(
    data = train_data1_msnf,
    name = "Data1",
    domain = "Data1",
    type = "continuous"
  ),
  list(
    data = train_data2_msnf,
    name = "Data2",
    domain = "Data2",
    type = "continuous"
  ),
  uid = "sample_ids"
)

# Generate settings_matrix: we only choose one set of parameters here
  # but you can test many to assess the relative stability of these choices
settings_matrix <- generate_settings_matrix(
  my_data_list,
  nrow = 1,
  max_k = 20,
  seed = 42
)

data_list_subsamples <- subsample_data_list(
  my_data_list,
  n_subsamples = 50, # calculate 30 subsamples
  subsample_fraction = 0.8 # for each subsample, use random 80% of patients
)

pairwise_aris <- subsample_pairwise_aris(
  data_list_subsamples,
  settings_matrix
)

## Adjusted Rand Index over subsampled clustering: 
pairwise_aris

# Run SNF and clustering
solutions_matrix <- batch_snf(
  my_data_list,
  settings_matrix
)

fraction_together <- fraction_clustered_together(
  data_list_subsamples,
  settings_matrix,
  solutions_matrix
)

## Fraction of pairs which cluster together 
fraction_together

## Co-cluster heatmap

# Extract co-cluster data
cocluster_data <- generate_cocluster_data(
  data_list = my_data_list,
  data_list_subsamples = data_list_subsamples,
  settings_matrix_row = settings_matrix[1, ]
)

# Plot co-clustering 
cocluster_heatmap(cocluster_data = cocluster_data)

  ###
  ###    STEP 11: (OPTIONAL) metaSNF
  ###

# Add target list
target_list <- generate_data_list(
  list(
    data = train_label_msnf,
    name = "Label",
    domain = "Label",
    type = "ordinal"
  ),
  list(
    data = train_age_msnf,
    name = "Age",
    domain = "Age",
    type = "continuous"
  ),
  uid = "sample_ids"
)

summarize_dl(target_list)
summarize_dl(my_data_list)

# Create settings matrix with many options 
settings_matrix <- generate_settings_matrix(
  my_data_list,
  nrow = 20,
  min_k = 20,
  max_k = 50,
  seed = 42
)

# Generate SNF solutions with all possible settings
solutions_matrix <- batch_snf(my_data_list, settings_matrix)

cluster_solutions <- get_cluster_solutions(solutions_matrix)

solutions_matrix_aris <- calc_aris(solutions_matrix)

meta_cluster_order <- get_matrix_order(solutions_matrix_aris)

ari_hm <- adjusted_rand_index_heatmap(
  solutions_matrix_aris,
  order = meta_cluster_order
)

# To show settings: 
settings_matrix_heatmap(settings_matrix, order = meta_cluster_order)

# If needed: 
# BiocManager::install("InteractiveComplexHeatmap")
## Interactive heatmap of Adjusted Rand Index
shiny_annotator(ari_hm)

  ###
  ###    STEP 12: (OPTIONAL) LABEL PROPAGATION
  ###

  ### USING SNFtool

train_list = list(train_data1, train_data2)
test_list = list(test_data1, test_data2)

new_label = groupPredict(train = train_list, 
                         test = test_list, 
                         groups = clusters_snf,
                         K = K,
                         alpha = alpha,
                         t = T, method = TRUE)

  ### USING metasnf

# A data list with just training subjects
train_data_list <- generate_data_list(
  list(
    data = train_data1_msnf,
    name = "Data1",
    domain = "Data1",
    type = "continuous"
  ),
  list(
    data = train_data2_msnf,
    name = "Data2",
    domain = "Data2",
    type = "continuous"
  ),
  uid = "sample_ids"
)

# A data list with training and testing subjects
full_data_list <- generate_data_list(
  list(
    data = data1_msnf,
    name = "Data1",
    domain = "Data1",
    type = "continuous"
  ),
  list(
    data = data2_msnf,
    name = "Data2",
    domain = "Data2",
    type = "continuous"
  ),
  uid = "sample_ids"
)

# Construct the target lists
train_target_list <- generate_data_list(
  list(
    data = train_label_msnf,
    name = "Label",
    domain = "Label",
    type = "ordinal"
  ),
  list(
    data = train_age_msnf,
    name = "Age",
    domain = "Age",
    type = "continuous"
  ),
  uid = "sample_ids"
)

train_solutions_matrix <- batch_snf(
  train_data_list,
  settings_matrix
)

# Generated predicted clustering for all 
extended_solutions_matrix <- extend_solutions(
  train_solutions_matrix,
  train_target_list
)

# choose clustering solution, here we choose a random one but 
  # users can select solutions in any manner they choose 
  # i.e. based on ARI, target p-values, etc
top_row <- extended_solutions_matrix[4, ]

# propagate labels
propagated_labels <- lp_solutions_matrix(top_row, full_data_list)

  ###
  ###    STEP 13: TEST INTEGRATIVE FEATURES
  ###

### Assessing target and confounder associations with clustering: 
extended_solutions_matrix <- extend_solutions( 
  solutions_matrix,
  data_list = my_data_list,
  target_list = target_list
)

target_pvals <- get_pvals(extended_solutions_matrix)

pval_heatmap(target_pvals, order = meta_cluster_order)

rep_solutions <- get_representative_solutions(
  solutions_matrix_aris,
  split_vector = c(18), ## one split means 2 representative solutions
  order = meta_cluster_order,
  extended_solutions_matrix
)

mc_manhattan_plot(
  rep_solutions,
  data_list = my_data_list,
  target_list = target_list,
  hide_x_labels = TRUE,
  point_size = 2,
  text_size = 12,
  threshold = 0.05,
  neg_log_pval_thresh = 5
)

  ###
  ###    STEP 13: TEST TEST TARGET AND CONFOUNDING FEATURES
  ###

pval_heatmap(target_pvals, order = meta_cluster_order)

rep_solutions <- get_representative_solutions(
  solutions_matrix_aris,
  split_vector = c(18), ## one split means 2 representative solutions
  order = meta_cluster_order,
  extended_solutions_matrix
)

mc_manhattan_plot(
  rep_solutions,
  data_list = my_data_list,
  target_list = target_list,
  hide_x_labels = TRUE,
  point_size = 2,
  text_size = 12,
  threshold = 0.05,
  neg_log_pval_thresh = 5
)

