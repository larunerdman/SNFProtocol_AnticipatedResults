
## set working directory (file where your data is)
setwd("C:/Users/ERDP5G/Desktop/Manuscripts/SNFProtocol/")

## read in example data
in_dat = readRDS("IBD_Example1and2.rds")

### install metasnf if needed
# devtools::install_github("BRANCHlab/metasnf")
library(metasnf)
library(SNFtool)
library(ggplot2)

## generate list of data that will be integrated
my_data_list = generate_data_list(
  list(
    data = in_dat[["Baseline"]],
    name = "BaselineData",
    domain = "Baseline",
    type = "continuous"
  ),
  list(
    data = in_dat[["Medication_History"]],
    name = "MedicationData",
    domain = "Medication",
    type = "mixed"
  ),
  uid = "sample_ids"
)

## summarize the data list you just created
summarize_dl(data_list = my_data_list)


      ###
      ### DEFINING TARGET VARIABLES AND CONFOUNDER VARIABLES
      ###

## Define your targets (outcomes)
str(in_dat[["Outcomes"]])
names(in_dat[["Outcomes"]])


cat_outcomes = c("sample_ids","IBD_dx","CD_UC_dx","PGA_dx", "extent_dx", "CD_loc_dx",
                 "L4a_dx", "L4b_dx", "disease_behavior_dx","B1_disease",
                 "Perianal_dx", "post_CS", "endo_healing",
                 "SFR_52_remission", "PGA_52_remission")
cont_outcomes = c("sample_ids","Age_dx")


## Define your confounders
str(in_dat[["Confounders"]])
names(in_dat[["Confounders"]])
cat_confounders = c("sample_ids","race", "sex", "insurance_provider")
cont_confounders = c("sample_ids","dep_index")

## create list of data you will use to compare with your SNF results
target_list <- generate_target_list(
  list(in_dat[["Outcomes"]][,cat_outcomes], "ClinicalCatTargets", "categorical"),
  list(in_dat[["Outcomes"]][,cont_outcomes], "ClinicalContTargets", "numeric"),
  list(in_dat[["Confounders" ]][,cat_confounders], "ConfoundersCat", "categorical"),
  list(in_dat[["Confounders" ]][,cont_confounders], "ConfoundersCont", "numeric"),
  uid = "sample_ids"
)

## summarize your target list
summarize_target_list(target_list)


      ###
      ###  PERFORMING META-SNF
      ###

## generate a matrix of settings which you will use for each SNF run
settings_matrix <- generate_settings_matrix(
  my_data_list,
  nrow = 20,
  min_k = 20,
  max_k = 50,
  seed = 42
)

## run SNF with each of the settings
solutions_matrix = batch_snf(my_data_list,
                             settings_matrix)

## cluster each of your SNF runs
cluster_solutions = get_cluster_solutions(solutions_matrix)

## create a matrix describine each of the SNF runs
solutions_matrix_aris <- calc_om_aris(solutions_matrix)

## plot a heatmap showing how similar each of the SNF runs are
adjusted_rand_index_heatmap(solutions_matrix_aris)

## extract order from SNF similarity heatmap (similarity of clustering solutions)
meta_cluster_order <- get_heatmap_order(solutions_matrix_aris)

## plot settings matrix with same order at SNF similarity
settings_matrix_heatmap(settings_matrix, order = meta_cluster_order)

## run SNF getting each of the SNF similarity matrices using each setting
batch_snf_results <- batch_snf(
  my_data_list,
  settings_matrix,
  return_similarity_matrices = TRUE
)

## extract the solutions matrix (same as what you have above)
solutions_matrix <- batch_snf_results$"solutions_matrix"

## extract all the similarity matrices -- in case you'd like to plot a heatmap or
  ## analyze the similarity matrix itself
similarity_matrices <- batch_snf_results$"similarity_matrices"

## compute silhouette scores on the SNF solutions
silhouette_scores <- calculate_silhouettes(solutions_matrix, similarity_matrices)

## If you get an eError: Package "clv" must be installed to use this function.
## so just go ahead and run install.packages("clv")
dunn_indices <- calculate_dunn_indices(solutions_matrix, similarity_matrices)

## Davies-Bouldin index
db_indices <- calculate_db_indices(solutions_matrix, similarity_matrices)

## create subsamples of your data to see how SNF changes with subsets of the data
data_list_subsamples <- subsample_data_list(
  my_data_list,
  n_subsamples = 30, # calculate 30 subsamples
  subsample_fraction = 0.8 # for each subsample, use random 80% of patients
)

## Find adjusted Rand index similarity of subsampled SNF solutions
pairwise_aris <- subsample_pairwise_aris(
  data_list_subsamples,
  settings_matrix
)


fraction_together <- fraction_clustered_together(
  data_list_subsamples,
  settings_matrix,
  solutions_matrix
)

## find p-values of all your SNF cluster solutions and their associations
## with all the target variables -- outcomes and confounders here
extended_solutions_matrix <- extend_solutions(solutions_matrix, target_list)

## extract the p-values
target_pvals <- p_val_select(extended_solutions_matrix)

## create a heatmap of the p-values
pvals_heatmap(target_pvals, order = meta_cluster_order)


#####
####  CHOOSE FINAL CLUSTERING SOLUTION, PLOT HEATMAP SOLUTION + MANHATTAN
#####

## choose your final clustering solution, here we choose the 9th
my_sim = similarity_matrices[[2]]

## extract your cluster assignments
my_clusts = unlist(solutions_matrix[2,grep(pattern = "subject",x = names(solutions_matrix))])

## plot clusters using SNFtool package function
displayClustersWithHeatmap(my_sim,my_clusts)


####
##    metasnf HEATMAP
####

## set figure margins
par(mar=c(3,3,3,3))

## create heatmap of selected similarity matrix
similarity_matrix_heatmap(
  similarity_matrix = my_sim,
  cluster_solution = my_clusts,
  scale_diag = "mean",
  log_graph = TRUE,
  data_list = target_list,
  left_hm = list(
    "IBD Diagnosis" = "IBD_dx",
    "IBD Severity" = "PGA_dx",
    "CD vs UC" = "CD_UC_dx"
  ),
  top_bar = list(
    "Age at Diagnosis" = "Age_dx"
  ),
  # top_hm = list(
  #   "Race" = "race"
  # ),
  heatmap_height = grid::unit(12, "cm"),
  heatmap_width = grid::unit(12, "cm"),
)

####
###   MANHATTAN PLOT
####

## extract p-values to plot (p-values of cluster solutions vs. target variables)
target_pvals <- p_val_select(extended_solutions_matrix)

## plot with bonferroni line
manhattan_plot(target_pvals[2,], threshold = 0.05, bonferroni_line = FALSE)

## plot without bonferroni line
manhattan_plot(target_pvals[2,], threshold = 0.05, bonferroni_line = TRUE)


####
###   Age at Dx Plot
####

## create vector of clusters
my_clusts2 = my_clusts
## change names to match sample IDs
names(my_clusts2) = sub(pattern = "subject_",replacement = "",x = names(my_clusts))

## match order of sample IDs in outcomes data frame
my_clusts2_ordered = my_clusts2[match(outcomes$sample_ids,names(my_clusts2))]

## create a data frame of cluster assignments and age at diagnosis
age_plt_df2 = data.frame(clusters = my_clusts2_ordered, Age_dx = outcomes$Age_dx)
head(age_plt_df2)

## make clusters into a factor variable
age_plt_df2$clusters = factor(age_plt_df2$clusters)

theme_set(
  theme_bw(base_size = 20)
)

## plot density of age at diangosis vs. cluster assignment
ggplot(age_plt_df2, aes(x = Age_dx, col = clusters)) +
  geom_density(lwd = 1.5)


####
###   CD vs UC plot
####

## create vector of clusters
my_clusts2 = my_clusts
## change names to match sample IDs
names(my_clusts2) = sub(pattern = "subject_",replacement = "",x = names(my_clusts))

## match order of sample IDs in outcomes data frame
my_clusts2_ordered = my_clusts2[match(outcomes$sample_ids,names(my_clusts2))]

## create a data frame of cluster assignments and age at diagnosis
cduc_plt_df2 = data.frame(clusters = my_clusts2_ordered, CDvsUC = outcomes$CD_UC_dx)
head(cduc_plt_df2)

## make clusters into a factor variable
cduc_plt_df2$clusters = factor(cduc_plt_df2$clusters)
cduc_plt_df2$CDvsUC = factor(cduc_plt_df2$CDvsUC, levels = 1:2,
                             labels = c("Crohn's Disease","Ulcerative Colitis/IBD-U"))
str(cduc_plt_df2$CDvsUC)

theme_set(
  theme_bw(base_size = 20)
)

## plot density of age at diangosis vs. cluster assignment
ggplot(cduc_plt_df2, aes(x = clusters, fill = CDvsUC)) +
  geom_bar(stat="count",position="dodge")

ggplot(cduc_plt_df2, aes(x = CDvsUC, fill = clusters)) +
  geom_bar(stat="count",position="dodge")


####
###   PGA plot
####

## create vector of clusters
my_clusts2 = my_clusts
## change names to match sample IDs
names(my_clusts2) = sub(pattern = "subject_",replacement = "",x = names(my_clusts))

## match order of sample IDs in outcomes data frame
my_clusts2_ordered = my_clusts2[match(outcomes$sample_ids,names(my_clusts2))]

## create a data frame of cluster assignments and age at diagnosis
pga_plt_df2 = data.frame(clusters = my_clusts2_ordered, PGA = outcomes$PGA_dx)
head(pga_plt_df2)

## make clusters into a factor variable
pga_plt_df2$clusters = factor(pga_plt_df2$clusters)
pga_plt_df2$IBDSeverity = factor(pga_plt_df2$PGA, levels = 2:4,
                             labels = c("mild","moderate","severe"))
str(pga_plt_df2$IBDSeverity)

theme_set(
  theme_bw(base_size = 20)
)

## plot density of age at diangosis vs. cluster assignment
ggplot(pga_plt_df2, aes(x = clusters, fill = IBDSeverity)) +
  geom_bar(stat="count",position="dodge")

ggplot(pga_plt_df2, aes(x = IBDSeverity, fill = clusters)) +
  geom_bar(stat="count",position="dodge")

