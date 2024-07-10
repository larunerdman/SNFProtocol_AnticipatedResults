
## set working directory (file where your data is)
setwd("C:/Users/ERDP5G/Desktop/Manuscripts/SNFProtocol/")

### install metasnf if needed
# devtools::install_github("BRANCHlab/metasnf")
library(metasnf)
library(SNFtool)
library(ggplot2)

## read in example data 
in_dat = readRDS("IBD_Example1and2.rds")

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
target_list <- generate_data_list(
  list(in_dat[["Outcomes"]][,cat_outcomes], "ClinicalCatTargets", "ClinicalCatTargets", "categorical"),
  list(in_dat[["Outcomes"]][,cont_outcomes],"ClinicalContTargets", "ClinicalContTargets", "numeric"),
  list(in_dat[["Confounders" ]][,cat_confounders],"ConfoundersCat", "ConfoundersCat", "categorical"),
  list(in_dat[["Confounders" ]][,cont_confounders],"ConfoundersCont", "ConfoundersCont", "numeric"),
  uid = "sample_ids", 
  remove_missing = FALSE
)

## summarize your target list
summarize_dl(target_list)


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

## create a matrix describing the adjusted Rand index for each of the SNF runs
solutions_matrix_aris <- calc_aris(solutions_matrix)

## get order of solutions matrix to apply going forward
meta_cluster_order <- get_matrix_order(solutions_matrix_aris)

heatmap_colours <- colorRampPalette(
  rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu"))
)(100)

## plot a heatmap showing how similar each of the SNF runs are
adjusted_rand_index_heatmap(
  solutions_matrix_aris,
  order = meta_cluster_order,
  col = heatmap_colours,
  show_row_names = TRUE,
  show_column_names = TRUE,
  rect_gp = grid::gpar(col = "gray", lwd = 1)
)


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

## find p-values of all your SNF cluster solutions and their associations
## with all the target variables -- outcomes and confounders here
extended_solutions_matrix <- extend_solutions(solutions_matrix, 
                                              target_list)


## extract the p-values
# get_pvals instead of pval_select
target_pvals <- get_pvals(extended_solutions_matrix)

## create a heatmap of the p-values
pval_heatmap(target_pvals, order = meta_cluster_order)


#####
####  CHOOSE FINAL CLUSTERING SOLUTION, PLOT HEATMAP SOLUTION + MANHATTAN
#####

solution_no = 2

## choose your final clustering solution, here we choose the 6th
my_sim = similarity_matrices[[solution_no]]

## extract your cluster assignments
my_clusts = unlist(solutions_matrix[solution_no, grep(pattern = "subject",x = names(solutions_matrix))])

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
target_pvals <- get_pvals(extended_solutions_matrix)

target_pvals[solution_no, ]

## plot with bonferroni line
esm_manhattan_plot(extended_solutions_matrix[solution_no,], threshold = 0.05, bonferroni_line = TRUE)

## plot without bonferroni line
esm_manhattan_plot(extended_solutions_matrix[solution_no,], threshold = 0.05, bonferroni_line = FALSE)


####
###   Age at Dx Plot
####

## create data frame of outcomes
outcomes = in_dat[["Outcomes"]]

# ## create vector of clusters
# my_clusts6 = my_clusts
## change names to match sample IDs
names(my_clusts) = sub(pattern = "subject_",replacement = "",x = names(my_clusts))

## match order of sample IDs in outcomes data frame
my_clusts_ordered = my_clusts[match(outcomes$sample_ids,names(my_clusts))]

## create a data frame of cluster assignments and age at diagnosis
age_plt_df = data.frame(clusters = my_clusts_ordered, Age_dx = outcomes$Age_dx)
head(age_plt_df)

## make clusters into a factor variable
age_plt_df$clusters = factor(age_plt_df$clusters, levels = 1:max(age_plt_df$clusters))

theme_set(
  theme_bw(base_size = 25)
)

## plot density of age at diangosis vs. cluster assignment
ggplot(age_plt_df, aes(x = Age_dx, col = clusters)) + 
  geom_density(lwd = 1.5)


####
###   CD vs UC plot
####

## create vector of clusters
# my_clusts = my_clusts
## change names to match sample IDs
# names(my_clusts) = sub(pattern = "subject_",replacement = "",x = names(my_clusts))

## match order of sample IDs in outcomes data frame
my_clusts_ordered = my_clusts[match(outcomes$sample_ids,names(my_clusts))]

## create a data frame of cluster assignments and age at diagnosis
cduc_plt_df = data.frame(clusters = my_clusts_ordered, CDvsUC = outcomes$CD_UC_dx)
head(cduc_plt_df)

## make clusters into a factor variable
cduc_plt_df$clusters = factor(cduc_plt_df$clusters, levels = 1:max(cduc_plt_df$clusters))
cduc_plt_df$CDvsUC = factor(cduc_plt_df$CDvsUC, levels = 1:2, 
                             labels = c("Crohn's Disease","Ulcerative Colitis/IBD-U"))
str(cduc_plt_df$CDvsUC)

theme_set(
  theme_bw(base_size = 25)
)

## plot density of age at diangosis vs. cluster assignment
ggplot(cduc_plt_df, aes(x = clusters, fill = CDvsUC)) + 
  geom_bar(stat="count",position="dodge")

ggplot(cduc_plt_df, aes(x = CDvsUC, fill = clusters)) + 
  geom_bar(stat="count",position="dodge")


####
###   PGA plot
####

## create vector of clusters
# my_clusts2 = my_clusts
## change names to match sample IDs
# names(my_clusts) = sub(pattern = "subject_",replacement = "",x = names(my_clusts))

## match order of sample IDs in outcomes data frame
my_clusts_ordered = my_clusts[match(outcomes$sample_ids,names(my_clusts))]

## create a data frame of cluster assignments and age at diagnosis
pga_plt_df = data.frame(clusters = my_clusts_ordered, PGA = outcomes$PGA_dx)
head(pga_plt_df)

## make clusters into a factor variable
pga_plt_df$clusters = factor(pga_plt_df$clusters)
pga_plt_df$IBDSeverity = factor(pga_plt_df$PGA, levels = 2:4, 
                             labels = c("mild","moderate","severe"))
str(pga_plt_df$IBDSeverity)

theme_set(
  theme_bw(base_size = 20)
)

## plot density of age at diangosis vs. cluster assignment
ggplot(pga_plt_df, aes(x = clusters, fill = IBDSeverity)) + 
  geom_bar(stat="count",position="dodge")

ggplot(pga_plt_df, aes(x = IBDSeverity, fill = clusters)) + 
  geom_bar(stat="count",position="dodge")

