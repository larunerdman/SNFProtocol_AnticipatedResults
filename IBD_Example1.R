###
###   IMPORTING NEEDED PACKAGES
###

library(SNFtool)
library(cluster)
library(RColorBrewer)
library(reshape2)
library(ggplot2)

### install metasnf if needed
# devtools::install_github("BRANCHlab/metasnf")
library(metasnf)

###
###   IMPORTING DATA
###

setwd("C:/Users/ERDP5G/Desktop/Manuscripts/SNFProtocol/")

## This data is a list with separate data frames for each set of 
  ## features to be integrated or compared to the SNF solution
ibd_list = readRDS("IBD_Example1and2.rds")

names(ibd_list)

###
###  CREATE CORRELATION MATRIX USING metasnf PACKAGE
###

## separate data types in data list
datatype_specific_data_list = generate_data_list(
  list(
    data = ibd_list[["Baseline"]],
    name = "BaselineData",
    domain = "Baseline",
    type = "continuous"
  ),
  list(
    data = ibd_list[["Medication_History"]][,c("sample_ids","antiTNF_exp","antiTNF_1st",
                                             "bio_naive.IFX","maint_IFX_lvl","ADA_maint")],
    name = "MedicationData",
    domain = "Medication",
    type = "categorical"
  ),
  list(
    data = ibd_list[["Medication_History"]][,c("sample_ids","IFX_startdose_before_age5",
                                             "IFX_startdose_before_age10",
                                             "ADA_startdose_before_age10",
                                             "VEDO_startdose_before_age5",
                                             "VEDO_startdose_before_age10")],
    name = "MedicationData",
    domain = "Medication",
    type = "continuous"
  ),
  uid = "sample_ids"
)

## calculate associations between features that will be integrated
my_assoc_matrix <- calculate_associations(datatype_specific_data_list)

## plot heatmap of correlations between features (p-values)
heatmap <- correlation_pval_heatmap(
  my_assoc_matrix
)


###
###   RUN SNF using SNFtool package 
###     (SNF will be run using metasnf in Example 2)
###

### divide list data to extract data sources to be integrated
  ## each of these is an R data.frame object
baseline = ibd_list[["Baseline"]]
meds = ibd_list[["Medication_History"]]

### check data contained in each of the data sources
names(baseline)
names(meds)

### confirm each data source contains the same samples 
identical(baseline$sample_ids, meds$sample_ids)

### Now extract target data (outcomes) and confounders
outcomes = ibd_list[["Outcomes"]]
### inspect contents of outcomes data frame
names(outcomes)
### check that the samples in the outcome data is equivalent to the baseline data
identical(baseline$sample_ids, outcomes$sample_ids)

### extract confounders
confounders = ibd_list[["Confounders"]]
### inspect contents of confounder data
names(confounders)
### check that the sample IDs of the confound data frame are equivalent 
  ### to those in the baseline data
identical(baseline$sample_ids, confounders$sample_ids)

##Set the hyper-parameters defaults for your SNF run : 
K=20 ##number of neighbors,must be greater than 1. usually(10~30) 
alpha=0.5 ##hyperparameter, usually (0.3~0.8) 
T=20 ###Number of Iterations, usually (10~50)

### Normalize baseline data 
baseline_norm = standardNormalization(baseline[,-1])
str(baseline_norm)
baseline_dist = dist2(X = baseline_norm, C = baseline_norm)
baseline_aff = affinityMatrix(baseline_dist,K = K,sigma = alpha)
estimateNumberOfClustersGivenGraph(baseline_aff,NUMC = 2:5)
baseline_clust3 = spectralClustering(baseline_aff,K = 3)

displayClustersWithHeatmap(baseline_dist, group = baseline_clust3,
                           col = brewer.pal(n = 9, name = "Blues"))

displayClustersWithHeatmap(baseline_aff, group = baseline_clust3,
                           col = brewer.pal(n = 9, name = "Blues"))

meds_dist = as.matrix(daisy(meds[,-1], metric = c("gower")))
meds_aff = affinityMatrix(meds_dist,K = K,sigma = alpha)
estimateNumberOfClustersGivenGraph(meds_aff,NUMC = 2:5)
meds_clust5 = spectralClustering(meds_aff,K = 4)

displayClustersWithHeatmap(meds_dist, group = meds_clust5,
                           col = brewer.pal(n = 9, name = "PuRd"))

displayClustersWithHeatmap(meds_aff, group = meds_clust5,
                           col = brewer.pal(n = 9, name = "PuRd"))


###
###   Plotting input data 
###

  ## plot baseline
dim(baseline_norm)
rownames(baseline_norm) = baseline[,1]
bl.mlt = melt(baseline_norm)
str(bl.mlt)

theme_set(
  theme_minimal(base_size = 5)
)

ggplot(bl.mlt, aes(x = Var2, y = Var1, fill = value)) + 
  geom_tile(show.legend = FALSE) +
  scale_fill_gradientn(colours = c("blue","white","red"))

  ## plot medication
rownames(meds) = meds$sample_ids
meds.mlt = melt(meds,id.vars = "sample_ids")
str(meds.mlt)

theme_set(
  theme_minimal(base_size = 5)
)

ggplot(meds_bin.mlt, aes(x = variable, y = sample_ids, fill = value)) + 
  geom_tile(show.legend = FALSE) +
  scale_fill_gradientn(colours = c("blue","white","red"))


###
###   SNF integration
###

## double check that you still have the same number of patients in your 
  ## similarity matrices (311 and 311)
dim(baseline_aff)
dim(meds_aff)

## SNF step 
ibd_snf = SNF(list(baseline_aff,meds_aff),K = K,t = T)

estimateNumberOfClustersGivenGraph(ibd_snf,NUMC = 2:5)

clusters2 = spectralClustering(ibd_snf,K = 2)
clusters3 = spectralClustering(ibd_snf,K = 3)
clusters4 = spectralClustering(ibd_snf,K = 4)
clusters5 = spectralClustering(ibd_snf,K = 5)

displayClustersWithHeatmap(ibd_snf, group = clusters2, 
                           col = brewer.pal(name = "OrRd",n =9)[c(1,2,5,6,7,8,9)])
displayClustersWithHeatmap(ibd_snf, group = clusters3, 
                           col = brewer.pal(name = "OrRd",n =9)[c(1,2,5,6,7,8,9)])

###
### ALLUVIAL PLOT OF CLUSTER OUTCOMES
###

alluv_clust_df = data.frame(clusters2,clusters3,clusters4,clusters5)
## get frequencies of each combination 

## make below into a data frame with frequencies of each combo: 
# table(paste0(clusters2,":",clusters3,":",clusters4,":",clusters5))


make_all_df = function(clust_df){
  
  count_tbl = table(apply(clust_df,1,function(x){paste(x, collapse = ":")}))
  
  out_df = data.frame(matrix(nrow = length(count_tbl),ncol = ncol(clust_df)))
  names(out_df) = paste0("Cluster",2:(ncol(clust_df)+1))
  out_df$Freq = NA
  
  for(i in 1:length(count_tbl)){
    # browser()
    out_df[i,] = c(unlist(strsplit(names(count_tbl)[i],split = ":")),count_tbl[i])
  }

  return(out_df)
}

alluv_df = make_all_df(clust_df = alluv_clust_df)
alluv_df$Freq = as.numeric(alluv_df$Freq)
head(alluv_df)

library(ggalluvial)

ggplot(data = alluv_df, 
       aes(axis1 = Cluster2, axis2 = Cluster3, 
           axis3 = Cluster4, axis4 = Cluster5, y = Freq))+
  scale_x_discrete(limits = c("Cluster2", "Cluster3", 
                              "Cluster4","Cluster5"), expand = c(.2, .05)) + 
  xlab("Cluster Assignment") +
  geom_alluvium(aes(fill = Cluster5)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal(base_size = 15)

### testing confounders
chisq.test(confounders$race,clusters2)
chisq.test(confounders$sex,clusters2)
chisq.test(confounders$insurance_provider,clusters2)
t.test(confounders$dep_index,clusters2)

### testing outcomes
chisq.test(outcomes$IBD_dx,clusters2)
t.test(outcomes$Age_dx,clusters2)
chisq.test(outcomes$CD_UC_dx,clusters2)
chisq.test(outcomes$PGA_dx,clusters2)


names(outcomes)

cat_outcomes = c("IBD_dx","CD_UC_dx","PGA_dx", "extent_dx", "CD_loc_dx",
                 "L4a_dx", "L4b_dx", "disease_behavior_dx","B1_disease",
                 "Perianal_dx", "post_CS", "endo_healing",
                 "SFR_52_remission", "PGA_52_remission")
cont_outcomes = c("Age_dx")


cat_confounders = c("race", "sex", "insurance_provider")

cont_confounders = c("dep_index")

cat_vars = c(cat_outcomes, cat_confounders)
cont_vars = c(cont_outcomes, cont_confounders)

all_vars = c(cat_outcomes, cont_outcomes, cat_confounders, cont_confounders)

p_df = data.frame(matrix(ncol = 5, nrow = length(all_vars)))
names(p_df) = c("Features","p.value2","p.value3","p.value4","p.value5")

p_df$Features = all_vars

for(i in 1:nrow(p_df)){
  
  feature = p_df$Features[i]
  
  if(feature %in% cat_vars){
    if(feature %in% cat_outcomes){
      p_df$p.value2[i] = chisq.test(outcomes[,feature],clusters2)$p.val
      p_df$p.value3[i] = chisq.test(outcomes[,feature],clusters3)$p.val
      p_df$p.value4[i] = chisq.test(outcomes[,feature],clusters4)$p.val
      p_df$p.value5[i] = chisq.test(outcomes[,feature],clusters5)$p.val
      
      # p_df$p.value2[i] = fisher.test(outcomes[,feature],clusters2)$p.val
      # p_df$p.value3[i] = fisher.test(outcomes[,feature],clusters3)$p.val
      # p_df$p.value4[i] = fisher.test(outcomes[,feature],clusters4)$p.val
    } else{
      p_df$p.value2[i] = chisq.test(confounders[,feature],clusters2)$p.val
      p_df$p.value3[i] = chisq.test(confounders[,feature],clusters3)$p.val
      p_df$p.value4[i] = chisq.test(confounders[,feature],clusters4)$p.val
      p_df$p.value5[i] = chisq.test(confounders[,feature],clusters5)$p.val
      
      # p_df$p.value2[i] = fisher.test(confounders[,feature],clusters2)$p.val
      # p_df$p.value3[i] = fisher.test(confounders[,feature],clusters3)$p.val
      # p_df$p.value4[i] = fisher.test(confounders[,feature],clusters4)$p.val
    }    
  } else{
    if(feature %in% cont_outcomes){
      p_df$p.value2[i] = t.test(outcomes[,feature],clusters2)$p.val
      p_df$p.value3[i] = t.test(outcomes[,feature],clusters3)$p.val
      p_df$p.value4[i] = t.test(outcomes[,feature],clusters4)$p.val
      p_df$p.value5[i] = t.test(outcomes[,feature],clusters5)$p.val
    } else{
      p_df$p.value2[i] = t.test(confounders[,feature],clusters2)$p.val
      p_df$p.value3[i] = t.test(confounders[,feature],clusters3)$p.val
      p_df$p.value4[i] = t.test(confounders[,feature],clusters4)$p.val
      p_df$p.value5[i] = t.test(confounders[,feature],clusters5)$p.val
    }
  }
}

p_df$Variable = NA
p_df$Variable[p_df$Features %in% names(confounders)] = "Confounder"
p_df$Variable[p_df$Features %in% names(outcomes)] = "Outcomes"
p_df$Features = factor(p_df$Features, levels = p_df$Features)

p_df$Clusters2 = -log10(p_df$p.value2)
p_df$Clusters3 = -log10(p_df$p.value3)
p_df$Clusters4 = -log10(p_df$p.value4)
p_df$Clusters5 = -log10(p_df$p.value5)

p_df.mlt = melt(p_df[,c("Features","Variable","Clusters2",
                        "Clusters3","Clusters4","Clusters5")])
str(p_df.mlt)

theme_set(
  theme_bw(base_size = 15)
)

ggplot(p_df.mlt, aes(x = Features, y = value, col = variable)) + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("blue","darkorchid","forestgreen","deeppink2")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_hline(yintercept = -log10(0.05), col = "red", lty = 2) + 
  geom_hline(yintercept = -log10(0.05/19), col = "black", lty = 2) + 
  ylab("-log10(p)")


####
###   GROUP PLOTS  -- ADD IBD DX OR CD/UC DX
####

## age

age_plt_df = data.frame(clusters2,outcomes$Age_dx)
names(age_plt_df) = c("Cluster","AgeatDiagnosis")
age_plt_df$Cluster = factor(age_plt_df$Cluster)

ggplot(age_plt_df, aes(x = AgeatDiagnosis, col = Cluster)) +
  geom_density(lwd = 1.5)

age_plt_df = data.frame(clusters5,outcomes$Age_dx)
names(age_plt_df) = c("Cluster","AgeatDiagnosis")
age_plt_df$Cluster = factor(age_plt_df$Cluster)

ggplot(age_plt_df, aes(x = AgeatDiagnosis, col = Cluster)) +
  geom_density(lwd = 1.5)

## PGA

pga_plt_df = data.frame(clusters2, outcomes$PGA_dx)
names(pga_plt_df) = c("Cluster","PGA")
pga_plt_df$Cluster = factor(pga_plt_df$Cluster)
pga_plt_df$PGA = factor(pga_plt_df$PGA)

ggplot(pga_plt_df, aes(x = PGA, fill = Cluster)) + 
  geom_bar(position = "dodge")

ggplot(pga_plt_df, aes(fill = PGA, x = Cluster)) + 
  geom_bar(position = "dodge")

pga_plt_df = data.frame(clusters5, outcomes$PGA_dx)
names(pga_plt_df) = c("Cluster","PGA")
pga_plt_df$Cluster = factor(pga_plt_df$Cluster)
pga_plt_df$PGA = factor(pga_plt_df$PGA)

ggplot(pga_plt_df, aes(x = PGA, fill = Cluster)) + 
  geom_bar(position = "dodge")

ggplot(pga_plt_df, aes(fill = PGA, x = Cluster)) + 
  geom_bar(position = "dodge")


## IBD
ibd_plt_df = data.frame(clusters2, outcomes$IBD_dx)
names(ibd_plt_df) = c("Cluster","IBDdiagnosis")
ibd_plt_df$Cluster = factor(ibd_plt_df$Cluster)
ibd_plt_df$IBDdiagnosis = factor(ibd_plt_df$IBDdiagnosis)

table(ibd_plt_df$Cluster,ibd_plt_df$IBDdiagnosis)

ggplot(ibd_plt_df, aes(x = IBDdiagnosis, fill = Cluster)) + 
  geom_bar(position = "dodge")

ggplot(ibd_plt_df, aes(fill = IBDdiagnosis, x = Cluster)) + 
  geom_bar(position = "dodge")

ibd_plt_df = data.frame(clusters5, outcomes$IBD_dx)
names(ibd_plt_df) = c("Cluster","IBDdiagnosis")
ibd_plt_df$Cluster = factor(ibd_plt_df$Cluster)
ibd_plt_df$IBDdiagnosis = factor(ibd_plt_df$IBDdiagnosis)

ggplot(ibd_plt_df, aes(x = IBDdiagnosis, fill = Cluster)) + 
  geom_bar(position = "dodge")

ggplot(ibd_plt_df, aes(fill = IBDdiagnosis, x = Cluster)) + 
  geom_bar(position = "dodge")

## CD/UC
cduc_plt_df = data.frame(clusters2, outcomes$CD_UC_dx)
names(cduc_plt_df) = c("Cluster","CDvsUC")
cduc_plt_df$Cluster = factor(cduc_plt_df$Cluster)
cduc_plt_df$CDvsUC = factor(cduc_plt_df$CDvsUC)

table(cduc_plt_df$Cluster,cduc_plt_df$CDvsUC)

ggplot(cduc_plt_df, aes(x = CDvsUC, fill = Cluster)) + 
  geom_bar(position = "dodge")

ggplot(cduc_plt_df, aes(fill = CDvsUC, x = Cluster)) + 
  geom_bar(position = "dodge")

cduc_plt_df = data.frame(clusters5, outcomes$CD_UC_dx)
names(cduc_plt_df) = c("Cluster","CDvsUC")
cduc_plt_df$Cluster = factor(cduc_plt_df$Cluster)
cduc_plt_df$CDvsUC = factor(cduc_plt_df$CDvsUC)

table(cduc_plt_df$Cluster,cduc_plt_df$CDvsUC)

ggplot(cduc_plt_df, aes(x = CDvsUC, fill = Cluster)) + 
  geom_bar(position = "dodge")

ggplot(cduc_plt_df, aes(fill = CDvsUC, x = Cluster)) + 
  geom_bar(position = "dodge")


## Dep index

dep_plt_df = data.frame(clusters2, confounders$dep_index)
names(dep_plt_df) = c("Cluster","DeprivationIndex")
dep_plt_df$Cluster = factor(dep_plt_df$Cluster)
# dep_plt_df$DeprivationIndex = factor(dep_plt_df$DeprivationIndex)

ggplot(dep_plt_df, aes(x = DeprivationIndex, col = Cluster)) +
  geom_density(lwd = 1.5)

dep_plt_df = data.frame(clusters5, confounders$dep_index)
names(dep_plt_df) = c("Cluster","DeprivationIndex")
dep_plt_df$Cluster = factor(dep_plt_df$Cluster)
# dep_plt_df$DeprivationIndex = factor(dep_plt_df$DeprivationIndex)

ggplot(dep_plt_df, aes(x = DeprivationIndex, col = Cluster)) +
  geom_density(lwd = 1.5)

