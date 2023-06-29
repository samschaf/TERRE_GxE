### Samantha Schaffner
### Aug 10, 2022

brain_blood_expr_perm <- function(DMR_genes, bg_genes, permutation_number){

ngenes <- length(DMR_genes)  
  
#### calculate the proportion of input DMR genes co-expressed in brain and blood at >= 0.5 TPM
#subset GTEx data to DMR genes
GTEx_sub <- GTEx[GTEx$ext_gene %in% DMR_genes,]
#check how many of the DMR genes are expressed in blood
nblood <- nrow(blood_expr <- GTEx_sub[GTEx_sub$Whole.Blood>0.5,]) #15
#check how many of the DMR genes are expressed in brain
nbrain <- nrow(brain_expr <- GTEx_sub[GTEx_sub$Brain...Amygdala>0.5 | GTEx_sub$Brain...Anterior.cingulate.cortex..BA24.>0.5 | GTEx_sub$Brain...Caudate..basal.ganglia.>0.5 | GTEx_sub$Brain...Cerebellar.Hemisphere>0.5 | GTEx_sub$Brain...Cerebellum>0.5 | GTEx_sub$Brain...Cortex>0.5 | GTEx_sub$Brain...Frontal.Cortex..BA9.>0.5 | GTEx_sub$Brain...Hippocampus>0.5 | GTEx_sub$Brain...Hypothalamus>0.5 | GTEx_sub$Brain...Nucleus.accumbens..basal.ganglia.>0.5 | GTEx_sub$Brain...Putamen..basal.ganglia.>0.5 | GTEx_sub$Brain...Spinal.cord..cervical.c.1.>0.5 | GTEx_sub$Brain...Substantia.nigra>0.5,]) #26
nbrainblood <- nrow(brain_blood_expr <- blood_expr[blood_expr$ens_gene %in% brain_expr$ens_gene,]) #14
propbrainblood <- nbrainblood/ngenes

#### calculate the proportion of randomly sampled input genes co-expressed in brain and blood at >= 0.5 TPM
# permuted
bootstrap_exp <-lapply(1:permutation_number, function(x){
set.seed(x)
random_genes <- sample(bg_genes, ngenes)
#subset GTEx data to DMR genes
GTEx_sub <- GTEx[GTEx$ext_gene %in% random_genes,]
#check how many of the DMR genes are expressed in blood
nblood <- nrow(blood_expr <- GTEx_sub[GTEx_sub$Whole.Blood>0.5,])
#check how many of the DMR genes are expressed in brain
nbrain <- nrow(brain_expr <- GTEx_sub[GTEx_sub$Brain...Amygdala>0.5 | GTEx_sub$Brain...Anterior.cingulate.cortex..BA24.>0.5 | GTEx_sub$Brain...Caudate..basal.ganglia.>0.5 | GTEx_sub$Brain...Cerebellar.Hemisphere>0.5 | GTEx_sub$Brain...Cerebellum>0.5 | GTEx_sub$Brain...Cortex>0.5 | GTEx_sub$Brain...Frontal.Cortex..BA9.>0.5 | GTEx_sub$Brain...Hippocampus>0.5 | GTEx_sub$Brain...Hypothalamus>0.5 | GTEx_sub$Brain...Nucleus.accumbens..basal.ganglia.>0.5 | GTEx_sub$Brain...Putamen..basal.ganglia.>0.5 | GTEx_sub$Brain...Spinal.cord..cervical.c.1.>0.5 | GTEx_sub$Brain...Substantia.nigra>0.5,])
nbrainblood_rnd <- nrow(brain_blood_expr <- blood_expr[blood_expr$ens_gene %in% brain_expr$ens_gene,])
propbrainblood_rnd <- nbrainblood_rnd/ngenes
propbrainblood_rnd })

bootstrap_exp <-do.call(rbind, bootstrap_exp)

print("Permutation P values for enrichment and depletion")
  #how many iterations of bootstrapping found MORE genes coexpressed than the real data? Divide this by the number of permutations to get a p-value
  enrich_p<- length(which(bootstrap_exp>=propbrainblood))/permutation_number
  #how many iterations of bootstrapping found MORE genes coexpressed than the real data? Divide this by the number of permutations to get a p-value
  depletion_p<- length(which(bootstrap_exp<=propbrainblood))/permutation_number
  print(paste("Enrichment: ", enrich_p, "; Depletion ", depletion_p, sep="")) }
