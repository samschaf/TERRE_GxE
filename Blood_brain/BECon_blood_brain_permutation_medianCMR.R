### Samantha Schaffner
### Sept 26, 2022

brain_blood_meth_perm <- function(nbrain_real, nCMRs, bg_CMRs, becon_ref, permutation_number=1000, verbose=FALSE){

###### Permuted proportion of CMRs co-methylated in blood and brain
  
# permuted
bootstrap_meth <-lapply(1:permutation_number, function(x){
set.seed(x)
if(verbose==TRUE){print(x)}
random_CMRs <- bg_CMRs[sample(1:nrow(bg_CMRs), nCMRs),]

#make a data frame for the random CMRs
random_CMRs$cmr <- as.character(random_CMRs$cmr)
cpgs <- unlist(strsplit(random_CMRs$cmr, split=","))
rnd_cmr_df <- lapply(1:length(cpgs), function(x){
  CMR_row <- random_CMRs[grep(cpgs[x], random_CMRs$cmr),]
  CMR_row$cpg <- cpgs[x]
  return(CMR_row)
})
rnd_cmr_df <- do.call(rbind, rnd_cmr_df)

### calculate correlation statistics across each random CMR
becon_ref_sub <- becon_ref[becon_ref$CpG %in% rnd_cmr_df$cpg,]
rnd_cmr_df_sub <- rnd_cmr_df[match(as.character(becon_ref_sub$CpG), rnd_cmr_df$cpg),]
becon_ref_sub <- cbind(becon_ref_sub, rnd_cmr_df_sub[,57:64])

#summarize correlation over the CMR
random_CMRs$median_BRAIN7 <- NA
random_CMRs$median_BRAIN10 <- NA
random_CMRs$median_BRAIN20 <- NA
random_CMRs$num_probes_in_becon <- NA

for (i in 1:nrow(random_CMRs)){
  cmr_becon <- becon_ref_sub[becon_ref_sub$cmr==random_CMRs$cmr[i],]
  random_CMRs$num_probes_in_becon[i] <- nrow(cmr_becon)
  if(random_CMRs$num_probes_in_becon[i]>0){
    random_CMRs$median_BRAIN7[i] <- median(cmr_becon$BRAIN7)
    random_CMRs$median_BRAIN10[i] <- median(cmr_becon$BRAIN10)
    random_CMRs$median_BRAIN20[i] <- median(cmr_becon$BRAIN20) }
}

random_CMRs$cor_becon <- (abs(random_CMRs$median_BRAIN7)>=0.3|abs(random_CMRs$median_BRAIN10)>=0.3|abs(random_CMRs$median_BRAIN20)>=0.3)
nbrain <- nrow(random_CMRs[random_CMRs$cor_becon==TRUE,])

nbrain })

bootstrap_meth <-do.call(rbind, bootstrap_meth)

print("Permutation P values for enrichment and depletion")
  #how many iterations of bootstrapping found MORE genes coexpressed than the real data? Divide this by the number of permutations to get a p-value
  enrich_p<- length(which(bootstrap_meth>=nbrain_real))/permutation_number
  #how many iterations of bootstrapping found MORE genes coexpressed than the real data? Divide this by the number of permutations to get a p-value
  depletion_p<- length(which(bootstrap_meth<=nbrain_real))/permutation_number
  print(paste("Enrichment: ", enrich_p, "; Depletion ", depletion_p, sep="")) }
