#includes full model used for TERRE EWAS, minus the weights

CMR_rlm_perm <- function(CMR_medB, meta, nDM_CMRs, nperm){
  
  rndDM_CMRs <- sapply(1:nperm, function(x){
  
  set.seed(x)
  meta$PD_null <- sample(meta$PD, nrow(meta))
  
  # run EWAS and create data frame to hold results
  print(paste("Running EWAS: permutation", x, "of", nperm, sep=" "))
  RLM <- lapply(1:ncol(CMR_medB), function(x) {
     rlm(CMR_medB[,x] ~ meta$PD_null + meta$age + meta$smoking + meta$alcohol1 + meta$head_trauma + meta$plate + meta$SentrixPosition_A + meta$CTP_PC1 + meta$CTP_PC2 + meta$CTP_PC3 + meta$CTP_PC4 + meta$CTP_PC5 + meta$CTP_PC6 + meta$ancestryPC1 + meta$ancestryPC2 + meta$ancestryPC3, psi = psi.huber, maxit=500)
  })
  pvals <- sapply(1:ncol(CMR_medB), function(x) {
    f.robftest(RLM[[x]], var = names(RLM[[x]]$coefficients)[2])$p.value })
  coefs <- sapply(1:ncol(CMR_medB), function(x) {
    RLM[[x]]$coefficients[2] })
  results <- data.frame(TargetID=colnames(CMR_medB), pval=pvals, adjP_BH=p.adjust(pvals, method="BH"), adjDB=coefs)
  
  # calculate number of CMRs passing thresholds
  DM_CMRs <- nrow(results[results$adjP_BH<=0.05 & abs(results$adjDB)>=0.03,])
  return(DM_CMRs) })
  
  # permutation p-value (enrichment/depletion)
  print("Permutation P values")
  #how many iterations of bootstrapping found MORE DM-CMRs than the real data? Divide this by the number of permutations to get a p-value
  enrich_p<- length(which(rndDM_CMRs>=nDM_CMRs))/nperm
  #how many iterations of bootstrapping found LESS DM-CMRs than the real data? Divide this by the number of permutations to get a p-value
  depletion_p<- length(which(rndDM_CMRs<=nDM_CMRs))/nperm
  print(paste("Enrichment: ", enrich_p, "; Depletion ", depletion_p, sep=""))
}
  


