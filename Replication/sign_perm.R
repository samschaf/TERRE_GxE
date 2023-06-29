sign_perm <- function(sig.cpgs, discovery_betas, discovery_DBs, validation_betas, validation_DBs, n_validation, permutation_number=1000){
   
   #sig.cpgs = data frame of CMR EWAS results
   #discovery_betas = data frame or matrix of betas tested for differential methylation discovery (CpGs as rows, samples as columns)
   #discovery_DBs = named vector of delta betas, corresponding to rownames of discovery_betas
   #validation_betas = data frame or matrix of betas tested for differential methylation validation (CpGs as rows, samples as columns)
   #validation_DBs = named vector of delta betas, corresponding to rownames of validation_betas
   #n_validation = number of CMRs that had the same effect direction in validation dataset
 
   bootstrap_effect <-lapply(1:permutation_number, function(x){
    set.seed(x)
   
    rnd_betas <- validation_betas[sample(1:nrow(validation_betas), nrow(sig.cpgs)),]
    rnd_DBs_disc <- discovery_DBs[match(rownames(rnd_betas),names(discovery_DBs))]
    rnd_DBs_val <- validation_DBs[match(rownames(rnd_betas),names(validation_DBs))]
    rnd_results <- data.frame(cpg=rownames(rnd_betas), DB_disc=rnd_DBs_disc, DB_val=rnd_DBs_val)
    rnd_hits <- nrow(rnd_results[sign(rnd_results$DB_disc)==sign(rnd_results$DB_val),]) })
   
   bootstrap_effect <- do.call(rbind, bootstrap_effect)
   
   print("Permutation P values for enrichment and depletion")
   #how many iterations of bootstrapping found MORE genes coexpressed than the real data? Divide this by the number of permutations to get a p-value
   enrich_p<- length(which(bootstrap_effect>=n_validation))/permutation_number
   #how many iterations of bootstrapping found MORE genes coexpressed than the real data? Divide this by the number of permutations to get a p-value
   depletion_p<- length(which(bootstrap_effect<=n_validation))/permutation_number
   print(paste("Enrichment: ", enrich_p, "; Depletion ", depletion_p, sep="")) }
