## ----setup, include=FALSE-----------------------------------------------------
library(tidyverse)
library(data.table)
knitr::opts_chunk$set(echo = TRUE)


## -----------------------------------------------------------------------------
terre_ancestry <- fread("../genotype_qc/TERRE_QC/raw_data.geno.maf.mind.sex_check.het_filter.ibd_filter.eigenvec")
terre_meta <- fread("/home1/NEURO/SHARE_DECIPHER/terre_meta_master.csv")
terre_meta$IID <- gsub(".*_(PAE.*)","\\1",terre_meta$IID)
terre_meta <- terre_meta %>% right_join(terre_ancestry,by=c("IID"="V2"))
load("/home1/NEURO/SHARE_DECIPHER/processed_DNAm_data/TERRE/TERRE_betas_combat.RData")#betas_sub
genotype <- fread("terre_data/raw_data.imputed.r2_30.maf_mismatch.traw")


## -----------------------------------------------------------------------------
colnames(genotype) <- gsub(".*_(PAE.*)","\\1",colnames(genotype))
colnames(genotype)
setcolorder(genotype,neworder = c(colnames(genotype)[1:6],terre_meta$IID))
betas_sub <- betas_sub[,colnames(betas_sub) %in% terre_meta$patient]
colnames(betas_sub) <- terre_meta$IID[match(colnames(betas_sub),terre_meta$patient)]
betas_sub <- betas_sub[,match(colnames(betas_sub),terre_meta$IID)]


## -----------------------------------------------------------------------------
all(colnames(betas_sub) == terre_meta$IID)
all(colnames(betas_sub) == colnames(genotype)[-c(1:6)])


## -----------------------------------------------------------------------------
methy_PC <- prcomp(t(betas_sub), center=T,rank.= 20)


## -----------------------------------------------------------------------------
cat_vars <- model.matrix(~0+plate, data=terre_meta %>% mutate(plate= as.factor(plate)))
for(i in 0:20){
  if(i == 0){
    covar<- cbind(cat_vars,terre_meta[,c("V3","V4","V5","age","men")])
  }
  else{
    pcs <- methy_PC$x[,1:i]
    covar<- cbind(pcs,cat_vars,terre_meta[,c("V3","V4","V5","age","men")])
  }
  write_delim(t(covar) %>% as.data.frame() %>% rownames_to_column("id"),sprintf("terre_data/covariates_%d_methy_PC.txt",i))
}


## -----------------------------------------------------------------------------
methy_annot <- fread("~/MethylationEPIC_v-1-0_B4.csv", skip = 7)


## -----------------------------------------------------------------------------
#SNP POS
write_delim(genotype[,.(SNP,CHR=paste0("chr",CHR),POS)],"terre_data/snp_pos.txt")
#SNPs
geno_vars <- colnames(genotype)[-c(1,3:6)]
fwrite(genotype[,..geno_vars],"terre_data/all_imputed_matrixeQTL.txt",sep = "\t",quote = F)
#Methy POS
fwrite(methy_annot[Name %in% rownames(betas_sub),.(geneid=Name, chr=paste0("chr",CHR),s1=MAPINFO,s2=MAPINFO)], "terre_data/probe_pos.txt",sep = "\t",quote=F)
#methy
fwrite(betas_sub %>% as.data.frame() %>%rownames_to_column("cpg"),"terre_data/methylation_combat.txt",sep="\t",quote=F)

#SNP POS
write_delim(genotype[CHR==21,.(SNP,CHR=paste0("chr",CHR),POS)],"terre_data/snp_pos_chr21.txt")
#SNPs
geno_vars <- colnames(genotype)[-c(1,3:6)]
write_delim(genotype[CHR==21,..geno_vars],"terre_data/all_imputed_matrixeQTL_chr21.txt")
#Methy POS
write_delim(methy_annot[Name %in% rownames(betas_sub),.(geneid=Name, chr=paste0("chr",CHR),s1=MAPINFO,s2=MAPINFO)][chr=="chr21"], "terre_data/probe_pos_chr21.txt")
chr21_cpg <- methy_annot[Name %in% rownames(betas_sub) & CHR == 21,]$Name
#methy
write_delim(betas_sub %>% as.data.frame() %>%rownames_to_column("cpg") %>% filter(cpg %in% chr21_cpg),"terre_data/methylation_combat_chr21.txt")
