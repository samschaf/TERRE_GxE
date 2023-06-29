library(data.table)
library(parallel)
probes <- fread("~/cis_mQTL_analyses/terre_data/probe_pos.txt")
snps <- fread("~/cis_mQTL_analyses/terre_data/snp_pos.txt")
# terre_cpgs <- fread("~/analysis_of_pd_dnam/gcta_analysis/terre_sig_DMR_sex_stratified.txt")
# terre_cpgs$chr <- gsub("chr","",terre_cpgs$chr)
all_cpgs <- colnames(fread("~/analysis_of_pd_dnam/gcta_analysis/terre_methylation_mvalue.phen", nrows = 1))[-c(1, 2)]
grm_cmd <- paste(
  "~/plink2",
  "--bfile", "~/analysis_of_pd_dnam/terre.hard_call",
  "--extract", "~/analysis_of_pd_dnam/gcta_analysis/terre_snps/%s_snps.txt",
  "--thread-num", "16",
  "--make-grm-list",
  "--out", "~/analysis_of_pd_dnam/gcta_analysis/terre_grm/%s",
  "&& gzip ~/analysis_of_pd_dnam/gcta_analysis/terre_grm/%s.grm"
)
cmd <- paste(
  "~/reacta",
  "--grm", "~/analysis_of_pd_dnam/gcta_analysis/terre_grm/%s",
  "--reml",
  "--pheno", "~/analysis_of_pd_dnam/gcta_analysis/terre_methylation_mvalue.phen",
  "--qcovar", "~/analysis_of_pd_dnam/gcta_analysis/terre_covariates.cov",
  "--mpheno", "%d",
  "--thread-num", "16",
  "--out", "~/analysis_of_pd_dnam/gcta_analysis/terre_output/%s"
)
run_reml <- function(i) {
  cpg <- probes$geneid[i]
  s1 <- probes$s1[i]
  s2 <- probes$s2[i]
  chr <- probes$chr[i]
  mpheno <- which(all_cpgs == cpg)
  if (!file.exists(sprintf("gcta_analysis/terre_snps/%s_snps.txt", cpg))) {
    fwrite(
      snps[CHR == chr & (POS >= s1 - 75000 & POS <= s2 + 75000), .(SNP)],
      sprintf("gcta_analysis/terre_snps/%s_snps.txt", cpg),
      quote = F,
      col.names = F,
      row.names = F
    )
  }
  if (!file.exists(sprintf("~/analysis_of_pd_dnam/gcta_analysis/terre_grm/%s.grm.gz", cpg))) {
    tmp <- sprintf(grm_cmd, cpg, cpg, cpg)
    system(tmp, intern = FALSE)
  }
  if (!file.exists(sprintf("~/analysis_of_pd_dnam/gcta_analysis/terre_output/%s.hsq", cpg))) {
    tmp <- sprintf(cmd, cpg, mpheno, cpg)
    system(tmp, intern = FALSE)
  }
  system(sprintf("rm ~/analysis_of_pd_dnam/gcta_analysis/terre_grm/%s*", cpg))
  system(sprintf("rm ~/analysis_of_pd_dnam/gcta_analysis/terre_snps/%s_snps.txt", cpg))
}
out <- mclapply(1:nrow(probes), run_reml, mc.cores = 64)
