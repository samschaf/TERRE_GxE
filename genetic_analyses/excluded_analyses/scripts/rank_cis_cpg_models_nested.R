library(data.table)
library(parallel)
library(tidyverse)
library(mice)
# Global data
num_cores <- 64L
setDTthreads(num_cores)
probe_pos <- fread("~/prs_ewas_integration/cis_mQTL_analyses/terre_data/probe_pos.txt")
snp_pos <- fread("~/prs_ewas_integration/cis_mQTL_analyses/terre_data/snp_pos.txt", key = "CHR,POS")
methy <- fread(
  "~/prs_ewas_integration/cis_mQTL_analyses/terre_data/methylation_combat.txt",
  key = "cpg"
)
probe_pos <- probe_pos[geneid %chin% methy$cpg]
geno <- fread(
  "~/prs_ewas_integration/cis_mQTL_analyses/terre_data/all_imputed_matrixeQTL.txt",
  key = "SNP"
)

shared_covariates <- data.table::transpose(
  na.omit(
    fread(
      "~/prs_ewas_integration/cis_mQTL_analyses/terre_data/covariates_CTP_PD.txt",
      fill = TRUE
    )
  ),
  make.names = "id"
)
env_data <- fread("/home1/NEURO/SHARE_DECIPHER/TERRE_pesticides/pesticides.csv")
mapping <- fread("/home1/NEURO/SHARE_DECIPHER/terre_meta_master.csv")

mapping$IID <- gsub("_PAE.*", "", mapping$IID)
env_data$num <- mapping$IID[match(env_data$num, mapping$patient)]
env_data <- env_data[!is.na(num)]
top_10_cases <- sort(colSums(env_data[, -c(1)], na.rm = T), decreasing = TRUE, index.return = TRUE)$ix[1:5]
envs <- colnames(env_data)[-c(1)][top_10_cases] # top represented
envs <- unique(c(envs, "i_ochl", "naclo4", "as", "h_triazine", "h_uree", "f_amide", "f_dithiocarb")) # environments of interest
pest_missing <- read.csv("/home1/NEURO/SHARE_DECIPHER/TERRE_pesticides/pesticides.csv")
pest_imputed <- read.csv("/home1/NEURO/SHARE_DECIPHER/TERRE_pesticides/pesticides_imputed.csv")

pre_mids <- rbind(cbind(data.frame(X_imputation_ = 0), pest_missing)[, colnames(pest_imputed)], pest_imputed)
pre_mids$num <- as.character(pre_mids$num)
pre_mids <- pre_mids[pre_mids$num %in% mapping$patient & !is.na(mapping$IID[match(pre_mids$num, mapping$patient)]), ]
pest_mids <- suppressWarnings(as.mids(pre_mids, .imp = "X_imputation_", .id = "num"))

# Utilities
cis_snps <- function(cpg) {
  i <- which(probe_pos$geneid == cpg)
  s1 <- probe_pos$s1[i]
  s2 <- probe_pos$s2[i]
  chr <- probe_pos$chr[i]
  snp_pos[CHR == chr & (POS >= s1 - 75000 & POS <= s2 + 75000), SNP]
}
generate_inputs <- function(row, sex) {
  set.seed(NULL)
  df <- data.frame(
    y = unlist(methy[row$cpg, -c(1), on = "cpg"])
  )
  covar <- shared_covariates
  if (!is.na(row$env)) {
    env <- row$env
    tmp_df <- complete(pest_mids, action = "long", include = TRUE)[, c(".imp", ".id", env)]
    colnames(tmp_df) <- c(".imp", ".id", "E")
    ix <- match(mapping$IID[match(pre_mids$num, mapping$patient)], colnames(geno)[-c(1)])
    tmp_df <- cbind(tmp_df, covar[ix, ], y = df$y[ix])
    if (grepl("G", row$model)) {
      G <- unlist(geno[row$SNP, -c(1), on = "SNP"])
      tmp_df$G <- G[ix]
      if (grepl("^male", sex)) {
        tmp_df <- tmp_df[tmp_df$men == 1, ]
      } else {
        tmp_df <- tmp_df[tmp_df$men == 0, ]
      }
      maf <- sum(tmp_df$G) / (2 * nrow(tmp_df))
      if (maf > 0.99 | maf < 0.01) {
        return(data.table())
      }
      exposures <- sum(tmp_df$E, na.rm = T) / nrow(tmp_df)
      if (exposures > 0.99 | exposures < 0.01) {
        return(data.table())
      }
      if (row$model == "GxE") { # GxE
        tmp_df$GxE <- tmp_df$G * tmp_df$E
        df_mids <- suppressWarnings(as.mids(tmp_df))
      } else { # G+E
        df_mids <- suppressWarnings(as.mids(tmp_df))
      }
    } else { # E only
      if (grepl("^male", sex)) {
        tmp_df <- tmp_df[tmp_df$men == 1, ]
      } else {
        tmp_df <- tmp_df[tmp_df$men == 0, ]
      }
      exposures <- sum(tmp_df$E, na.rm = T) / nrow(tmp_df)
      if (exposures > 0.99 | exposures < 0.01) {
        return(data.table())
      }
      df_mids <- suppressWarnings(as.mids(tmp_df))
    }
    return(df_mids)
  } else { # G only
    df$G <- unlist(geno[row$SNP, -c(1), on = "SNP"])
    df <- na.omit(cbind(df, covar))
    if (grepl("^male", sex)) {
      df <- df[df$men == 1, ]
    } else {
      df <- df[df$men == 0, ]
    }
    df$men <- NULL
    maf <- sum(df$G) / (2 * nrow(df))
    if (maf > 0.99 | maf < 0.01) {
      return(data.table())
    }
    return(df)
  }
}


fit_model <- function(rows) { # row of snp cpg model env data.frame
  output <- lapply(
    1:nrow(rows),
    function(i) {
      row <- rows[i, ]
      set.seed(NULL)
      df <- generate_inputs(row, argv[[2]])
      if (is_empty(df)) {
        return(data.table())
      }
      base_model <- with(df, lm(y ~ PD + V3 + V4 + V5 + age + head_trauma + CTP_PC1 + CTP_PC2 + CTP_PC3 + CTP_PC4 + CTP_PC5))
      fit_gxe <- with(df, lm(y ~ PD + G + E + GxE + V3 + V4 + V5 + age + head_trauma + CTP_PC1 + CTP_PC2 + CTP_PC3 + CTP_PC4 + CTP_PC5))
      fit_ge <- with(df, lm(y ~ PD + G + E + V3 + V4 + V5 + age + head_trauma + CTP_PC1 + CTP_PC2 + CTP_PC3 + CTP_PC4 + CTP_PC5))
      fit_e <- with(df, lm(y ~ PD + E + V3 + V4 + V5 + age + head_trauma + CTP_PC1 + CTP_PC2 + CTP_PC3 + CTP_PC4 + CTP_PC5))
      fit_g <- with(df, lm(y ~ PD + G + V3 + V4 + V5 + age + head_trauma + CTP_PC1 + CTP_PC2 + CTP_PC3 + CTP_PC4 + CTP_PC5))
      f_test_gxe <- anova(fit_gxe, fit_ge)
      f_test_gee <- anova(fit_ge, fit_e)
      f_test_geg <- anova(fit_ge, fit_g)
      f_test_g <- anova(fit_g, base_model)
      f_test_e <- anova(fit_e, base_model)
      f_stat_gxe <- f_test_gxe$out$`1 ~~ 2`$result[1]
      f_p_gxe <- f_test_gxe$out$`1 ~~ 2`$result[4]
      f_stat_gee <- f_test_gee$out$`1 ~~ 2`$result[1]
      f_p_gee <- f_test_gee$out$`1 ~~ 2`$result[4]
      f_stat_geg <- f_test_geg$out$`1 ~~ 2`$result[1]
      f_p_geg <- f_test_geg$out$`1 ~~ 2`$result[4]
      f_stat_g <- f_test_g$out$`1 ~~ 2`$result[1]
      f_p_g <- f_test_g$out$`1 ~~ 2`$result[4]
      f_stat_e <- f_test_e$out$`1 ~~ 2`$result[1]
      f_p_e <- f_test_e$out$`1 ~~ 2`$result[4]
      res <- data.table(
        snp = row$SNP,
        cpg = row$cpg,
        env = row$env,
        f_stat_gxe,
        f_p_gxe,
        f_stat_gee,
        f_p_gee,
        f_stat_geg,
        f_p_geg,
        f_stat_g,
        f_p_g,
        f_stat_e,
        f_p_e
      )
      return(res)
    }
  )
  return(rbindlist(output))
}

manifest_from_cpg <- function(cpg) {
  snps <- cis_snps(cpg)
  if (length(snps) > 0) {
    ge_models <- expand_grid(cpg, SNP = snps, env = envs, model = c("GxE"))
    return(data.table(ge_models))
  } else {
    return(data.table())
  }
}
argv <- commandArgs(trailingOnly = TRUE)
probes <- fread(argv[[1]])$cpg


# Generate manifest of all regressions
manifest_list <- mclapply(probes, manifest_from_cpg, mc.cores = num_cores)
manifest <- rbindlist(manifest_list)
head(manifest)

methy <- methy[unique(probes), on = "cpg"]
geno <- geno[unique(manifest$SNP), on = "SNP"]
gc()
# Run all regressions in parallel
setDTthreads(1)
out_path <- paste0(argv[[2]], "nested_dnam_breakdown.txt")
if (!file.exists(out_path)) {
  gc()
  out_list <- mclapply(
    splitIndices(nrow(manifest), num_cores),
    function(i) {
      fit_model(manifest[i, ])
    },
    mc.cores = num_cores
  )
}
setDTthreads(num_cores)
fwrite(rbindlist(out_list), out_path, sep = "\t", row.names = F, quote = F)
