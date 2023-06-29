library(data.table)
library(parallel)
library(tidyr)
library(dplyr)
library(rlang)
library(tibble)
library(mice)
library(glue)
# Global data
argv <- commandArgs(trailingOnly = TRUE) # cpgs |  male or female

num_cores <- 64L
setDTthreads(num_cores)
manifest <- fread(argv[[1]])
probe_pos <- fread("~/prs_ewas_integration/cis_mQTL_analyses/terre_data/probe_pos.txt")
snp_pos <- fread("~/prs_ewas_integration/cis_mQTL_analyses/terre_data/snp_pos.txt", key = "CHR,POS")
methy <- fread(
  glue("~/prs_ewas_integration/cis_mQTL_analyses/{argv[[2]]}_methy_2022.txt.gz"),
  key = "cpg"
)
probe_pos <- probe_pos[geneid %chin% methy$cpg]
geno <- fread(
  glue("~/prs_ewas_integration/cis_mQTL_analyses/{argv[[2]]}_geno_2022.txt.gz"),
  key = "SNP"
)

shared_covariates <- fread(
  glue("~/prs_ewas_integration/cis_mQTL_analyses/{argv[[2]]}_meta_2022_w_weights.txt.gz"),
)
shared_covariates$plate <- as.character(shared_covariates$plate)
shared_covariates$plate <- as.character(shared_covariates$alcohol1)
shared_covariates$plate <- as.character(shared_covariates$smoking)
env_data <- fread("/home1/NEURO/SHARE_DECIPHER/TERRE_pesticides/pesticides.csv")

env_data$num <- colnames(geno)[-c(1)][match(env_data$num, colnames(geno)[-c(1)])]
env_data <- env_data[!is.na(num)]
pest_missing <- read.csv("/home1/NEURO/SHARE_DECIPHER/TERRE_pesticides/pesticides.csv")
pest_imputed <- read.csv("/home1/NEURO/SHARE_DECIPHER/TERRE_pesticides/pesticides_imputed.csv")

pre_mids <- rbind(cbind(data.frame(X_imputation_ = 0), pest_missing)[, colnames(pest_imputed)], pest_imputed)
pre_mids$num <- as.character(pre_mids$num)
pre_mids <- pre_mids[pre_mids$num %in% colnames(geno), ]
pest_mids <- suppressWarnings(as.mids(pre_mids, .imp = "X_imputation_", .id = "num"))
envs <- colnames(pest_mids$data)
# Utilities
cis_snps <- function(cpg) {
  i <- which(probe_pos$geneid == cpg)
  s1 <- probe_pos$s1[i]
  s2 <- probe_pos$s2[i]
  chr <- probe_pos$chr[i]
  snp_pos[CHR == chr & (POS >= s1 - 75000 & POS <= s2 + 75000), SNP]
}

generate_inputs <- function(row, methy, geno) {
  set.seed(NULL)
  print(row)
  df <- data.frame(
    y = unlist(methy[row$cpg, -c(1), on = "cpg"])
  )
  covar <- shared_covariates
  if (!is.na(row$env) & row$env != "") {
    env <- row$env
    tmp_df <- complete(pest_mids, action = "long", include = TRUE)[, c(".imp", ".id", env)]
    colnames(tmp_df) <- c(".imp", ".id", "E")
    ix <- match(pre_mids$num, colnames(geno)[-c(1)])
    tmp_df <- cbind(tmp_df, covar[ix, ], y = df$y[ix])
    if (grepl("G", row$model)) {
      G <- unlist(geno[row$SNP, -c(1), on = "SNP"])
      tmp_df$G <- G[ix]
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
    maf <- sum(df$G) / (2 * nrow(df))
    if (maf > 0.99 | maf < 0.01) {
      return(data.table())
    }
    return(df)
  }
}


fit_model <- function(row, methy, geno) { # row of snp cpg model env data.frame
  base_formula <- y ~ PD + age + smoking + alcohol1 + head_trauma + CTP_PC1 + CTP_PC2 + CTP_PC3 + CTP_PC4 + ancestryPC1 + ancestryPC2 + ancestryPC3 + plate + SentrixPosition_A
  set.seed(NULL)
  df <- generate_inputs(row, methy, geno)
  if (is_empty(df)) {
    return(data.table())
  }
  if (!is.na(row$env) & row$env != "") {
    if (grepl("G", row$model)) {
      if (row$model == "GxE") { # GxE
        fit <- with(df, lm(formula(format(update(base_formula, ~ . + G + E + GxE))), weights = weights))
      } else { # G+E
        fit <- with(df, lm(formula(format(update(base_formula, ~ . + G + E))), weights = weights))
      }
      base_model <- with(df, lm(formula(format(base_formula)), weights = weights))
    } else { # E only
      base_model <- with(df, lm(formula(format(base_formula)), weights = weights))
      fit <- with(df, lm(formula(format(update(base_formula, ~ . + E))), weights = weights))
    }
    stats <- column_to_rownames(summary(pool(fit)) %>% select(-df), var = "term")
    f_test <- anova(fit, base_model)
    f_stat <- f_test$out$`1 ~~ 2`$result[1]
    f_p <- f_test$out$`1 ~~ 2`$result[4]
    aic <- mean(sapply(fit$analyses, AIC))
    delta_aic <- aic - mean(sapply(base_model$analyses, AIC))
    pd_base <- (summary(pool(base_model)) %>% filter(term == "PD"))$estimate
    pd_alt <- (summary(pool(fit)) %>% filter(term == "PD"))$estimate
    R2 <- pool.r.squared(fit)[1]
  } else { # G only
    fit <- lm(update(base_formula, ~ . + G), data = df, weights = weights)
    base_model <- lm(base_formula, data = df, weights = weights)
    f_test <- anova(fit, base_model)
    f_stat <- f_test$F[2]
    f_p <- f_test$`Pr(>F)`[2]
    stats <- coef(summary(fit))
    R2 <- summary(fit)$r.squared
    aic <- AIC(fit)
    delta_aic <- aic - AIC(base_model)
    pd_base <- coef(base_model)["PD"]
    pd_alt <- coef(fit)["PD"]
  }
  if (!is.na(row$SNP) & row$SNP != "") {
    G <- stats["G", ]
  } else {
    G <- c(NA, NA, NA, NA)
  }
  if (!is.na(row$env) & row$env != "") {
    E <- stats["E", ]
  } else {
    E <- c(NA, NA, NA, NA)
  }
  if (row$model == "GxE") {
    GxE <- stats["GxE", ]
  } else {
    GxE <- c(NA, NA, NA, NA)
  }
  names(G) <- paste0("G", c("est", "se", "t", "p"))
  names(E) <- paste0("E", c("est", "se", "t", "p"))
  names(GxE) <- paste0("GxE", c("est", "se", "t", "p"))
  res <- c(row, G, E, GxE)
  res$R2 <- R2
  res$f <- f_stat
  res$f_p <- f_p
  res$aic <- aic
  res$delta_aic <- delta_aic
  res$pd_alt <- pd_alt
  res$pd_base <- pd_base
  res$pd_change <- pd_alt - pd_base
  res$pd_percent_change <- ((pd_alt - pd_base) / pd_alt) * 100
  return(as.data.table(res))
}


head(manifest)
dim(manifest)
methy <- methy[unique(manifest$cpg), on = "cpg"]
geno <- geno[unique(manifest$SNP), on = "SNP"]
geno <- geno[!duplicated(SNP)]
gc()
# Run all regressions in parallel
setDTthreads(1)
gc()
out_list <- mclapply(
  1:nrow(manifest),
  function(i) {
    fit_model(manifest[i, ], methy = methy[cpg == manifest[i, ]$cpg], geno = geno[SNP == manifest[i, ]$SNP])
  },
  mc.cores = num_cores
)
setDTthreads(num_cores)

out_path <- paste0(argv[[2]], "_top_models_pd_change.txt")
fwrite(rbindlist(out_list), out_path, sep = "\t", row.names = F, quote = F)
