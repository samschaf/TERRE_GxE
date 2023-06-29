library(data.table)
library(parallel)
library(tidyr)
library(dplyr)
library(rlang)
library(tibble)

# Global data
num_cores <- 8
setDTthreads(num_cores)
probe_pos <- fread("~/prs_ewas_integration/cis_mQTL_analyses/digpd_data/probe_pos_CTP.txt")
snp_pos <- fread("~/prs_ewas_integration/cis_mQTL_analyses/digpd_data/snp_pos_CTP.txt", key = "CHR,POS")
methy <- fread(
  "~/prs_ewas_integration/cis_mQTL_analyses/digpd_data/methylation_combat_CTP.txt",
  key = "cpg"
)
probe_pos <- probe_pos[geneid %chin% methy$cpg]
geno <- fread(
  "~/prs_ewas_integration/cis_mQTL_analyses/digpd_data/all_imputed_matrixeQTL_CTP.txt",
  key = "SNP",
  fill=TRUE
)
geno <- na.omit(geno)

shared_covariates <- data.table::transpose(
  na.omit(
    fread(
      "~/prs_ewas_integration/cis_mQTL_analyses/digpd_data/covariates_CTP.txt",
      fill = TRUE
    )
  ),
  make.names = "id"
)

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
        tmp_df <- tmp_df[tmp_df$sex == 1, ]
      } else {
        tmp_df <- tmp_df[tmp_df$sex == 0, ]
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
        tmp_df <- tmp_df[tmp_df$sex == 1, ]
      } else {
        tmp_df <- tmp_df[tmp_df$sex == 0, ]
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
      df <- df[df$sex == 1, ]
    } else {
      df <- df[df$sex == 0, ]
    }
    df$sex <- NULL
    maf <- sum(df$G) / (2 * nrow(df))
    if (maf > 0.99 | maf < 0.01) {
      return(data.table())
    }
    return(df)
  }
}


fit_model <- function(rows) { # row of snp cpg model env data.frame
  output <- lapply(1:nrow(rows),
    function(i){
      row <- rows[i,]
      set.seed(NULL)
      df <- generate_inputs(row,argv[[2]])
      if (is_empty(df)) {
        return(data.table())
      }
      if (!is.na(row$env)) {
        if (grepl("G", row$model)) {
          if (row$model == "GxE") { # GxE
            fit <- with(df, lm(y ~ PD + G + E + GxE + V3 + V4 + V5 + age + head_trauma + CTP_PC1 + CTP_PC2 + CTP_PC3 + CTP_PC4 + CTP_PC5))
          } else { # G+E
            fit <- with(df, lm(y ~ PD + G + E + V3 + V4 + V5 + age + head_trauma + CTP_PC1 + CTP_PC2 + CTP_PC3 + CTP_PC4 + CTP_PC5))
          }
          base_model <- with(df, lm(y ~ PD + V3 + V4 + V5 + age + head_trauma + CTP_PC1 + CTP_PC2 + CTP_PC3 + CTP_PC4 + CTP_PC5))
        } else { # E only
          base_model <- with(df, lm(y ~ PD + V3 + V4 + V5 + age + head_trauma + CTP_PC1 + CTP_PC2 + CTP_PC3 + CTP_PC4 + CTP_PC5))
          fit <- with(df, lm(y ~ PD + E + V3 + V4 + V5 + age + head_trauma + CTP_PC1 + CTP_PC2 + CTP_PC3 + CTP_PC4 + CTP_PC5))
        }
        stats <- column_to_rownames(summary(pool(fit)) %>% select(-df), var = "term")
        f_test <- anova(fit, base_model)
        f_stat <- f_test$out$`1 ~~ 2`$result[1]
        f_p <- f_test$out$`1 ~~ 2`$result[4]
        aic <- mean(sapply(fit$analyses, AIC))
        R2 <- pool.r.squared(fit)[1]
      } else { # G only
        fit <- lm(y ~ ., data = df)
        base_model <- lm(y ~ . - G, data = df)
        f_test <- anova(fit, base_model)
        f_stat <- f_test$F[2]
        f_p <- f_test$`Pr(>F)`[2]
        stats <- coef(summary(fit))
        R2 <- summary(fit)$r.squared
        aic <- AIC(fit)
      }
      if (!is.na(row$SNP)) {
        G <- stats["G", ]
      } else {
        G <- c(NA, NA, NA, NA)
      }
      if (!is.na(row$env)) {
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
      return(as.data.table(res))
    }
  )
  return(rbindlist(output))
}

manifest_from_cpg <- function(cpg) {
  snps <- cis_snps(cpg)
  if (length(snps) > 0) {
    g_models <- tibble(cpg, SNP = snps, env = NA, model = "G")
  #  e_models <- tibble(cpg, SNP = NA, env = envs, model = "E")
  #  ge_models <- expand_grid(cpg, SNP = snps, env = envs, model = c("G+E", "GxE"))
    return(data.table(g_models))
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
manifest <- manifest[SNP %in% geno$SNP]
methy <- methy[unique(probes), on = "cpg"]
geno <- geno[unique(manifest$SNP), on = "SNP"]
gc()
# Run all regressions in parallel
setDTthreads(1)
for (mod in c("G")){#, "G+E", "GxE", "E")) {
  out_path <- paste0(argv[[2]], "_", mod, "_dnam_breakdown.txt")
  manifest_cur <- manifest[model == mod]
  gc()
  out_list <- mclapply(
      splitIndices(nrow(manifest_cur),num_cores),
      function(i) {
        fit_model(manifest_cur[i,])
      },
      mc.cores = num_cores
    )
    setDTthreads(num_cores)
    fwrite(rbindlist(out_list), out_path, sep = "\t", row.names = F, quote = F)
}
