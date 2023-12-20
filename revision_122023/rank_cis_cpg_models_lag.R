library(data.table)
library(parallel)
library(tidyr)
library(dplyr)
library(sfsmisc)
library(rlang)
library(tibble)
library(mice)
library(MASS)
library(glue)
# Global data
#argv <- commandArgs(trailingOnly = TRUE) # cpgs |  male or female
set.seed(Sys.time())
argv <- list("~/analysis_of_pd_dnam/sex_combined_pd_cmr_cpgs.txt", "female")
num_cores <- 64L
setDTthreads(num_cores)

probe_pos <- fread("~/pd_grs_ewas_integration/cis_mQTL_analyses/terre_data/probe_pos.txt")
snp_pos <- fread("~/pd_grs_ewas_integration/cis_mQTL_analyses/terre_data/snp_pos.txt", key = "CHR,POS")
methy <- fread(
  glue("~/pd_grs_ewas_integration/cis_mQTL_analyses/{argv[[2]]}_methy_2022.txt.gz"),
  key = "cpg"
)
probe_pos <- probe_pos[geneid %chin% methy$cpg]
geno <- fread(
  glue("~/pd_grs_ewas_integration/cis_mQTL_analyses/{argv[[2]]}_geno_2022.txt.gz"),
  key = "SNP"
)

shared_covariates <- fread(
  glue("~/pd_grs_ewas_integration/cis_mQTL_analyses/{argv[[2]]}_meta_2022_w_weights.txt.gz"),
)
shared_covariates$plate <- as.character(shared_covariates$plate)
shared_covariates$alcohol1 <- as.character(shared_covariates$alcohol1)
shared_covariates$smoking <- as.character(shared_covariates$smoking)
env_data <- fread("/home1/NEURO/SHARE_DECIPHER/TERRE_pesticides/lag_not_imputed.csv")# fread("/home1/NEURO/SHARE_DECIPHER/TERRE_pesticides/pesticides.csv")
# mapping <- fread("/home1/NEURO/SHARE_DECIPHER/terre_meta_master.csv")

# mapping$IID <- gsub("_PAE.*", "", mapping$IID)
env_data$num <- colnames(geno)[-c(1)][match(env_data$num, colnames(geno)[-c(1)])]
env_data <- env_data[!is.na(num)]
pest_missing <- read.csv("/home1/NEURO/SHARE_DECIPHER/TERRE_pesticides/lag_not_imputed.csv")#read.csv("/home1/NEURO/SHARE_DECIPHER/TERRE_pesticides/pesticides.csv")
pest_imputed <- read.csv("/home1/NEURO/SHARE_DECIPHER/TERRE_pesticides/lag_imputed.csv")#read.csv("/home1/NEURO/SHARE_DECIPHER/TERRE_pesticides/pesticides_imputed.csv")
pre_mids <- rbind(cbind(data.frame(X_imputation_ = 0), pest_missing)[, colnames(pest_imputed)], pest_imputed)
pre_mids$num <- pre_mids$NUM # NOTE ADDED FOR LAG DATA
pre_mids$num <- as.character(pre_mids$num)
pre_mids <- pre_mids %>% dplyr::select(last_col(),everything(),-NUM)
pest_lag_cols <- which(grepl("_lag",colnames(pre_mids)) & unlist(apply(pre_mids[pre_mids$num %in% colnames(geno),],2,function(x) !all(is.na(x)) & length(unique(x)) > 2)))
pre_mids <- pre_mids[pre_mids$num %in% colnames(geno), c(1,2,pest_lag_cols) ]
pre_mids <- pre_mids %>% dplyr::mutate(across(contains("lag"),~as.character(cut(.,breaks=c(0.1,median(.,na.rm=T),max(.,na.rm=T)),labels=c("Recent","Past")))))
pre_mids <- pre_mids %>% dplyr::mutate(across(contains("lag"),~replace(.,is.na(.),"None")))
pest_mids <- suppressWarnings(as.mids(pre_mids, .imp = "X_imputation_", .id = "num"))
envs <- colnames(pest_mids$data)
# Utilities
cis_snps <- function(cpg) {
  i <- which(probe_pos$geneid == cpg)
  s1 <- probe_pos$s1[i]
  s2 <- probe_pos$s2[i]
  chr <- probe_pos$chr[i]
  snp_pos[CHR == chr & (POS >= s1 - 75000 & POS <= s2)]$SNP
}

### @TODO ADD LAG TIME VARIABLE IN HERE!
generate_inputs <- function(row, methy, geno) {
  df <- data.frame(
    y = unlist(methy[row$cpg, -c(1), on = "cpg"])
  )
  covar <- shared_covariates
  if (!is.na(row$env)) {
    env <- row$env
    tmp_df <- complete(pest_mids, action = "long", include = TRUE)[, c(".imp", ".id", env)]
    tmp_df$EPast <- as.numeric(tmp_df[,c(env)] == "Past")
    tmp_df$ERecent <- as.numeric(tmp_df[,c(env)] == "Recent")
    tmp_df[,c(env)] <- NULL
    # colnames(tmp_df) <- c(".imp", ".id","EPast","ERecent")
    ix <- match(pre_mids$num, colnames(geno)[-c(1)])
    tmp_df <- cbind(tmp_df, covar[ix, ], y = df$y[ix])
    if (grepl("G", row$model)) {
      G <- unlist(geno[row$SNP, -c(1), on = "SNP"])
      tmp_df$G <- G[ix]
      maf <- sum(tmp_df$G) / (2 * nrow(tmp_df))
      if (maf > 0.95 | maf < 0.05) {
        return(data.table())
      }
      exposures <- sum(tmp_df$EPast + tmp_df$ERecent, na.rm = T) / nrow(tmp_df)
      if (exposures > 0.95 | exposures < 0.05) {
        return(data.table())
      }
      if (row$model == "GxE") { # GxE
        tmp_df$GxEPast <- tmp_df$G * tmp_df$EPast
        tmp_df$GxERecent <- tmp_df$G * tmp_df$ERecent
        df_mids <- suppressWarnings(as.mids(tmp_df))
      } else { # G+E
        df_mids <- suppressWarnings(as.mids(tmp_df))
      }
    } else { # E only
      exposures <- sum(tmp_df$EPast + tmp_df$ERecent, na.rm = T) / nrow(tmp_df)
      if (exposures > 0.95 | exposures < 0.05) {
        return(data.table())
      }
      df_mids <- suppressWarnings(as.mids(tmp_df))
    }
    return(df_mids)
  } else { # G only
    ix <- match(covar, colnames(geno)[-c(1)])
    df$G <- unlist(geno[row$SNP, -c(1), on = "SNP"])
    df <- na.omit(cbind(df, covar))
    maf <- sum(df$G) / (2 * nrow(df))
    if (maf > 0.95 | maf < 0.05) {
      return(data.table())
    }
    return(df)
  }
}


fit_model <- function(rows, methy, geno) { # row of snp cpg model env data.frame
  base_formula <- y ~ PD + age + smoking + alcohol1 + head_trauma + CTP_PC1 + CTP_PC2 + CTP_PC3 + CTP_PC4 + CTP_PC5 + CTP_PC6 + ancestryPC1 + ancestryPC2 + ancestryPC3 + plate+  SentrixPosition_A
  output <- lapply(
    1:nrow(rows),
    function(i) {
      tryCatch({
          row <- rows[i, ]
          set.seed(NULL)
          df <- generate_inputs(row, methy, geno)
          if (!class(df) == "mids") {
            return(data.table())
          }
          if (!is.na(row$env)) {
            if (grepl("G", row$model)) {
              if (row$model == "GxE") { # GxE
                fit <- with(df, rlm(formula(format(update(base_formula, ~ . + G + EPast + ERecent+ GxEPast + GxERecent))), weights = weights, psi = psi.huber,maxit = 100))
              } else { # G+E
                fit <- with(df, rlm(formula(format(update(base_formula, ~ . + G + EPast + ERecent))), weights = weights, psi = psi.huber,maxit = 100))
              }
              base_model <- with(df, rlm(formula(format(base_formula)), weights = weights, psi = psi.huber,maxit = 100))
            } else { # E only
              base_model <- with(df, rlm(formula(format(base_formula)), weights = weights, psi = psi.huber,maxit = 100))
              fit <- with(df, rlm(formula(format(update(base_formula, ~ . + EPast + ERecent))), weights = weights, psi = psi.huber,maxit = 100))
            }
            dfcom <- nrow(shared_covariates)- length(fit$analyses[[1]]$coefficients - 1)  # degrees of freedom
            stats <- column_to_rownames(summary(pool(fit,dfcom = dfcom)) %>% dplyr::select(-df), var = "term")
            f_test <- D1(fit, base_model,dfcom = dfcom)$result 
            f_stat <- f_test[1]
            f_p <- f_test[4]
            aic <- mean(sapply(fit$analyses, AIC))
            delta_aic <- aic - mean(sapply(base_model$analyses, AIC))
            pd_base <- (summary(pool(base_model)) %>% filter(term == "PD"))$estimate
            pd_alt <- (summary(pool(fit)) %>% filter(term == "PD"))$estimate
            #R2 <- pool.r.squared(fit)[1]
          } else { # G only
            fit <- rlm(update(base_formula, ~ . + G), data = df, weights = weights, psi = psi.huber,maxit = 100)
            base_model <- rlm(base_formula, data = df, weights = weights, psi = psi.huber,maxit = 100)
            f_test <- f.robftest(fit, "G")
            f_stat <- f_test$statistic
            f_p <- f_test$p.value
            stats <- coef(summary(fit))
            stats <- cbind(stats,p=NA)
            stats["G","p"] <- f_p
            #R2 <- summary(fit)$r.squared
            aic <- AIC(fit)
            delta_aic <- aic - AIC(base_model)
            pd_base <- coef(base_model)["PD"]
            pd_alt <- coef(fit)["PD"]
          }
          if (!is.na(row$SNP)) {
            G <- stats["G", ]
          } else {
            G <- c(NA, NA, NA, NA)
          }
          if (!is.na(row$env)) {
            EPast <- stats["EPast", ]
            ERecent <- stats["ERecent", ]
          } else {
            EPast <- c(NA, NA, NA, NA)
            ERecent <- c(NA, NA, NA, NA)
          }
          if (row$model == "GxE") {
            GxEPast <- stats["GxEPast", ]
            GxERecent <- stats["GxERecent", ]
          } else {
            GxEPast <- c(NA, NA, NA, NA)
            GxERecent <- c(NA, NA, NA, NA)
          }
          names(G) <- paste0("G", c("est", "se", "t", "p"))
          names(EPast) <- paste0("EPast", c("est", "se", "t", "p"))
          names(ERecent) <- paste0("ERecent", c("est", "se", "t", "p"))
          names(GxEPast) <- paste0("GxEPast", c("est", "se", "t", "p"))
          names(GxERecent) <- paste0("GxERecent", c("est", "se", "t", "p"))
          res <- c(row, G, EPast,ERecent, GxEPast,GxERecent)
         # res$R2 <- R2
          res$f <- f_stat
          res$f_p <- f_p
          res$aic <- aic
          res$delta_aic <- delta_aic
          res$pd_alt <- pd_alt
          res$pd_base <- pd_base
          res$pd_change <- pd_alt - pd_base
          res$pd_percent_change <- ((pd_alt - pd_base) / pd_alt) * 100
          return(as.data.table(res))
        },
        error = function(e){
          print(e)
          data.table()
          }
        )
      }
    )
  return(rbindlist(output))
}

manifest_from_cpg <- function(cpg) {
  snps <- cis_snps(cpg)
  if (length(snps) > 0) {
    g_models <- tibble(cpg, SNP = snps, env = NA, model = "G")
    e_models <- tibble(cpg, SNP = NA, env = envs, model = "E")
    ge_models <- expand_grid(cpg, SNP = snps, env = envs, model = c("G+E", "GxE"))
    return(data.table(dplyr::bind_rows(g_models, e_models, ge_models)))
  } else {
    return(data.table())
  }
}
probes <- fread(argv[[1]])$cpg


# Generate manifest of all regressions
manifest_list <- mclapply(probes, manifest_from_cpg, mc.cores = num_cores)
manifest <- rbindlist(manifest_list)
head(manifest)
dim(manifest)
methy <- methy[unique(probes), on = "cpg"]
geno <- geno[unique(manifest$SNP), on = "SNP"]
geno <- geno[!duplicated(SNP)]
gc()
# Run all regressions in parallel
for (mod in c("GxE", "G+E","E", "G")) {
  setDTthreads(1)
  out_path <- paste0(argv[[2]], "_", mod, "_lag_dnam_breakdown.txt")
  if(file.exists(out_path)){
    print(glue("{out_path} already exists, skipping..."))
    next
  }
  manifest_cur <- manifest[model == mod][1:100]
  out_list <- mclapply(
    splitIndices(nrow(manifest_cur), num_cores),
    function(i) {
      fit_model(manifest_cur[i, ], methy = methy[cpg %in% manifest_cur[i, ]$cpg], geno = geno[SNP %in% manifest_cur[i, ]$SNP])
    },
    mc.cores = num_cores
  )
  print(head(out_list[[1]]))
  setDTthreads(num_cores)
  fwrite(rbindlist(out_list), out_path, sep = "\t", row.names = F, quote = F)
}
