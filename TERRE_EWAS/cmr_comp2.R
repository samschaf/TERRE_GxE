cmr_comp2 <- function (cmrs, Mdata, cmethod = c("pca", "mean", "median"), 
          verbose = T) 
{
  cmethod <- match.arg(cmethod)
  nrMm = nrow(Mdata)
  ncMm = ncol(Mdata)
  nrTst = min(1000, nrMm)
  ncTst = min(1000, ncMm)
  if ((nrMm > ncMm) && (length(grep("cg", colnames(Mdata)[1:ncTst])) < 
                        1) && (length(grep("cg", rownames(Mdata)[1:nrTst])) > 
                               0)) {
    if (verbose) 
      print("Note: Probes should be in colums, so slowing down now to gently transpose the data for you.\nPlease be patient and if this operation is not over by the time you finish reading this, then it would appear that your computer is rather slow, so you may want to go for a stroll, or just chill and comeback later.")
    Mdata = t(Mdata)
  }
  if (is.null(cmrs)) {
    if (verbose) 
      print("CMRs not provided, estimating now. This will take awhile..")
    cmrs = unlist(cmr(Mdata, verbose = verbose), recursive = F)
  }
  else {
    #if (length(cmrs) < 24) 
      if (length(cmrs[[1]][[1]]) != 1) 
        cmrs = unlist(cmrs, recursive = F)
  }
  ncmr = length(cmrs)
  cmr_frstPrb_nam = sapply(cmrs, function(x) {
    x[[1]]
  })
  if (ncmr > 1000) 
    messs = " CMRs, this may take awhile."
  else messs = " CMRs."
  if (verbose) {
    print(paste0("Found ", ncmr, messs))
    print("Starting CMRs composite estimation.")
  }
  if (cmethod == "pca") {
    res = as.matrix(do.call(cbind, lapply(cmrs, function(x) {
      #print(x) #added in for troubleshooting
      pc1 = prcomp(Mdata[, x], retx = FALSE, center = TRUE, 
                   scale. = FALSE, rank. = 1)
      pcld = pc1$rotation
      return(abs(Mdata[, x] %*% pcld/sum(pcld)))
    })))
  }
  else if (cmethod == "mean") {
    res = as.matrix(do.call(cbind, lapply(cmrs, function(x) {
      pc1 = prcomp(Mdata[, x], retx = FALSE, center = TRUE, 
                   scale. = FALSE, rank. = 1)
      return(apply(Mdata[, x], 1, mean))
    })))
  }
  else {
    res = as.matrix(do.call(cbind, lapply(cmrs, function(x) {
      #print(x) #troubleshooting
      #changed below line to keep complete cases of Mdata (was giving errors with NAs)
      #this removes individuals with NAs for all the probes within a CMR
      pc1 = prcomp(Mdata[, x][complete.cases(Mdata[,x]),], retx = FALSE, center = TRUE, 
                   scale. = FALSE, rank. = 1)
      return(apply(Mdata[, x], 1, median))
    })))
  }
  rownames(res) = rownames(Mdata)
  colnames(res) = cmr_frstPrb_nam
  if (verbose) 
    print("Finished.")
  return(res)
}