cmrplot <- function(sig.cpgs, metadata, fdat, betas, cmr, path, sex_strat=TRUE, verbose=TRUE, width=480, height=480){

  require(reshape2)
  require(ggplot2)
  require(gridExtra)
  require(dplyr)
  
  #get the probe names of the cmr for plotting
  cmrs <- unique(sig.cpgs[sig.cpgs$threshold_adjDB=="TRUE","cmrGR"])
  if (verbose==TRUE){print(as.character(cmrs[cmr]))} #print cmr coordinates
  chr <- paste("chr", sig.cpgs[sig.cpgs$cmrGR==cmrs[cmr],"CHR"], sep="")
  cmr_probes <- unlist(strsplit(as.character(sig.cpgs[sig.cpgs$cmrGR==cmrs[cmr],"cmr"]), split=","))

  #get annotation data
  #fdat_sub <- fdat[fdat$TargetID %in% cmr_probes,]
  gene <- unique(unlist(strsplit(as.character(sig.cpgs[sig.cpgs$cmrGR==cmrs[cmr],"geneSymbol"]),split=",")))
  group <- unique(unlist(strsplit(as.character(sig.cpgs[sig.cpgs$cmrGR==cmrs[cmr],"UCSC_REFGENE_GROUP"]),split=";")))
  island <- unique(unlist(strsplit(as.character(sig.cpgs[sig.cpgs$cmrGR==cmrs[cmr],"RELATION_TO_UCSC_CPG_ISLAND"]),split=";")))
  lab <- paste("CMR",cmr,gene,group,island,sep="_")
  if (verbose==TRUE){ print(lab) } #print gene and context

  #get betas
  betas_sub <- as.data.frame(betas[rownames(betas) %in% cmr_probes,])
  betas_sub$cpg <- rownames(betas_sub)
  betas_sub <- reshape2::melt(betas_sub)

  #add coordinates
  if (verbose==TRUE){print("Annotating coordinates")}
  betas_sub$coord <- NA
  for (i in 1:nrow(betas_sub)){
    betas_sub$coord[i] <- fdat[fdat$TargetID==betas_sub$cpg[i],"MAPINFO"]
  }  

  #add sex and PD information
  betas_sub$PD <- NA
  betas_sub$Sex <- NA
  betas_sub$sex_PD <- NA

  metadata$sex_PD <- NULL
  metadata[metadata$reportedSex=="F" & metadata$PD==0,"sex_PD"] <- "F_ctrl"
  metadata[metadata$reportedSex=="F" & metadata$PD==1,"sex_PD"] <- "F_case"
  metadata[metadata$reportedSex=="M" & metadata$PD==0,"sex_PD"] <- "M_ctrl"
  metadata[metadata$reportedSex=="M" & metadata$PD==1,"sex_PD"] <- "M_case"

  metadata$patient <- as.character(metadata$patient)
  betas_sub$variable <- as.character(betas_sub$variable)
  
  for (i in 1:nrow(betas_sub)){
    #print(i)
    betas_sub$PD[i] <- as.character(metadata[metadata$patient==betas_sub$variable[i],"PD"])
    betas_sub$Sex[i] <- as.character(metadata[metadata$patient==betas_sub$variable[i],"reportedSex"])
    betas_sub$sex_PD[i] <- metadata[metadata$patient==betas_sub$variable[i],"sex_PD"]
  }

  if (sex_strat==TRUE){
    #calculate group means for drawing lines (sex stratified)
    F_ctrl_means <- betas_sub[betas_sub$sex_PD=="F_ctrl",] %>%
      group_by(coord) %>%
      summarise(mean = mean(value, na.rm=TRUE))
    F_ctrl_means$sex_PD <- "F_ctrl"
    F_case_means <- betas_sub[betas_sub$sex_PD=="F_case",] %>%
      group_by(coord) %>%
      summarise(mean = mean(value, na.rm=TRUE))
    F_case_means$sex_PD <- "F_case"
    M_ctrl_means <- betas_sub[betas_sub$sex_PD=="M_ctrl",] %>%
      group_by(coord) %>%
      summarise(mean = mean(value, na.rm=TRUE))
    M_ctrl_means$sex_PD <- "M_ctrl"
    M_case_means <- betas_sub[betas_sub$sex_PD=="M_case",] %>%
      group_by(coord) %>%
      summarise(mean = mean(value, na.rm=TRUE))
    M_case_means$sex_PD <- "M_case"
    means <- rbind(F_ctrl_means, F_case_means, M_ctrl_means, M_case_means)
    colnames(means)[2] <- "value" }

#calculate group means for drawing lines (sexes combined)
  ctrl_means <- betas_sub[betas_sub$PD==0,] %>%
    group_by(coord) %>%
    summarise(mean = mean(value, na.rm=TRUE))
  ctrl_means$PD <- "Control"
  case_means <- betas_sub[betas_sub$PD==1,] %>%
    group_by(coord) %>%
    summarise(mean = mean(value, na.rm=TRUE))
  case_means$PD <- "Case"
  all.equal(ctrl_means$coord, case_means$coord)
  means_combined <- rbind(ctrl_means, case_means)
  colnames(means_combined)[2] <- "value"

  #make categorical variables factors
  betas_sub$PD <- gsub(1, "Case", gsub(0, "Control", betas_sub$PD))
  betas_sub$PD <- as.factor(betas_sub$PD)
  betas_sub$Sex <- as.factor(betas_sub$Sex)
  means_combined$PD <- as.factor(means_combined$PD)

  all <- ggplot(betas_sub, aes(as.numeric(coord), value)) +
    geom_point(aes(colour=PD, group=PD), data=betas_sub) +
    scale_y_continuous("Beta value",limits=c(0,1)) + theme_classic() + geom_line(data=means_combined, aes(group=PD, colour=PD)) + scale_colour_manual(values=c("black","grey48")) + scale_x_continuous(paste(chr, "coordinate", sep=" ")) + ggtitle("All individuals") 

  if (sex_strat==FALSE){ 
    
    png(paste(path, lab, ".png", sep="", width=width, height=height))
    return(all)
    dev.off()
    
    } else{
  
  f <- ggplot(betas_sub[betas_sub$Sex==0,], aes(as.numeric(coord), value)) +
  geom_point(aes(colour=sex_PD, group=sex_PD), data=betas_sub[betas_sub$Sex=="F",]) +
  scale_y_continuous("Beta value", limits=c(0,1)) + theme_classic() + geom_line(data=means[means$sex_PD %in% c("F_ctrl", "F_case"),], aes(group=sex_PD, colour=sex_PD)) + scale_colour_manual(values=c("pink4","pink2"), labels=c("Case", "Control"), name="PD") + scale_x_continuous(paste(chr, "coordinate", sep=" ")) + ggtitle("Females") 

  m <- ggplot(betas_sub[betas_sub$Sex==1,], aes(as.numeric(coord), value)) +
  geom_point(aes(colour=sex_PD, group=sex_PD), data=betas_sub[betas_sub$Sex=="M",]) +
  scale_y_continuous("Beta value", limits=c(0,1)) + theme_classic() + geom_line(data=means[means$sex_PD %in% c("M_ctrl", "M_case"),], aes(group=sex_PD, colour=sex_PD)) + scale_colour_manual(values=c("lightsteelblue4","lightsteelblue2"), labels=c("Case", "Control"), name="PD") + scale_x_continuous(paste(chr, "coordinate", sep=" ")) + ggtitle("Males")

  png(paste(path, lab, ".png", sep=""), width=width, height=height)
  grid.arrange(all, f, m, nrow=3)
  dev.off()
  
  }}