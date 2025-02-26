#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script for null modelling (beta NTI, Raup-Crick Beta-Diversity, and Elements of Metacommunity)
#Reference: http://uu.diva-portal.org/smash/get/diva2:1373632/DATASET09.txt
#v1.0 (All metrics are now being saved)
#v1.1 (Found a serious bug as Raup-Crick Beta-Diversity can be used in both incidence and Bray-Curtis model)

# library(phyloseq)
# library(vegan)
# library(ape)
# library(picante)
# library(ecodist)
# library(metacom)

qpe <- function(physeq) {
  # TODO: remove
  grouping_column <- "Groups"

  # Process input data
  abund_table <- phyloseq::otu_table(physeq)
  meta_table <- phyloseq::sample_data(physeq)
  OTU_taxonomy <- phyloseq::tax_table(physeq)
  OTU_tree <- phyloseq::phy_tree(physeq)

  #Convert the data to phyloseq format
  OTU = otu_table(as.matrix(abund_table), taxa_are_rows = FALSE)
  TAX = tax_table(as.matrix(OTU_taxonomy))
  SAM = sample_data(meta_table)

  physeq<-merge_phyloseq(phyloseq(OTU, TAX),SAM,OTU_tree)

  #Pruning and subsampling
  physeq<-prune_taxa(taxa_sums(physeq)>10, physeq)

  #Rarefy to minimum sample size
  physeq_rel = rarefy_even_depth(physeq, sample.size = min(sample_sums(physeq)))

  summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
      if (na.rm)
        sum(!is.na(x))
      else
        length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                     c(N    = length2(xx[[col]], na.rm=na.rm),
                       mean = mean   (xx[[col]], na.rm=na.rm),
                       sd   = sd     (xx[[col]], na.rm=na.rm)
                     )
                   },
                   measurevar
    )

    # Rename the "mean" column
    datac <- rename(datac, measurevar = "mean")

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
  }


  abund_table_ems<-otu_table(physeq_rel)
  meta_table_ems<-sample_data(physeq_rel)

  collated_coherence<-NULL
  collated_turnover<-NULL
  collated_boundary<-NULL
  collated_sitescores<-NULL
  collated_pairwise_RC_abundance<-NULL
  collated_pairwise_RC_incidence<-NULL
  collated_pairwise_bNTI<-NULL

  for(i in 1:nlevels(meta_table_ems$Groups)){
    abund_table_ems_group<-abund_table_ems[meta_table_ems$Groups==levels(meta_table_ems$Groups)[i],,drop=F]
    abund_table_ems_group = abund_table_ems_group[,which(colSums(abund_table_ems_group) != 0)]

    #Now calculate Elements of Metacommunity
    met_ems_group=Metacommunity( abund_table_ems_group, scores = 1, method = "r1", sims = 999, order = T, binary = F, verbose=T, allowEmpty=T)
    om_ems_group=OrderMatrix(abund_table_ems_group, outputScores=T, binary=F)

    coherence<-t(data.frame(row.names=met_ems_group$Coherence[-nrow(met_ems_group$Coherence),1],stats=met_ems_group$Coherence[-nrow(met_ems_group$Coherence),2]))
    rownames(coherence)<-levels(meta_table$Groups)[i]

    turnover<-t(data.frame(row.names=met_ems_group$Turnover[-nrow(met_ems_group$Turnover),1],stats=met_ems_group$Turnover[-nrow(met_ems_group$Turnover),2]))
    rownames(turnover)<-levels(meta_table$Groups)[i]

    boundary<-t(data.frame(row.names=met_ems_group$Boundary[,1],stats=met_ems_group$Boundary[,2]))
    rownames(boundary)<-levels(meta_table$Groups)[i]

    sitescores<-data.frame(om_ems_group$sitescores)
    colnames(sitescores)<-c("sitescores")
    sitescores$Groups<-levels(meta_table_ems$Groups)[i]

    #Now calculate Raup-Crick abundance based
    results=raup_crick_abundance(abund_table_ems_group, set_all_species_equal = F, plot_names_in_col1 = F, reps=999)
    pairwise_RC_abundance<-reshape2::melt(as.matrix(results))
    pairwise_RC_abundance$Groups<-levels(meta_table_ems$Groups)[i]


    #Now calculate Raup-Crick incidence based
    results=raup_crick_incidence(abund_table_ems_group,plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
    pairwise_RC_incidence<-reshape2::melt(as.matrix(results))
    pairwise_RC_incidence$Groups<-levels(meta_table_ems$Groups)[i]


    #Now calculate betaMNTD
    m=match.phylo.data(OTU_tree, t(abund_table_ems_group)) #Extract subtree for each group
    OTU_tree_ems_group=m$phy

    abund_table_ems_group<-t(abund_table_ems_group)

    beta.mntd.weighted = as.matrix(comdistnt(t(abund_table_ems_group),cophenetic(OTU_tree_ems_group),abundance.weighted=T));


    beta.reps = 999;
    rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(abund_table_ems_group),ncol(abund_table_ems_group),beta.reps));
    dim(rand.weighted.bMNTD.comp);

    for (rep in 1:beta.reps) {

      rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(abund_table_ems_group),taxaShuffle(cophenetic(OTU_tree_ems_group)),abundance.weighted=T,exclude.conspecifics = F));

      print(c(date(),rep));

    }

    weighted.bNTI = matrix(c(NA),nrow=ncol(abund_table_ems_group),ncol=ncol(abund_table_ems_group));
    dim(weighted.bNTI);

    for (columns in 1:(ncol(abund_table_ems_group)-1)) {
      for (rows in (columns+1):ncol(abund_table_ems_group)) {

        rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
        weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
        rm("rand.vals");

      };
    };
    rownames(weighted.bNTI) = colnames(abund_table_ems_group);
    colnames(weighted.bNTI) = colnames(abund_table_ems_group);
    pairwise_bNTI<-reshape2::melt(as.matrix(weighted.bNTI))
    pairwise_bNTI$Groups<-levels(meta_table_ems$Groups)[i]

    #Now collate all the statistics together
    if(is.null(collated_coherence)){collated_coherence<-coherence} else {collated_coherence<-rbind(collated_coherence,coherence)}
    if(is.null(collated_boundary)){collated_boundary<-boundary} else {collated_boundary<-rbind(collated_boundary,boundary)}
    if(is.null(collated_turnover)){collated_turnover<-turnover} else {collated_turnover<-rbind(collated_turnover,turnover)}
    if(is.null(collated_sitescores)){collated_sitescores<-sitescores} else {collated_sitescores<-rbind(collated_sitescores,sitescores)}
    if(is.null(collated_pairwise_RC_abundance)){collated_pairwise_RC_abundance<-pairwise_RC_abundance} else {collated_pairwise_RC_abundance<-rbind(collated_pairwise_RC_abundance,pairwise_RC_abundance)}
    if(is.null(collated_pairwise_RC_incidence)){collated_pairwise_RC_incidence<-pairwise_RC_incidence} else {collated_pairwise_RC_incidence<-rbind(collated_pairwise_RC_incidence,pairwise_RC_incidence)}
    if(is.null(collated_pairwise_bNTI)){collated_pairwise_bNTI<-pairwise_bNTI} else {collated_pairwise_bNTI<-rbind(collated_pairwise_bNTI,pairwise_bNTI)}


  }

  collated_RC<-summarySE(collated_pairwise_RC_incidence, measurevar = "value", groupvars = "Groups")
  rownames(collated_RC)<-collated_RC[,1]
  collated_RC<-collated_RC[,-1]
  collated_bNTI<-summarySE(collated_pairwise_bNTI, measurevar = "value", groupvars = "Groups")
  rownames(collated_bNTI)<-collated_bNTI[,1]
  collated_bNTI<-collated_bNTI[,-1]

  return(structure(list(
    coherence = collated_coherence,
    boundary = collated_boundary,
    turnover = collated_turnover,
    sitescores = collated_sitescores,
    pairwise_RC = collated_pairwise_RC_abundance,
    pairwise_bNTI = collated_pairwise_bNTI,
    RC = collated_RC,
    bNTI = collated_bNTI
  ), className = "ECQuantitativeProcessEstimate"))
}

# write.csv(collated_coherence,file=paste("Coherence_",label,".csv",sep=""))
# write.csv(collated_boundary,file=paste("Boundary_",label,".csv",sep=""))
# write.csv(collated_turnover,file=paste("Turnover_",label,".csv",sep=""))
# write.csv(collated_sitescores,file=paste("Sitescores_",label,".csv",sep=""))
# write.csv(collated_pairwise_RC_abundance,file=paste("PairwiseRC_",label,".csv",sep=""))
# write.csv(collated_pairwise_bNTI,file=paste("PairwisebNTI_",label,".csv",sep=""))
# write.csv(collated_RC,file=paste("RC_",label,".csv",sep=""))
# write.csv(collated_bNTI,file=paste("bNTI_",label,".csv",sep=""))
