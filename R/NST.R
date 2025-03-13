#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script for null modelling (ST/NST/MST)
#Authors: Umer, Anna
#v1.0 (All metrics are now being saved)

# library(phyloseq)
# library(vegan)
# library(ape)
# library(NST)
# library(ggplot2)

nst_distance_measure = list("manhattan", "mManhattan", "euclidean", "mEudclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "mGower", "morisita", "horn", "binomial", "chao", "cao")

nst <- function(
  abund_table,
  meta_table,
  number_of_randomizations = 1000,
  distance_measure= "cao", #"manhattan" "mManhattan" "euclidean" "mEuclidean"  "canberra" "bray" "kulczynski" "jaccard" "gower" "altGower" "mGower" "morisita" "horn" "binomial" "chao" "cao"
  abundance.weighted = FALSE, #Logic, consider abundances or not (just presence/absence). default is TRUE.
  #Jaccard with abundance.weighted=TRUE is called Ruzicka
  null_model = "PF", #"EE" "EP" "EF" "PE" "PP" "PF" "FE" "FP" "FF" with details given below:
  #Abbreviation|Ways_to_constrain_taxa_occurrence_frequency|Ways_to_constrain_richness_in_each_sample
  #============|===========================================|=========================================
  # EE           Equiprobable                                Equiprobable
  # EP           Equiprobable                                Proportional
  # EF           Equiprobable                                Fixed
  # PE           Proportional                                Equiprobable
  # PP           Proportional                                Proportional
  # PF           Proportional                                Fixed
  # FE           Fixed                                       Equiprobable
  # FP           Fixed                                       Proportional
  # FF           Fixed                                       Fixed

  # As to occurrence frequency:
  #   "Equiprobable" means that all taxa have equal probability to occur;
  #   "Proportional" means that the occurrence probability of a taxon is proportional to its observed occurrence frequency;
  #   "Fixed" means that the occurrence frequency of a taxon is fixed as observed.
  # As to species richness in each sample:
  #   "Equiprobable" means that all samples have equal probability to contain a taxon;
  #   "Proportional" means the occurrence probability in a sample is proportional to the observed richness in this sample;
  #   "Fixed" means the occurrence frequency of a taxon is fixed as observed
  SES = TRUE, #Logic, whether to calculate standardized effect size, which is (observed dissimilarity - mean of null dissimilarity)/standard deviation of null dissimilarity. default is FALSE.
  RC = FALSE # Logic, whether to calculate modified Raup-Crick metric, which is percentage of null dissimilarity lower than observed dissimilarity x 2 - 1. default is FALSE.
) {
  #Bug in tNST (abund_table should be of type "matrix")

  comm <- as(abund_table, "matrix")
  group <- as(meta_table[,"Groups",drop=F], "matrix")

  tnst = NST::tNST(
    comm = comm,
    group = group,
    rand = number_of_randomizations,
    dist.method = distance_measure,
    null.model = null_model,
    output.rand = TRUE,
    nworker = parallel::detectCores() - 1,
    SES = SES,
    RC = RC
  )

  #Extract mean NST/ST/MST values of groups
  df<-NULL
  measure_name<-names(tnst$index.grp)[grepl("^NST",names(tnst$index.grp))]
  tmp<-tnst$index.grp[,c("group",measure_name)]
  tmp$measure="NST"
  names(tmp)<-c("Groups","value","measure")
  df<-tmp
  measure_name<-names(tnst$index.grp)[grepl("^ST",names(tnst$index.grp))]
  tmp<-tnst$index.grp[,c("group",measure_name)]
  tmp$measure="ST"
  names(tmp)<-c("Groups","value","measure")
  df<-rbind(df,tmp)
  measure_name<-names(tnst$index.grp)[grepl("^MST",names(tnst$index.grp))]
  tmp<-tnst$index.grp[,c("group",measure_name)]
  tmp$measure="MST"
  names(tmp)<-c("Groups","value","measure")
  df<-rbind(df,tmp)

  return(structure(list(df = df), className="ECStochasticityRatio"))
}

plot.ECStochasticityRatio <- function(value, measures = c("NST", "ST", "MST")) {
  # NST_width=4
  # NST_height=8
  
  value <- value[value$measure %in% measures,]

  p <- ggplot(value, aes(x = Groups, y = value, fill = Groups)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%0.2f", round(value, digits = 2))), vjust = -0.3, size = 3.5) +
    facet_wrap(~ measure, scales = "free_y") +
    ylim(0, 1.1) +
    ylab("Stochasticity Ratios (scaled to 1)") +
    theme(panel.spacing = unit(2, "lines"))

  return(p)
}

save_data <- function() {
  #Perform PANOVA
  nst.pova=nst.panova(nst.result=tnst, rand=number_of_randomizations)
  #Convert group numbers to actual names
  for(i in 1:nlevels(meta_table$Groups)){
    nst.pova$group1<-gsub(paste("^",i,"$",sep=""),levels(meta_table$Groups)[i],nst.pova$group1)
    nst.pova$group2<-gsub(paste("^",i,"$",sep=""),levels(meta_table$Groups)[i],nst.pova$group2)
  }

  write.csv(nst.pova,file=paste("Stochasticity-Ratios_",null_model,"_",distance_measure,"_",as.character(number_of_randomizations),"_",as.character(SES),"_",as.character(RC),"_",as.character(abundance.weighted),"_",label,"_PANOVA",".csv",sep=""))
  write.csv(tnst$index.pair.grp,file=paste("Stochasticity-Ratios_",null_model,"_",distance_measure,"_",as.character(number_of_randomizations),"_",as.character(SES),"_",as.character(RC),"_",as.character(abundance.weighted),"_",label,"_PAIRWISE",".csv",sep=""))
}
