#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script for null modelling (beta NTI, Raup-Crick Beta-Diversity, and Elements of Metacommunity)
#Reference: http://uu.diva-portal.org/smash/get/diva2:1373632/DATASET09.txt
#v1.0 (All metrics are now being saved)
#v1.1 (Found a serious bug as Raup-Crick Beta-Diversity can be used in both incidence and Bray-Curtis model)

#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script for visualisation of null modelling (beta NTI, Raup-Crick Beta-Diversity, and Elements of Metacommunity)
#Reference: http://uu.diva-portal.org/smash/get/diva2:1373632/DATASET09.txt
#Authors: Umer, Anna, and Simon
#Versions: 1.2 (fixed ordering issue)

# library(phyloseq)
# library(vegan)
# library(ape)
# library(picante)
# library(ecodist)
# library(metacom)
# library(ggplot2)
# library(yarrr)

qpe <- function(physeq) {
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
    met_ems_group = metacom::Metacommunity(abund_table_ems_group, scores = 1, method = "r1", sims = 999, order = T, binary = F, verbose = T, allowEmpty = T)
    om_ems_group = metacom::OrderMatrix(abund_table_ems_group, outputScores = T, binary = F)

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
    m=picante::match.phylo.data(OTU_tree, t(abund_table_ems_group)) #Extract subtree for each group
    OTU_tree_ems_group=m$phy

    abund_table_ems_group<-t(abund_table_ems_group)

    beta.mntd.weighted = as.matrix(picante::comdistnt(t(abund_table_ems_group), stats::cophenetic(OTU_tree_ems_group), abundance.weighted = TRUE))


    beta.reps = 999;
    rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(abund_table_ems_group),ncol(abund_table_ems_group),beta.reps));
    dim(rand.weighted.bMNTD.comp);

    for (rep in 1:beta.reps) {

      rand.weighted.bMNTD.comp[,,rep] = as.matrix(picante::comdistnt(t(abund_table_ems_group),taxaShuffle(cophenetic(OTU_tree_ems_group)),abundance.weighted=T,exclude.conspecifics = F));

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

plot.ECQuantitativeProcessEstimate <- function(value) {
  PairwisebNTI <- value$pairwise_bNTI
  PairwiseRC <- value$pairwise_RC
  RC <- value$RC
  Coherence <- value$coherence
  BoundaryClump <- value$boundary
  Turnover <- value$turnover

  label="Hypothesis1"
  ordering<-c(
    "COV",
    "CH"
  )
  colours <- c(
    "red",
    "blue",
    #Next colors are for lines mainly used in the PCoA script
    "#000080","#4876FF","#CAE1FF","#9FB6CD","#1E90FF","#00F5FF","#00C957",grey.colors(1000));
  RC_width=8
  RC_height=8
  QPE_width=8
  QPE_height=8
  EMS_width=8
  EMS_height=8
  #PARAMETERS ###########################

  QPE_table<-cbind(Var1=as.character(PairwiseRC$Var1), Var2=as.character(PairwiseRC$Var2),bNTI=PairwisebNTI[,"value",drop=F], RC=PairwiseRC[,"value",drop=F], Groups=PairwiseRC[,"Groups",drop=F])
  #Now we want to get rid of any rows that have NA there to get the comparisons from
  # N x N to N (N-1)/2. A way to do this to use complete.cases
  QPE_table<-QPE_table[complete.cases(QPE_table),]
  names(QPE_table)<-c("Var1","Var2","bNTI","RC","Groups")
  QPE_table$Groups<-factor(as.character(QPE_table$Groups))

  QPE_df<-NULL

  for(i in levels(QPE_table$Groups)){
    tmp<-QPE_table[QPE_table$Groups==i,]
    sp_mask<-abs(tmp$bNTI)>2
    hs_count<-sum(tmp$bNTI[sp_mask]<2)
    vs_count<-sum(tmp$bNTI[sp_mask]>2)
    sig_count<-sum(sp_mask)
    total_count<-nrow(tmp)
    nonsig_count<-nrow(tmp[!sp_mask,])
    dl_count<-sum(tmp[!sp_mask,"RC"]>0.95)
    hd_count<-sum(tmp[!sp_mask,"RC"]<(-0.95))
    ed_count<-nonsig_count-dl_count-hd_count
    tmp2<-data.frame(measure="Homogeneous Selection",value=(hs_count/total_count*100),Groups=i)
    tmp2<-rbind(tmp2,data.frame(measure="Variable Selection",value=(vs_count/total_count*100),Groups=i))
    tmp2<-rbind(tmp2,data.frame(measure="Dispersal Limitation",value=(dl_count/total_count*100),Groups=i))
    tmp2<-rbind(tmp2,data.frame(measure="Undominated",value=(ed_count/total_count*100),Groups=i))
    tmp2<-rbind(tmp2,data.frame(measure="Homogenizing Dispersal",value=(hd_count/total_count*100),Groups=i))
    if(is.null(QPE_df)){QPE_df<-tmp2} else {QPE_df<-rbind(QPE_df,tmp2)}
  }

  #Change the ordering of the Groups
  QPE_df$Groups<-factor(as.character(QPE_df$Groups),levels=ordering)
  #pdf(paste("QPE_",label,".pdf",sep=""),width=QPE_width,height=QPE_height)

  p <- ggplot(QPE_df, aes(x = Groups, y = value, fill = Groups)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%0.2f", round(value, digits = 2))), vjust = -0.3, size = 3.5) +
    facet_wrap(~ measure, strip.position = "left", ncol = 1, scales = "free_y") +
    scale_fill_manual(values = colours) +
    ylim(0, 110) +
    ylab("% Assembly Processes") +
    theme(panel.spacing = unit(2, "lines"), axis.text.x = element_text(angle = 90, hjust = 1))

  return(p)
}

plot_ems <- function() {
  #We move on to EMS (Elements of Metacommunity Structure)
  #UMER: Not confident about interpretation, MAY BE BUGGY!
  collated_community_types<-NULL
  for(i in rownames(Coherence)){
    community_type="Random"
    if(Coherence[i,"p"]<0.05){
      #Top left figure
      if(Coherence[i,"z"]<(-1.96)){
        community_type="Checkerboard"
      } else if(Coherence[i,"z"]>1.96) {
        #Top middle left
        if(Turnover[i,"z"]<(-1.96)){
          if(BoundaryClump[i,"p"]<0.05){
            if(BoundaryClump[i,"index"]<0){
              community_type="Nested Hyperdispersed species loss"
            } else {
              community_type="Nested Clumped species loss"
            }
          }
          else{
            community_type="Nested Random species loss"
          }
        }
        #Top middle right
        else if(Turnover[i,"z"]>1.96)
        {
          if(BoundaryClump[i,"p"]<0.05){
            if(BoundaryClump[i,"index"]<0){
              community_type="Evenly spaced"
            } else {
              community_type="Clementsian"
            }
          }
          else{
            community_type="Gleasonian"
          }
        }
        if(Turnover[i,"p"]>0.05){
          community_type=paste("Quasi-structure",community_type)
        }
      }

    }
    #Collate all the information together
    if(is.null(collated_community_types)){collated_community_types<-community_type} else {collated_community_types<-c(collated_community_types,community_type)}
  }

  EMS<-cbind(Coherence,Turnover,BoundaryClump,Metacommunity=collated_community_types)
  names(EMS)<-c("Coherence_Abs","Coherence_z","Coherence_p","Coherence_simMean","Coherence_simVariance",
                "Turnover_turn","Turnover_z","Turnover_p","Turnover_simMean","Turnover_simVariance",
                "Clumping_index","Clumping_p","Clumping_df","Metacommunity")
  write.csv(EMS,file=paste("EMS_",label,".csv",sep=""))

  EMS<-cbind(EMS,Groups=rownames(EMS))
  #Change the ordering of the Groups
  EMS$Groups<-factor(as.character(EMS$Groups),levels=ordering)
  pdf(paste("EMS_",label,"_DONOTUSE.pdf",sep=""),width=EMS_width,height=EMS_height)
  p<-ggplot(EMS,aes(Groups,Coherence_z,color=Turnover_z,size=Clumping_index,shape=Metacommunity))
  p<-p+geom_point()
  p <- p + geom_hline(yintercept = 1.96,linetype="dotted")
  p <-p +geom_text(aes(x=1,y=1.96, label="1.96\n"), colour="blue",size=2)

  p <- p + geom_hline(yintercept = -1.96,linetype="dotted")
  p <-p +geom_text(aes(x=1,y=-1.96, label="\n-1.96"), colour="blue",size=2)

  p <- p + ylab("Coherence (z-value)")
  p <- p + scale_color_continuous("Turnover (z-value)")
  p <- p + scale_size_continuous("Boundary clumping (Morisita's index)")
  p<-p+theme(strip.background = element_rect(fill="white"))+theme(panel.spacing = unit(2, "lines"),
                                                                  axis.text.x = element_text(angle = 90, hjust = 1))

  p<-p+theme_minimal()

  return(p)
}

# pdf(paste("RC_",label,".pdf",sep=""),width=RC_width,height=RC_height)
plot_rc <- function(RC) {
  RC<-cbind(RC,Groups=rownames(RC))
  #Change the ordering of the Groups
  RC$Groups<-factor(as.character(RC$Groups),levels=ordering)

  p<-ggplot(RC,aes(Groups,value,colour=Groups))
  p<-p+geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1, lty=1)+
    geom_point(size=5)

  p <- p + geom_hline(yintercept = 0,linetype="dotted")
  p <-p +geom_text(aes(x=1,y=0, label="0\n"), colour="blue",size=2)

  p <- p + geom_hline(yintercept = 1,linetype="dotted")
  p <-p +geom_text(aes(x=1,y=1, label="+1\n"), colour="blue",size=2)

  p <- p + geom_hline(yintercept = -1,linetype="dotted")
  p <-p +geom_text(aes(x=1,y=-1, label="\n-1"), colour="blue",size=2)

  p<-p+theme_minimal()
  p<-p+scale_colour_manual(values=colours)
  p<-p+ylab("Incidence-based beta-diversity (±SE)")
  p<-p+theme(strip.background = element_rect(fill="white"))+theme(panel.spacing = unit(2, "lines"),
                                                                  axis.text.x = element_text(angle = 90, hjust = 1))

  return(p)
}
