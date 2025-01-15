#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script for null modelling (beta NTI, Raup-Crick Beta-Diversity, and Elements of Metacommunity)
#Reference: http://uu.diva-portal.org/smash/get/diva2:1373632/DATASET09.txt
#v1.0 (All metrics are now being saved)
#v1.1 (Found a serious bug as Raup-Crick Beta-Diversity can be used in both incidence and Bray-Curtis model)

library(phyloseq)
library(vegan)
library(ape)
library(picante)
library(ecodist)
library(metacom)

#PARAMETERS ###########################
physeq<-import_biom("../../data/VSEARCH/feature_w_tax.biom")
meta_table<-read.csv("../../data/meta_table.csv",header=T,row.names=1)
#Load the tree using ape package
OTU_tree <- read.tree("../../data/VSEARCH/tree.nwk")
#/PARAMETERS ###########################

abund_table<-otu_table(physeq)
abund_table<-t(abund_table)
#Uncomment if you'd like to get rid of samples below a certain library size
abund_table<-abund_table[rowSums(abund_table)>=5000,]


OTU_taxonomy<-as.data.frame(tax_table(physeq))
colnames(OTU_taxonomy)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Otus")

#Ensure that all columns of OTU_taxonomy are character and not factors
OTU_taxonomy[] <- lapply(OTU_taxonomy, function(x) as.character(x))
OTU_taxonomy[is.na(OTU_taxonomy)]<-""
OTU_taxonomy$Otus<-gsub("D_6__|s__","",OTU_taxonomy$Otus)
OTU_taxonomy$Genus<-gsub("D_5__|g__","",OTU_taxonomy$Genus)
OTU_taxonomy$Family<-gsub("D_4__|f__","",OTU_taxonomy$Family)
OTU_taxonomy$Order<-gsub("D_3__|o__","",OTU_taxonomy$Order)
OTU_taxonomy$Class<-gsub("D_2__|c__","",OTU_taxonomy$Class)
OTU_taxonomy$Phylum<-gsub("D_1__|p__","",OTU_taxonomy$Phylum)
OTU_taxonomy$Kingdom<-gsub("D_0__|d__","",OTU_taxonomy$Kingdom)

#Remove singletons and adjust OTU_taxonomy
abund_table<-abund_table[,colSums(abund_table)>1]
OTU_taxonomy<-OTU_taxonomy[colnames(abund_table),]

#get rid of contaminants with "Unassigned", "Chloroplast" and "Mitochondria" assignment", and "non classified" at Phylum level
abund_table<-abund_table[,!(OTU_taxonomy$Kingdom %in% c("Unassigned") | OTU_taxonomy$Phylum=="" | OTU_taxonomy$Order %in% c("Chloroplast") | OTU_taxonomy$Family %in% c("Mitochondria"))]


#extract subset of abund_table for which samples also exists in meta_table
abund_table<-abund_table[rownames(abund_table) %in% rownames(meta_table),]
#when reducing the abund_table, there is a high likelihood that an OTU was only present in a sample that is removed, so we shrink
#the abund_table to get rid of empty columns
abund_table<-abund_table[,colSums(abund_table)>0]
#make your meta_table smaller by only considering samples that appear in abund_table
meta_table<-meta_table[rownames(abund_table),]
#make OTU_taxonomy smaller by only considering OTUs that appear in abund_table
OTU_taxonomy<-OTU_taxonomy[colnames(abund_table),]
#At this point we have abund_table, meta_table, and OTU_taxonomy are ready and their dimensions should match
#/DATA IMPORT############################################################

#PARAMETERS CHANGE THE GROUPING COLUMN AS YOU DESIRE############################
#In the hypothesis space, all you need is to select the rows in meta_table you are interested in
#and then allocate a column to meta_table$Groups that you want to use.

label="Hypothesis1"
meta_table<-meta_table[!meta_table$Within_Dam_Sample %in% c("Piezometer"),]
meta_table$Groups<-as.factor(as.character(meta_table$Within_Dam_Sample))

#PARAMETERS CHANGE THE GROUPING COLUMN AS YOU DESIRE############################

#Recent version of R puts apostrophe in the OTU tip labels so we will just remove them if that exist
OTU_tree$tip.label<-gsub("'","",OTU_tree$tip.label)

#Adjust abund_table to contain only those rows that got selected in the Hypothesis space
abund_table<-abund_table[rownames(meta_table),]
#After adjustment, get rid of OTUs that are all empty
abund_table<-abund_table[,colSums(abund_table)>0]
#Adjust OTU taxonomy
OTU_taxonomy<-OTU_taxonomy[colnames(abund_table),]

#Convert the data to phyloseq format
OTU = otu_table(as.matrix(abund_table), taxa_are_rows = FALSE)
TAX = tax_table(as.matrix(OTU_taxonomy))
SAM = sample_data(meta_table)

physeq<-merge_phyloseq(phyloseq(OTU, TAX),SAM,OTU_tree)

#Pruning and subsampling
physeq<-prune_taxa(taxa_sums(physeq)>10, physeq)

#Rarefy to minimum sample size
physeq_rel = rarefy_even_depth(physeq, sample.size = min(sample_sums(physeq)))


#incidence_based Raup_Crick: https://github.com/jfq3/QsNullModels/blob/master/R/raup_crick.R
raup_crick=function(spXsite, plot_names_in_col1=FALSE, classic_metric=FALSE, 
                    split_ties=TRUE, reps=999, set_all_species_equal=FALSE, 
                    as.distance.matrix=TRUE, report_similarity=FALSE){
  
  # expects a species by site matrix for spXsite, with row names for plots, 
  # or optionally plots named in column 1.  By default calculates a modification 
  # of the Raup-Crick metric (standardizing the metric to range from -1 to 1 
  # instead of 0 to 1). Specifying classic_metric=TRUE instead calculates the 
  # original Raup-Crick metric that ranges from 0 to 1. The option split_ties 
  # (defaults to TRUE) adds half of the number of null observations that are 
  # equal to the observed number of shared species to the calculation- this 
  # is highly recommended.  The argument report_similarity defaults to FALSE 
  # so the function reports a dissimilarity (which is appropriate as a measure 
  # of beta diversity).  Setting report_similarity=TRUE returns a measure of 
  # similarity, as Raup and Crick originally specified.  If ties are split 
  # (as we recommend) the dissimilarity (default) and similarity (set 
  # report_similarity=TRUE) calculations can be flipped by multiplying by -1 
  # (for our modification, which ranges from -1 to 1) or by subtracting the 
  # metric from 1 (for the classic metric which ranges from 0 to 1). If ties 
  # are not split (and there are ties between the observed and expected shared 
  # number of species) this conversion will not work. The argument reps specifies
  # the number of randomizations (a minimum of 999 is recommended- default was 
  # 9999).  set_all_species_equal weights all species equally in the null model 
  # instead of weighting species by frequency of occupancy.  
  
  # Note that the choice of how many plots (rows) to include has a real impact 
  # on the metric, as species and their occurrence frequencies across the set 
  # of plots is used to determine gamma and the frequency with which each 
  # species is drawn from the null model	
  
  # this section moves plot names in column 1 (if specified as being present) 
  # into the row names of the matrix and drops the column of names
  if(plot_names_in_col1){
    row.names(spXsite)<-spXsite[,1]
    spXsite<-spXsite[,-1]
  }
  
  # count number of sites and total species richness across all plots (gamma)
  n_sites<-nrow(spXsite)
  gamma<-ncol(spXsite)
  
  # make the spXsite matrix into a pres/abs. (overwrites initial spXsite matrix):
  ceiling(spXsite/max(spXsite))->spXsite
  
  # create an occurrence vector- used to give more weight to widely distributed 
  # species in the null model:
  occur<-apply(spXsite, MARGIN=2, FUN=sum)
  
  # NOT recommended- this is a non-trivial change to the metric:
  # sets all species to occur with equal frequency in the null model
  # e.g.- discards any occupancy frequency information
  if(set_all_species_equal){
    occur<-rep(1,gamma)
  }
  
  # determine how many unique species richness values are in the dataset
  # this is used to limit the number of null communities that have to be 
  # calculated
  alpha_levels<-sort(unique(apply(spXsite, MARGIN=1, FUN=sum)))
  
  # make_null:
  
  # alpha_table is used as a lookup to help identify which null distribution
  # to use for the tests later.  It contains one row for each combination of 
  # alpha richness levels. 
  
  alpha_table<-data.frame(c(NA), c(NA))
  names(alpha_table)<-c("smaller_alpha", "bigger_alpha")
  col_count<-1
  
  # null_array will hold the actual null distribution values.  Each element
  # of the array corresponds to a null distribution for each combination of 
  # alpha values.  The alpha_table is used to point to the correct null 
  # distribution- the row numbers of alpha_table correspond to the [[x]] 
  # indices of the null_array.  Later the function will find the row of 
  # alpha_table with the right combination of alpha values.  That row number 
  # is used to identify the element of null_array that contains the correct 
  # null distribution for that combination of alpha levels. 
  
  null_array<-list()
  
  # looping over each combination of alpha levels:
  
  for(a1 in 1:length(alpha_levels)){
    for(a2 in a1:length(alpha_levels)){
      
      # build a null distribution of the number of shared species for a 
      # pair of alpha values:
      null_shared_spp<-NULL
      for(i in 1:reps){
        
        # two empty null communities of size gamma:
        com1<-rep(0,gamma)
        com2<-rep(0,gamma)
        
        # add alpha1 number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma, alpha_levels[a1], replace=FALSE, prob=occur)]<-1
        
        # same for com2:
        com2[sample(1:gamma, alpha_levels[a2], replace=FALSE, prob=occur)]<-1
        
        # how many species are shared in common?
        null_shared_spp[i]<-sum((com1+com2)>1)
        
      }
      
      # store null distribution, record values for alpha 1 and 2 in the alpha_table to 
      # help find the correct null distribution later:
      null_array[[col_count]]<-null_shared_spp
      
      alpha_table[col_count, which(names(alpha_table)=="smaller_alpha")]<-alpha_levels[a1]
      alpha_table[col_count, which(names(alpha_table)=="bigger_alpha")]<-alpha_levels[a2]
      
      # increment the counter for the columns of the alpha table/ elements of the null array
      col_count<-col_count+1	
      
    }
    
  }
  
  # create a new column with both alpha levels to match on:
  alpha_table$matching<-paste(alpha_table[,1], alpha_table[,2], sep="_")
  
  #####################
  # do the test:
  
  # build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))
  
  # for each pair of sites (duplicates effort now to make a full matrix instead 
  # of a half one- but this part should be minimal time as compared to the null 
  # model building)
  for(i in 1:n_sites){
    for(j in 1:n_sites){
      
      # how many species are shared between the two sites:
      n_shared_obs<-sum((spXsite[i,]+spXsite[j,])>1)
      
      # what was the observed richness of each site?
      obs_a1<-sum(spXsite[i,])
      obs_a2<-sum(spXsite[j,])
      
      # place these alphas into an object to match against alpha_table (sort so 
      # smaller alpha is first)
      obs_a_pair<-sort(c(obs_a1, obs_a2))
      
      # match against the alpha table- row index identifies which element of the 
      # null array contains the correct null distribution for the observed 
      # combination of alpha values:
      null_index<-which(alpha_table$matching==paste(obs_a_pair[1], obs_a_pair[2], sep="_"))
      
      # how many null observations is the observed value tied with?
      num_exact_matching_in_null<-sum(null_array[[null_index]]==n_shared_obs)
      
      # how many null values are bigger than the observed value?
      num_greater_in_null<-sum(null_array[[null_index]]>n_shared_obs)
      
      rc<-(num_greater_in_null)/reps
      
      if(split_ties){
        
        rc<-((num_greater_in_null+(num_exact_matching_in_null)/2)/reps)
      }
      
      if(!classic_metric){
        
        # our modification of raup crick standardizes the metric to range 
        # from -1 to 1 instead of 0 to 1
        
        rc<-(rc-.5)*2
      }
      
      #  at this point rc represents an index of dissimilarity- multiply by -1 
      # to convert to a similarity as specified in the original 1979 Raup Crick paper
      if(report_similarity & !classic_metric){
        rc<- rc*-1
      }
      
      # the switch to similarity is done differently if the original 0 to 1 range 
      # of the metric is used:
      if(report_similarity & classic_metric){
        rc<- 1-rc
      }
      
      # store the metric in the results matrix:
      results[i,j]<-round(rc, digits=2)
      
    }
  }
  
  if(as.distance.matrix){
    results<-as.dist(results)
  }	
  
  return(results)
  
}





# Second step of the QPE approach (abundance-based Raup-Crick beta-diversity)
raup_crick_abundance = function(spXsite, plot_names_in_col1=TRUE, classic_metric=FALSE, split_ties=TRUE, reps=9999, set_all_species_equal=FALSE, as.distance.matrix=TRUE, report_similarity=FALSE){
  
  ##expects a species by site matrix for spXsite, with row names for plots, or optionally plots named in column 1.  By default calculates a modification of the Raup-Crick metric (standardizing the metric to range from -1 to 1 instead of 0 to 1). Specifying classic_metric=TRUE instead calculates the original Raup-Crick metric that ranges from 0 to 1. The option split_ties (defaults to TRUE) adds half of the number of null observations that are equal to the observed number of shared species to the calculation- this is highly recommended.  The argument report_similarity defaults to FALSE so the function reports a dissimilarity (which is appropriate as a measure of beta diversity).  Setting report_similarity=TRUE returns a measure of similarity, as Raup and Crick originally specified.  If ties are split (as we recommend) the dissimilarity (default) and similarity (set report_similarity=TRUE) calculations can be flipped by multiplying by -1 (for our modification, which ranges from -1 to 1) or by subtracting the metric from 1 (for the classic metric which ranges from 0 to 1). If ties are not split (and there are ties between the observed and expected shared number of species) this conversion will not work. The argument reps specifies the number of randomizations (a minimum of 999 is recommended- default is 9999).  set_all_species_equal weights all species equally in the null model instead of weighting species by frequency of occupancy.
  
  
  ##Note that the choice of how many plots (rows) to include has a real impact on the metric, as species and their occurrence frequencies across the set of plots is used to determine gamma and the frequency with which each species is drawn from the null model
  
  
  ##this section moves plot names in column 1 (if specified as being present) into the row names of the matrix and drops the column of names
  if(plot_names_in_col1){
    row.names(spXsite)<-spXsite[,1]
    spXsite<-spXsite[,-1]
  }
  
  
  ## count number of sites and total species richness across all plots (gamma)
  n_sites<-nrow(spXsite)
  gamma<-ncol(spXsite)
  
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))
  
  ##make the spXsite matrix into a new, pres/abs. matrix:
  ceiling(spXsite/max(spXsite))->spXsite.inc
  
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur<-apply(spXsite.inc, MARGIN=2, FUN=sum)
  
  ##create an abundance vector- used to give more weight to abundant species in the second step of the null model:
  abundance<-apply(spXsite, MARGIN=2, FUN=sum)
  
  ##make_null:
  
  ##looping over each pairwise community combination:
  
  for(null.one in 1:(nrow(spXsite)-1)){
    for(null.two in (null.one+1):nrow(spXsite)){
      
      null_bray_curtis<-NULL
      for(i in 1:reps){
        
        ##two empty null communities of size gamma:
        com1<-rep(0,gamma)
        com2<-rep(0,gamma)
        
        ##add observed number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma, sum(spXsite.inc[null.one,]), replace=FALSE, prob=occur)]<-1
        com1.samp.sp = sample(which(com1>0),(sum(spXsite[null.one,])-sum(com1)),replace=TRUE,prob=abundance[which(com1>0)]);
        com1.samp.sp = cbind(com1.samp.sp,1); # head(com1.samp.sp);
        com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2],com1.samp.sp[,1],FUN=sum)); colnames(com1.sp.counts) = 'counts'; # head(com1.sp.counts);
        com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts)); # head(com1.sp.counts);
        com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts; # com1;
        #sum(com1) - sum(spXsite[null.one,]); ## this should be zero if everything work properly
        rm('com1.samp.sp','com1.sp.counts');
        
        ##same for com2:
        com2[sample(1:gamma, sum(spXsite.inc[null.two,]), replace=FALSE, prob=occur)]<-1
        com2.samp.sp = sample(which(com2>0),(sum(spXsite[null.two,])-sum(com2)),replace=TRUE,prob=abundance[which(com2>0)]);
        com2.samp.sp = cbind(com2.samp.sp,1); # head(com2.samp.sp);
        com2.sp.counts = as.data.frame(tapply(com2.samp.sp[,2],com2.samp.sp[,1],FUN=sum)); colnames(com2.sp.counts) = 'counts'; # head(com2.sp.counts);
        com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts)); # head(com2.sp.counts);
        com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts; # com2;
        # sum(com2) - sum(spXsite[null.two,]); ## this should be zero if everything work properly
        rm('com2.samp.sp','com2.sp.counts');
        
        null.spXsite = rbind(com1,com2); # null.spXsite;
        
        ##calculate null bray curtis
        null_bray_curtis[i] = ecodist::distance(null.spXsite,method='bray-curtis');
        
      }; # end reps loop
      
      ## empirically observed bray curtis
      obs.bray = ecodist::distance(spXsite[c(null.one,null.two),],method='bray-curtis');
      
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null = sum(null_bray_curtis==obs.bray);
      
      ##how many null values are smaller than the observed *dissimilarity*?
      num_less_than_in_null = sum(null_bray_curtis<obs.bray);
      
      rc = (num_less_than_in_null )/reps; # rc;
      
      if(split_ties){
        
        rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
      };
      
      
      if(!classic_metric){
        
        ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
        
        rc = (rc-.5)*2
      };
      
      results[null.two,null.one] = round(rc,digits=2); ##store the metric in the results matrix
      
      print(c(null.one,null.two,date()));
      
    }; ## end null.two loop
    
  }; ## end null.one loop
  
  if(as.distance.matrix){ ## return as distance matrix if so desired
    results<-as.dist(results)
  }	
  
  return(results)
  
}; ## end function


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
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
  results=raup_crick(abund_table_ems_group,plot_names_in_col1 = F, reps = 999, as.distance.matrix = T, set_all_species_equal = F)
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

write.csv(collated_coherence,file=paste("Coherence_",label,".csv",sep=""))
write.csv(collated_boundary,file=paste("Boundary_",label,".csv",sep=""))
write.csv(collated_turnover,file=paste("Turnover_",label,".csv",sep=""))
write.csv(collated_sitescores,file=paste("Sitescores_",label,".csv",sep=""))
write.csv(collated_pairwise_RC_abundance,file=paste("PairwiseRC_",label,".csv",sep=""))
write.csv(collated_pairwise_bNTI,file=paste("PairwisebNTI_",label,".csv",sep=""))
write.csv(collated_RC,file=paste("RC_",label,".csv",sep=""))
write.csv(collated_bNTI,file=paste("bNTI_",label,".csv",sep=""))
