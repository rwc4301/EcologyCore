#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script for calculation of phylogenetic alpha diversity metrics such as NRI/NTI
#to give an account of stochastic versus deterministic nature of microbial communities
#v1.2 (All metrics are now being saved)

environmental_filtering <- function(physeq, runs = 999, iterations = 1000, top_n = 2000, abundance.weighted = FALSE, null.model = c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap")) {
  null.model <- match.arg(null.model)

  # Process input data
  abund_table <- phyloseq::otu_table(physeq)
  meta_table <- phyloseq::sample_data(physeq)
  OTU_taxonomy <- phyloseq::tax_table(physeq)
  OTU_tree <- phyloseq::phy_tree(physeq)

  #We extract N most abundant OTUs
  abund_table<-abund_table[,order(colSums(abund_table),decreasing=TRUE)][,1:min(top_n,dim(abund_table)[2])]

  #Adjust OTU_tree
  OTU_tree<-ape::drop.tip(OTU_tree,OTU_tree$tip.label[!OTU_tree$tip.label %in% colnames(abund_table)])


  #Sometimes the filtering will remove the OTUs and these won't exist in the OTU_tree
  abund_table<-abund_table[,OTU_tree$tip.label]

  #There is a bug in adespatial as it doesn't like phyloseq's class and won't calculate
  #the values, so I am forcing it to become matrix
  comm <- as(abund_table,"matrix")
  dist <- cophenetic(OTU_tree)

  #Reference http://kembellab.ca/r-workshop/biodivR/SK_Biodiversity_R.html
  abund_table.sesmpd <- ses.mpd(comm, dist, null.model = null.model, abundance.weighted = abundance.weighted, runs = runs, iterations = iterations)
  abund_table.sesmntd <- ses.mntd(comm, dist, null.model = null.model, abundance.weighted = abundance.weighted, runs = runs, iterations = iterations)

  #Write all the metrics in a file for further analyses elsewhere
  data_to_write<-data.frame(abund_table.sesmpd[,"mpd.obs.z",drop=F],abund_table.sesmntd[,"mntd.obs.z",drop=F])
  #Convert SES scores to NRI/NTI
  data_to_write<-data_to_write*-1
  colnames(data_to_write)<-c("NRI","NTI")
  #write.csv(data_to_write,file=paste("Environmental_Filtering","_",second_label,"_",label,"_",null.model,".csv",sep=""))
  #/Write all the metrics in a file for further analyses elsewhere

  #Positive SES values (abund_table.sesmpd$mpd.obs.z > 0) and high quantiles (abund_table.sesmpd$mpd.obs.p > 0.95)
  #indicate phylogenetic evenness, while negative SES values and low quantiles (abund_table.sesmpd$mpd.obs.p < 0.05)
  #indicate phylogenetic clustering, relative to the null model. MPD is generally thought to be more sensitive to
  #tree-wide patterns of phylogenetic clustering and eveness, while MNTD is more sensitive to patterns of evenness
  #and clustering closer to the tips of the phylogeny.
  df<-rbind(data.frame(value=-1*abund_table.sesmpd$mpd.obs.z,meta_table,measure="NRI"),
            data.frame(value=-1*abund_table.sesmntd$mntd.obs.z,meta_table,measure="NTI"))

  #To do anova, we will convert our data.frame to data.table
  library(data.table)
  grouping_column="Groups"
  #Since we can't pass a formula to data.table, I am creating
  #a dummy column .group. so that I don't change names in the formula
  dt<-data.table(data.frame(df,.group.=df[,grouping_column]))

  #I am also specifying a p-value cutoff for the ggplot2 strips
  pValueCutoff<-0.05
  pval<-dt[, list(pvalue = sprintf("%.2g",
                                   tryCatch(summary(aov(value ~ .group.))[[1]][["Pr(>F)"]][1],error=function(e) NULL))),
           by=list(measure)]

  #Filter out pvals that we don't want
  pval<-pval[!pval$pvalue=="",]
  pval<-pval[as.numeric(pval$pvalue)<=pValueCutoff,]

  #I am using sapply to generate significances for pval$pvalue using the cut function.
  pval$pvalue<-sapply(as.numeric(pval$pvalue),function(x){as.character(cut(x,breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))})

  #Update df$measure to change the measure names if the grouping_column has more than three classes
  if(length(unique(as.character(df[,grouping_column])))>2){
    df$measure<-as.character(df$measure)
    if(dim(pval)[1]>0){
      for(i in seq(1:dim(pval)[1])){
        df[df$measure==as.character(pval[i,measure]),"measure"]=paste(as.character(pval[i,measure]),as.character(pval[i,pvalue]))
      }
    }
    df$measure<-as.factor(df$measure)
  }

  #Get all possible combination of values in the grouping_column
  s<-NULL
  if(provide_your_own_pvalue_combinations){
    s<-provided_combination
  } else{
    s<-combn(unique(as.character(df[,grouping_column])),2)
  }

  #df_pw will store the pair-wise p-values
  df_pw<-NULL
  for(k in unique(as.character(df$measure))){
    #We need to calculate the coordinate to draw pair-wise significance lines
    #for this we calculate bas as the maximum value
    bas<-max(df[(df$measure==k),"value"])

    #Calculate increments as % of the maximum values
    inc<-0.05*(bas-min(df[(df$measure==k),"value"]))

    #Give an initial increment
    bas<-bas+inc
    for(l in 1:dim(s)[2]){

      #Do a pair-wise anova
      tmp<-c(k,s[1,l],s[2,l],bas,paste(sprintf("%.2g",tryCatch(summary(aov(as.formula(paste("value ~",grouping_column)),data=df[(df$measure==k) & (df[,grouping_column]==s[1,l] | df[,grouping_column]==s[2,l]),] ))[[1]][["Pr(>F)"]][1],error=function(e) NULL)),"",sep=""))

      #Ignore if anova fails
      if(!is.na(as.numeric(tmp[length(tmp)]))){

        #Only retain those pairs where the p-values are significant
        if(as.numeric(tmp[length(tmp)])<0.05){
          if(is.null(df_pw)){df_pw<-tmp}else{df_pw<-rbind(df_pw,tmp)}

          #Generate the next position
          bas<-bas+inc
        }
      }
    }
  }

  if(!is.null(df_pw)){
    if(sum(class(df_pw) %in% c("character"))>0){
      df_pw<-t(as.matrix(df_pw))
    }
    df_pw<-data.frame(row.names=NULL,df_pw)
    names(df_pw)<-c("measure","from","to","y","p")
  }

  return(structure(list(df = df, df_pw = df_pw, pval = pval, meta_table = meta_table, data_to_write = data_to_write), className = "ECEnvironmentalFiltering"))
}

#' @import picante
ses.mpd <- function(
    samp,
    dis,
    null.model = c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"),
    abundance.weighted = FALSE,
    runs = 999,
    iterations = 1000,
    ncores = parallel::detectCores() - 1
) {
  dis <- as.matrix(dis)
  mpd.obs <- mpd(samp, dis, abundance.weighted = abundance.weighted)
  null.model <- match.arg(null.model)

  # Setup parallel backend
  cl <- parallel::makeCluster(ncores)
  parallel::clusterExport(cl, c("mpd", "taxaShuffle", "randomizeMatrix", "iterations"))

  # Parallel computation of null model
  mpd.rand <- switch(null.model,
    taxa.labels = t(parallel::parSapply(cl, 1:runs, function(x) mpd(samp, taxaShuffle(dis), abundance.weighted = abundance.weighted))),
    richness = t(parallel::parSapply(cl, 1:runs, function(x) mpd(randomizeMatrix(samp, null.model = "richness"), dis, abundance.weighted))),
    frequency = t(parallel::parSapply(cl, 1:runs, function(x) mpd(randomizeMatrix(samp, null.model = "frequency"), dis, abundance.weighted))),
    sample.pool = t(parallel::parSapply(cl, 1:runs, function(x) mpd(randomizeMatrix(samp, null.model = "richness"), dis, abundance.weighted))),
    phylogeny.pool = t(parallel::parSapply(cl, 1:runs, function(x) mpd(randomizeMatrix(samp, null.model = "richness"), taxaShuffle(dis), abundance.weighted))),
    independentswap = t(parallel::parSapply(cl, 1:runs, function(x) mpd(randomizeMatrix(samp, null.model = "independentswap", iterations), dis, abundance.weighted))),
    trialswap = t(parallel::parSapply(cl, 1:runs, function(x) mpd(randomizeMatrix(samp, null.model = "trialswap", iterations), dis, abundance.weighted)))
  )

  parallel::stopCluster(cl)

  # Calculate statistics
  mpd.rand.mean <- apply(X = mpd.rand, MARGIN = 2, FUN = mean, na.rm = TRUE)
  mpd.rand.sd <- apply(X = mpd.rand, MARGIN = 2, FUN = sd, na.rm = TRUE)
  mpd.obs.z <- (mpd.obs - mpd.rand.mean)/mpd.rand.sd
  mpd.obs.rank <- apply(X = rbind(mpd.obs, mpd.rand), MARGIN = 2, FUN = rank)[1, ]
  mpd.obs.rank <- ifelse(is.na(mpd.rand.mean), NA, mpd.obs.rank)

  return (data.frame(
    ntaxa = specnumber(samp),
    mpd.obs,
    mpd.rand.mean,
    mpd.rand.sd,
    mpd.obs.rank,
    mpd.obs.z,
    mpd.obs.p = mpd.obs.rank/(runs + 1),
    runs = runs,
    row.names = row.names(samp)
  ))
}

#' @import picante
ses.mntd <- function(
  samp,
  dis,
  null.model = c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"),
  abundance.weighted = FALSE,
  runs = 999,
  iterations = 1000,
  ncores = parallel::detectCores() - 1
) {
  dis <- as.matrix(dis)
  mntd.obs <- mntd(samp, dis, abundance.weighted)
  null.model <- match.arg(null.model)

  # Setup parallel backend
  cl <- parallel::makeCluster(ncores)
  parallel::clusterExport(cl, c("mntd", "taxaShuffle", "randomizeMatrix", "iterations"))

  # Parallel computation of null model
  mntd.rand <- switch(null.model,
    taxa.labels = t(parallel::parSapply(cl, 1:runs, function(x) mntd(samp, taxaShuffle(dis), abundance.weighted))),
    richness = t(parallel::parSapply(cl, 1:runs, function(x) mntd(randomizeMatrix(samp, null.model = "richness"), dis, abundance.weighted))),
    frequency = t(parallel::parSapply(cl, 1:runs, function(x) mntd(randomizeMatrix(samp, null.model = "frequency"), dis, abundance.weighted))),
    sample.pool = t(parallel::parSapply(cl, 1:runs, function(x) mntd(randomizeMatrix(samp, null.model = "richness"), dis, abundance.weighted))),
    phylogeny.pool = t(parallel::parSapply(cl, 1:runs, function(x) mntd(randomizeMatrix(samp, null.model = "richness"), taxaShuffle(dis), abundance.weighted))),
    independentswap = t(parallel::parSapply(cl, 1:runs, function(x) mntd(randomizeMatrix(samp, null.model = "independentswap", iterations), dis, abundance.weighted))),
    trialswap = t(parallel::parSapply(cl, 1:runs, function(x) mntd(randomizeMatrix(samp, null.model = "trialswap", iterations), dis, abundance.weighted)))
  )

  parallel::stopCluster(cl)

  # Calculate statistics
  mntd.rand.mean <- apply(X = mntd.rand, MARGIN = 2, FUN = mean, na.rm = TRUE)
  mntd.rand.sd <- apply(X = mntd.rand, MARGIN = 2, FUN = sd, na.rm = TRUE)
  mntd.obs.z <- (mntd.obs - mntd.rand.mean)/mntd.rand.sd
  mntd.obs.rank <- apply(X = rbind(mntd.obs, mntd.rand), MARGIN = 2, FUN = rank)[1, ]
  mntd.obs.rank <- ifelse(is.na(mntd.rand.mean), NA, mntd.obs.rank)

  return (data.frame(
    ntaxa = specnumber(samp),
    mntd.obs,
    mntd.rand.mean,
    mntd.rand.sd,
    mntd.obs.rank,
    mntd.obs.z,
    mntd.obs.p = mntd.obs.rank/(runs + 1),
    runs = runs,
    row.names = row.names(samp)
  ))
}
