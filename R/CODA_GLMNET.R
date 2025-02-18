#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script to perform variable selection through penalized regression on the set of all pairwise log-ratios.
#Version: 1.2 #Added expression levels for binary outcomes

# library(coda4microbiome)
# library(phyloseq)
# library(ggplot2)
# library(gplots)
# library(mixOmics)

occ_threshold<-function (m, threshold, max_absent = 0)
{
  occs <- colSums(m > max_absent)
  goodspecies <- occs >= threshold
  return(m[, goodspecies])
}

coda_glmnet_analysis <- function(physeq) {
  res<-NULL

  # TODO: remove
  grouping_column <- "Groups"

  # Process input data
  abund_table <- phyloseq::otu_table(physeq)
  meta_table <- phyloseq::sample_data(physeq)
  OTU_taxonomy <- phyloseq::tax_table(physeq)
  OTU_tree <- phyloseq::phy_tree(physeq)

  #First apply occupancy threshold to remove low-occupancy taxa
  abund_table<-occ_threshold(abund_table,occupancy_threshold)

  #Now choose top_N abundant taxa
  abund_table<-abund_table[,order(colSums(abund_table),decreasing=TRUE)][,1:min(top_N,dim(abund_table)[2])]

  #Apply normalisation just for visualisation of discriminating transcripts
  normalised_table<-abund_table

  #Apply normalisation methods
  if(normalisation_method=="logrelative"){
    normalised_table<-log((normalised_table+1)/(rowSums(normalised_table)+dim(normalised_table)[2]))
  } else if(normalisation_method=="CSS"){
    #From metagenomeSeq website
    data.metagenomeSeq = newMRexperiment(t(normalised_table),
                                         featureData=NULL, libSize=NULL, normFactors=NULL) #using filtered data
    p = cumNormStat(data.metagenomeSeq) #default is 0.5
    data.cumnorm = cumNorm(data.metagenomeSeq, p=p)
    #data.cumnorm
    normalised_table = t(MRcounts(data.cumnorm, norm=TRUE, log=TRUE))
  } else if (normalisation_method %in% c("TSS+ILR","TSS+CLR")){
    TSS.divide = function(x){
      x/sum(x)
    }
    if(normalisation_method=="TSS+ILR"){
      normalised_table<-logratio.transfo(t(apply(normalised_table+1, 1, TSS.divide)),logratio="ILR")
    } else if (normalisation_method=="TSS+CLR"){
      normalised_table<-logratio.transfo(t(apply(normalised_table+1, 1, TSS.divide)),logratio="CLR")
    }
  }
  normalised_table<-as(normalised_table,"matrix")


  #Now loop through all groups in meta_table$Groups
  for (i in levels(meta_table$Groups)){
    #Extract subset of abund_table and meta_table for each group
    mt<-meta_table[meta_table$Groups==i,,drop=F]
    at<-abund_table[rownames(mt),,drop=F]
    #Adjust at to get rid of empty feature column
    at<-at[,colSums(at)>0]
    #Now through all the environmental covariates

  tryCatch({
    for (j in environmental_covariates){
      #Now apply only when the environmental_covariate doesn't have same value and the algorithm doesn't pick it up as binary outcome variable
      if(length(unique(mt[,j]))>1){
        #Get rid of missing data
        mt2<-mt[complete.cases(mt[,j]),,drop=FALSE]
        at2<-at[rownames(mt2),,drop=FALSE]
        #The choice of choosing the penalized parameter as lambda="lambda.min" instead of "lambda.1se" is purely
        #to increase the number of variables selected
        #Check if there are only two values
        if(length(unique(mt2[,j]))==2){
          if(class(mt2[,j])!="factor"){
            mt2[,j]<-as.factor(as.character(mt2[,j]))
          }
          res<-coda4microbiome::coda_glmnet(x=at2,y=mt2[,j])
        } else {
          res<-coda4microbiome::coda_glmnet(x=at2,y=mt2[,j],lambda="lambda.min",showPlots=FALSE)
        }
      }
    }
  }, error = function(e) {
    message(conditionMessage(e))
  })
  }

  class(res) <- "CompositionalRegression"
  return (res)
}

plot.CompositionalRegression <- function(res) {

  #Now draw the expression of the selected taxa
  df<-reshape2::melt(normalised_table[,res$taxa.name])
  colnames(df)<-c("Sample","Feature","Value")
  df<-data.frame(df,Groups=mt2[as.character(df$Sample),j])

  p<-ggplot(df,aes(Groups,Value,colour=Groups))+ylab(normalisation_method)
  p<-p+geom_boxplot(outlier.size=0,show.legend=FALSE,position="identity")+geom_jitter(position = position_jitter(height = 0, width=0), size=2)
  p<-p+facet_wrap( ~ Feature , scales="free_x",nrow=1)
  p<-p+theme_bw()
  p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+theme(strip.text.x = element_text(size = 16, colour = "black", angle = 90))
  p<-p+scale_color_manual("Groups",values=c("#ff7e24","#7967ed"))
  p<-p+guides(colour=FALSE) #FALSE
  pdf(paste("Expressions_plot_",label,"_",i,"_",j,".pdf",sep=""),width=ceiling((length(res$taxa.name)*80/200)+2.6),height=15)
  print(p)
  dev.off()

  pdf(paste("Signature_plot_",label,"_",i,"_",j,".pdf",sep=""),height=max(3,ceiling(length(res$taxa.num)/4)*height_adjustment),width=20)
  plot(res$`signature plot`)
  dev.off()
  pdf(paste("Predictions_plot_",label,"_",i,"_",j,".pdf",sep=""),height=5,width=10)
  plot(res$`predictions plot`)
  dev.off()
}
