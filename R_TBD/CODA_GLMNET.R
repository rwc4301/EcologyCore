#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script to perform variable selection through penalized regression on the set of all pairwise log-ratios. 
#Version: 1.2 #Added expression levels for binary outcomes

library(coda4microbiome)
library(phyloseq)
library(ggplot2)
library(gplots)
library(mixOmics)


#PARAMETERS ###########################
which_level="Otus" #Otus Genus Family Order Class Phylum
normalisation_method<-"TSS+CLR"
occupancy_threshold<-1 #Retain only those taxa that are observed in X number of samples
top_N=100 #How many most abundant taxa you want to use
height_adjustment=0.9 #height adjustment for the signature_plot
physeq<-import_biom("../../data/VSEARCH/feature_w_tax.biom")
meta_table<-read.csv("../../data/meta_table2.csv",header=T,row.names=1,check.names=FALSE)
#/PARAMETERS ###########################

abund_table<-otu_table(physeq)
abund_table<-t(abund_table)
#Uncomment if you'd like to get rid of samples below a certain library size
abund_table<-abund_table[rowSums(abund_table)>=5000,]


OTU_taxonomy<-data.frame(as(tax_table(physeq),"matrix"))
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
#and then allocate a column to meta_table$Groups that you want to use to capture all the environments
#to regress against the environmental covariates (columns of meta_table) listed in environmental_covariates


label="Hypothesis1"
meta_table<-meta_table[!meta_table$Within_Dam_Sample %in% c("Piezometer"),]
meta_table$Groups<-paste(meta_table$Within_Dam_Sample,meta_table$Sample_Time)
meta_table$Groups<-factor(meta_table$Groups,levels=c(
  "Community Scoop Hole May",
  "Community Scoop Hole July",
  "Downstream May",
  "Downstream July",
  "Excavated Scoop Hole  May",
  "Excavated Scoop Hole  July",
  "Hand Pump May",
  "Hand Pump July",
  "Open Well May",
  "Open Well July"
))

environmental_covariates<-c(
  "pH",
  "Temp",
  "conductivity",
  "turbidity",
  "Phosphate_mg_by_L",
  "Phosphors_mg_by_L",
  "D50_grain_size",
  "Mean_Gradient_of_Catchment"
)


#PARAMETERS CHANGE THE GROUPING COLUMN AS YOU DESIRE############################

#Adjust abund_table to contain only those rows that got selected in the Hypothesis space
abund_table<-abund_table[rownames(meta_table),]
#After adjustment, get rid of OTUs that are all empty
abund_table<-abund_table[,colSums(abund_table)>0]
#Adjust OTU taxonomy
OTU_taxonomy<-OTU_taxonomy[colnames(abund_table),]

#COLLATE OTUS AT A PARTICULAR LEVEL#######################################
new_abund_table<-NULL
if(which_level=="Otus"){
  new_abund_table<-abund_table
} else {
  list<-unique(OTU_taxonomy[,which_level])
  new_abund_table<-NULL
  for(i in list){
    tmp<-data.frame(rowSums(abund_table[,rownames(OTU_taxonomy)[OTU_taxonomy[,which_level]==i],drop=FALSE]))
    if(i==""){colnames(tmp)<-c("__Unknowns__")} else {
      #colnames(tmp)<-paste("",i,sep="")
      colnames(tmp)<-gsub(";+$","",paste(sapply(OTU_taxonomy[OTU_taxonomy[,which_level]==i,][1,1:which(colnames(OTU_taxonomy)==which_level)],as.character),collapse=";"))
    }
    if(is.null(new_abund_table)){new_abund_table<-tmp} else {new_abund_table<-cbind(tmp,new_abund_table)}
  }
}

new_abund_table<-as.data.frame(as(new_abund_table,"matrix"))
abund_table<-new_abund_table
#/COLLATE OTUS AT A PARTICULAR LEVEL#######################################

occ_threshold<-function (m, threshold, max_absent = 0) 
{
  occs <- colSums(m > max_absent)
  goodspecies <- occs >= threshold
  return(m[, goodspecies])
}

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
      res<-NULL
      if(length(unique(mt2[,j]))==2){
        if(class(mt2[,j])!="factor"){
          mt2[,j]<-as.factor(as.character(mt2[,j]))
        }
        res<-coda_glmnet(x=at2,y=mt2[,j])
        
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
        
        
      } else {
        res<-coda_glmnet(x=at2,y=mt2[,j],lambda="lambda.min",showPlots=FALSE)
      }
      pdf(paste("Signature_plot_",label,"_",i,"_",j,".pdf",sep=""),height=max(3,ceiling(length(res$taxa.num)/4)*height_adjustment),width=20)
      plot(res$`signature plot`)
      dev.off()
      pdf(paste("Predictions_plot_",label,"_",i,"_",j,".pdf",sep=""),height=5,width=10)
      plot(res$`predictions plot`)
      dev.off()
    }
  }
}, error = function(e) {
  message(conditionMessage(e))
})
}

