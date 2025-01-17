#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script for calculation of beta dispersion

library(phyloseq)
library(vegan)
library(ggplot2)
library(ape)
library(phangorn)

#PARAMETERS ###########################
which_distance<-"bray" #bray unifrac wunifrac
physeq<-import_biom("../../data/VSEARCH/feature_w_tax.biom")
meta_table<-read.csv("../../data/meta_table.csv",header=T,row.names=1)
#Load the tree using ape package
OTU_tree <- read.tree("../../data/VSEARCH/tree.nwk")
OTU_tree$tip.label<-gsub("'","",OTU_tree$tip.label)
height_image=12
width_image=25
use_provided_colors=TRUE
colours<-c("#F8766D","#F8766D","#A3A500","#A3A500","#00BF7D","#00BF7D","#00B0F6","#00B0F6","#E76BF3","#E76BF3")
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

# Hypothesis 1
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

#PARAMETERS CHANGE THE GROUPING COLUMN AS YOU DESIRE############################

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

#I am using a function called midpoint() from phangorn package to root the tree in the middle
physeq<-merge_phyloseq(phyloseq(OTU, TAX),SAM,OTU_tree)

#Get all pair-wise combinations
s<-combn(unique(as.character(meta_table$Groups)),2)

#Go through all pair-wise combinations
df<-NULL
for (i in 1:dim(s)[2]){
  #Get a pruned subset of samples and taxa
  physeq_subset<-prune_samples(sample_data(physeq)$Groups %in% c(s[1,i],s[2,i]),physeq)
  physeq_subset<-prune_taxa(taxa_sums(physeq_subset)>0,physeq_subset)
  
  mod<-NULL
  if(which_distance=="bray"){
    mod<-betadisper(phyloseq::distance(physeq_subset,method="bray"),sample_data(physeq_subset)$Groups,type="centroid")
  } else if(which_distance=="wunifrac") {
    mod<-betadisper(phyloseq::distance(physeq_subset,method="wunifrac"),sample_data(physeq_subset)$Groups,type="centroid")
  } else if(which_distance=="unifrac"){
    mod<-betadisper(phyloseq::distance(physeq_subset,method="unifrac"),sample_data(physeq_subset)$Groups,type="centroid")
  }
  
  dist<-mod[["distances"]]
  df2<-data.frame(row.names=names(dist),Value=dist,meta_table[names(dist),"Groups",drop=F],Comparison=paste(s[1,i],"-",s[2,i]))
  p.value<-summary(aov(Value~ Groups,data=df2))[[1]][["Pr(>F)"]][1]
  title_string<-paste("p =",sprintf("%.5g",p.value),cut(p.value,breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")))
  df2$p.value<-p.value
  df2$Comparison<-paste(df2$Comparison," (",title_string,")",sep="")
  if(is.null(df)){df<-df2}else{df<-rbind(df,df2)}
}

#Prune df to those pair-wise comparisons that are significant
df<-df[df$p.value<=0.05,]
if(dim(df)[1]>0){
  q<-ggplot(df,aes(Groups,Value,colour=Groups))+ylab("Distance to Centroid")
  q<-q+geom_boxplot()+geom_jitter(position = position_jitter(height = 0, width=0))
  #q<-q+geom_point(size=5,alpha=0.2)
  q<-q+theme_bw()
  q<-q+facet_wrap(. ~ Comparison, drop=TRUE,scales="free",nrow=1)
  q<-q+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+theme(strip.text.x = element_text(size = 16, colour = "black", angle = 90))

  if(use_provided_colors){
    q<-q+scale_color_manual("Groups",values=colours)
  }
  pdf(paste("betadisper_",which_distance,"_",label,".pdf",sep=""),width=width_image,height=height_image)
  print(q)
  dev.off()
}
