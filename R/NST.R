#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script for null modelling (ST/NST/MST)
#Authors: Umer, Anna
#v1.0 (All metrics are now being saved)

library(phyloseq)
library(vegan)
library(ape)
library(NST)
library(ggplot2)

#PARAMETERS ###########################
physeq<-import_biom("../../data/VSEARCH/feature_w_tax.biom")
meta_table<-read.csv("../../data/meta_table.csv",header=T,row.names=1)
number_of_randomizations=1000
distance_measure= "cao" #"manhattan" "mManhattan" "euclidean" "mEuclidean"  "canberra" "bray" "kulczynski" "jaccard" "gower" "altGower" "mGower" "morisita" "horn" "binomial" "chao" "cao"
abundance.weighted= FALSE #Logic, consider abundances or not (just presence/absence). default is TRUE.
#Jaccard with abundance.weighted=TRUE is called Ruzicka 
null_model="PF" #"EE" "EP" "EF" "PE" "PP" "PF" "FE" "FP" "FF" with details given below:
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
SES = TRUE #Logic, whether to calculate standardized effect size, which is (observed dissimilarity - mean of null dissimilarity)/standard deviation of null dissimilarity. default is FALSE.
RC = FALSE # Logic, whether to calculate modified Raup-Crick metric, which is percentage of null dissimilarity lower than observed dissimilarity x 2 - 1. default is FALSE.
NST_width=4
NST_height=8
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

colours <- c(
  "red",
  "blue",
  #Next colors are for lines mainly used in the PCoA script
  "#000080","#4876FF","#CAE1FF","#9FB6CD","#1E90FF","#00F5FF","#00C957",grey.colors(1000));
#PARAMETERS CHANGE THE GROUPING COLUMN AS YOU DESIRE############################

#Adjust abund_table to contain only those rows that got selected in the Hypothesis space
abund_table<-abund_table[rownames(meta_table),]
#After adjustment, get rid of OTUs that are all empty
abund_table<-abund_table[,colSums(abund_table)>0]
#Adjust OTU taxonomy
OTU_taxonomy<-OTU_taxonomy[colnames(abund_table),]

#Bug in tNST (abund_table should be of type "matrix")
abund_table<-as(abund_table,"matrix")

tnst=tNST(comm=abund_table, group=meta_table[,"Groups",drop=F], 
          rand=number_of_randomizations,
          dist.method = distance_measure,
          null.model=null_model,
          output.rand=TRUE, nworker=1,
          SES=SES, RC=RC)

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


pdf(paste("Stochasticity-Ratios_",null_model,"_",distance_measure,"_",as.character(number_of_randomizations),"_",as.character(SES),"_",as.character(RC),"_",as.character(abundance.weighted),"_",label,"_FIGURE",".pdf",sep=""),width=NST_width,height=NST_height)
p<-ggplot(df, aes(x=Groups, y=value, fill=Groups)) 
p<-p+geom_bar(stat="identity")+theme_minimal()
p<-p+geom_text(aes(label=sprintf("%0.2f", round(value, digits = 2))), vjust=-0.3, size=3.5)
p<-p+facet_wrap(~measure, strip.position="left", ncol=1,scales="free_y")
p<-p+scale_fill_manual(values=colours)
p<-p+ylim(0,1.1)
p<-p+ylab("Stochasticity Ratios (scaled to 1)")
p<-p+theme(strip.background = element_rect(fill="white"))+theme(panel.spacing = unit(2, "lines"),
                                                                axis.text.x = element_text(angle = 90, hjust = 1))
print(p)
dev.off()

#Perform PANOVA
nst.pova=nst.panova(nst.result=tnst, rand=number_of_randomizations)
#Convert group numbers to actual names
for(i in 1:nlevels(meta_table$Groups)){
  nst.pova$group1<-gsub(paste("^",i,"$",sep=""),levels(meta_table$Groups)[i],nst.pova$group1)
  nst.pova$group2<-gsub(paste("^",i,"$",sep=""),levels(meta_table$Groups)[i],nst.pova$group2)
}

write.csv(nst.pova,file=paste("Stochasticity-Ratios_",null_model,"_",distance_measure,"_",as.character(number_of_randomizations),"_",as.character(SES),"_",as.character(RC),"_",as.character(abundance.weighted),"_",label,"_PANOVA",".csv",sep=""))
write.csv(tnst$index.pair.grp,file=paste("Stochasticity-Ratios_",null_model,"_",distance_measure,"_",as.character(number_of_randomizations),"_",as.character(SES),"_",as.character(RC),"_",as.character(abundance.weighted),"_",label,"_PAIRWISE",".csv",sep=""))



