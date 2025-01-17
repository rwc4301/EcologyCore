#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script for finding log2 fold different species using DESeq2

library(phyloseq)
library(vegan)
library(ggplot2)
library(plyr)
library(DESeq2)
library(stringr)

#PARAMETERS ###########################
which_level<-"Otus" #Phylum Class Order Family Genus Otus
sig = 0.05
fold = 2
physeq<-import_biom("../../data/VSEARCH/feature_w_tax.biom")
meta_table<-read.csv("../../data/meta_table2.csv",header=T,row.names=1)
height_image=15
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
#and then allocate a column to meta_table$Groups that you want to use, and PERMANOVA_variables to give variables

meta_table<-meta_table[!meta_table$Within_Dam_Sample %in% c("Piezometer"),]
meta_table$Groups<-paste(meta_table$Within_Dam_Sample,meta_table$Sample_Time)


#Hypothesis 1
labels="Excavated_scoop_hole"
meta_table<-meta_table[meta_table$Groups %in% c("Excavated Scoop Hole  May","Excavated Scoop Hole  July"),]
meta_table$Groups<-factor(meta_table$Groups,levels=c("Excavated Scoop Hole  May","Excavated Scoop Hole  July"))

# #Hypothesis 2
# labels="Community_scoop_hole"
# meta_table<-meta_table[meta_table$Groups %in% c("Community Scoop Hole May","Community Scoop Hole July"),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Community Scoop Hole May","Community Scoop Hole July"))

# #Hypothesis 3
# labels="Hand_pump"
# meta_table<-meta_table[meta_table$Groups %in% c("Hand Pump May","Hand Pump July"),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Hand Pump May","Hand Pump July"))

# #Hypothesis 4
# labels="Open_well"
# meta_table<-meta_table[meta_table$Groups %in% c("Open Well May","Open Well July"),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Open Well May","Open Well July"))

# #Hypothesis 5
# labels="Downstream"
# meta_table<-meta_table[meta_table$Groups %in% c("Downstream May","Downstream July"),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Downstream May","Downstream July"))

# #Hypothesis 6
# labels="Excavated Scoop Hole  May_Community Scoop Hole May"
# meta_table<-meta_table[meta_table$Groups %in% c("Excavated Scoop Hole  May","Community Scoop Hole May"),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Excavated Scoop Hole  May","Community Scoop Hole May"))

# #Hypothesis 7
# labels="Excavated Scoop Hole  May_Open Well May"
# meta_table<-meta_table[meta_table$Groups %in% c("Excavated Scoop Hole  May","Open Well May"  ),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Excavated Scoop Hole  May","Open Well May"  ))

# #Hypothesis 8
# labels="Excavated Scoop Hole  May_Hand Pump May"
# meta_table<-meta_table[meta_table$Groups %in% c("Excavated Scoop Hole  May","Hand Pump May" ),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Excavated Scoop Hole  May","Hand Pump May" ))

# #Hypothesis 9
# labels="Excavated Scoop Hole  May_Downstream May"
# meta_table<-meta_table[meta_table$Groups %in% c("Excavated Scoop Hole  May","Downstream May" ),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Excavated Scoop Hole  May","Downstream May" ))

# #Hypothesis 10
# labels="Community Scoop Hole May_Open Well May"
# meta_table<-meta_table[meta_table$Groups %in% c("Community Scoop Hole May","Open Well May" ),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Community Scoop Hole May","Open Well May" ))

# #Hypothesis 11
# labels="Community Scoop Hole May_Hand Pump May"
# meta_table<-meta_table[meta_table$Groups %in% c("Community Scoop Hole May","Hand Pump May"),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Community Scoop Hole May","Hand Pump May"))

# #Hypothesis 12
# labels="Community Scoop Hole May_Downstream May"
# meta_table<-meta_table[meta_table$Groups %in% c("Community Scoop Hole May","Downstream May"),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Community Scoop Hole May","Downstream May") )

# #Hypothesis 13
# labels="Open Well May_Hand Pump May"
# meta_table<-meta_table[meta_table$Groups %in% c("Open Well May","Hand Pump May"),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Open Well May","Hand Pump May"))

# #Hypothesis 14
# labels="Open Well May_Downstream May"
# meta_table<-meta_table[meta_table$Groups %in% c("Open Well May","Downstream May"),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Open Well May","Downstream May"))

# #Hypothesis 15
# labels="Hand Pump May_Downstream May"
# meta_table<-meta_table[meta_table$Groups %in% c("Hand Pump May","Downstream May"),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Hand Pump May","Downstream May"))

# #Hypothesis 16
# labels="Community Scoop Hole July_Excavated Scoop Hole  July"
# meta_table<-meta_table[meta_table$Groups %in% c("Community Scoop Hole July","Excavated Scoop Hole  July"),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Community Scoop Hole July","Excavated Scoop Hole  July"))

# #Hypothesis 17
# labels="Community Scoop Hole July_Hand Pump July"
# meta_table<-meta_table[meta_table$Groups %in% c("Community Scoop Hole July","Hand Pump July"),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Community Scoop Hole July","Hand Pump July"))

# #Hypothesis 18
# labels="Community Scoop Hole July_Open Well July"
# meta_table<-meta_table[meta_table$Groups %in% c("Community Scoop Hole July","Open Well July"),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Community Scoop Hole July","Open Well July"))

# #Hypothesis 19
# labels="Community Scoop Hole July_Downstream July"
# meta_table<-meta_table[meta_table$Groups %in% c("Community Scoop Hole July","Downstream July"),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Community Scoop Hole July","Downstream July"))

# #Hypothesis 20
# labels="Excavated Scoop Hole  July_Hand Pump July"
# meta_table<-meta_table[meta_table$Groups %in% c("Excavated Scoop Hole  July", "Hand Pump July"),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Excavated Scoop Hole  July", "Hand Pump July"))

# #Hypothesis 21
# labels="Excavated Scoop Hole  July_Open Well July"
# meta_table<-meta_table[meta_table$Groups %in% c("Excavated Scoop Hole  July","Open Well July" ),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Excavated Scoop Hole  July","Open Well July" ))

# #Hypothesis 22
# labels="Excavated Scoop Hole  July_Downstream July"
# meta_table<-meta_table[meta_table$Groups %in% c("Excavated Scoop Hole  July", "Downstream July"),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Excavated Scoop Hole  July", "Downstream July"))

# #Hypothesis 23
# labels="Hand Pump July_Open Well July"
# meta_table<-meta_table[meta_table$Groups %in% c("Hand Pump July","Open Well July"),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Hand Pump July","Open Well July"))

# #Hypothesis 24
# labels="Hand Pump July_Downstream July"
# meta_table<-meta_table[meta_table$Groups %in% c("Hand Pump July","Downstream July"),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Hand Pump July","Downstream July"))

# #Hypothesis 25
# labels="Open Well July_Downstream July"
# meta_table<-meta_table[meta_table$Groups %in% c("Open Well July","Downstream July"),]
# meta_table$Groups<-factor(meta_table$Groups,levels=c("Open Well July","Downstream July") )




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
    tmp<-data.frame(rowSums(abund_table[,rownames(OTU_taxonomy)[OTU_taxonomy[,which_level]==i],drop=F]))
    if(i==""){colnames(tmp)<-c("__Unknowns__")} else {colnames(tmp)<-paste("",i,sep="")}
    if(is.null(new_abund_table)){new_abund_table<-tmp} else {new_abund_table<-cbind(tmp,new_abund_table)}
  }
}

new_abund_table<-as.data.frame(as(new_abund_table,"matrix"))
abund_table<-new_abund_table
#/COLLATE OTUS AT A PARTICULAR LEVEL#######################################


#We will convert our table to DESeqDataSet object
countData = round(as(abund_table, "matrix"), digits = 0)
# We will add 1 to the countData otherwise DESeq will fail with the error:
# estimating size factors
# Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  :
# every gene contains at least one zero, cannot compute log geometric means
countData<-(t(countData+1)) 

dds <- DESeqDataSetFromMatrix(countData, meta_table, as.formula(~ Groups))

#Reference:https://github.com/MadsAlbertsen/ampvis/blob/master/R/amp_test_species.R
#Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
#Some reason this doesn't work: data_deseq_test = DESeq(dds, test="wald", fitType="parametric")
data_deseq_test = DESeq(dds)

## Extract the results
res = results(data_deseq_test, cooksCutoff = FALSE)
res_tax = cbind(as.data.frame(res), as.matrix(countData[rownames(res), ]), OTU = rownames(res))

plot.point.size = 2
label=F
tax.display = NULL
tax.aggregate = "OTU"

res_tax_sig = subset(res_tax, padj < sig & fold < abs(log2FoldChange))

res_tax_sig <- res_tax_sig[order(res_tax_sig$padj),]

## Plot the data
### MA plot
res_tax$Significant <- ifelse(rownames(res_tax) %in% rownames(res_tax_sig) , "Yes", "No")
res_tax$Significant[is.na(res_tax$Significant)] <- "No"
p1 <- ggplot(data = res_tax, aes(x = baseMean, y = log2FoldChange, color = Significant)) +
  geom_point(size = plot.point.size) +
  scale_x_log10() +
  scale_color_manual(values=c("black", "red")) +
  labs(x = "Mean abundance", y = "Log2 fold change")+theme_bw()
if(label == T){
  if (!is.null(tax.display)){
    rlab <- data.frame(res_tax, Display = apply(res_tax[,c(tax.display, tax.aggregate)], 1, paste, collapse="; "))
  } else {
    rlab <- data.frame(res_tax, Display = res_tax[,tax.aggregate])
  }
  p1 <- p1 + geom_text(data = subset(rlab, Significant == "Yes"), aes(label = Display), size = 4, vjust = 1)
}
pdf(paste("NB_MA_",which_level,"_",labels,".pdf",sep=""))
print(p1)
dev.off()

res_tax_sig_abund = cbind(as.data.frame(countData[rownames(res_tax_sig), ]), OTU = rownames(res_tax_sig), padj = res_tax[rownames(res_tax_sig),"padj"]) 

#Apply normalisation (either use relative or log-relative transformation)
data<-log((abund_table+1)/(rowSums(abund_table)+dim(abund_table)[2]))
data<-as.data.frame(data)

#Now we plot taxa significantly different between the categories
df<-NULL
sig_otus<-res_tax[rownames(res_tax_sig),"OTU"]

for(i in sig_otus){
  tmp<-NULL
  if(which_level=="Otus"){
    tmp<-data.frame(data[,i],meta_table$Groups,rep(paste(paste(i,gsub(";+$","",paste(sapply(OTU_taxonomy[i,],as.character),collapse=";")))," padj = ",sprintf("%.5g",res_tax[i,"padj"]),sep=""),dim(data)[1]))
  } else {
    tmp<-data.frame(data[,i],meta_table$Groups,rep(paste(i," padj = ",sprintf("%.5g",res_tax[i,"padj"]),sep=""),dim(data)[1]))
  }
  
  if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)} 
}
colnames(df)<-c("Value","Groups","Taxa")

library(ggplot2)

p<-ggplot(df,aes(Groups,Value,colour=Groups))+ylab("Log-relative normalised")
p<-p+geom_boxplot(outlier.size = 0)+geom_jitter(position = position_jitter(height = 0, width=0),alpha=0.5,outlier.colour = NULL)+theme_bw()+
  facet_wrap( ~ Taxa , scales="free_x",nrow=1)
p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+theme(strip.text.x = element_text(size = 16, colour = "black", angle = 90))
pdf(paste("NB_significant_",which_level,"_",labels,".pdf",sep=""),width=ceiling((length(sig_otus)*80/200)+2.6),height=height_image)
print(p)
dev.off()

data_to_write<-res_tax_sig[,c("baseMean","log2FoldChange","pvalue","padj")]
tmp<-sapply(as.character(rownames(res_tax_sig)),function(x){aggregate(data[,x],by=list(meta_table$Groups),FUN=mean)[,2]})
tmp<-as.data.frame(t(tmp))
colnames(tmp)<-levels(meta_table$Groups)
data_to_write<-cbind(data_to_write,Upregulated=as.character(sapply(rownames(tmp),function(x){levels(meta_table$Groups)[which.max(tmp[x,])]})))

if(which_level=="Otus"){
  rownames(data_to_write)<-as.character(sapply(rownames(data_to_write),function(x) rep(paste(paste(x,gsub(";+$","",paste(sapply(OTU_taxonomy[x,],as.character),collapse=";")))))))
}
write.csv(data_to_write,paste("NB_significant_",which_level,"_",labels,".csv",sep=""))
