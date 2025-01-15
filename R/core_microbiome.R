#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script for core analysis
#Reference: https://microbiome.github.io/tutorials/Core.html
#v2.0 fixed reordering bug

rm(list=ls())

library(phyloseq)
library(ggplot2)
library(viridis)
library(microbiome)
library(RColorBrewer)
library(cowplot)


#PARAMETERS ###########################
physeq<-import_biom("../../data/VSEARCH/feature_w_tax.biom")
meta_table<-read.csv("../../data/meta_table.csv",header=T,row.names=1)
which_level<-"Phylum" #Phylum Class Order Family Genus Otus
height_image_heatmap=16
legend_text_size=6
legend_title_size=8
what_detection="absolute" #absolute relative
minimum_prevalence=0.85
text_size=12
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

label="Hypothesis1"
meta_table<-meta_table[!meta_table$Within_Dam_Sample %in% c("Piezometer"),]
meta_table$Groups<-as.factor(as.character(meta_table$Within_Dam_Sample))
meta_table$Type<-as.factor(as.character(meta_table$Sand_Dam))
meta_table$Type2<-as.character(meta_table$Sample_Time)
meta_table$Type2<-factor(meta_table$Type2,levels=c("May","July"))

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




#Convert the data to phyloseq format
OTU = otu_table(as.matrix(abund_table), taxa_are_rows = FALSE)
TAX = tax_table(as.matrix(OTU_taxonomy))
SAM = sample_data(meta_table)

physeq<-NULL
if(which_level=="Otus"){
  physeq<-merge_phyloseq(phyloseq(OTU, TAX),SAM)
} else {
  physeq<-merge_phyloseq(phyloseq(OTU),SAM)
}

# keep only taxa with positive sums
pseq.2 <- prune_taxa(taxa_sums(physeq) > 0, physeq)

# Calculate compositional version of the data
# (relative abundances)
pseq.rel <- microbiome::transform(pseq.2, "compositional")

prevalences <- seq(.05, 1, .05)
detections<-NULL
pseq_to_plot<-NULL
if(what_detection=="relative"){
  #Detection with Relative Abundances
  detections <- 10^seq(log10(1e-3), log10(.2), length = 20)
  pseq_to_plot<-pseq.rel
  
} else if(what_detection=="absolute"){
  #Detection with Absolute Count
  detections <- 10^seq(log10(1), log10(max(abundances(pseq.2))/10), length = 20)
  detections <- round(detections)
  pseq_to_plot<-pseq.2
}



datacore <- plot_core(pseq_to_plot, plot.type = "heatmap", 
                         prevalences = prevalences,
                         detections = detections,
                         #colours = gray(seq(1,0,length=20)),
                         colours=rev(brewer.pal(5, "Spectral")),
                          min.prevalence = minimum_prevalence, horizontal = TRUE)

if(which_level=="Otus"){
  # get the data used for plotting 
  df <- datacore$data 


  # get the list of OTUs
  list <- df$Taxa 

  # get the taxonomy data
  tax <- tax_table(pseq.2)
  tax <- as.data.frame(tax)

  # add the OTus to last column
  tax$OTU <- rownames(tax)#OTU

  # select taxonomy of only 
  # those OTUs that are used in the plot
  tax2 <- dplyr::filter(tax, rownames(tax) %in% list) 

  # We will merege all the column into one ##c="OTU"
  tax.unit <- tidyr::unite(tax2, Taxa_level,c("OTU",names(tax2)[-length(names(tax2))]), sep = "_;", remove = FALSE)
  tax.unit$Taxa_level <- gsub(pattern="_;",replacement=";", tax.unit$Taxa_level)
  tax.unit$Taxa_level <- gsub(pattern=";+$",replacement="", tax.unit$Taxa_level)
  # add this new information into the plot data df
  df$Taxa <- tax.unit$Taxa_level

  #Refactorize df$Taxa with respect to list
  df$Taxa<-factor(df$Taxa,levels=as.character(sapply(levels(list),function(x){unique(df$Taxa)[grep(paste(x,";",sep=""),unique(df$Taxa))]})))


  # replace the data in the plot object
  datacore$data <- df
} else {
  datacore$data<-datacore$data[datacore$data$Taxa!="__Unknowns__",]}

#Now plot it
pdf(paste("Core_Heatmap_",label,"_",which_level,"_",what_detection,".pdf",sep=""),width=ceiling(nlevels(datacore$data$Taxa)/5)+2,height=height_image_heatmap)
p<-datacore 
p<-p+theme_bw()  + theme_cowplot(font_size = text_size) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x  = element_text(angle=90,vjust=0.5,hjust = 1),legend.title=element_text(size=legend_title_size),
        legend.text=element_text(size=legend_text_size))



if(what_detection=="relative"){
  p<-p + xlab("Relative Abundance (%)") 
}  
print(p)
dev.off()



