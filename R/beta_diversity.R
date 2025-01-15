#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script for calculation of beta diversity (Principle Coordinate Analysis)

#Sometimes bray-curtis distance doesn't work with the newer version of phyloseq
#install it using:
#library("devtools")
#install_github("joey711/phyloseq")

rm(list=ls())
library(phyloseq)
library(vegan)
library(ggplot2)
library(ape)
library(phangorn)
library(stringr)
library(grid)

#PARAMETERS ###########################
which_level<-"Otus" #Phylum Class Order Family Genus Otus
kind <- "se" #sd se (sd is for drawing ellipse based on sd, se is for drawing ellipse based on standard errors)
which_distance<-"bray" #bray unifrac wunifrac
physeq<-import_biom("../../data/VSEARCH/feature_w_tax.biom")
meta_table<-read.csv("../../data/meta_table.csv",header=T,row.names=1)
#Load the tree using ape package
OTU_tree <- read.tree("../../data/VSEARCH/tree.nwk")
exclude_legends=FALSE
draw_mean_values_text=FALSE
point_size=6
point_opacity=1
draw_glow=FALSE
point_glow_opacity=0.2
point_glow_differential=3
mean_values_text_size=8
mean_values_text_opacity=1
draw_confidence_intervals=TRUE
draw_ellipses_and_not_polygons=TRUE
pairwise_connections_and_not_longitudinal_connections=TRUE
should_connections_end_in_arrows=FALSE
opacity_ellipses_polygons=0.8
linesize_ellipses_polygons=0.8
linetype_ellipses_polygons="solid" #blank solid dashed dotted dotdash longdash twodash
linking_samples_line_size=1
linking_samples_line_opacity=0.5
linking_samples_linetype="solid" #blank solid dashed dotted dotdash longdash twodash
legend_text_size=20
legend_title_size=20
axis_title_size=30
text_size=30
axis_text_size=30
height_image=15
width_image=15
use_provided_colors=FALSE
colours <- c("#ffa172", "#81fc76", "#68aeff","#c320d8","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00",grey.colors(1000));
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
#and then allocate a column to meta_table$Groups that you want to use, meta_table$Connections to connect them
#additionally if you provide a second meta_table$Subconnections, it will connect the group averages
#You can use meta_table$Type to assign shape, meta_table$Type2 to have multiple ellipses within the meta_Table$Groups variable, and PERMANOVA_variables to give variables

# Hypothesis 1
# label="Hypothesis1"
# meta_table<-meta_table[!meta_table$Within_Dam_Sample %in% c("Piezometer"),]
# meta_table$Groups<-as.factor(as.character(meta_table$Within_Dam_Sample))
# meta_table$Type<-as.character(meta_table$Sample_Time)
# meta_table$Type<-factor(meta_table$Type,levels=c("May","July"))
# meta_table$Type2<-as.character(meta_table$Sample_Time)
# meta_table$Type2<-factor(meta_table$Type2,levels=c("May","July"))
# meta_table$Connections<-as.character(meta_table$Groups)
# meta_table$Subconnections<-meta_table$Type2
# PERMANOVA_variables<-c("Sand_Dam","Sample_Time","Within_Dam_Sample")


# Hypothesis 1
label="Hypothesis2"
meta_table<-meta_table[!meta_table$Within_Dam_Sample %in% c("Piezometer"),]
meta_table$Groups<-factor(as.character(meta_table$Sample_Time))
meta_table$Type<-as.factor(as.character(meta_table$Within_Dam_Sample))
meta_table$Type2<-as.factor(as.character(meta_table$Within_Dam_Sample))
#meta_table$Connections<-as.character(meta_table$Groups)
#meta_table$Subconnections<-meta_table$Type2
PERMANOVA_variables<-c("Sand_Dam","Sample_Time","Within_Dam_Sample")



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

#Convert the data to phyloseq format
OTU = otu_table(as.matrix(abund_table), taxa_are_rows = FALSE)
TAX = tax_table(as.matrix(OTU_taxonomy))
SAM = sample_data(meta_table)

physeq<-NULL
if(which_level=="Otus"){
  #physeq<-merge_phyloseq(phyloseq(OTU, TAX),SAM,midpoint(OTU_tree))
  physeq<-merge_phyloseq(phyloseq(OTU, TAX),SAM,OTU_tree)
} else {
  physeq<-merge_phyloseq(phyloseq(OTU),SAM)
}



#Reference: http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
#Data frame df_ell contains values to show ellipses. It is calculated with function veganCovEllipse which is hidden in vegan package. This function is applied to each level of NMDS (group) and it uses also function cov.wt to calculate covariance matrix.
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#coloring function  
gg_color_hue<-function(n){
  hues=seq(15,375,length=n+1)
  hcl(h=hues,l=65,c=100)[1:n]
}

sol<-NULL
if(which_distance=="bray"){
  sol<-cmdscale(phyloseq::distance(physeq,"bray"),eig=T)
} else if(which_distance=="wunifrac" & which_level=="Otus") {
  sol<-cmdscale(phyloseq::distance(physeq,"wunifrac"),eig=T)
} else if(which_distance=="unifrac" & which_level=="Otus"){
  sol<-cmdscale(phyloseq::distance(physeq,"unifrac"),eig=T)
}

#Check to see if the meta_table$Type2 is assigned or not
if(is.null(meta_table$Type2)){meta_table$Type2<-meta_table$Groups}

if(!is.null(sol)){
  PCOA=data.frame(x=sol$points[,1],y=sol$points[,2],meta_table)

  plot.new()
  ord<-ordiellipse(sol, interaction(meta_table$Groups,meta_table$Type2),display = "sites", kind = kind, conf = 0.95, label = T)
  dev.off()


  #Generate ellipse points
  df_ell <- data.frame()
  for(h in levels(PCOA$Groups)){
    for(g in levels(PCOA$Type2)){
      if(paste(h,".",g,sep="")!="" && (paste(h,".",g,sep="") %in% names(ord))){
        tryCatch(df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCOA[PCOA$Groups==h & PCOA$Type2==g,],
                                                      veganCovEllipse(ord[[paste(h,".",g,sep="")]]$cov,ord[[paste(h,".",g,sep="")]]$center,ord[[paste(h,".",g,sep="")]]$scale)))
                                   ,Groups=h,Type2=g)),error=function(e) NULL)
   }
  }
  }
  
  if (sum(dim(df_ell))>0){
    colnames(df_ell)<-c("x","y","Groups","Type2")
  }
  df_ell$Groups<-factor(df_ell$Groups,levels=levels(PCOA$Groups))
  df_ell$Type2<-factor(df_ell$Type2,levels=levels(PCOA$Type2))
  
  #Generate mean values from PCOA plot grouped on
  PCOA.mean=aggregate(PCOA[,1:2],list(group=PCOA$Groups),mean)

  #Connecting samples based on lines with meta_table$Connections and meta_table$Subconnections ###########
  PCOA_lines<-NULL
  
  #Check if meta_table$Subconnections exists
  if(!is.null(meta_table$Subconnections)){
    #Step 1, populate Connections_IDs and get rid of singletons using PCOA$Connections
    Connections_IDs<-as.character(PCOA$Connections)
    Connections_IDs<-names(table(Connections_IDs)[table(Connections_IDs)>1])
    Subconnections_IDs_mask<-as.character(sapply(Connections_IDs,function(x){if(length(unique(PCOA[PCOA$Connections==x,"Subconnections"]))<2) x else "___APPROVED__"}))
    #Step 2, loop through each Connection_IDs and see we can find multiple PCOA$Subconnections and then filter Connections_IDs further
    Connections_IDs<-Connections_IDs[!Connections_IDs %in% Subconnections_IDs_mask]
    Connections_IDs<-sort(Connections_IDs)
    if(length(Connections_IDs)>0){
      for(p in Connections_IDs){
        Subconnections_IDs<-unique(PCOA[PCOA$Connections==p,"Subconnections"])
        if(is.factor(Subconnections_IDs)){
          Subconnections_IDs<-as.character(Subconnections_IDs)
        } else {
            Subconnections_IDs<-sort(Subconnections_IDs)
        }
        S<-NULL
        if(pairwise_connections_and_not_longitudinal_connections){
        #Get pair-wise combinations from Subconnections_IDs
         S<-combn(Subconnections_IDs,2)
        } else{
        #Get longitudinal connections from Subconnections_IDs
          S<-t(cbind(Subconnections_IDs[-length(Subconnections_IDs)],Subconnections_IDs[-1]))
        }
        for(ii in 1:ncol(S)){
          tmp<-data.frame(t(colMeans(PCOA[PCOA$Connections==p & PCOA$Subconnections==S[1,ii],c("x","y"),drop=F])),t(colMeans(PCOA[PCOA$Connections==p & PCOA$Subconnections==S[2,ii],c("x","y"),drop=F])),p)
          colnames(tmp)<-c("xfrom","yfrom","xto","yto","ID")
          if(is.null(PCOA_lines)){PCOA_lines<-tmp} else {PCOA_lines<-rbind(PCOA_lines,tmp)}
        }
      }
    }
  } else {
    #To connect lines between samples we need to extract connections
    Connections_IDs<-as.character(PCOA$Connections)
    #Next we filter out Connections_IDs that are singletons and also uniquify them
    Connections_IDs<-names(table(Connections_IDs)[table(Connections_IDs)>1])
    Connections_IDs<-sort(Connections_IDs)
    if(length(Connections_IDs)>0){
      #We iterate through the IDs one at a time
      for(p in Connections_IDs){
        rownames_list<-rownames(PCOA[PCOA$Connections %in% p,,drop=F])
        S<-NULL
        if(pairwise_connections_and_not_longitudinal_connections){
          #Get pair-wise combinations from rownames_list
          S<-combn(rownames_list,2)
        } else{
          #Get longitudinal connections from rownames_list
          S<-t(cbind(rownames_list[-length(rownames_list)],rownames_list[-1]))
        }
        for(ii in 1:ncol(S)){
          tmp<-cbind(PCOA[S[1,ii],c("x","y")],PCOA[S[2,ii],c("x","y")],p)
          colnames(tmp)<-c("xfrom","yfrom","xto","yto","ID")
          if(is.null(PCOA_lines)){PCOA_lines<-tmp} else {PCOA_lines<-rbind(PCOA_lines,tmp)}
        }
      }
    }
  }
  #/#Connecting samples based on lines with meta_table$Connections and meta_table$Subconnections ###########
  
  
  cols=gg_color_hue(length(unique(PCOA$Groups)))
  
  p<-ggplot(data=PCOA,aes(x,y,colour=Groups))
  if(!"Type" %in% colnames(meta_table)){
    p<-p + geom_point(aes(PCOA$x,PCOA$y,colour=PCOA$Groups),inherit.aes=F,alpha=point_opacity,size=point_size)
  } else{
    p<-p + geom_point(aes(PCOA$x,PCOA$y,colour=PCOA$Groups, shape=PCOA$Type),inherit.aes=F,alpha=point_opacity,size=point_size)
    p<-p+scale_shape("Type")
  }
  
  if(draw_glow){
    p<-p + geom_point(alpha=point_glow_opacity,size = point_size+point_glow_differential,show.legend=FALSE)
  }
  
  p<-p+theme_bw()
  if(draw_mean_values_text){
    p<-p+ annotate("text",x=PCOA.mean$x,y=PCOA.mean$y,label=PCOA.mean$group,size=mean_values_text_size,colour=cols,alpha=mean_values_text_opacity,vjust=0.3)
  }
  if (sum(dim(df_ell))>0){
  if(draw_confidence_intervals){
    if(identical(levels(meta_table$Groups), levels(meta_table$Type2))){
      if(draw_ellipses_and_not_polygons){
        p<-p+ geom_path(data=df_ell, aes(x=x, y=y), size=linesize_ellipses_polygons, linetype=linetype_ellipses_polygons,alpha=opacity_ellipses_polygons, show.legend = FALSE)
      } else {
        p<-p+ geom_polygon(data=df_ell, aes(x=x, y=y,fill=Groups), size=linesize_ellipses_polygons, linetype=linetype_ellipses_polygons,alpha=opacity_ellipses_polygons, show.legend = FALSE)
      }
    } else {
      for (q in unique(as.numeric(meta_table$Groups))){
        for(r in unique(as.numeric(meta_table$Type2))){
          if(draw_ellipses_and_not_polygons){
            p<-p+ geom_path(data=df_ell[as.numeric(df_ell$Groups)==q & as.numeric(df_ell$Type2)==r,], aes(x=x, y=y), size=linesize_ellipses_polygons, colour=cols[q],linetype=r,alpha=opacity_ellipses_polygons, show.legend = FALSE)
          } else {
            p<-p+ geom_polygon(data=df_ell[as.numeric(df_ell$Groups)==q & as.numeric(df_ell$Type2)==r,], aes(x=x, y=y,fill=Groups), size=linesize_ellipses_polygons, colour=cols[q],linetype=r,alpha=opacity_ellipses_polygons, show.legend = FALSE)
          }
          }
      }
    }
  }
  }
  if(exclude_legends){
    p<-p+guides(colour=FALSE)
  }
  
  #only draw lines connecting dots if the lines are available
  if(!is.null(PCOA_lines)){
    arrow<-NULL
    if(should_connections_end_in_arrows){
      arrow=arrow(length=unit(0.2,"inches"))
    }
    p<-p+geom_segment(data=PCOA_lines,inherit.aes=FALSE,aes(x=xfrom,y=yfrom,xend=xto,yend=yto),colour="grey20",size=linking_samples_line_size,alpha=linking_samples_line_opacity,linetype=linking_samples_linetype, show.legend = FALSE,arrow=arrow)
  }
  
  
  p<-p+xlab(paste("Dim1 (",sprintf("%.4g",(sol$eig[1]/sum(sol$eig))*100),"%)",sep=""))+ylab(paste("Dim2 (",sprintf("%.4g",(sol$eig[2]/sum(sol$eig))*100),"%)",sep=""))
  
  p<-p+theme(legend.title=element_text(size=legend_title_size),
             legend.text=element_text(size=legend_text_size),
             text = element_text(size=text_size),
             axis.text=element_text(size=axis_text_size),
             axis.title=element_text(size=axis_title_size))
  
 
  if(use_provided_colors){
    p<-p+scale_color_manual("Groups",values=colours)
    p<-p+scale_fill_manual("Groups",values=colours)
  }
  
  
   pdf(paste("PCOA_",which_distance,"_",which_level,"_",label,".pdf",sep=""),width=width_image,height=height_image)
  print(p)
  dev.off()
  
  dist<-phyloseq::distance(physeq,which_distance)
  capture.output(adonis(as.formula(paste("dist ~",paste(PERMANOVA_variables,collapse="+"))), data=meta_table[rownames(otu_table(physeq)),]),file=paste("ADONIS_",which_distance,"_",which_level,"_",label,".txt",sep=""))
  
  
}