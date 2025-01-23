#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script for calculation of beta dispersion

# library(phyloseq)
# library(vegan)
# library(ggplot2)
# library(ape)
# library(phangorn)

beta_dispersion_analysis <- function(physeq) {
  # TODO: remove
  grouping_column <- "Groups"

  # Process input data
  abund_table <- phyloseq::otu_table(physeq)
  meta_table <- phyloseq::sample_data(physeq)
  OTU_taxonomy <- phyloseq::tax_table(physeq)
  OTU_tree <- phyloseq::phy_tree(physeq)

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
}
