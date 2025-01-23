#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script for null modelling (ST/NST/MST)
#Authors: Umer, Anna
#v1.0 (All metrics are now being saved)

# library(phyloseq)
# library(vegan)
# library(ape)
# library(NST)
# library(ggplot2)

nst <- function(physeq) {
  # TODO: remove
  grouping_column <- "Groups"

  # Process input data
  abund_table <- phyloseq::otu_table(physeq)
  meta_table <- phyloseq::sample_data(physeq)
  OTU_taxonomy <- phyloseq::tax_table(physeq)
  OTU_tree <- phyloseq::phy_tree(physeq)

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
}

nst_plot <- function() {
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
}
