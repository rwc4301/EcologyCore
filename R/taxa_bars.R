#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script for proportional representation of N most abundant species

# library(phyloseq)
# library(vegan)
# library(plyr)
# library(ape)
# library(grid)

taxa_bars_analysis <- function(abund_table, meta_table) {
  #Apply proportion normalisation
  x<-abund_table/rowSums(abund_table)
  x<-x[,order(colSums(x),decreasing=TRUE)]


  taxa_list<-colnames(x)[1:min(dim(x)[2],N)]
  if (c("__Unknowns__") %in% taxa_list){
    taxa_list<-colnames(x)[1:min(dim(x)[2],N+1)]
    taxa_list<-taxa_list[!grepl("__Unknowns__",taxa_list)]
  }
  N<-length(taxa_list)

  #Generate a new table with everything added to Others
  new_x<-NULL
  if(N==dim(x)[2]){
    new_x<-data.frame(x[,colnames(x) %in% taxa_list])
  } else {
    new_x<-data.frame(x[,colnames(x) %in% taxa_list],Others=rowSums(x[,!colnames(x) %in% taxa_list]))
  }
  if(which_level=="Otus"){
    if(N==dim(x)[2]){
      colnames(new_x)<-c(paste(colnames(new_x),sapply(colnames(new_x),function(x) gsub(";+$","",paste(sapply(OTU_taxonomy[x,],as.character),collapse=";")))))
    } else {
      colnames(new_x)<-c(paste(colnames(new_x)[-(N+1)],sapply(colnames(new_x)[-(N+1)],function(x) gsub(";+$","",paste(sapply(OTU_taxonomy[x,],as.character),collapse=";")))),"Others")
    }
  }

  df<-NULL
  for (i in 1:dim(new_x)[2]){
    tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Groups=meta_table$Groups)
    if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
  }
  colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00",rainbow(300)[seq(1,300,6)]);

  #The first step to get levels is to find all the possible values in a given column in df$sample by outputing the following command on terminal:
  # cat(paste("levels=c(",paste(paste("\"",unique(as.character(df$Sample)),"\"",sep=""),collapse=","),")",sep=""))
  # then use df$Sample<-factor(as.character(df$Sample),levels=c()) list
}

taxa_bars_plot <- function() {
  p<-ggplot(df,aes(Sample,Value,fill=Taxa))+geom_bar(stat="identity")+facet_grid(. ~ Groups, drop=TRUE,scale="free",space="free_x",switch=switch_strip)
  p<-p+scale_fill_manual(values=colours[1:(N+1)],guide=guide_legend(ncol = how_many_columns_for_legend))
  p<-p+theme_bw()+ylab("Proportions")
  p<-p+ scale_y_continuous(expand = c(0.02,0))+theme(strip.background = element_rect(fill="gray85"),
                                                     panel.spacing = unit(0.3, "lines"),
                                                     legend.position = "bottom",
                                                     strip.text = element_text(size=strip_text_size,angle=90),
                                                     legend.text=element_text(size=legend_text_size),
                                                     text = element_text(size=text_size),
                                                     axis.text=element_text(size=axis_text_size),
                                                     axis.title=element_text(size=axis_title_size))

  if(reveal_sample_names){
    p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  } else {
    p<-p+theme(axis.text.x=element_blank(),axis.ticks=element_blank())
  }



  pdf(paste("TAXAplot_",which_level,"_",label,".pdf",sep=""),height=height_image,width=width_image)
  print(p)
  dev.off()
}
