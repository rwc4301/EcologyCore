collate_taxonomy <- function(abund_table, taxa_table, rank = "Otus") {
  new_abund_table<-NULL
  if(which_level=="Otus"){
    new_abund_table<-abund_table
  } else {
    list<-unique(OTU_taxonomy[,which_level])
    new_abund_table<-NULL
    for(i in list){
      tmp<-data.frame(rowSums(abund_table[,rownames(OTU_taxonomy)[OTU_taxonomy[,which_level]==i]]))
      if(i==""){colnames(tmp)<-c("__Unknowns__")} else {colnames(tmp)<-paste("",i,sep="")}
      if(is.null(new_abund_table)){new_abund_table<-tmp} else {new_abund_table<-cbind(tmp,new_abund_table)}
    }
  }
  new_abund_table<-as.data.frame(as(new_abund_table,"matrix"))
  abund_table<-new_abund_table

  return(abund_table)
}
