# library(phyloseq)
# library(stringr)
# library(data.table)
# library(vegan)
# library(ggplot2)
# library(grid) #We need grid to draw the arrows

## Also works for picrust data - need to round abund_table

alpha_diversity <- function(abund_table, meta_table, grouping_column) {
  #Calculate Richness
  R<-vegan::rarefy(abund_table,min(rowSums(abund_table)))
  df_R<-data.frame(sample=names(R),value=R,measure=rep("Richness",length(R)))

  #Calculate Shannon entropy
  H<-vegan::diversity(abund_table)
  df_H<-data.frame(sample=names(H),value=H,measure=rep("Shannon",length(H)))

  #Calculate Simpson diversity index
  simp <- vegan::diversity(abund_table, "simpson")
  df_simp<-data.frame(sample=names(simp),value=simp,measure=rep("Simpson",length(simp)))

  #Calculate Fisher alpha
  alpha <- vegan::fisher.alpha(abund_table)
  df_alpha<-data.frame(sample=names(alpha),value=alpha,measure=rep("Fisher alpha",length(alpha)))

  #Calculate Pielou's evenness
  S <- vegan::specnumber(abund_table)
  J <- H/log(S)
  df_J<-data.frame(sample=names(J),value=J,measure=rep("Pielou's evenness",length(J)))

  #Uncomment to retain everything
  df<-rbind(df_R,df_H,df_simp,df_alpha,df_J)

  rownames(df)<-NULL

  #Incorporate categorical data in df
  df<-data.frame(df,meta_table[as.character(df$sample),])

  #To do anova, we will convert our data.frame to data.table

  #Since we can't pass a formula to data.table, I am creating
  #a dummy column .group. so that I don't change names in the formula
  dt<-data.table::as.data.table(data.frame(df,.group.=df[,grouping_column]))

  # pval<-data.table::as.data.table(dt)[, list(pvalue = sprintf("%.2g", tryCatch({
  #   pv <- summary(stats::aov(value ~ .group.))[[1]][["Pr(>F)"]][1]
  #   print(pv)
  #   pv
  # }, error = function(e) NULL))),
  # by = list(measure)
  # ]

  return(dt)
}

alpha_diversity_2 <- function(dt, pval, meta_table, grouping_column) {
  #I am also specifying a p-value cutoff for the ggplot2 strips
  pValueCutoff<-0.05

  #Filter out pvals that we don't want
  pval<-pval[!pval$pvalue=="",]
  pval<-pval[as.numeric(pval$pvalue)<=pValueCutoff,]

  #I am using sapply to generate significances for pval$pvalue using the cut function.
  pval$pvalue<-sapply(as.numeric(pval$pvalue),function(x){as.character(cut(x,breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))})

  df <- as.data.frame(dt)

  #### IT'S FUCKED!!!!!

  #Update df$measure to change the measure names if the grouping_column has more than three classes
  # if(length(unique(as.character(meta_table[,grouping_column])))>2){
  #   df$measure<-as.character(df$measure)
  #   if(dim(pval)[1]>0){
  #     for(i in seq(1:dim(pval)[1])){
  #       df[df$measure==as.character(pval[i,measure]),"measure"]=paste(as.character(pval[i,measure]),as.character(pval[i,pvalue]))
  #     }
  #   }
  #   df$measure<-as.factor(df$measure)
  # }

  #Get all possible combination of values in the grouping_column
  s<-combn(unique(as.character(df[,grouping_column])),2)

  #df_pw will store the pair-wise p-values
  df_pw<-NULL
  for(k in unique(as.character(df$measure))){
    #We need to calculate the coordinate to draw pair-wise significance lines
    #for this we calculate bas as the maximum value
    bas<-max(df[(df$measure==k),"value"])

    #Calculate increments as % of the maximum values
    inc<-0.05*(bas-min(df[(df$measure==k),"value"]))

    #Give an initial increment
    bas<-bas+inc
    for(l in 1:dim(s)[2]){

      tmp<-NULL
      #If it is paired-data we are interested in
      if(is.null(meta_table$Connections)){
        #Do a normal anova
        tmp<-c(k,s[1,l],s[2,l],bas,paste(sprintf("%.2g",tryCatch(summary(aov(as.formula(paste("value ~",grouping_column)),data=df[(df$measure==k) & (df[,grouping_column]==s[1,l] | df[,grouping_column]==s[2,l]),] ))[[1]][["Pr(>F)"]][1],error=function(e) NULL)),"",sep=""))
      } else {
        #Do a paired anova with Error(Connections/Groups)
        tmp<-c(k,s[1,l],s[2,l],bas,paste(sprintf("%.2g",tryCatch(summary(aov(as.formula(paste("value ~",grouping_column,"+",paste("Error(","Connections/",grouping_column,")",sep=""))),data=df[(df$measure==k) & (df[,grouping_column]==s[1,l] | df[,grouping_column]==s[2,l]),] ))[[1]][[1]][["Pr(>F)"]][1],error=function(e) NULL)),"",sep=""))
      }


      #Ignore if anova fails
      if(!is.na(as.numeric(tmp[length(tmp)]))){

        #Only retain those pairs where the p-values are significant
        if(as.numeric(tmp[length(tmp)])<0.05){
          if(is.null(df_pw)){df_pw<-tmp}else{df_pw<-rbind(df_pw,tmp)}

          #Generate the next position
          bas<-bas+inc
        }
      }
    }
  }

  if(!is.null(df_pw)){
    if(sum(class(df_pw) %in% c("character"))>0){
      df_pw<-t(as.matrix(df_pw))
    }
    df_pw<-data.frame(row.names=NULL,df_pw)
    names(df_pw)<-c("measure","from","to","y","p")
  }

  return(list(df, df_pw))
}

#### Plotting Function ####

#' @import ggplot2
plot_alpha_diversity <- function(df, df_pw) {
  p<-NULL
  if(!"Type" %in% colnames(meta_table)){
    if(!"Type2" %in% colnames(meta_table)){
      p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column,group=grouping_column),data=df)
    } else {
      p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column,group=grouping_column,shape="Type2"),data=df)
    }
  } else {
    p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column,group=grouping_column,shape="interaction(Type,Type2)"),data=df)
  }

  p<-p+geom_boxplot(outlier.size=0,show.legend=FALSE)#+geom_jitter(position = position_jitter(height = 0, width=0),show.legend=FALSE)
  p<-p+theme_bw()
  p<-p+geom_point(size=point_size,alpha=point_opacity)
  p<-p+facet_wrap(~measure,scales="free_y",nrow=number_of_rows)+ylab("Observed Values")+xlab("Samples")

  if(!is.null(df_pw)){
    #This loop will generate the lines and signficances
    for(i in 1:dim(df_pw)[1]){
      p<-p+geom_path(inherit.aes=F,aes(x,y),data = data.frame(x = c(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"])),which(levels(df[,grouping_column])==as.character(df_pw[i,"to"]))), y = c(as.numeric(as.character(df_pw[i,"y"])),as.numeric(as.character(df_pw[i,"y"]))), measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), color="black",lineend = "butt",arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
      p<-p+geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(x=(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"]))+which(levels(df[,grouping_column])==as.character(df_pw[i,"to"])))/2,y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))),size=pairwise_text_size)
      if(exclude_pvalues_text_from_drawing){
        p<-p+geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(x=(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"]))+which(levels(df[,grouping_column])==as.character(df_pw[i,"to"])))/2,y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),label=paste("p=",as.character(as.numeric(as.character(df_pw[i,"p"])))),sep=""),size=pairwise_text_size,vjust=-1)
      }
    }
  }
  if(use_provided_colors){
    p<-p+scale_color_manual(grouping_column,values=colours)
  }
  #if crashes at panel.margin change it to panel.spacing, if crashes at panel.spacing, change it to panel.margin
  p<-p+theme(
    strip.background = element_rect(fill="white"),
    panel.spacing = unit(2, "lines"),
    strip.text = element_text(size=strip_text_size),
    legend.text = element_text(size=legend_text_size),
    text = element_text(size=text_size),
    axis.text = element_text(size=axis_text_size),
    axis.title = element_text(size=axis_title_size),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

  if(legends_position_bottom){
    p<-p+theme(legend.key = element_blank(),  #removes the box around each legend item
               legend.position = "bottom", #legend at the bottom
               legend.direction = "horizontal",
               legend.box = "horizontal",
               legend.box.just = "centre")
  }
  if(exclude_legends){
    p<-p+guides(colour=FALSE)
  }

  return (p)
}

theme_plot <- function(p) {


  return(p)
}

#### Write Data Out ####

write_data <- function() {
  pdf(paste("ANOVA_diversity_",which_level,"_",label,".pdf",sep=""),height=height_image,width=width_image)
  print(p)
  dev.off()

  #Write all the metrics in a file for further analyses elsewhere
  data_to_write<-data.frame(df_R[,"value",drop=F],df_H[,"value",drop=F],df_simp[,"value",drop=F],df_alpha[,"value",drop=F],df_J[,"value",drop=F])
  colnames(data_to_write)<-c("Richness","Shannon","Simpson","FisherAlpha","PielouEvenness")
  write.csv(data_to_write,paste("Diversity_",which_level,"_",label,".csv",sep=""))
  #/Write all the metrics in a file for further analyses else where
}
