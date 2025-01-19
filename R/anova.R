# Helper function to perform anova across groups in a data frame

# TODO: not working - calculating pval table weirdly fails only when called inside the package
anova <- function(df, grouping_column) {
  #To do anova, we will convert our data.frame to data.table

  #Since we can't pass a formula to data.table, I am creating
  #a dummy column .group. so that I don't change names in the formula
  dt<-data.table::as.data.table(data.frame(df,.group.=df[,grouping_column]))

  pval<-data.table::as.data.table(dt)[, list(pvalue = sprintf("%.2g", tryCatch({
    pv <- summary(stats::aov(value ~ .group.))[[1]][["Pr(>F)"]][1]
    print(pv)
    pv
  }, error = function(e) NULL))),
  by = list(measure)
  ]

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
