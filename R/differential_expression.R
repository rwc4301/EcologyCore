#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script for finding log2 fold different species using DESeq2

# library(phyloseq)
# library(vegan)
# library(ggplot2)
# library(plyr)
# library(DESeq2)
# library(stringr)

# The 'design' argument should be set to a formula https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/formula
# By default, just regress against the Groups column, but you may want to include batch here too if accounting for batch effects
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
de_create <- function(abund_table, meta_table, OTU_taxonomy, design = as.formula("~ Groups")) {
  which_level<-"Genus" #Phylum Class Order Family Genus Otus
  sig = 0.05
  fold = 2

  # Process input data
  # abund_table <- phyloseq::otu_table(physeq)
  # meta_table <- phyloseq::sample_data(physeq)
  # OTU_taxonomy <- phyloseq::tax_table(physeq)

  #We will convert our table to DESeqDataSet object
  countData = round(as(abund_table, "matrix"), digits = 0)
  # We will add 1 to the countData otherwise DESeq will fail with the error:
  # estimating size factors
  # Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  :
  # every gene contains at least one zero, cannot compute log geometric means
  countData<-(t(countData+1))

  dds <- DESeq2::DESeqDataSetFromMatrix(countData, meta_table, design)

  #Reference:https://github.com/MadsAlbertsen/ampvis/blob/master/R/amp_test_species.R
  #Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
  #Some reason this doesn't work: dds = DESeq(dds, test="wald", fitType="parametric")
  dds = DESeq2::DESeq(dds)

  res_names <- DESeq2::resultsNames(dds)

  return(list(res_names, dds))
}

de_summarise_results <- function(dds, name = NULL) {
  res = DESeq2::results(dds, cooksCutoff = FALSE, name = name)
  return(summary(res))
}

de_get_results <- function(dds, name = NULL) {
  ## Extract the results
  res = DESeq2::results(dds, cooksCutoff = FALSE, name = name)
  res_tax = cbind(as.data.frame(res), as.matrix(countData[rownames(res), ]), OTU = rownames(res))

  # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#log-fold-change-shrinkage-for-visualization-and-ranking
  # or to shrink log fold changes association with condition:
  # res <- DESeq2::lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm") # type = apeglm | ashr | normal

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

  return(list(res_tax, df))
}

ma_plot_differential_expression <- function(res_tax) {
  plot.point.size = 2
  label=F
  tax.display = NULL
  tax.aggregate = "OTU"

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

  return(p1)
  # pdf(paste("NB_MA_",which_level,"_",labels,".pdf",sep=""))
  # print(p1)
  # dev.off()
}

plot_differential_expression <- function(df) {
  height_image=15

  p <- ggplot(df, aes(Groups, Value, colour = Groups)) +
    geom_boxplot(outlier.size = 0) +
    geom_jitter(position = position_jitter(height = 0, width=0),alpha=0.5,outlier.colour = NULL) +
    ylab("Log-relative normalised") +
    facet_wrap( ~ Taxa , scales="free_x",nrow=1) +
    theme(strip.text.x = element_text(size = 16, colour = "black", angle = 90))

  return (p)
}

differential_expression_write <- function() {
  data_to_write<-res_tax_sig[,c("baseMean","log2FoldChange","pvalue","padj")]
  tmp<-sapply(as.character(rownames(res_tax_sig)),function(x){aggregate(data[,x],by=list(meta_table$Groups),FUN=mean)[,2]})
  tmp<-as.data.frame(t(tmp))
  colnames(tmp)<-levels(meta_table$Groups)
  data_to_write<-cbind(data_to_write,Upregulated=as.character(sapply(rownames(tmp),function(x){levels(meta_table$Groups)[which.max(tmp[x,])]})))

  if(which_level=="Otus"){
    rownames(data_to_write)<-as.character(sapply(rownames(data_to_write),function(x) rep(paste(paste(x,gsub(";+$","",paste(sapply(OTU_taxonomy[x,],as.character),collapse=";")))))))
  }
  write.csv(data_to_write,paste("NB_significant_",which_level,"_",labels,".csv",sep=""))
}
