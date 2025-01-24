#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script for core analysis
#Reference: https://microbiome.github.io/tutorials/Core.html
#v2.0 fixed reordering bug

# library(phyloseq)
# library(ggplot2)
# library(viridis)
# library(microbiome)
# library(RColorBrewer)
# library(cowplot)

core_microbiome <- function(physeq) {
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
    datacore$data<-datacore$data[datacore$data$Taxa!="__Unknowns__",]
  }

  return(datacore)
}
