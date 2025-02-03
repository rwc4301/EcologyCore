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

core_microbiome <- function(abund_table, meta_table, OTU_taxonomy, OTU_tree, what_detection = "relative", taxa_rank = "Species", minimum_prevalence = 0.85, short_names = FALSE) {
  # Process input data
  # abund_table <- phyloseq::otu_table(physeq)
  # meta_table <- phyloseq::sample_data(physeq)
  # OTU_taxonomy <- phyloseq::tax_table(physeq)
  # OTU_tree <- phyloseq::phy_tree(physeq)

  physeq <- to_phyloseq(abund_table, taxa_table, meta_table, taxa_tree, taxa_rank)

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
    detections <- 10^seq(log10(1), log10(max(microbiome::abundances(pseq.2))/10), length = 20)
    detections <- round(detections)
    pseq_to_plot<-pseq.2
  }

  datacore <- microbiome::plot_core(
    pseq_to_plot,
    plot.type = "heatmap",
    prevalences = prevalences,
    detections = detections,
    # colours = gray(seq(1,0,length=20)),
    colours = rev(RColorBrewer::brewer.pal(5, "Spectral")),
    min.prevalence = minimum_prevalence,
    horizontal = TRUE
  )

  # get the data used for plotting
  df <- datacore$data

  if(which_level=="Otus") {
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



  } else {
    df <- datacore$data[datacore$data$Taxa!="__Unknowns__",]
  }

  # Function to shorten taxa names
  shorten_levels <- sapply(levels(df$Taxa), function(x) {
    # Split each level at the ';' and retain only the last element, unless that is non-descriptive in which case take a longer name
    elements <- strsplit(x, ";")[[1]]
    if (grepl("uncultured", x) | grepl("Subgroup", x)) {
      #return(paste(tail(elements, n = 2), sep = ";"))
      paste(elements[length(elements) - 1], elements[length(elements)], sep = ";")
    } else {
      tail(elements, n = 1)
    }
  })

  if (short_names) {
    levels(df$Taxa) <- shorten_levels
  }

  # Add a group dummy column for facet title
  if (length(unique(meta_table$Groups)) == 1) {
    df$.group. <- meta_table$Groups[1]
  }

  # replace the data in the plot object
  datacore$data <- df

  datacore <- datacore + ggplot2::facet_wrap(~ .group.)

  return(datacore)
}
