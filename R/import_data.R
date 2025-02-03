import_data <- function(biom_path, meta_path, tree_path = NULL, round_abund = FALSE) {
  # Wrap everything in the physeq class and return
  physeq <- NULL
  if (!is.null(tree_path)) {
    physeq <- phyloseq::import_biom(biom_path, treefilename = tree_path)
    taxa_tree <- phyloseq::phy_tree(physeq)

    # Recent version of R puts apostrophe in the OTU tip labels so we will just remove them if that exist
    taxa_tree$tip.label<-gsub("'","",taxa_tree$tip.label)
  } else {
    physeq <- phyloseq::import_biom(biom_path)
  }
  meta_table <- read.csv(meta_path, header = TRUE, row.names = 1)

  abund_table <- phyloseq::otu_table(physeq)
  abund_table <- t(abund_table)

  # Some cases where you want to round up abundances, e.g. for picrust pathways
  if (round_abund) {
    abund_table<-round(abund_table)
  }

  taxa_tree <- NULL
  if (!is.null(tree_path)) {
    taxa_tree <- ape::read.tree(tree_path)

    # Recent version of R puts apostrophe in the OTU tip labels so we will just remove them if that exist
    taxa_tree$tip.label<-gsub("'","",taxa_tree$tip.label)
  }

  # Uncomment if you'd like to get rid of samples below a certain library size
  abund_table<-abund_table[rowSums(abund_table)>=5000,]

  colnames(phyloseq::tax_table(physeq)) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

  taxa_table <- as.data.frame(phyloseq::tax_table(physeq))

  # Ensure that all columns of taxa_table are character and not factors
  #taxa_table[] <- lapply(taxa_table, function(x) as.character(x))
  taxa_table[is.na(taxa_table)]<-""

  # Remove prefixes from taxonomy names
  taxa_table$Species<-gsub("D_6__|s__","",taxa_table$Species)
  taxa_table$Genus<-gsub("D_5__|g__","",taxa_table$Genus)
  taxa_table$Family<-gsub("D_4__|f__","",taxa_table$Family)
  taxa_table$Order<-gsub("D_3__|o__","",taxa_table$Order)
  taxa_table$Class<-gsub("D_2__|c__","",taxa_table$Class)
  taxa_table$Phylum<-gsub("D_1__|p__","",taxa_table$Phylum)
  taxa_table$Kingdom<-gsub("D_0__|d__|k__","",taxa_table$Kingdom)

  # Remove singletons and adjust taxa_table
  abund_table<-abund_table[,colSums(abund_table)>1]
  taxa_table<-taxa_table[colnames(abund_table),]

  # Get rid of contaminants with "Unassigned", "Chloroplast" and "Mitochondria" assignment", and "non classified" at Phylum level
  abund_table<-abund_table[,!(taxa_table$Kingdom %in% c("Unassigned") | taxa_table$Phylum=="" | taxa_table$Order %in% c("Chloroplast") | taxa_table$Family %in% c("Mitochondria"))]

  # Extract subset of abund_table for which samples also exists in meta_table
  abund_table<-abund_table[rownames(abund_table) %in% rownames(meta_table),]

  # There is a high likelihood that there are OTUs only present in a sample that was removed, so shrink abund_table to get rid of empty columns
  abund_table<-abund_table[,colSums(abund_table)>0]

  # Shrink meta table by only considering samples that appear in abund_table
  meta_table<-meta_table[rownames(abund_table),]

  # Shrink taxonomy table by only considering OTUs that appear in abund_table
  taxa_table<-taxa_table[colnames(abund_table),]

  # Coerce into phyloseq objects
  meta_table <- phyloseq::sample_data(meta_table)
  taxa_table <- phyloseq::tax_table(as.matrix(taxa_table))

  # Merge the edited tables into a new phyloseq object and return
  return (phyloseq::merge_phyloseq(phyloseq::phyloseq(abund_table, taxa_table), meta_table, taxa_tree))
}

# Use to combine multiple phyloseq objects together for example when working with multiple batches
merge_data <- function(phyloseqs) {

}

collate_taxonomy_2 <- function(physeq, rank = "Species") {
  if (!rank %in% phyloseq::rank_names(physeq)) {
    stop(sprintf("Tried to collate taxonomy around unknown rank %s", rank))
  }
  return (phyloseq::tax_glom(physeq, taxrank = rank))
}

collate_taxonomy <- function(abund_table, taxa_table, rank = "Species") {
  new_abund_table <- NULL

  if (rank == "Species") {
    new_abund_table <- abund_table
  } else {
    list <- unique(taxa_table[, rank])

    for (i in list) {
      tmp <- data.frame(rowSums(abund_table[, rownames(taxa_table)[taxa_table[, rank] == i]]))
      if (i == "") {
        colnames(tmp) <- c("__Unknowns__")
      } else {
        colnames(tmp) <- paste("", i, sep="")
      }
      if (is.null(new_abund_table)) {
        new_abund_table <- tmp
      } else {
        new_abund_table <- cbind(tmp, new_abund_table)
      }
    }
  }

  new_abund_table <- as.data.frame(as(new_abund_table, "matrix"))
  return(new_abund_table)
}

to_phyloseq <- function(abund_table, OTU_taxonomy, meta_table, taxa_tree, taxa_rank) {
  #Convert the data back to phyloseq format
  OTU = phyloseq::otu_table(as.matrix(abund_table), taxa_are_rows = FALSE)
  TAX = phyloseq::tax_table(as.matrix(OTU_taxonomy))
  SAM = phyloseq::sample_data(meta_table)

  # Something weird seems to happen with all the transposing abundance tables where sample names get set as taxa names
  # Can't quite figure out where it happens so just bodging it with this conditional
  if (!identical(phyloseq::taxa_names(TAX), phyloseq::taxa_names(OTU))) {
    OTU@taxa_are_rows = FALSE
    warning("Taxa names didn't match first time, applied workaround. This is a bug in EcologyCore.")
    #stop("Taxa names are not equal between taxa and abundance tables.")
  }

  physeq<-NULL
  if (taxa_rank == "Species") {
    #physeq<-merge_phyloseq(phyloseq(OTU, TAX),SAM,midpoint(taxa_tree))
    physeq<-phyloseq::merge_phyloseq(phyloseq::phyloseq(OTU, TAX),SAM,taxa_tree)
  } else {
    physeq<-phyloseq::merge_phyloseq(phyloseq::phyloseq(OTU),SAM)
  }

  return(physeq)
}
