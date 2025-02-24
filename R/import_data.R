import_rds <- function(physeq_path) {
  return (readRDS(physeq_path))
}

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

  # abund_table = read.delim("Combined/feature-table.tsv", header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")

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

  # Seriously weird bug here
  if (!identical(phyloseq::taxa_names(phyloseq::otu_table(abund_table)), phyloseq::taxa_names(phyloseq::tax_table(as.matrix(taxa_table))))) {
    abund_table@taxa_are_rows = !abund_table@taxa_are_rows
  }

  # Merge the edited tables into a new phyloseq object and return
  return (join_phyloseq(abund_table, taxa_table, meta_table, taxa_tree))
}

join_phyloseq <- function(abund_table, taxa_table, meta_table, taxa_tree = NULL) {
  return (phyloseq::merge_phyloseq(phyloseq::otu_table(abund_table, taxa_are_rows = FALSE), phyloseq::tax_table(as.matrix(taxa_table)), phyloseq::sample_data(meta_table), taxa_tree))
}

# Use to combine multiple phyloseq objects together for example when working with multiple batches
merge_data <- function(phyloseqs) {
  batch1 <- phyloseq::import_biom("data/Batch 1/asv_even_taxon.biom", treefilename = "data/Batch 1/tree_rooted.nwk")
  batch2 <- phyloseq::import_biom("data/Batch 2/asv_even_taxon2.biom", treefilename = "data/Batch 2/rooted_tree2.nwk")
  batch3 <- phyloseq::import_biom("data/Batch 3/asv_even_taxon3.biom", treefilename = "data/Batch 3/rooted_tree3.nwk")

  # Unclear whether ASV nums are consistent across batches, so edit taxa names
  b1_otu <- as.data.frame(otu_table(batch1))
  b1_tax <- as.data.frame(tax_table(batch1))
  b1_tree <- phy_tree(batch1)
  b1_names <- paste(taxa_names(batch1), "1", sep = "_")

  # Append batch num to taxa names
  rownames(b1_otu) <- b1_names
  rownames(b1_tax) <- b1_names
  b1_tree$tip.label <- b1_names

  # Check valid
  if (intersect(rownames(b1_otu), rownames(b1_tax), b1_tree$tip.label) > 0) {
    # Rebuild phyloseq
    b1_otu <- otu_table(as.matrix(b1_otu), taxa_are_rows = TRUE)
    b1_tax <- tax_table(as.matrix(b1_tax))

    batch1 <- merge_phyloseq(b1_otu, b1_tax, b1_tree)
  }
}

metadata_frame <- function(physeq) {
  ma <- as.data.frame(physeq@sam_data@.Data)
  rownames(ma) <- physeq@sam_data@row.names
  colnames(ma) <- physeq@sam_data@names

  return(ma)
}

expand_otu_names <- function(otu_names, taxa_table, use_short_names = FALSE) {
  if (!use_short_names) {
    return(paste(otu_names, sapply(otu_names, function(x) {
      # For each OTU, replace NA with "" for each taxonomic level.
      tax_vec <- sapply(taxa_table[x, ], function(val) {
        if (is.na(val)) "" else as.character(val)
      })
      cleaned_tax <- gsub(";+$", "", paste(tax_vec, collapse = ";"))
      return(cleaned_tax)
    })))
  } else {
    # Generate short names
    names_out <- sapply(otu_names, function(x) {
      # Get taxonomic levels for OTU x, replacing NA with empty strings.
      taxa_levels <- sapply(taxa_table[x, ], function(val) ifelse(is.na(val), "", as.character(val)))

      # Attempt to use the "Species" level if available.
      if (!is.null(taxa_levels["Species"]) && taxa_levels["Species"] != "") {
        return(paste(taxa_levels["Genus"], taxa_levels["Species"]))
      } else if (!is.null(taxa_levels["Genus"]) && taxa_levels["Genus"] != "") {
        # If species is missing but genus is available, just return the genus.
        return(taxa_levels["Genus"])
      } else {
        # If collating at a higher taxonomy leads to missing lower-level info,
        # try to find the highest non-empty taxonomic level.
        nonempty_indices <- which(taxa_levels != "" & taxa_levels != "uncultured")
        if (length(nonempty_indices) > 0) {
          last_meaningful <- max(nonempty_indices)
          return(paste0("Unknown ", taxa_levels[last_meaningful]))
        } else {
          return("Unknown taxonomy")
        }
      }
    })

    # Handle duplicates by adding an index
    duplicated_names <- table(names_out)
    for (name in names(duplicated_names)) {
      if (duplicated_names[name] > 1) {
        indices <- which(names_out == name)
        names_out[indices] <- paste0(names_out[indices], "_", seq_along(indices))
      }
    }

    return(names_out)
  }
}

  # if(which_level=="Otus"){
  #   if(N==dim(x)[2]){
  #     colnames(new_x)<-c(paste(colnames(new_x),sapply(colnames(new_x),function(x) gsub(";+$","",paste(sapply(OTU_taxonomy[x,],as.character),collapse=";")))))
  #   } else {
  #     colnames(new_x)<-c(paste(colnames(new_x)[-(N+1)],sapply(colnames(new_x)[-(N+1)],function(x) gsub(";+$","",paste(sapply(OTU_taxonomy[x,],as.character),collapse=";")))),"Others")
  #   }
  # }

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

collate_taxonomy_fast <- function(physeq, rank = "Species") {
  if (!rank %in% phyloseq::rank_names(physeq)) {
    stop(sprintf("Tried to collate taxonomy around unknown rank %s", rank))
  }

  # Convert to data.table format
  require(data.table)
  otu_tab <- as.data.table(as.data.frame(phyloseq::otu_table(physeq)))
  tax_tab <- as.data.table(as.data.frame(phyloseq::tax_table(physeq)))

  # Add taxonomy column to abundance data
  otu_tab[, tax_group := tax_tab[[rank]]]

  # Aggregate by taxonomic group
  agg_tab <- otu_tab[, lapply(.SD, sum),
                     by = tax_group,
                     .SDcols = names(otu_tab)[1:(ncol(otu_tab)-1)]]

  # Convert back to phyloseq format
  otu_mat <- as.matrix(agg_tab[, -"tax_group"])
  rownames(otu_mat) <- agg_tab$tax_group

  # Create new taxonomy table
  tax_mat <- matrix(agg_tab$tax_group,
                    ncol = ncol(tax_tab),
                    dimnames = list(agg_tab$tax_group, colnames(tax_tab)))

  # Return new phyloseq object
  return(phyloseq::phyloseq(
    phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE),
    phyloseq::tax_table(tax_mat),
    phyloseq::sample_data(physeq)
  ))
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
