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
      # Get taxonomic levels for OTU x; replace any NA with empty strings.
      taxa_levels <- sapply(taxa_table[x, ], function(val) {
        if (is.na(val)) "" else as.character(val)
      })
      
      # Extract cleaned species and genus values with explicit NA checking.
      species_val <- if (!is.null(taxa_levels["Species"]) && !is.na(taxa_levels["Species"]) &&
                          taxa_levels["Species"] != "") {
        taxa_levels["Species"]
      } else {
        ""
      }
      genus_val <- if (!is.null(taxa_levels["Genus"]) && !is.na(taxa_levels["Genus"]) &&
                        taxa_levels["Genus"] != "") {
        taxa_levels["Genus"]
      } else {
        ""
      }
      
      # Use the species level if available.
      if (species_val != "") {
        return(paste(genus_val, species_val))
      } else if (genus_val != "") {
        # If species is missing but genus is available, just return the genus.
        return(genus_val)
      } else {
        # If both species and genus are missing, look for any non-empty taxonomic level.
        nonempty_indices <- which(taxa_levels != "" & taxa_levels != "uncultured")
        if (length(nonempty_indices) > 0) {
          last_meaningful <- max(nonempty_indices)
          return(paste0("Unknown ", taxa_levels[last_meaningful]))
        } else {
          return("Unknown taxonomy")
        }
      }
    })
    
    # Handle duplicates by appending a numeric suffix.
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

collate_taxonomy_fast <- function(physeq, rank = "Species") {
  if (!rank %in% phyloseq::rank_names(physeq)) {
    stop(sprintf("Tried to collate taxonomy around unknown rank %s", rank))
  }
  
  # Ensure data.table is loaded
  require(data.table)
  
  # Convert to data.table format with proper scoping
  otu_tab <- data.table::setDT(as.data.frame(phyloseq::otu_table(physeq)))
  tax_tab <- data.table::setDT(as.data.frame(phyloseq::tax_table(physeq)))
  
  # Add taxonomy column to abundance data (fixed syntax)
  otu_tab[, "tax_group" := tax_tab[[rank]]]
  
  # Aggregate by taxonomic group
  agg_tab <- otu_tab[, lapply(.SD, sum), 
                     by = "tax_group",
                     .SDcols = names(otu_tab)[names(otu_tab) != "tax_group"]]
  
  # Convert back to phyloseq format
  otu_mat <- as.matrix(agg_tab[, .SD, .SDcols = names(agg_tab)[names(agg_tab) != "tax_group"]])
  rownames(otu_mat) <- agg_tab$tax_group
  
  # Create new taxonomy table - fill with NA except for the aggregated rank
  tax_mat <- matrix(NA, 
                   nrow = nrow(agg_tab),
                   ncol = ncol(tax_tab),
                   dimnames = list(agg_tab$tax_group, colnames(tax_tab)))
  tax_mat[, rank] <- agg_tab$tax_group
  
  # Return new phyloseq object
  return(phyloseq::phyloseq(
    phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE),
    phyloseq::tax_table(tax_mat),
    phyloseq::sample_data(physeq)
  ))
} 