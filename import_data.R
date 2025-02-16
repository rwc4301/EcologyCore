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