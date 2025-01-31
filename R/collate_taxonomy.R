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
