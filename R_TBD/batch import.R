batch1 <- phyloseq::import_biom("data/Batch 1/asv_even_taxon.biom")#, treefilename = "data/Batch 1/tree_rooted.nwk")
batch2 <- phyloseq::import_biom("data/Batch 2/asv_even_taxon2.biom")#, treefilename = "data/Batch 2/rooted_tree2.nwk")
batch3 <- phyloseq::import_biom("data/Batch 3/asv_even_taxon3.biom")#, treefilename = "data/Batch 3/rooted_tree3.nwk")

append_batch_id <- function(batch1, i) {
  # Unclear whether ASV nums are consistent across batches, so edit taxa names
  b1_otu <- as.data.frame(otu_table(batch1))
  b1_tax <- as.data.frame(tax_table(batch1))
  # b1_tree <- phy_tree(batch1)
  b1_names <- paste(taxa_names(batch1), i, sep = "_")

  # Append batch num to taxa names
  rownames(b1_otu) <- b1_names
  rownames(b1_tax) <- b1_names
  # b1_tree$tip.label <- b1_names

  # Check valid
  if (length(intersect(rownames(b1_otu), rownames(b1_tax))) > 0) {
    # Rebuild phyloseq
    b1_otu <- otu_table(as.matrix(b1_otu), taxa_are_rows = TRUE)
    b1_tax <- tax_table(as.matrix(b1_tax))

    batch1 <- merge_phyloseq(b1_otu, b1_tax)#, b1_tree)

    return (batch1)
  } else {
    stop ("Didn't work :(")
  }
}

batch1 <- append_batch_id(batch1, "1")
batch2 <- append_batch_id(batch2, "2")
batch3 <- append_batch_id(batch3, "3")

