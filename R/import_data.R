import_data <- function(physeq_path, meta_path, tree_path = NULL, round_abund = FALSE) {
  physeq <- phyloseq::import_biom(physeq_path)
  meta_table <- read.csv(meta_path, header = TRUE, row.names = 1)

  abund_table <- phyloseq::otu_table(physeq)
  abund_table <- t(abund_table)

  # Some cases where you want to round up abundances, e.g. for picrust pathways
  if (round_abund) {
    abund_table<-round(abund_table)
  }

  OTU_tree <- NULL
  if (!is.null(tree_path)) {
    OTU_tree <- ape::read.tree(tree_path)

    # Recent version of R puts apostrophe in the OTU tip labels so we will just remove them if that exist
    OTU_tree$tip.label<-gsub("'","",OTU_tree$tip.label)
  }

  #Uncomment if you'd like to get rid of samples below a certain library size
  abund_table<-abund_table[rowSums(abund_table)>=5000,]
  OTU_taxonomy<-as.data.frame(phyloseq::tax_table(physeq))
  colnames(OTU_taxonomy)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Otus")

  #Ensure that all columns of OTU_taxonomy are character and not factors
  OTU_taxonomy[] <- lapply(OTU_taxonomy, function(x) as.character(x))
  OTU_taxonomy[is.na(OTU_taxonomy)]<-""
  OTU_taxonomy$Otus<-gsub("D_6__|s__","",OTU_taxonomy$Otus)
  OTU_taxonomy$Genus<-gsub("D_5__|g__","",OTU_taxonomy$Genus)
  OTU_taxonomy$Family<-gsub("D_4__|f__","",OTU_taxonomy$Family)
  OTU_taxonomy$Order<-gsub("D_3__|o__","",OTU_taxonomy$Order)
  OTU_taxonomy$Class<-gsub("D_2__|c__","",OTU_taxonomy$Class)
  OTU_taxonomy$Phylum<-gsub("D_1__|p__","",OTU_taxonomy$Phylum)
  OTU_taxonomy$Kingdom<-gsub("D_0__|d__","",OTU_taxonomy$Kingdom)

  #Remove singletons and adjust OTU_taxonomy
  abund_table<-abund_table[,colSums(abund_table)>1]
  OTU_taxonomy<-OTU_taxonomy[colnames(abund_table),]

  #get rid of contaminants with "Unassigned", "Chloroplast" and "Mitochondria" assignment", and "non classified" at Phylum level
  abund_table<-abund_table[,!(OTU_taxonomy$Kingdom %in% c("Unassigned") | OTU_taxonomy$Phylum=="" | OTU_taxonomy$Order %in% c("Chloroplast") | OTU_taxonomy$Family %in% c("Mitochondria"))]

  #extract subset of abund_table for which samples also exists in meta_table
  abund_table<-abund_table[rownames(abund_table) %in% rownames(meta_table),]
  #when reducing the abund_table, there is a high likelihood that an OTU was only present in a sample that is removed, so we shrink
  #the abund_table to get rid of empty columns
  abund_table<-abund_table[,colSums(abund_table)>0]
  #make your meta_table smaller by only considering samples that appear in abund_table
  meta_table<-meta_table[rownames(abund_table),]
  #make OTU_taxonomy smaller by only considering OTUs that appear in abund_table
  OTU_taxonomy<-OTU_taxonomy[colnames(abund_table),]

  #At this point we have abund_table, meta_table, and OTU_taxonomy are ready and their dimensions should match
  return(list(abund_table, OTU_taxonomy, meta_table, OTU_tree))
}
