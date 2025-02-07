#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script for proportional representation of N most abundant species

# library(phyloseq)
# library(vegan)
# library(plyr)
# library(ape)
# library(grid)
# library(tidyverse)
# library(reshape2)
# library(vegan)
# library(ggsci)

# Adapted from https://github.com/ShadeLab/PAPER_Shade_CurrOpinMicro
# input dataset needs to be rarefied and the rarefaction depth included
#

# Class Descriptions for EcologyCore Core (ECC):
## ECCOccupancyAbundance
## ECCRankedSimilarity

# set threshold_model to 'elbow' to use first-order difference
core_occupancy_abundance_model <- function(physeq, threshold_model = 2) {
  # Add input validation
  if (is.null(physeq)) {
    stop("Input phyloseq object cannot be NULL")
  }

  abund_table <- as.data.frame(phyloseq::otu_table(physeq))

  # Ensure consistent orientation of abundance table
  if (!phyloseq::taxa_are_rows(physeq)) {
    abund_table <- t(abund_table)
  }

  # Add check for zero row sums before decostand
  if (any(rowSums(abund_table) == 0)) {
    warning("Removing samples with zero counts")
    abund_table <- abund_table[rowSums(abund_table) > 0,]
  }

  # Add check for zero column sums
  if (any(colSums(abund_table) == 0)) {
    warning("Removing OTUs with zero counts")
    abund_table <- abund_table[,colSums(abund_table) > 0]
  }

  # If a tax table is supplied we'll use to generate readable names for ASVs
  taxa_table <- NULL
  if (!is.null(physeq@tax_table)) {
    taxa_table <- as.data.frame(phyloseq::tax_table(physeq))
  }

  meta_table <- as.data.frame(phyloseq::sample_data(physeq))

  #TODO: check SampleID column doesn't exist first
  meta_table$SampleID <- rownames(meta_table)

  # TODO: clean up hypothesis testing
  label="Abstraction_Method"
  #meta_table$Groups <- as.factor(as.character(paste("Run ", meta_table$run, sep = "")))
  meta_table$Groups <- "Test Group"
  meta_table$Type <- NULL
  meta_table$Type2 <- NULL
  meta_table$Connections <- NULL

  #meta_table <- meta_table[meta_table$Groups == "Run 1",]

  # Validate that sample names exist in abundance table
  valid_samples <- intersect(rownames(meta_table), colnames(abund_table))
  if (length(valid_samples) == 0) {
    stop("No matching sample names between metadata and abundance table")
  }

  # Subset both tables to matching samples
  meta_table <- meta_table[valid_samples,]
  abund_table <- abund_table[, valid_samples]

  # Remove empty OTUs after subsetting
  abund_table <- abund_table[rowSums(abund_table) > 0,]

  #abund_table <- t(abund_table)

  # What proportion of samples does each ASV occupy
  otu_PA <- 1*((abund_table>0)==1)                                              # presence-absence data
  otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                       # occupancy calculation

  # Modify the relative abundance calculation to avoid NaN
  otu_rel <- apply(decostand(abund_table, method="total", MARGIN=2), 1, function(x) {
    if (all(is.na(x))) return(0)
    mean(x, na.rm = TRUE)
  })

  occ_abund_table <- data.frame(otu_occ=otu_occ, otu_rel=otu_rel) %>%           # combining occupancy and abundance data frame
    rownames_to_column('otu')

  # # Collate response and set class name
  # res <- list(abund_table = abund_table, occ_abund_table = occ_abund_table)
  # class(res) <- "ECCOccupancyAbundance"

  otu <- abund_table

  PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>%
    gather(SampleID, abun, -otu) %>%
    left_join(meta_table, by = 'SampleID') %>%
    group_by(otu, Groups) %>%   # maybe group by type instead
    summarise(time_freq=sum(abun>0)/length(abun),            # frequency of detection between time points
              coreTime=ifelse(time_freq == 1, 1, 0)) %>%     # 1 only if occupancy 1 with specific time, 0 if not
    group_by(otu) %>%
    summarise(sumF=sum(time_freq),
              sumG=sum(coreTime),
              nS=length(Groups),
              Index=(sumF+sumG)/nS)                 # calculating weighting Index based on number of time points detected and

  otu_ranked <- occ_abund_table %>%
    left_join(PresenceSum, by='otu') %>%
    transmute(otu=otu,                           # check what transmute function does
              rank=Index) %>%
    arrange(desc(rank))

  # Calculating the contribution of ranked OTUs to the BC similarity
  BCaddition <- NULL

  nReads=1000

  # calculating BC dissimilarity based on the 1st ranked OTU
  otu_start=otu_ranked$otu[1]
  start_matrix <- as.matrix(otu[otu_start,])
  #start_matrix <- t(start_matrix)
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/(2*nReads))
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_s <- data.frame(x_names,x)
  names(df_s)[2] <- 1
  BCaddition <- rbind(BCaddition,df_s)
  # calculating BC dissimilarity based on addition of ranked OTUs from 2nd to 500th.
  # Can be set to the entire length of OTUs in the dataset, however it might take
  # some time if more than 5000 OTUs are included.
  for(i in 2:500){
    otu_add=otu_ranked$otu[i]
    add_matrix <- as.matrix(otu[otu_add,])
    #add_matrix <- t(add_matrix)
    start_matrix <- rbind(start_matrix, add_matrix)
    x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads))
    x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
    df_a <- data.frame(x_names,x)
    names(df_a)[2] <- i
    BCaddition <- left_join(BCaddition, df_a, by=c('x_names'))
  }
  # calculating the BC dissimilarity of the whole dataset (not needed if the second loop is already including all OTUs)
  x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
  x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
  df_full <- data.frame(x_names,x)
  names(df_full)[2] <- length(rownames(otu))
  BCfull <- left_join(BCaddition,df_full, by='x_names')

  rownames(BCfull) <- BCfull$x_names
  temp_BC <- BCfull
  temp_BC$x_names <- NULL
  temp_BC_matrix <- as.matrix(temp_BC)

  BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>%
    gather(comparison, BC, -rank) %>%
    group_by(rank) %>%
    summarise(MeanBC=mean(BC)) %>%            # mean Bray-Curtis dissimilarity
    arrange(desc(-MeanBC)) %>%
    mutate(proportionBC=MeanBC/max(MeanBC))   # proportion of the dissimilarity explained by the n number of ranked OTUs
  Increase=BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
  increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
  BC_ranked <- left_join(BC_ranked, increaseDF)
  BC_ranked <- BC_ranked[-nrow(BC_ranked),]

  #Creating thresholds for core inclusion

  # elbow method (first order difference) (script modified from https://pommevilla.github.io/random/elbows.html)
  fo_difference <- function(pos){
    left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos
    right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)
    return(left - right)
  }
  BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)

  threshold <- 0
  if (threshold_model == "elbow") {
    # get highest first-order difference in LCBDs
    threshold <- which.max(BC_ranked$fo_diffs)
  } else {
    # get last taxa where increase in explanatory value >= 2%
    threshold <- last(as.numeric(as.character(BC_ranked$rank[(BC_ranked$IncreaseBC >= 1 + (threshold_model / 100))])))
  }

  # Use Sloan neutral model to prioritize OTUs
  # Fitting neutral model (Burns et al., 2016 (ISME J) - functions are in the sncm.fit.R)
  spp=t(otu)
  taxon=as.vector(rownames(otu))

  #Models for the whole community
  obs.np=sncm.fit(spp, taxon, stats=FALSE, pool=NULL)
  sta.np=sncm.fit(spp, taxon, stats=TRUE, pool=NULL)

  above.pred=sum(obs.np$freq > (obs.np$pred.upr), na.rm=TRUE)/sta.np$Richness  # fraction of OTUs above prediction
  below.pred=sum(obs.np$freq < (obs.np$pred.lwr), na.rm=TRUE)/sta.np$Richness  # fraction of OTUs below prediction

  #Create a column defining "core" OTUs
  occ_abund_table$fill <- 'no'
  occ_abund_table$fill[occ_abund_table$otu %in% otu_ranked$otu[1:threshold]] <- 'core'

  # Collate response and set class name
  res <- list(
    abund_table = abund_table,
    taxa_table = taxa_table,
    meta_table = meta_table,
    occ_abund_table = occ_abund_table,
    BC_ranked = BC_ranked,
    threshold = threshold,
    threshold_model = threshold_model,
    obs.np = obs.np,
    sta.np = sta.np,
    above.pred = above.pred,
    below.pred = below.pred,
    otu_ranked = otu_ranked
  )
  class(res) <- "ECCOccupancyAbundance"

  return (res)
}

plot.ECCOccupancyAbundance <- function(value) {
  # Occupancy abundance plot:
  # p <- ggplot(data = value$occ_abund_table, aes(x = log10(otu_rel), y = otu_occ)) +
  #   geom_point(pch = 21, fill = 'white') +
  #   labs(x = "log10(mean relative abundance)", y = "Occupancy")

  p <- ggplot() +
    geom_point(data = value$occ_abund_table[value$occ_abund_table$fill == 'no',], aes(x = log10(otu_rel), y = otu_occ), pch=21, fill='white', alpha=0.2) +
    geom_point(data = value$occ_abund_table[value$occ_abund_table$fill != 'no',], aes(x = log10(otu_rel), y = otu_occ), pch=21, fill='blue', size=1.8) +
    geom_line(color='black', data= value$obs.np, size=1, aes(y= value$obs.np$freq.pred, x=log10(value$obs.np$p)), alpha=.25) +
    geom_line(color='black', lty='twodash', size=1, data= value$obs.np, aes(y= value$obs.np$pred.upr, x=log10(value$obs.np$p)), alpha=.25)+
    geom_line(color='black', lty='twodash', size=1, data= value$obs.np, aes(y= value$obs.np$pred.lwr, x=log10(value$obs.np$p)), alpha=.25)+
    labs(x="log10(mean relative abundance)", y="Occupancy")

  #TODO: highlight on the same occupancy-abundance plot OTUs that are above and below the neutral model prediction
  #TODO: add to the plot the above.pred and below.pred values

  print(p)
}

plot.ECCRankedSimilarity <- function(value) {
  BC_ranked <- value$BC_ranked

  #Creating plot of Bray-Curtis similarity
  p <- ggplot(BC_ranked[1:100,], aes(x=factor(BC_ranked$rank[1:100], levels=BC_ranked$rank[1:100]))) +
    geom_point(aes(y=proportionBC)) +
    theme_classic() + theme(strip.background = element_blank(),axis.text.x = element_text(size=7, angle=45)) +
    geom_vline(xintercept = value$threshold, lty=3, col='red', cex=.5) +
    # geom_vline(xintercept=lastCall, lty=3, col='blue', cex=.5) +
    labs(x='Ranked Taxa',y='Contribution to Bray-Curtis Dissimilarity')

  if (value$threshold_model == "elbow") {
    p <- p + annotate(geom = "text", x = value$threshold + 14, y = 0.1, label = sprintf("Elbow method (%d)", value$threshold), color = "red")
  } else {
    p <- p + annotate(geom = "text", x = value$threshold, y = 0.5, label = sprintf("Last %d%% increase (%d)", value$threshold_model, value$threshold), color = "blue")
  }

  print(p)
}

plot.ECCHeatmap <- function(value, what_detection = "absolute") {
  prevalences <- seq(.05, 1, .05)

  #Creating a plot of core taxa occupancy by time point
  core <- value$occ_abund_table$otu[value$occ_abund_table$fill == 'core']
  core_otu <- value$abund_table[rownames(value$abund_table) %in% core,]

  # # keep only taxa with positive sums
  # pseq.2 <- prune_taxa(taxa_sums(physeq) > 0, physeq)
  #
  # # Calculate compositional version of the data
  # # (relative abundances)
  # pseq.rel <- microbiome::transform(pseq.2, "compositional")

  detections <- NULL
  if(what_detection=="relative"){
    #Detection with Relative Abundances
    detections <- 10^seq(log10(1e-3), log10(.2), length = 20)
  } else if(what_detection=="absolute") {
    #Detection with Absolute Count
    detections <- 10^seq(log10(1), log10(max(microbiome::abundances(core_otu))), length = 20)
    detections <- round(detections)
  }

  if (!is.null(value$taxa_table)) {
    rownames(core_otu) <- expand_otu_names(rownames(core_otu), value$taxa_table)
  }

  p <- microbiome::plot_core(
    core_otu,
    prevalences = prevalences,
    detections = detections,
    plot.type = "heatmap",
    colours = rev(RColorBrewer::brewer.pal(5, "Spectral")),
    # taxa.order = value$otu_ranked$otu,
    horizontal = TRUE,
  )

  return(p)


  # otu_relabun <- decostand(otu, method="total", MARGIN=2)
  #
  # plotDF <- data.frame(otu = as.factor(row.names(otu_relabun)), otu_relabun) %>%
  #   gather(SampleID, relabun, -otu) %>%
  #   left_join(meta_table, by = 'SampleID') %>%
  #   left_join(otu_ranked, by='otu') %>%
  #   filter(otu %in% core) %>%
  #   group_by(otu, Groups) %>%
  #   summarise(time_freq=sum(relabun>0)/length(relabun),
  #             coreTime=ifelse(time_freq == 1, 1, 0),
  #             detect=ifelse(time_freq > 0, 1, 0))
  #
  # plotDF$otu <- factor(plotDF$otu, levels=otu_ranked$otu[1:34])
  #
  # p <- ggplot(plotDF,aes(x=otu, time_freq,fill=factor(Groups))) +
  #   geom_bar(stat = 'identity', position = 'dodge') +
  #   coord_flip() +
  #   scale_x_discrete(limits = rev(levels(plotDF$otu))) +
  #   theme(axis.text = element_text(size=6)) +
  #   labs(x='Ranked OTUs', y='Occupancy by site', fill="Sampling date") +
  #   scale_fill_npg()
}

core_taxa <- function(abund_table, meta_table, OTU_taxonomy, OTU_tree, N = 25) {
  # TODO: remove
  grouping_column <- "Groups"

  # Process input data
  # abund_table <- phyloseq::otu_table(physeq)
  # meta_table <- phyloseq::sample_data(physeq)
  # OTU_taxonomy <- phyloseq::tax_table(physeq)
  # OTU_tree <- phyloseq::phy_tree(physeq)

  #Apply proportion normalisation
  x<-abund_table/rowSums(abund_table)
  x<-x[,order(colSums(x),decreasing=TRUE)]


  taxa_list<-colnames(x)[1:min(dim(x)[2],N)]
  if (c("__Unknowns__") %in% taxa_list){
    taxa_list<-colnames(x)[1:min(dim(x)[2],N+1)]
    taxa_list<-taxa_list[!grepl("__Unknowns__",taxa_list)]
  }
  N<-length(taxa_list)

  #Generate a new table with everything added to Others
  new_x<-NULL
  if(N==dim(x)[2]){
    new_x<-data.frame(x[,colnames(x) %in% taxa_list])
  } else {
    new_x<-data.frame(x[,colnames(x) %in% taxa_list],Others=rowSums(x[,!colnames(x) %in% taxa_list]))
  }
  if(which_level=="Otus"){
    if(N==dim(x)[2]){
      colnames(new_x)<-c(paste(colnames(new_x),sapply(colnames(new_x),function(x) gsub(";+$","",paste(sapply(OTU_taxonomy[x,],as.character),collapse=";")))))
    } else {
      colnames(new_x)<-c(paste(colnames(new_x)[-(N+1)],sapply(colnames(new_x)[-(N+1)],function(x) gsub(";+$","",paste(sapply(OTU_taxonomy[x,],as.character),collapse=";")))),"Others")
    }
  }

  df<-NULL
  for (i in 1:dim(new_x)[2]){
    tmp<-data.frame(row.names=NULL,Sample=rownames(new_x),Taxa=rep(colnames(new_x)[i],dim(new_x)[1]),Value=new_x[,i],Groups=meta_table$Groups)
    if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
  }
  colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00",rainbow(300)[seq(1,300,6)]);

  #The first step to get levels is to find all the possible values in a given column in df$sample by outputing the following command on terminal:
  # cat(paste("levels=c(",paste(paste("\"",unique(as.character(df$Sample)),"\"",sep=""),collapse=","),")",sep=""))
  # then use df$Sample<-factor(as.character(df$Sample),levels=c()) list

  return(df)
}

plot_core_taxa <- function(df, N = 25, switch_strip = "x", legend_columns = 1) {
  reveal_sample_names=TRUE
  legend_text_size=20
  axis_title_size=30
  text_size=30
  axis_text_size=30
  strip_text_size=30
  height_image=25
  width_image=60

  p <- ggplot(df, aes(Sample, Value, fill = Taxa)) +
    geom_bar(stat = "identity") +
    facet_grid(. ~ Groups, drop = TRUE, scale = "free", space = "free_x", switch = switch_strip) +
    scale_fill_manual(values = colours[1:(N+1)], guide = guide_legend(ncol = legend_columns)) +
    ylab("Proportions") +
    scale_y_continuous(expand = c(0.02, 0)) +
    theme(
      strip.background = element_rect(fill="gray85"),
      panel.spacing = unit(0.3, "lines"),
      strip.text = element_text(size=strip_text_size,angle=90),
    )

  if(reveal_sample_names) {
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  } else {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks = element_blank())
  }

  return(p)
}
