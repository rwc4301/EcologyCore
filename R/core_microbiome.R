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

# set model_threshold to 'elbow' to use first-order difference
# model = occ_abund for new occupancy model; set model_threshold for last % increase
# model = occupancy for old model based on prevalence and set model_threshold for min prevalence
#' @import tibble
#' @import dplyr
#' @import tidyr
core_microbiome <- function(physeq, group_by = NULL, group = "All", model = "occ_abund", model_threshold = 2) {
  library(dplyr)

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

  # TODO: cycle multiple groups in hypothesis testing
  if (is.null(group_by)) {
    meta_table["Groups"] <- "All"
  } else {
    meta_table["Groups"] <- meta_table[, group_by]
  }

  if (!(group %in% meta_table$Groups)) {
    stop(sprintf("Couldn't subset samples, no group found with name %s in column %s", group, group_by))
  } else {
    meta_table <- meta_table[meta_table$Groups == group,]
  }

  if (length(unique(meta_table$Groups)) == 1) {
    label <- meta_table$Groups[1]
  }

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
  otu_rel <- apply(vegan::decostand(abund_table, method="total", MARGIN=2), 1, function(x) {
    if (all(is.na(x))) return(0)
    mean(x, na.rm = TRUE)
  })

  occ_abund_table <- data.frame(otu_occ=otu_occ, otu_rel=otu_rel) %>%           # combining occupancy and abundance data frame
    tibble::rownames_to_column('otu')

  # # Collate response and set class name
  # res <- list(abund_table = abund_table, occ_abund_table = occ_abund_table)
  # class(res) <- "ECCOccupancyAbundance"

  otu <- abund_table

  PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>%
    tidyr::gather(SampleID, abun, -otu) %>%
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

  # Modified original script to accept non-rarefied data

  # Pre-calculate the total counts for each sample from the full OTU table
  sample_totals <- colSums(otu)

  # --- 1. Start with the Top OTU ---

  otu_start <- otu_ranked$otu[1]
  start_matrix <- as.matrix(otu[otu_start, ])
  #start_matrix <- t(start_matrix)  # Now rows correspond to the first OTU's abundances per sample

  # Calculate Bray-Curtis dissimilarity for each pair of samples,
  # but normalize by the full sample totals
  x <- apply(combn(ncol(start_matrix), 2), 2, function(pair) {
    i1 <- pair[1]
    i2 <- pair[2]
    numerator <- sum(abs(start_matrix[, i1] - start_matrix[, i2]))
    denominator <- sample_totals[i1] + sample_totals[i2]
    numerator / denominator
  })
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(pair) {
    paste(colnames(start_matrix)[pair], collapse = ' - ')
  })
  df_s <- data.frame(x_names, x)
  names(df_s)[2] <- 1  # Label this column with the number of OTUs used (here, "1")
  BCaddition <- df_s  # Initialize our results data frame

  # --- 2. Loop Over Additional OTUs (from the 2nd to the 500th) ---

  for(i in 2:min(500, nrow(otu_ranked))) {
    # Ensure the OTU identifier is treated as a character string
    otu_add <- as.character(otu_ranked$otu[i])
    add_matrix <- as.matrix(otu[otu_add, ])
    #add_matrix <- t(add_matrix)

    # Append the new OTU's counts (as a row) to our cumulative matrix.
    start_matrix <- rbind(start_matrix, add_matrix)

    # Recalculate the Bray-Curtis dissimilarity over the cumulative subset,
    # but use the full sample totals for normalization.
    x <- apply(combn(ncol(start_matrix), 2), 2, function(pair) {
      i1 <- pair[1]
      i2 <- pair[2]
      numerator <- sum(abs(start_matrix[, i1] - start_matrix[, i2]))
      denominator <- sample_totals[i1] + sample_totals[i2]
      numerator / denominator
    })
    x_names <- apply(combn(ncol(start_matrix), 2), 2, function(pair) {
      paste(colnames(start_matrix)[pair], collapse = ' - ')
    })

    df_a <- data.frame(x_names, x)
    names(df_a)[2] <- i  # Label this column with the number of OTUs used (i)

    # Merge the new results into our cumulative results data frame
    BCaddition <- left_join(BCaddition, df_a, by = "x_names")
  }

  # --- 3. Calculate the Bray-Curtis Dissimilarity Using the Full OTU Table ---

  x <- apply(combn(ncol(otu), 2), 2, function(pair) {
    i1 <- pair[1]
    i2 <- pair[2]
    numerator <- sum(abs(otu[, i1] - otu[, i2]))
    denominator <- sample_totals[i1] + sample_totals[i2]
    numerator / denominator
  })
  x_names <- apply(combn(ncol(otu), 2), 2, function(pair) {
    paste(colnames(otu)[pair], collapse = ' - ')
  })
  df_full <- data.frame(x_names, x)
  names(df_full)[2] <- length(rownames(otu))  # Label with the total number of OTUs

  # Merge the full data results with the cumulative results
  BCfull <- left_join(BCaddition, df_full, by = "x_names")

  # # calculating BC dissimilarity based on the 1st ranked OTU
  # nReads=1000
  # otu_start=otu_ranked$otu[1]
  # start_matrix <- as.matrix(otu[otu_start,])
  # start_matrix <- t(start_matrix)
  # x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/(2*nReads))
  # x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  # df_s <- data.frame(x_names,x)
  # names(df_s)[2] <- 1
  # BCaddition <- rbind(BCaddition,df_s)
  # # calculating BC dissimilarity based on addition of ranked OTUs from 2nd to 500th.
  # # Can be set to the entire length of OTUs in the dataset, however it might take
  # # some time if more than 5000 OTUs are included.
  # for(i in 2:500){
  #   otu_add=otu_ranked$otu[i]
  #   add_matrix <- as.matrix(otu[otu_add,])
  #   add_matrix <- t(add_matrix)
  #   start_matrix <- rbind(start_matrix, add_matrix)
  #   x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads))
  #   x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  #   df_a <- data.frame(x_names,x)
  #   names(df_a)[2] <- i
  #   BCaddition <- left_join(BCaddition, df_a, by=c('x_names'))
  # }
  # # calculating the BC dissimilarity of the whole dataset (not needed if the second loop is already including all OTUs)
  # x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
  # x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))
  # df_full <- data.frame(x_names,x)
  # names(df_full)[2] <- length(rownames(otu))
  # BCfull <- left_join(BCaddition,df_full, by='x_names')

  rownames(BCfull) <- BCfull$x_names
  temp_BC <- BCfull
  temp_BC$x_names <- NULL
  temp_BC_matrix <- as.matrix(temp_BC)

  BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>%
    tidyr::gather(comparison, BC, -rank) %>%
    group_by(rank) %>%
    summarise(MeanBC=mean(BC)) %>%            # mean Bray-Curtis dissimilarity
    arrange(desc(-MeanBC)) %>%
    mutate(proportionBC=MeanBC/max(MeanBC))   # proportion of the dissimilarity explained by the n number of ranked OTUs
  Increase=BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
  increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))
  BC_ranked <- left_join(BC_ranked, increaseDF)
  BC_ranked <- BC_ranked[-nrow(BC_ranked),]

  # Creating thresholds for core inclusion
  # First (conservative) method is the elbow method for first-order difference - https://pommevilla.github.io/random/elbows.html
  # Second finds all taxa until increase in explanatory value reaches < 2%

  # elbow method (first order difference)
  fo_difference <- function(pos){
    left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos
    right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)
    return(left - right)
  }
  BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)

  fo_threshold <- which.max(BC_ranked$fo_diffs)
  # Get candidate ranks that meet the threshold for IncreaseBC
  pi_candidates <- as.numeric(as.character(BC_ranked$rank[BC_ranked$IncreaseBC >= 1 + (model_threshold / 100)]))
  if (length(pi_candidates) == 0) {
    warning("No OTUs meet the minimum increase threshold. Setting pi_threshold to the maximum available rank.")
    pi_threshold <- max(as.numeric(as.character(BC_ranked$rank)))
  } else {
    pi_threshold <- last(pi_candidates)
  }

  print(pi_candidates)
  print(fo_threshold)

  # We will also select ranked OTUs which meet the threshold criteria for MinCore designation.
  ranks <- vctrs::vec_c(factor("1"),
              BC_ranked$rank[BC_ranked$IncreaseBC >= 1 + (model_threshold / 100)])

  occ_abund_table$Type <- 'None'
  # If you want to include the first set of taxa as MinCore, you can uncomment this line:
  occ_abund_table$Type[occ_abund_table$otu %in% otu_ranked$otu[1:fo_threshold]] <- 'MinCore'

  # Note the use of parentheses to correctly create the sequence from the lower to the upper index.
  occ_abund_table$Type[occ_abund_table$otu %in% otu_ranked$otu[(fo_threshold + 1):pi_threshold]] <- 'Core'
  #occ_abund_table$Type[occ_abund_table$otu %in% otu_ranked$otu[ranks]] <- 'MinCore'

  # Give a title
  occ_abund_table$Title <- label

  # Use Sloan neutral model to prioritize OTUs
  # Fitting neutral model (Burns et al., 2016 (ISME J) - functions are in the sncm.fit.R)
  spp=t(otu)
  taxon=as.vector(rownames(otu))

  print(sprintf("Sample group has %d taxa", length(taxon)))
  if (length(taxon) > 4000) {
    warning("Binomial model may encounter issues when working with large numbers of taxa")
  }

  #Models for the whole community
  # obs.np=sncm.fit(spp, taxon, stats=FALSE, pool=NULL)
  # sta.np=sncm.fit(spp, taxon, stats=TRUE, pool=NULL)
  #
  # above.pred=sum(obs.np$freq > (obs.np$pred.upr), na.rm=TRUE)/sta.np$Richness  # fraction of OTUs above prediction
  # below.pred=sum(obs.np$freq < (obs.np$pred.lwr), na.rm=TRUE)/sta.np$Richness  # fraction of OTUs below prediction

  # Collate response and set class name
  res <- list(
    label = label,
    abund_table = abund_table,
    taxa_table = taxa_table,
    meta_table = meta_table,
    occ_abund_table = occ_abund_table,
    BC_ranked = BC_ranked,
    threshold1 = fo_threshold,
    threshold2 = pi_threshold,
    model_threshold = model_threshold,
    # obs.np = obs.np,
    # sta.np = sta.np,
    # above.pred = above.pred,
    # below.pred = below.pred,
    obs.np = NULL,
    sta.np = NULL,
    above.pred = NULL,
    below.pred = NULL,
    otu_ranked = otu_ranked
  )
  class(res) <- "ECCoreMicrobiome"

  return (res)
}

plot.ECCOccupancyAbundance <- function(value) {
  # Occupancy abundance plot:
  # p <- ggplot(data = value$occ_abund_table, aes(x = log10(otu_rel), y = otu_occ)) +
  #   geom_point(pch = 21, fill = 'white') +
  #   labs(x = "log10(mean relative abundance)", y = "Occupancy")

  p <- ggplot() +
    geom_point(data = value$occ_abund_table, aes(x = log10(otu_rel), y = otu_occ, fill = Type), pch=21, alpha=0.2) +
    scale_colour_manual(values = c("blue", "blue", "white"), aesthetics = c("fill")) +
    labs(x="log10(mean relative abundance)", y="Occupancy")

  if (!is.null(value$obs.np)) {
    p <- p +
      geom_line(color='black', data= value$obs.np, size=1, aes(y= value$obs.np$freq.pred, x=log10(value$obs.np$p)), alpha=.25) +
      geom_line(color='black', lty='twodash', size=1, data= value$obs.np, aes(y= value$obs.np$pred.upr, x=log10(value$obs.np$p)), alpha=.25)+
      geom_line(color='black', lty='twodash', size=1, data= value$obs.np, aes(y= value$obs.np$pred.lwr, x=log10(value$obs.np$p)), alpha=.25)
  }

  p <- p + facet_wrap(~ Title)

  p <- p + guides(fill="none")

    # geom_point(data = value$occ_abund_table[value$occ_abund_table$Type == 'None',], aes(x = log10(otu_rel), y = otu_occ), pch=21, fill='white', alpha=0.2) +
    # geom_point(data = value$occ_abund_table[value$occ_abund_table$Type != 'MinCore',], aes(x = log10(otu_rel), y = otu_occ), pch=21, fill='blue', size=1.8) +
    # geom_point(data = value$occ_abund_table[value$occ_abund_table$Type != 'Core',], aes(x = log10(otu_rel), y = otu_occ), pch=21, fill='lightblue', size=1.8) +

  #TODO: highlight on the same occupancy-abundance plot OTUs that are above and below the neutral model prediction
  #TODO: add to the plot the above.pred and below.pred values

  return(p)
}

plot.ECCRankedSimilarity <- function(value) {
  BC_ranked <- value$BC_ranked
  BC_ranked$Title <- value$label

  #Creating plot of Bray-Curtis similarity
  p <- ggplot(BC_ranked[1:100,], aes(x=factor(BC_ranked$rank[1:100], levels=BC_ranked$rank[1:100]))) +
    geom_point(aes(y=proportionBC)) +
    theme_classic() + theme(strip.background = element_blank(),axis.text.x = element_text(size=7, angle=45)) +
    geom_vline(xintercept = value$threshold, lty=3, col='red', cex=.5) +
    # geom_vline(xintercept=lastCall, lty=3, col='blue', cex=.5) +
    labs(x='Ranked Taxa',y='Contribution to Bray-Curtis Dissimilarity')

    # p <- p + annotate(geom = "text", x = value$threshold1 + 14, y = 0.1, label = sprintf("Elbow method (%d)", value$threshold1), color = "red")
    p <- p + annotate(geom = "text", x = value$threshold2, y = 0.5, label = sprintf("Last %d%% increase (%d)", value$model_threshold, value$threshold2), color = "blue")

  p <- p + facet_wrap(~ Title)

  return(p)
}

plot.ECCoreMicrobiome <- function(value, what_detection = "absolute", horizontal = TRUE) {
  prevalences <- seq(.05, 1, .05)

  #Creating a plot of core taxa occupancy by time point
  core <- value$occ_abund_table$otu[value$occ_abund_table$Type == 'MinCore']
  core_otu <- value$abund_table[rownames(value$abund_table) %in% core,]

  if (length(core) == 0) {
    stop("No taxa included in core")
  }

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
    detections <- 10^seq(log10(1), log10(max(microbiome::abundances(core_otu))/10), length = 20)
    detections <- round(detections)
  }

  if (!is.null(value$taxa_table)) {
    rownames(core_otu) <- expand_otu_names(rownames(core_otu), value$taxa_table, use_short_names = TRUE)
  }

  p <- microbiome::plot_core(
    core_otu,
    prevalences = prevalences,
    detections = detections,
    plot.type = "heatmap",
    colours = rev(RColorBrewer::brewer.pal(5, "Spectral")),
    # taxa.order = value$otu_ranked$otu,
    horizontal = horizontal,
  )

  p$data$Title <- value$label

  p <- p + facet_wrap(~ Title)
  #p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


  return(p)


  # otu_relabun <- vegan::decostand(otu, method="total", MARGIN=2)
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

taxa_bars <- function(physeq, group_by = NULL, N = 25, which_level = "Genus") {
  # Process input data
  abund_table <- phyloseq::otu_table(physeq)
  meta_table <- phyloseq::sample_data(physeq)
  OTU_taxonomy <- phyloseq::tax_table(physeq)
  OTU_tree <- phyloseq::phy_tree(physeq)

  # taxa should be columns
  if (phyloseq::taxa_are_rows(physeq)) {
    abund_table <- t(abund_table)
  }

  if (is.null(group_by)) {
    meta_table["Groups"] <- "All"
  } else {
    meta_table["Groups"] <- meta_table[, group_by]
  }

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

  df <- NULL
  for (i in 1:dim(new_x)[2]) {
    tmp <- data.frame(
      row.names = NULL,
      Sample = rownames(new_x),
      Taxa = rep(colnames(new_x)[i], dim(new_x)[1]),
      Value = new_x[,i],
      Groups = meta_table$Groups
    )

    if (i == 1) {
      df <- tmp
    } else {
      df <- rbind(df, tmp)
    }
  }
  colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00",rainbow(300)[seq(1,300,6)]);

  #The first step to get levels is to find all the possible values in a given column in df$sample by outputing the following command on terminal:
  # cat(paste("levels=c(",paste(paste("\"",unique(as.character(df$Sample)),"\"",sep=""),collapse=","),")",sep=""))
  # then use df$Sample<-factor(as.character(df$Sample),levels=c()) list

  meta_table <- as(sample_data(physeq), "data.frame")
  meta_table$Sample <- rownames(meta_table)

  df <- merge(df, meta_table, by = "Sample")

  return(structure(df, className = "ECTaxaBars"))
}

plot.ECTaxaBars <- function(value, N = 25, switch_strip = "x", legend_columns = 1, shorten = TRUE) {
  reveal_sample_names=TRUE
  legend_text_size=20
  axis_title_size=30
  text_size=30
  axis_text_size=30
  strip_text_size=16
  height_image=25
  width_image=60

  taxa_ids <- unique(value$Taxa)
  taxa_names <- get_taxa_names(taxa_ids, physeq, TRUE)

  value$Taxa <- ifelse(value$Taxa %in% taxa_ids, taxa_names[match(value$Taxa, taxa_ids)], value$Taxa)

  factor_levels <- c("Others", sort(unique(value$Taxa[value$Taxa != "Others"])))
  value$Taxa <- factor(value$Taxa, levels = factor_levels, ordered = TRUE)

  #colour_scale <- c("lightgrey", randomcoloR::distinctColorPalette(length(levels(value$Taxa)) - 1, altCol = TRUE))
  colour_scale <- c("lightgrey", colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n = length(levels(value$Taxa)) - 1))

  p <- ggplot(value, aes(Sample, Value, fill = Taxa)) +
    geom_bar(stat = "identity") +
    facet_grid(. ~ Groups, drop = TRUE, scale = "free", space = "free_x", switch = switch_strip) +
    scale_fill_manual(values = colour_scale) + #, guide = guide_legend(ncol = legend_columns)) +
    ylab("Proportions") +
    scale_y_continuous(expand = c(0.02, 0)) +
    theme_light() +
    theme(
      strip.background = element_rect(fill = "#cdcdcd"),
      strip.text = element_text(color = "black"),
      legend.position = "bottom"
      # strip.background = element_rect(fill="gray85"),
      # panel.spacing = unit(0.3, "lines"),
      # strip.text = element_text(size=strip_text_size,angle=90),
    )

  if(reveal_sample_names) {
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  } else {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks = element_blank())
  }

  return(p)
}

#Reference: https://microbiome.github.io/tutorials/Core.html
core_microbiome_legacy <- function(physeq, what_detection = "relative", taxa_rank = "Species", minimum_prevalence = 0.85, short_names = FALSE) {
  # Process input data
  # abund_table <- phyloseq::otu_table(physeq)
  # meta_table <- phyloseq::sample_data(physeq)
  # OTU_taxonomy <- phyloseq::tax_table(physeq)
  # OTU_tree <- phyloseq::phy_tree(physeq)

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

