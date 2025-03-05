# library(phyloseq)
# library(vegan)
# library(ggplot2)
# library(ape)
# library(phangorn)
# library(stringr)
# library(grid)

beta_diversity <- function(physeq, PERMANOVA_variables, distance_metrics = c("bray", "unifrac", "wunifrac"), kind = "se", confidence_limit = 0.95) {
  #Reference: http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
  #Data frame df_ell contains values to show ellipses. It is calculated with function veganCovEllipse which is hidden in vegan package. This function is applied to each level of NMDS (group) and it uses also function cov.wt to calculate covariance matrix.
  veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100)
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }

  if(is.null(sample_data(physeq)$Groups)) {
    stop("Groups column not set in sample data")
  }

  # Check to see if the meta_table$Type2 is assigned or not
  if(is.null(sample_data(physeq)$Type2)) {
    sample_data(physeq)$Type2 <- sample_data(physeq)$Groups
  }

  meta_table <- as(sample_data(physeq), "data.frame") # this method is required to coerce sample_data to data.frame

  if (class(meta_table) != "data.frame") {
    stop(sprintf("Metadata could not be coerced to the correct format; adonis2 requires a 'data.frame' and was given a '%s'", class(meta_table)))
  }

  mds <- list()
  for (d in distance_metrics) {
    # Check we've been given a supported distance metric and error if not
    # if (!(d %in% names(metrics))) {
    #   warning(sprintf("Unknown distance metric %s, passing over.", d))
    #   next
    # }
    # if (taxa_rank != "Species" && d != "bray") {
    #   warning(sprintf("Unsupported taxa rank %s used with distance metric %s, passing over.", taxa_rank, d))
    #   next
    # }

    dist <- phyloseq::distance(physeq, d)
    sol <- cmdscale(dist, eig = TRUE)
    # aov <- vegan::adonis(as.formula(paste("dist ~", paste(PERMANOVA_variables, collapse = "+"))), data = meta_table[rownames(phyloseq::otu_table(physeq)), ])
    aov <- vegan::adonis(as.formula(paste("dist ~", paste(PERMANOVA_variables, collapse = "+"))), data = meta_table)

    if (!is.null(sol)) {
      sol <- append(sol, list(metric = d))
      sol <- append(sol, list(aov = aov))
      mds <- append(mds, list(sol))
    }
  }

  PCOA <- NULL
  PCOA_lines <- NULL
  df_ord <- NULL

  for (sol in mds) {
  # if(!is.null(sol)) {
    if (is.null(PCOA)) {
      PCOA <- data.frame(x = sol$points[,1], y = sol$points[,2], metric = sol$metric, meta_table)
    } else {
      PCOA <- rbind(PCOA, data.frame(x = sol$points[,1], y = sol$points[,2], metric = sol$metric, meta_table))
    }

    plot.new()
    ord<-vegan::ordiellipse(sol, interaction(meta_table$Groups,meta_table$Type2),display = "sites", kind = kind, conf = confidence_limit, label = T)
    dev.off()


    #Generate ellipse points
    df_ell <- data.frame()
    for(h in levels(PCOA$Groups)) {
      for(g in levels(PCOA$Type2)) {
        if(paste(h, ".", g, sep="") != "" && (paste(h, ".", g, sep="") %in% names(ord))) {
          tryCatch(
            {
              df_ell <- rbind(
                df_ell,
                cbind(
                  as.data.frame(
                    with(
                      PCOA[PCOA$Groups == h & PCOA$Type2 == g, ],
                      veganCovEllipse(ord[[paste(h, ".", g, sep="")]]$cov,
                                      ord[[paste(h, ".", g, sep="")]]$center,
                                      ord[[paste(h, ".", g, sep="")]]$scale)
                    )
                  ),
                  Groups = h,
                  Type2 = g,
                  metric = as.character(sol$metric)
                )
              )
            },
            error = function(e) warning(e)
          )
        }
      }
    }

    if (sum(dim(df_ell))>0){
      colnames(df_ell)<-c("x","y","Groups","Type2")
    }
    df_ell$Groups<-factor(df_ell$Groups,levels=levels(PCOA$Groups))
    df_ell$Type2<-factor(df_ell$Type2,levels=levels(PCOA$Type2))

    if (is.null(df_ord)) {
      df_ord <- df_ell
    } else {
      df_ord <- rbind(df_ord, df_ell)
    }

    #Generate mean values from PCOA plot grouped on
    PCOA.mean=aggregate(PCOA[,1:2],list(group=PCOA$Groups),mean)

    #Connecting samples based on lines with meta_table$Connections and meta_table$Subconnections ###########
    # PCOA_lines<-NULL

    #Check if meta_table$Subconnections exists
    if(!is.null(meta_table$Subconnections)){
      #Step 1, populate Connections_IDs and get rid of singletons using PCOA$Connections
      Connections_IDs<-as.character(PCOA$Connections)
      Connections_IDs<-names(table(Connections_IDs)[table(Connections_IDs)>1])
      Subconnections_IDs_mask<-as.character(sapply(Connections_IDs,function(x){if(length(unique(PCOA[PCOA$Connections==x,"Subconnections"]))<2) x else "___APPROVED__"}))
      #Step 2, loop through each Connection_IDs and see we can find multiple PCOA$Subconnections and then filter Connections_IDs further
      Connections_IDs<-Connections_IDs[!Connections_IDs %in% Subconnections_IDs_mask]
      Connections_IDs<-sort(Connections_IDs)
      if(length(Connections_IDs)>0){
        for(p in Connections_IDs){
          Subconnections_IDs<-unique(PCOA[PCOA$Connections==p,"Subconnections"])
          if(is.factor(Subconnections_IDs)){
            Subconnections_IDs<-as.character(Subconnections_IDs)
          } else {
              Subconnections_IDs<-sort(Subconnections_IDs)
          }
          S<-NULL
          if(pairwise_connections_and_not_longitudinal_connections){
          #Get pair-wise combinations from Subconnections_IDs
           S<-combn(Subconnections_IDs,2)
          } else{
          #Get longitudinal connections from Subconnections_IDs
            S<-t(cbind(Subconnections_IDs[-length(Subconnections_IDs)],Subconnections_IDs[-1]))
          }
          for(ii in 1:ncol(S)){
            tmp<-data.frame(t(colMeans(PCOA[PCOA$Connections==p & PCOA$Subconnections==S[1,ii],c("x","y"),drop=F])),t(colMeans(PCOA[PCOA$Connections==p & PCOA$Subconnections==S[2,ii],c("x","y"),drop=F])),p)
            colnames(tmp)<-c("xfrom","yfrom","xto","yto","ID")
            if(is.null(PCOA_lines)){PCOA_lines<-tmp} else {PCOA_lines<-rbind(PCOA_lines,tmp)}
          }
        }
      }
    } else {
      #To connect lines between samples we need to extract connections
      Connections_IDs<-as.character(PCOA$Connections)
      #Next we filter out Connections_IDs that are singletons and also uniquify them
      Connections_IDs<-names(table(Connections_IDs)[table(Connections_IDs)>1])
      Connections_IDs<-sort(Connections_IDs)
      if(length(Connections_IDs)>0){
        #We iterate through the IDs one at a time
        for(p in Connections_IDs){
          rownames_list<-rownames(PCOA[PCOA$Connections %in% p,,drop=F])
          S<-NULL
          if(pairwise_connections_and_not_longitudinal_connections){
            #Get pair-wise combinations from rownames_list
            S<-combn(rownames_list,2)
          } else{
            #Get longitudinal connections from rownames_list
            S<-t(cbind(rownames_list[-length(rownames_list)],rownames_list[-1]))
          }
          for(ii in 1:ncol(S)){
            tmp<-cbind(PCOA[S[1,ii],c("x","y")],PCOA[S[2,ii],c("x","y")],p)
            colnames(tmp)<-c("xfrom","yfrom","xto","yto","ID")
            if(is.null(PCOA_lines)){PCOA_lines<-tmp} else {PCOA_lines<-rbind(PCOA_lines,tmp)}
          }
        }
      }
    }
    #/#Connecting samples based on lines with meta_table$Connections and meta_table$Subconnections ###########
  }

  # dist<-phyloseq::distance(physeq,which_distance)
  # capture.output(adonis(as.formula(paste("dist ~",paste(PERMANOVA_variables,collapse="+"))), data=meta_table[rownames(otu_table(physeq)),]),file=paste("ADONIS_",which_distance,"_",which_level,"_",label,".txt",sep=""))

  colnames(df_ord)[5] <- "metric"    #weird bug here

  return(structure(list(
    df = PCOA,
    PCOA_lines = PCOA_lines,
    df_ord = df_ord,
    mds = mds,
    PERMANOVA_variables = PERMANOVA_variables,
    meta_table = meta_table
  ), className = "ECBetaDiversity"))
}

#' @import ggplot2
plot.ECBetaDiversity <- function(value) {
  df = value$df
  PCOA_lines = value$PCOA_lines
  df_ord = value$df_ord
  mds = value$mds
  PERMANOVA_variables = value$PERMANOVA_variables
  meta_table = value$meta_table

  point_size = 2
  point_opacity = 0.25
  # number_of_rows = 1
  # use_provided_colors = FALSE

  draw_glow=FALSE
  draw_mean_values_text=FALSE
  draw_confidence_intervals=TRUE
  draw_ellipses_and_not_polygons=FALSE

  opacity_ellipses_polygons=0.2
  linesize_ellipses_polygons=0.8
  linetype_ellipses_polygons="solid" #blank solid dashed dotted dotdash longdash twodash

  exclude_legends=FALSE

  #coloring function
  gg_color_hue<-function(n){
    hues=seq(15,375,length=n+1)
    hcl(h=hues,l=65,c=100)[1:n]
  }

  cols=gg_color_hue(length(unique(df$Groups)))

  plots <- list()
  for (metric in unique(df$metric)) {
    PCOA <- df[df$metric == metric, ]
    df_ell <- df_ord[df_ord$metric == metric, ]



    p <- ggplot(data = PCOA, aes(x, y, colour = Groups))
    p<-p+facet_grid(~ metric, scales = "free", labeller = labeller(metric = as_labeller(metrics)))

    if(!"Type" %in% colnames(meta_table)){
      p<-p + geom_point(aes(x, y, colour = Groups),inherit.aes=F,alpha=point_opacity,size=point_size)
    } else {
      p<-p + geom_point(aes(x, y, colour = Groups, shape = Type),inherit.aes=F,alpha=point_opacity,size=point_size)
      p<-p+scale_shape("Type")
    }

    if(draw_glow){
      p<-p + geom_point(alpha=point_glow_opacity,size = point_size+point_glow_differential,show.legend=FALSE)
    }

    p <- plot_theme_default(p)


    if(draw_mean_values_text){
      p<-p+ annotate("text",x=PCOA.mean$x,y=PCOA.mean$y,label=PCOA.mean$group,size=mean_values_text_size,colour=cols,alpha=mean_values_text_opacity,vjust=0.3)
    }
    if (sum(dim(df_ell))>0){
    if(draw_confidence_intervals){
      if(identical(levels(meta_table$Groups), levels(meta_table$Type2))){
        if(draw_ellipses_and_not_polygons){
          p<-p+ geom_path(data=df_ell, aes(x=x, y=y, group = metric), size=linesize_ellipses_polygons, linetype=linetype_ellipses_polygons,alpha=opacity_ellipses_polygons, show.legend = FALSE)
        } else {
          p<-p+ geom_polygon(data=df_ell, aes(x=x, y=y, group = Groups, fill=Groups), size=linesize_ellipses_polygons, linetype=linetype_ellipses_polygons,alpha=opacity_ellipses_polygons, show.legend = FALSE)
        }
      } else {
        for (q in unique(as.numeric(meta_table$Groups))){
          for(r in unique(as.numeric(meta_table$Type2))){
            if(draw_ellipses_and_not_polygons){
              p<-p+ geom_path(data=df_ell[as.numeric(df_ell$Groups)==q & as.numeric(df_ell$Type2)==r,], aes(x=x, y=y), size=linesize_ellipses_polygons, colour=cols[q],linetype=r,alpha=opacity_ellipses_polygons, show.legend = FALSE)
            } else {
              p<-p+ geom_polygon(data=df_ell[as.numeric(df_ell$Groups)==q & as.numeric(df_ell$Type2)==r,], aes(x=x, y=y,fill=Groups), size=linesize_ellipses_polygons, colour=cols[q],linetype=r,alpha=opacity_ellipses_polygons, show.legend = FALSE)
            }
            }
        }
      }
    }
    }
    if(exclude_legends){
      p<-p+guides(colour=FALSE)
    }

    #only draw lines connecting dots if the lines are available
    if(!is.null(PCOA_lines)){
      arrow<-NULL
      if(should_connections_end_in_arrows){
        arrow=arrow(length=unit(0.2,"inches"))
      }
      p<-p+geom_segment(data=PCOA_lines,inherit.aes=FALSE,aes(x=xfrom,y=yfrom,xend=xto,yend=yto),colour="grey20",size=linking_samples_line_size,alpha=linking_samples_line_opacity,linetype=linking_samples_linetype, show.legend = FALSE,arrow=arrow)
    }

    sol <- mds[[which(sapply(mds, function(x) x$metric == metric))]]
    if (!is.null(sol)) {
      p <- p +
        xlab(sprintf("Dim 1 (%.4g%%)", sol$eig[1] / sum(sol$eig) * 100)) +
        ylab(sprintf("Dim 2 (%.4g%%)", sol$eig[2] / sum(sol$eig) * 100))
    }

    # Display the PERMANOVA stats: pipe aov table to a string and display in a text grob
    str <- paste(capture.output(sol$aov$aov.tab[PERMANOVA_variables, c("R2", "Pr(>F)")]), collapse = "\n")
    text <- grid::textGrob(paste("PERMANOVA", str, sep = "\n"), gp = grid::gpar(fontsize = 12, fontfamily = "mono"))

    plots <- append(plots, list(p, text))
  }

  # Layout matrix set column-wise to display text grobs underneath plots
  layout_matrix <- matrix(seq_len(length(plots)), nrow = 2, byrow = FALSE)

  obj <- do.call(gridExtra::grid.arrange, c(plots, list(layout_matrix = layout_matrix, heights = c(2, 1))))

  structure(list(gtable = obj), class = "BetaDiversityPlot")
}

print.BetaDiversityPlot <- function(x) {
  grid::grid.draw(x$gtable)
}

beta_diversity_write <- function() {
  # pdf(paste("PCOA_",which_distance,"_",which_level,"_",label,".pdf",sep=""),width=width_image,height=height_image)
  # print(p)
  # dev.off()
  #
  # dist<-phyloseq::distance(physeq,which_distance)
  # capture.output(adonis(as.formula(paste("dist ~",paste(PERMANOVA_variables,collapse="+"))), data=meta_table[rownames(otu_table(physeq)),]),file=paste("ADONIS_",which_distance,"_",which_level,"_",label,".txt",sep=""))
  #
}
