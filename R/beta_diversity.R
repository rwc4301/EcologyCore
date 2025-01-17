# library(phyloseq)
# library(vegan)
# library(ggplot2)
# library(ape)
# library(phangorn)
# library(stringr)
# library(grid)

beta_diversity_analysis <- function(abund_table, meta_table, OTU_taxonomy, OTU_tree, grouping_column) {
  #Convert the data to phyloseq format
  OTU = otu_table(as.matrix(abund_table), taxa_are_rows = FALSE)
  TAX = tax_table(as.matrix(OTU_taxonomy))
  SAM = sample_data(meta_table)

  physeq<-NULL
  if(which_level=="Otus"){
    #physeq<-merge_phyloseq(phyloseq(OTU, TAX),SAM,midpoint(OTU_tree))
    physeq<-merge_phyloseq(phyloseq(OTU, TAX),SAM,OTU_tree)
  } else {
    physeq<-merge_phyloseq(phyloseq(OTU),SAM)
  }



  #Reference: http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
  #Data frame df_ell contains values to show ellipses. It is calculated with function veganCovEllipse which is hidden in vegan package. This function is applied to each level of NMDS (group) and it uses also function cov.wt to calculate covariance matrix.
  veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100)
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }

  #coloring function
  gg_color_hue<-function(n){
    hues=seq(15,375,length=n+1)
    hcl(h=hues,l=65,c=100)[1:n]
  }

  sol<-NULL
  if(which_distance=="bray"){
    sol<-cmdscale(phyloseq::distance(physeq,"bray"),eig=T)
  } else if(which_distance=="wunifrac" & which_level=="Otus") {
    sol<-cmdscale(phyloseq::distance(physeq,"wunifrac"),eig=T)
  } else if(which_distance=="unifrac" & which_level=="Otus"){
    sol<-cmdscale(phyloseq::distance(physeq,"unifrac"),eig=T)
  }

  #Check to see if the meta_table$Type2 is assigned or not
  if(is.null(meta_table$Type2)){meta_table$Type2<-meta_table$Groups}

  if(!is.null(sol)){
    PCOA=data.frame(x=sol$points[,1],y=sol$points[,2],meta_table)

    plot.new()
    ord<-ordiellipse(sol, interaction(meta_table$Groups,meta_table$Type2),display = "sites", kind = kind, conf = 0.95, label = T)
    dev.off()


    #Generate ellipse points
    df_ell <- data.frame()
    for(h in levels(PCOA$Groups)){
      for(g in levels(PCOA$Type2)){
        if(paste(h,".",g,sep="")!="" && (paste(h,".",g,sep="") %in% names(ord))){
          tryCatch(df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCOA[PCOA$Groups==h & PCOA$Type2==g,],
                                                        veganCovEllipse(ord[[paste(h,".",g,sep="")]]$cov,ord[[paste(h,".",g,sep="")]]$center,ord[[paste(h,".",g,sep="")]]$scale)))
                                     ,Groups=h,Type2=g)),error=function(e) NULL)
     }
    }
    }

    if (sum(dim(df_ell))>0){
      colnames(df_ell)<-c("x","y","Groups","Type2")
    }
    df_ell$Groups<-factor(df_ell$Groups,levels=levels(PCOA$Groups))
    df_ell$Type2<-factor(df_ell$Type2,levels=levels(PCOA$Type2))

    #Generate mean values from PCOA plot grouped on
    PCOA.mean=aggregate(PCOA[,1:2],list(group=PCOA$Groups),mean)

    #Connecting samples based on lines with meta_table$Connections and meta_table$Subconnections ###########
    PCOA_lines<-NULL

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

    return(PCOA)
  }

  return(NULL)
}

beta_diversity_plot <- function(PCOA) {
  cols=gg_color_hue(length(unique(PCOA$Groups)))

  p<-ggplot(data=PCOA,aes(x,y,colour=Groups))
  if(!"Type" %in% colnames(meta_table)){
    p<-p + geom_point(aes(PCOA$x,PCOA$y,colour=PCOA$Groups),inherit.aes=F,alpha=point_opacity,size=point_size)
  } else{
    p<-p + geom_point(aes(PCOA$x,PCOA$y,colour=PCOA$Groups, shape=PCOA$Type),inherit.aes=F,alpha=point_opacity,size=point_size)
    p<-p+scale_shape("Type")
  }

  if(draw_glow){
    p<-p + geom_point(alpha=point_glow_opacity,size = point_size+point_glow_differential,show.legend=FALSE)
  }

  p<-p+theme_bw()
  if(draw_mean_values_text){
    p<-p+ annotate("text",x=PCOA.mean$x,y=PCOA.mean$y,label=PCOA.mean$group,size=mean_values_text_size,colour=cols,alpha=mean_values_text_opacity,vjust=0.3)
  }
  if (sum(dim(df_ell))>0){
  if(draw_confidence_intervals){
    if(identical(levels(meta_table$Groups), levels(meta_table$Type2))){
      if(draw_ellipses_and_not_polygons){
        p<-p+ geom_path(data=df_ell, aes(x=x, y=y), size=linesize_ellipses_polygons, linetype=linetype_ellipses_polygons,alpha=opacity_ellipses_polygons, show.legend = FALSE)
      } else {
        p<-p+ geom_polygon(data=df_ell, aes(x=x, y=y,fill=Groups), size=linesize_ellipses_polygons, linetype=linetype_ellipses_polygons,alpha=opacity_ellipses_polygons, show.legend = FALSE)
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


  p<-p+xlab(paste("Dim1 (",sprintf("%.4g",(sol$eig[1]/sum(sol$eig))*100),"%)",sep=""))+ylab(paste("Dim2 (",sprintf("%.4g",(sol$eig[2]/sum(sol$eig))*100),"%)",sep=""))

  p<-p+theme(legend.title=element_text(size=legend_title_size),
             legend.text=element_text(size=legend_text_size),
             text = element_text(size=text_size),
             axis.text=element_text(size=axis_text_size),
             axis.title=element_text(size=axis_title_size))


  if(use_provided_colors){
    p<-p+scale_color_manual("Groups",values=colours)
    p<-p+scale_fill_manual("Groups",values=colours)
  }


  pdf(paste("PCOA_",which_distance,"_",which_level,"_",label,".pdf",sep=""),width=width_image,height=height_image)
  print(p)
  dev.off()

  dist<-phyloseq::distance(physeq,which_distance)
  capture.output(adonis(as.formula(paste("dist ~",paste(PERMANOVA_variables,collapse="+"))), data=meta_table[rownames(otu_table(physeq)),]),file=paste("ADONIS_",which_distance,"_",which_level,"_",label,".txt",sep=""))
}

beta_diversity_write <- function() {
  # TBD
}
