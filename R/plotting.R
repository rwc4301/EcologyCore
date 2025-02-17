# height_image=12
# width_image=25
# use_provided_colors=TRUE
# colours<-c("#F8766D","#F8766D","#A3A500","#A3A500","#00BF7D","#00BF7D","#00B0F6","#00B0F6","#E76BF3","#E76BF3")


# legend_text_size=6
# legend_title_size=8
# text_size=12

format <- list("svg", "pdf", "html", "png")

colours <- c("#ffa172", "#81fc76", "#68aeff","#c320d8","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00",grey.colors(1000))

plot_theme_default <- function(
  plot,                                 # ggplot object to theme
  text_size = 16,
  point_size = 5,
  point_opacity = 0.8,
  axis_text_size = 12,
  strip_text_size = 16,
  increment_divider = 2,
  exclude_pvalues_text_from_drawing = FALSE,
  legend.position = "left", # left | right | top | bottom
  exclude_legends = TRUE,
  pairwise_text_size = 7,
  number_of_rows = 1,
  legend_text_size = 12,
  legend_title_size = 30,
  axis_title_size = 16,
  height_image = 10,
  width_image = 25,
  use_provided_colors = FALSE
  #colours <- c("#ffa172", "#81fc76", "#68aeff","#c320d8","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00",grey.colors(1000))
) {
  #if crashes at panel.margin change it to panel.spacing, if crashes at panel.spacing, change it to panel.margin
  plot <- plot + theme_light() + theme(
    # Axis theme
    axis.text = element_text(size = axis_text_size),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_text(size = axis_title_size),
    # Legend theme
    legend.text = element_text(size = legend_text_size),
    legend.title = element_text(size = legend_title_size),
    # Panel theme
    panel.spacing = unit(2, "lines"),
    # Strip theme
    strip.background = element_rect(fill = "#cdcdcd"),
    strip.text = element_text(size = strip_text_size, color = "black"),
    # Text theme
    text = element_text(size = text_size),
  )

  if(legend.position == "bottom") {
    plot <- plot + theme(
      legend.key = element_blank(),  #removes the box around each legend item
      legend.position = "bottom", #legend at the bottom
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.box.just = "centre"
    )
  }

  if(exclude_legends) {
    plot <- plot + guides(colour = FALSE)
  }


  if(use_provided_colors){
    p<-p+scale_color_manual("Groups",values=colours)
    p<-p+scale_fill_manual("Groups",values=colours)
  }

  return(plot)
}

save_plot <- function(plot, filename, format = "pdf", height, width) {
  switch(format,
    "pdf" = {
      pdf(filename, height = height, width = width)
      print(plot)
      dev.off()
    },
    "svg" = {

    },
    "html" = {

    },
    "png" = {

    }
  )
}
