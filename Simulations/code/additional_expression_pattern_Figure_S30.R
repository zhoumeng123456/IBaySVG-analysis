library(magick)
library(cowplot)
library(patchwork)

#############################################################
# meta_process()
#
# Description:
#   Preprocess spatial expression data for visualization.
#   Combines spatial coordinates and gene expression matrix,
#   and rescales spatial coordinates to [0, 1].
#
# Arguments:
#   position1 : N × 2 matrix of spatial coordinates (x, y).
#   count1    : Gene expression matrix (N × G).
#
# Returns:
#   A data.frame containing:
#     - Normalized spatial coordinates (x, y)
#     - Gene expression values for all genes
#
############################################################
meta_process =function(position1,count1){
  meta=cbind(position1,t(count1))
  colnames(meta)[c(1:2)] = c('x','y')
  meta=as.data.frame(meta)
  meta$x <- (meta$x - min(meta$x))/(max(meta$x) - min(meta$x))
  meta$y <- (meta$y - min(meta$y))/(max(meta$y) - min(meta$y))
  return(meta)
}

############################################################
# pattern_plot2()
#
# Description:
#   Visualizes the spatial expression pattern of a selected gene
#   using a scatter plot over spatial coordinates.
#
#   The function generates a ggplot object where each point
#   represents a spatial location and color encodes gene expression.
#
# Arguments:
#   pltdat    : Data frame containing spatial coordinates (x, y)
#               and gene expression values.
#   igene     : Index of the gene to visualize.
#   xy        : Logical. If FALSE, spatial coordinates are parsed
#               from rownames formatted as "xXy".
#   main      : Logical. Whether to display a plot title.
#   titlesize : Relative size of the title text.
#   pointsize : Size of spatial points.
#   xpand     : Expansion parameter for x-axis.
#   ypand     : Expansion parameter for y-axis.
#   title     : Optional custom title.
#
# Returns:
#   A ggplot object representing the spatial expression pattern.
#
############################################################
pattern_plot2 <- function(pltdat, igene, xy = T, main = F, titlesize = 2, 
                          pointsize = 3, xpand = 0, ypand = 1, title = NULL) {
  if (!xy) {
    xy <- matrix(as.numeric(do.call(rbind, strsplit(as.character(pltdat[, 
                                                                        1]), split = "x"))), ncol = 2)
    rownames(xy) <- as.character(pltdat[, 1])
    colnames(xy) <- c("x", "y")
    pd <- cbind.data.frame(xy, pltdat[, 2:ncol(pltdat)])
  } else {
    pd <- pltdat
  }
  
  pal <- colorRampPalette(c("mediumseagreen", "lightyellow2", "deeppink"))
  gpt <- ggplot(pd, aes(x = x, y = y, color = pd[, igene + 2])) + geom_point(size = pointsize) + 
    # scale_color_gradientn(colours=pal(5))+
    scale_color_gradientn(colours = pal(5)) + scale_x_discrete(expand = c(xpand, 
                                                                          ypand)) + scale_y_discrete(expand = c(xpand, ypand)) + coord_equal() + 
    # labs(title = colnames(pd)[igene+2], x = NULL, y = NULL)+
    #theme_bw()
    theme_void()
  
  if (main) {
    if (is.null(title)) {
      title = colnames(pd)[igene + 2]
    }
    out = gpt + labs(title = title, x = NULL, y = NULL) + theme(legend.position = "none", 
                                                                plot.title = element_text(hjust = 0.5, size = rel(titlesize)))
  } else {
    out = gpt + labs(title = NULL, x = NULL, y = NULL) + theme(legend.position = "none")
  }
  return(out)
}



###reproduce the results

##generate the simulated data by the sim_create() or normal_create() illustrated in generate_simu.R
result=sim_create(xspace = "polynomial1",yspace ="polynomial1",gene_size = 20,sv_mark = c(1,1),svgene_size=0.5,inf_size = 0.01,
                   seed =2,domainnum = 2,use_covariate=FALSE)#result=normal_create(inf_size = 0.01,domainnum = 2,kernel="exp",svgene_size=0.5,seed=1,gene_size=20,tau1=9,tau2=1,num_cores =8)

##plot the spatial pattern and select the proper one
meta1=meta_process(result[[2]],t(result[[1]]))
for(genenum in c(1:10)){
  p1=pattern_plot2(meta1,genenum,xpand=0.01,ypand =0.01,main = FALSE)
  print(p1)
}
ggsave(paste0("select_pattern.png"), plot = p1, width = 4, height = 4, dpi = 300)


#load and combine the results after generating nine additional spatial structures.
name1=c("linearfocal","focalperiod","linearperiod","ZINNGP","sigmoid","polynomial1","polynomial2","polynomial3","polynomial4")
img_paths <-"Simulations/result_data/additional spatial pattern"
images <- lapply(1:9, function(genenum) {
  img <- image_read(here::here(img_paths,paste0(name1[genenum],".png")))
})

name2=c("Linear-Focal","Focal-Period","Linear-Period","ZINNGP","Sigmoid","Polynomial1","Polynomial2","Polynomial3","Polynomial4")
image_plots <- lapply(seq_along(images), function(i) {
  ggdraw() + 
    draw_image(images[[i]]) +  
    draw_label(
      name2[i],              
      x = 0.5, y = 1.06, 
      hjust = 0.5, vjust = 1,
      size = 9, 
      fontface = "bold"
    )+
    theme(
      plot.margin = margin(10, 0, 10, 0)  
    )
})

combined_plot <- 
  (wrap_elements(full = (image_plots[[4]] | image_plots[[1]] | image_plots[[2]])) /
     wrap_elements(full = (image_plots[[3]] |image_plots[[5]] | image_plots[[6]]))/
     wrap_elements(full = (image_plots[[7]] |image_plots[[8]] | image_plots[[9]])))+
  plot_layout(guides = "collect")+
  plot_annotation(
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    tag_levels = list(c("A", "B","C")), 
    tag_prefix = "",
    tag_sep = "",
    theme = theme(plot.tag = element_text(size = 12, face = "bold")
    )  
  )

print(combined_plot)















