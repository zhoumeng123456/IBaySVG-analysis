
#produce the plot of spatial expression pattern of genes by SRTsim, an independent, reproducible, and flexible SRT simulation framework.
#available at github library 'xzhoulab/SRTsim'
shinySRT1 <- SRTsim_shiny()
simSRT1 <- Shiny2SRT(shinySRT1)
count=simSRT1@simCounts
count=as.matrix(count)
meta1=meta_process(position,count)
g=pattern_plot2(meta1,1,xpand=-2,ypand = 2,main = FALSE)
ggsave(paste0("period.png"), plot = g, width = 6, height = 6, dpi = 150)



##read the plot
name1=c("linear","focal","period","complex","linear_gene","focal_gene","period_gene","complex_gene")
img_paths <-"RealData/result_data/plot of figures1/"
images <- lapply(1:8, function(genenum) {
  img <- image_read(here::here(img_paths,paste0(name1[genenum],".png")))
})


##layout and save
name2=c("Linear","Focal","Period","Complex","AQP4","COX6C","CAMK2N1","AGR2")
image_plots <- lapply(seq_along(images), function(i) {
  ggdraw() + 
    draw_image(images[[i]]) +  
    draw_label(
      name2[i],              
      x = 0.5, y = 0.95, 
      hjust = 0.5, vjust = 1,
      size = 10, 
      color = "black",
      fontface = "bold"
    )
})

combined_plot <- (
  wrap_elements(full = (image_plots[[1]] | image_plots[[2]] | image_plots[[3]]|image_plots[[4]])) /
    wrap_elements(full = (image_plots[[5]] |image_plots[[6]] | image_plots[[7]]|image_plots[[8]] ))
) +
  plot_layout(guides = "collect",heights = c(0.8,0.8))+
  plot_annotation(
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    tag_levels = list(c("A", "B")),  
    tag_prefix = "",
    tag_sep = "",
    theme = theme(plot.tag = element_text(size = 14, face = "bold"))  
  )

print(combined_plot)

ggsave(paste0("combined_fig.png"), plot = combined_plot, width = 8, height = 5.5, dpi = 300)


