library(scatterpie)


##plot for dlpfc dataset
#function for plot the scatterpie
plot_scatterpie <- function(position, calpha, pie_scale = 0.4) {
  
  col_manual <- c(
    "Inhib" = "#FFDDDD",
    "Oligo" = "#B0E0E6",
    "OPC" ="#FFFACD" ,
    "Excit" ="#32CD32" ,
    "MicroOligo" = "#FFEB99",
    "Astro" = "#8A2BE2",
    "EndoMural" = "#FF1493"
  )
  focus_types <- c("Astro", "EndoMural", "OPC")
  calpha[, focus_types] <- calpha[, focus_types] * 2 
  calpha <- calpha / rowSums(calpha)
  
  oligo_mask <- calpha[, "Oligo"] < 0.75
  oligo_mask1 <- calpha[, "Oligo"] > 0.75
  calpha[oligo_mask, "Oligo"] <- calpha[oligo_mask, "Oligo"] * 0.6  
  calpha[oligo_mask1, "MicroOligo"] <- calpha[oligo_mask1, "MicroOligo"] * 0.6  
  calpha[oligo_mask1, "Oligo"] <- calpha[oligo_mask1, "Oligo"] * 1.2 
  
  calpha <- calpha / rowSums(calpha)
  
  meta <- cbind(position, calpha)
  colnames(meta)[1:2] <- c('x', 'y')
  meta <- as.data.frame(meta)
  
  meta$x <- (meta$x - min(meta$x)) / (max(meta$x) - min(meta$x))
  meta$y <- (meta$y - min(meta$y)) / (max(meta$y) - min(meta$y))
  
  ggplot() +
    geom_scatterpie(
      data = meta, 
      aes(x = x, y = y),
      cols = colnames(meta)[3:9], 
      color = NA, 
      pie_scale = pie_scale
    ) +
    coord_fixed(ratio = 1) +
    scale_fill_manual(values = col_manual) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.border = element_blank(), 
      legend.text = element_text(size = 14),  
      plot.tag = element_text(size = 18,face="bold")
    ) +
    labs(x = NULL, y = NULL)  
}

#enhance the contrast
enhance_contrast <- function(calpha, power=3) {
  calpha_trans <- calpha^power
  calpha_trans / rowSums(calpha_trans)
}


#load the data for dlpfc with samedonor
dir="data/Realdataset/dlpfc samedonor"
position4=as.matrix(read.csv(here::here(dir,"matrix4_position_samedonor.csv"),row.names = 1))
calpha4=as.matrix(read.csv(here::here(dir,"matrix4_celltype_samedonor.csv"),row.names = 1))
position3=as.matrix(read.csv(here::here(dir,"matrix3_position_samedonor.csv"),row.names = 1))
calpha3=as.matrix(read.csv(here::here(dir,"matrix3_celltype_samedonor.csv"),row.names = 1))
position2=as.matrix(read.csv(here::here(dir,"matrix2_position_samedonor.csv"),row.names = 1))
calpha2=as.matrix(read.csv(here::here(dir,"matrix2_celltype_samedonor.csv"),row.names = 1))
position1=as.matrix(read.csv(here::here(dir,"matrix1_position_samedonor.csv"),row.names = 1))
calpha1=as.matrix(read.csv(here::here(dir,"matrix1_celltype_samedonor.csv"),row.names = 1))

#or load the data for dlpfc with acrossdonor(share the common enhance and produce function)
dir="data/Realdataset/dlpfc acrossdonor"
position4=as.matrix(read.csv(here::here(dir,"matrix4_position_acrossdonor.csv"),row.names = 1))
calpha4=as.matrix(read.csv(here::here(dir,"matrix4_celltype_acrossdonor.csv"),row.names = 1))
position3=as.matrix(read.csv(here::here(dir,"matrix3_position_acrossdonor.csv"),row.names = 1))
calpha3=as.matrix(read.csv(here::here(dir,"matrix3_celltype_acrossdonor.csv"),row.names = 1))
position2=as.matrix(read.csv(here::here(dir,"matrix2_position_acrossdonor.csv"),row.names = 1))
calpha2=as.matrix(read.csv(here::here(dir,"matrix2_celltype_acrossdonor.csv"),row.names = 1))
position1=as.matrix(read.csv(here::here(dir,"matrix1_position_acrossdonor.csv"),row.names = 1))
calpha1=as.matrix(read.csv(here::here(dir,"matrix1_celltype_acrossdonor.csv"),row.names = 1))


#enhance and produce the plot
calpha1 <- enhance_contrast(calpha1, power=4)  
calpha2 <- enhance_contrast(calpha2, power=4)  
calpha3 <- enhance_contrast(calpha3, power=4)  
calpha4 <- enhance_contrast(calpha4, power=4)  

p1 <- plot_scatterpie(position1, calpha1)
p2 <- plot_scatterpie(position2, calpha2)
p3 <- plot_scatterpie(position3, calpha3)
p4 <- plot_scatterpie(position4, calpha4)


#layout and save
combined_plot <- (p1 | p2 | p3 | p4) +
  plot_layout(guides = "collect") +  
  plot_annotation(
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    tag_levels = 'A',
    tag_prefix = "",
    tag_sep = "",
    theme = theme(
      plot.tag = element_text(size = 18),  
      legend.position = "bottom",  
      legend.box = "horizontal",
      legend.direction = "horizontal",  
      legend.text = element_text(size = 14),
      legend.key.size = unit(1.5, "cm")  
    )
  ) &
  guides(fill = guide_legend(nrow = 1, title = NULL, keywidth = unit(1.5, "cm")))

#print(combined_plot)

ggsave("cell_type_samedonor.png", plot = combined_plot, width = 16, height =5, dpi = 300)



##plot for scc
#load the data for scc
plot_scatterpie_scc <- function(position, calpha, pie_scale = 0.4) {
  # col_manual <- c(
  #   "Inhib" = "#FFDDDD",
  #   "Oligo" = "#B0E0E6",
  #   "OPC" ="#FFFACD" ,
  #   "Excit" ="#32CD32" ,
  #   "MicroOligo" = "#FFEB99",
  #   "Astro" = "#8A2BE2",
  #   "EndoMural" = "#FF1493"
  # )
  col_manual<-c(
    "Keratinocyte" = "#FFDDDD",
    "TSK" = "#32CD32",
    "Melanocyte"="#FF1493"
  )
  calpha <- calpha / rowSums(calpha)
  
  meta <- cbind(position, calpha)
  colnames(meta)[1:2] <- c('x', 'y')
  meta <- as.data.frame(meta)
  
  # 归一化坐标
  meta$x <- (meta$x - min(meta$x)) / (max(meta$x) - min(meta$x))
  meta$y <- (meta$y - min(meta$y)) / (max(meta$y) - min(meta$y))
  
  ggplot() +
    geom_scatterpie(
      data = meta, 
      aes(x = x, y = y),
      cols = colnames(meta)[3:5], 
      color = NA, 
      pie_scale = pie_scale
    ) +
    coord_fixed(ratio = 1) +
    scale_fill_manual(values = col_manual) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      # 移除坐标轴相关元素
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.border = element_blank(),  # 移除绘图区灰色边框
      legend.text = element_text(size = 14),  # 增加图例文字大小
      plot.tag = element_text(size = 18,face="bold")
    ) +
    labs(x = NULL, y = NULL)  # 移除坐标轴标签
}

dir="data/Realdataset/scc"
position3=as.matrix(read.csv(here::here(dir,"matrix3_position_scc.csv"),row.names = 1))
calpha3=as.matrix(read.csv(here::here(dir,"matrix3_celltype_scc.csv"),row.names = 1))
position2=as.matrix(read.csv(here::here(dir,"matrix2_position_scc.csv"),row.names = 1))
calpha2=as.matrix(read.csv(here::here(dir,"matrix2_celltype_scc.csv"),row.names = 1))
position1=as.matrix(read.csv(here::here(dir,"matrix1_position_scc.csv"),row.names = 1))
calpha1=as.matrix(read.csv(here::here(dir,"matrix1_celltype_scc.csv"),row.names = 1))

#produce the plot
p1 <- plot_scatterpie_scc(position1, calpha1, pie_scale = 0.6)
p2 <- plot_scatterpie_scc(position2, calpha2, pie_scale = 0.6)
p3 <- plot_scatterpie_scc(position3, calpha3, pie_scale = 0.6)

#layout and save
combined_plot <- (p1 | p2 | p3) +
  plot_layout(guides = "collect") + 
  plot_annotation(
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    tag_levels = 'A',
    tag_prefix = "",
    tag_sep = "",
    theme = theme(
      plot.tag = element_text(size = 18),  
      legend.position = "bottom",  
      legend.box = "horizontal",
      legend.direction = "horizontal",  
      legend.text = element_text(size = 14),  
      legend.key.size = unit(1.5, "cm")  
    )
  ) &
  guides(fill = guide_legend(nrow = 1, title = NULL, keywidth = unit(1.5, "cm")))

#print(combined_plot)

ggsave("cell_type_scc.png", plot = combined_plot, width = 11, height =4, dpi = 300)






