library(jpeg)
library(png)
library(cowplot)



############################################################################################################################################
                                                #plots of spot cluster identified by IBaySVG 

############################################################################################################################################
plot_total_list=list()#the total list to save
dir="RealData/result_data/spotcluster"

#plot of spot cluster identified by IBaySVG in dlpfc dataset with same donor
plot_list=list()
load(paste0(here(dir,"spot_IBaySVG_dlpfc_samedonor.RData")))
for(datanum in c(1:4)){
  position1=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc samedonor"),paste0("matrix",datanum,"_position_samedonor.csv")),row.names = 1))
  coordinates=as.data.frame(position1)
  coordinates$cluster <- result_total[[datanum]]
  if(datanum==1 |datanum==2){
    custom_colors <- c("#bdb5e1" ,"#f9d580","#54beaa" ,"#fccccb","#eca680","#99b9e9","#b0d992","#e3716e","#7ac7e2","#f3deb7")
  }else{
    custom_colors <- c( "#bdb5e1","#eca680","#54beaa","#fccccb","#b0d992","#99b9e9","#f9d580","#e3716e","#7ac7e2","#f3deb7")
  }
  plot <- ggplot(coordinates, aes(x = x, y = y, color = as.factor(cluster))) +
    geom_point(size = 2) +
    scale_color_manual(values = custom_colors, name = "Cluster") +
    guides(color = guide_legend(override.aes = list(size = 5),nrow = 1)) +  
    theme(
      axis.title = element_blank(),
      legend.position="bottom",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      legend.background = element_blank(),
      legend.key.size = unit(1.5, "cm"),  
      legend.key.height = unit(1, "cm"),  
      legend.key.width = unit(0.8, "cm"),  
      legend.text = element_text(size = 18),  
      legend.title = element_text(size = 20),
      plot.tag = element_text(size = 22, face = "bold")
    )
  plot_list[[datanum]]=plot
}
plot_total_list[["dlpfc_samedonor"]]=plot_list
  

#plot of spot cluster identified by IBaySVG in dlpfc dataset with across donor
load(paste0(here(dir,"spot_IBaySVG_dlpfc_acrossdonor.RData")))
plot_list=list()
for(datanum in c(1:4)){
  position1=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc acrossdonor"),paste0("matrix",datanum,"_position_acrossdonor.csv")),row.names = 1))
  coordinates=as.data.frame(position1)
  coordinates$cluster <- result_total[[datanum]]
  custom_colors <- c("#54beaa" ,"#fccccb","#eca680","#bdb5e1","#b0d992","#99b9e9","#f9d580","#e3716e","#7ac7e2","#f3deb7")
  plot <- ggplot(coordinates, aes(x = x, y = y, color = as.factor(cluster))) +
    geom_point(size = 2) +
    scale_color_manual(values = custom_colors, name = "Estimated") +
    guides(color = guide_legend(override.aes = list(size = 5),nrow = 1)) + 
    theme(
      axis.title = element_blank(),
      legend.position="bottom",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      legend.background = element_blank(),
      legend.key.size = unit(1.5, "cm"),  
      legend.key.height = unit(1, "cm"), 
      legend.key.width = unit(0.8, "cm"),   
      legend.text = element_text(size = 18),  
      legend.title = element_text(size = 20),
      plot.tag = element_text(size = 22, face = "bold")
    )
  plot_list[[datanum]]=plot
  
}
plot_total_list[["dlpfc_acrossdonor"]]=plot_list
  

#plot of spot cluster identified by IBaySVG in scc dataset 
load(here(dir,"spot_IBaySVG_scc.RData"))
plot_list=list()
for(datanum in c(1:3)){
  position1=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/scc"),paste0("matrix",datanum,"_position_scc.csv")),row.names = 1))
  coordinates=as.data.frame(position1)
  coordinates$cluster <- result_total[[datanum]]
  if(datanum==1){
    custom_colors <- c("#99b9e9","#bdb5e1" ,"#54beaa" ,"#f9d580","#eca680","#fccccb","#e3716e","#b0d992","#7ac7e2","#f3deb7")
  }else if(datanum==2){
    custom_colors <- c("#99b9e9","#bdb5e1" ,"#54beaa" ,"#f9d580","#eca680","#e3716e","#b0d992","#fccccb","#7ac7e2","#f3deb7")
  }else{
    custom_colors <- c("#99b9e9","#f9d580","#fccccb","#54beaa" ,"#bdb5e1","#eca680","#e3716e","#b0d992","#7ac7e2","#f3deb7")
  }
  plot <- ggplot(coordinates, aes(x = x, y = y, color = as.factor(cluster))) +
    geom_point(size = 2) +
    scale_color_manual(values = custom_colors, name = "Estimated") +
    guides(color = guide_legend(override.aes = list(size = 5),nrow = 1)) +  
    theme(
      axis.title = element_blank(),
      legend.position="bottom",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      legend.background = element_blank(),
      legend.key.size = unit(1.5, "cm"),  
      legend.key.height = unit(1, "cm"),  
      legend.key.width = unit(0.8, "cm"),   
      legend.text = element_text(size = 18),  
      legend.title = element_text(size = 20),
      plot.tag = element_text(size = 22, face = "bold")
    )
  plot_list[[datanum]]=plot
}
plot_total_list[["scc"]]=plot_list


############################################################################################################################################
                                                        #plots of real spot cluster or HE-stained picture

############################################################################################################################################
plot_total_list_real=list()#the total list to save
dir="RealData/result_data/spotcluster"

#plot of real spot cluster in dlpfc dataset with same donor
load(paste0(here(dir,"spot_real_dlpfc_samedonor.RData")))
plot_list=list()
for(datanum in c(1:4)){
  position1=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc samedonor"),paste0("matrix",datanum,"_position_samedonor.csv")),row.names = 1))
  coordinates=as.data.frame(position1)
  coordinates$cluster <- layer_barcode[[datanum]]
  custom_colors <- c("#54beaa", "#b0d992","#eca680",  "#f9d580", "#bdb5e1", "#fccccb","#e3716e","#99b9e9","#7ac7e2","#f3deb7")
  plot <- ggplot(coordinates, aes(x = x, y = y, color = as.factor(cluster))) +
    geom_point(size = 2) +
    scale_color_manual(values = custom_colors, name = "Cluster") +
    guides(color = guide_legend(override.aes = list(size = 5),nrow = 1)) +  
    theme(
      axis.title = element_blank(),
      legend.position="bottom",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      legend.background = element_blank(),
      legend.key.size = unit(1.5, "cm"),  
      legend.key.height = unit(1, "cm"),  
      legend.key.width = unit(0.8, "cm"),  
      legend.text = element_text(size = 18),  
      legend.title = element_text(size = 20),
      plot.tag = element_text(size = 22, face = "bold")
    )
  plot_list[[datanum]]=plot
  
}
plot_total_list_real[["dlpfc_samedonor"]]=plot_list


#plot of real spot cluster in dlpfc dataset with across donor
load(paste0(here(dir,"spot_real_dlpfc_acrossdonor.RData")))
plot_list=list()
for(datanum in c(1:4)){
  position1=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc acrossdonor"),paste0("matrix",datanum,"_position_acrossdonor.csv")),row.names = 1))
  coordinates=as.data.frame(position1)
  coordinates$cluster <- layer_barcode[[datanum]]
  custom_colors <- c("#54beaa", "#b0d992","#eca680",  "#f9d580", "#bdb5e1", "#fccccb","#e3716e","#99b9e9","#7ac7e2","#f3deb7")
  plot <- ggplot(coordinates, aes(x = x, y = y, color = as.factor(cluster))) +
    geom_point(size = 2) +
    scale_color_manual(values = custom_colors, name = "Real") +
    guides(color = guide_legend(override.aes = list(size = 5),nrow = 1)) +  
    theme(
      axis.title = element_blank(),
      legend.position="bottom",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      legend.background = element_blank(),
      legend.key.size = unit(1.5, "cm"),  
      legend.key.height = unit(1, "cm"),  
      legend.key.width = unit(0.8, "cm"),   
      legend.text = element_text(size = 18),  
      legend.title = element_text(size = 20),
      plot.tag = element_text(size = 22, face = "bold")
    )
  plot_list[[datanum]]=plot
  
}
plot_total_list_real[["dlpfc_acrossdonor"]]=plot_list

###plot of real HE picture in scc dataset
plot_list=list()
for(datanum in c(1:3)){
  img <- readPNG(here(dir,paste0("rep",datanum,"_scc.png")))
  img_plot <- ggdraw() +draw_image(img)
  plot_list[[datanum]]=img_plot
}
plot_total_list_real[["scc"]]=plot_list


  
############################################################################################################################################
                                                                  #layout and save:Figure3

############################################################################################################################################

load(here(dir,"plot_IBaySVG_spotcluster.RData"))
load(here(dir,"plot_real_spotcluster.RData"))

plot1 <- plot_total_list[["dlpfc_samedonor"]][[1]] + theme(legend.position = "none")
plot2 <- plot_total_list[["dlpfc_samedonor"]][[2]] + theme(legend.position = "none")
plot3 <- plot_total_list[["dlpfc_samedonor"]][[3]] + theme(legend.position = "none")
plot4 <- plot_total_list[["dlpfc_samedonor"]][[4]] + theme(legend.position = "none")
plot5 <- plot_total_list[["dlpfc_acrossdonor"]][[1]] + theme(legend.position = "none")
plot6 <- plot_total_list[["dlpfc_acrossdonor"]][[2]] + theme(legend.position = "none")
plot7 <- plot_total_list[["dlpfc_acrossdonor"]][[3]] + theme(legend.position = "none")
plot8 <- plot_total_list[["dlpfc_acrossdonor"]][[4]] + theme(legend.position = "none")
plot9 <- plot_total_list[["scc"]][[1]] + theme(legend.position = "none")
plot10 <- plot_total_list[["scc"]][[2]]+ theme(legend.position = "none")
plot11 <- plot_total_list[["scc"]][[3]] + theme(legend.position = "none")

plot12 <- plot_total_list_real[["dlpfc_samedonor"]][[1]] + theme(legend.position = "none")
plot13 <- plot_total_list_real[["dlpfc_samedonor"]][[2]]+ theme(legend.position = "none")
plot14 <- plot_total_list_real[["dlpfc_samedonor"]][[3]] + theme(legend.position = "none")
plot15 <- plot_total_list_real[["dlpfc_samedonor"]][[4]] + theme(legend.position = "none")
plot16 <- plot_total_list_real[["dlpfc_acrossdonor"]][[1]] + theme(legend.position = "none")
plot17 <- plot_total_list_real[["dlpfc_acrossdonor"]][[2]] + theme(legend.position = "none")
plot18 <- plot_total_list_real[["dlpfc_acrossdonor"]][[3]] + theme(legend.position = "none")
plot19 <- plot_total_list_real[["dlpfc_acrossdonor"]][[4]] + theme(legend.position = "none")
plot20=plot_total_list_real[["scc"]][[1]]
plot21=plot_total_list_real[["scc"]][[2]]
plot22=plot_total_list_real[["scc"]][[3]]

legend_A <- get_legend(plot_total_list_real[["dlpfc_acrossdonor"]][[1]] + 
                         theme(
                           legend.position = "bottom",
                           legend.title = element_text(size = 15), 
                           legend.text  = element_text(size = 15)  
                         ))
legend_B <- get_legend(plot_total_list[["scc"]][[1]] +theme(
  legend.position = "bottom",
  legend.title = element_text(size = 15),  
  legend.text  = element_text(size = 15)  
))
legend_C <- get_legend(plot_total_list[["dlpfc_acrossdonor"]][[2]] +     theme(
  legend.position = "bottom",
  legend.title = element_text(size = 15),  
  legend.text  = element_text(size = 15)   
))


combined_plot <- 
  ((wrap_elements((plot12 | plot13 |plot14 |plot15)) + 
      labs(tag = "A",title="Real") + 
      theme(plot.tag = element_text(size = 18, face = "bold", margin = margin(0)),
            plot.tag.position = c(0.01, 1), plot.title = element_text(hjust = 0.02,size = 14,face = "bold",margin = margin(t = 20)   
            ))) |
  (wrap_elements((plot16 | plot17 |plot18 |plot19)) + 
        labs(tag = "B",title="Real") + 
        theme(plot.tag = element_text(size = 18, face = "bold", margin = margin(0)),
              plot.tag.position = c(0.01, 1), plot.title = element_text(hjust = 0.02,size = 14,face = "bold", margin = margin(t = 20)))))/
  ((wrap_elements((plot1 |plot2|plot3|plot4)) + 
      labs(title="Estimated") + 
      theme(plot.title = element_text(hjust = 0.02,size = 14,face = "bold",margin = margin(t = 20)   # 👈 关键：往下推
      ))) |
  (wrap_elements((plot5 |plot6|plot7|plot8)) + 
        labs(title="Estimated") + 
        theme(plot.title = element_text(hjust = 0.02,size = 14,face = "bold",margin = margin(t = 20)))))/
  legend_A/legend_C/
  ((wrap_elements(plot20 | plot21| plot22) + 
      labs(tag = "C",title="Real") + 
      theme(plot.tag = element_text(size = 18, face = "bold", margin = margin(0)),
            plot.tag.position = c(0.01, 1),plot.title = element_text(hjust = 0.02,size = 14,face = "bold",margin = margin(t = 20))))|
     (wrap_elements(plot9|plot10| plot11) + 
        labs(title="Estimated") + 
        theme(plot.title = element_text(hjust = 0.02,size = 14,face = "bold",margin = margin(t = 20)))))/legend_B+plot_layout(heights = c(1,1,0.12,0.12,1,0.12))

ggsave(paste0("spotcluster_alldataset.png"), plot = combined_plot, width = 16, height = 10, dpi = 300)

