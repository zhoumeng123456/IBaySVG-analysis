#######################################################################################################################################
                                                    #comparative analysis in real dataset

#######################################################################################################################################
#Here we first take the results of dlpfc dataset with samedonor for an example.

dataset="dlpfc_samedonor"
load(here::here("RealData/result_data/realdata svgene",paste0("data_",dataset,"_svgene_list.RData")))

#1.genes are easily to identify
#compare to all seven single-sample methods
all_inter_sv=Reduce(intersect,list(data_dlpfc_samedonor_svgene_list[[1]][[6]],data_dlpfc_samedonor_svgene_list[[3]][[6]],
                                         data_dlpfc_samedonor_svgene_list[[2]][[6]],data_dlpfc_samedonor_svgene_list[[4]][[6]],
                                         data_dlpfc_samedonor_svgene_list[[5]][[6]],data_dlpfc_samedonor_svgene_list[[8]][[6]],
                                         data_dlpfc_samedonor_svgene_list[[9]][[6]],"IBaySVG"=data_dlpfc_samedonor_svgene_list[[6]]))#2
#compare to five conventional single-sample methods
all_inter_sv=Reduce(intersect,list(data_dlpfc_samedonor_svgene_list[[1]][[6]],data_dlpfc_samedonor_svgene_list[[3]][[6]],
                                   data_dlpfc_samedonor_svgene_list[[2]][[6]],data_dlpfc_samedonor_svgene_list[[4]][[6]],
                                   data_dlpfc_samedonor_svgene_list[[5]][[6]],"IBaySVG"=data_dlpfc_samedonor_svgene_list[[6]]))#26

#compare to GASTON and StarTrail
length(intersect(all_inter_sv,data_dlpfc_samedonor_svgene_list[[9]][[6]]))#11
length(intersect(all_inter_sv,data_dlpfc_samedonor_svgene_list[[8]][[6]]))#15


#2.identified genes limited by Intersect strategy
#genes detected by at least one single-sample method under intersection integration
all_union_inter_sv=Reduce(union,list(data_dlpfc_samedonor_svgene_list[[1]][[6]],data_dlpfc_samedonor_svgene_list[[3]][[6]],
                                   data_dlpfc_samedonor_svgene_list[[2]][[6]],data_dlpfc_samedonor_svgene_list[[4]][[6]],
                                   data_dlpfc_samedonor_svgene_list[[5]][[6]],data_dlpfc_samedonor_svgene_list[[8]][[6]],
                                   data_dlpfc_samedonor_svgene_list[[9]][[6]]))#1544

length(intersect(all_union_inter_sv,data_dlpfc_samedonor_svgene_list[[6]]))#561

#genes that are overlooked by all single-sample methods with intersection-based integration
interproblem=setdiff(data_dlpfc_samedonor_svgene_list[[6]],all_union_inter_sv)#531

#genes are detected by at least one single-sample approach Using the Cauchy method
print(length(intersect(Reduce(union,list(data_dlpfc_samedonor_svgene_list[[1]][[7]],data_dlpfc_samedonor_svgene_list[[2]][[7]]
                                         ,data_dlpfc_samedonor_svgene_list[[3]][[7]],data_dlpfc_samedonor_svgene_list[[4]][[7]]
                                         ,data_dlpfc_samedonor_svgene_list[[5]][[7]])),interproblem)))#237

#genes are detected by at least one single-sample approach Using the HMP method
print(length(intersect(Reduce(union,list(data_dlpfc_samedonor_svgene_list[[1]][[8]],data_dlpfc_samedonor_svgene_list[[2]][[8]]
                                         ,data_dlpfc_samedonor_svgene_list[[3]][[8]],data_dlpfc_samedonor_svgene_list[[4]][[8]]
                                         ,data_dlpfc_samedonor_svgene_list[[5]][[8]])),interproblem)))#420


#3.identified genes limited by Union strategy
#genes identified by three or more of the single-sample methods using the union-based integration
all_union_union_list=list(data_dlpfc_samedonor_svgene_list[[1]][[5]],data_dlpfc_samedonor_svgene_list[[3]][[5]],
                                     data_dlpfc_samedonor_svgene_list[[2]][[5]],data_dlpfc_samedonor_svgene_list[[4]][[5]],
                                     data_dlpfc_samedonor_svgene_list[[5]][[5]],data_dlpfc_samedonor_svgene_list[[8]][[5]],
                                     data_dlpfc_samedonor_svgene_list[[9]][[5]])

gene_counts <- table(unlist(all_union_union_list))
genes_in_three_or_more <- names(gene_counts[gene_counts >= 3])#1795

#genes undetected by IBaySVG
unionproblem=setdiff(genes_in_three_or_more,data_dlpfc_samedonor_svgene_list[[6]])#1010

cauchy_list=table(unlist(list(data_dlpfc_samedonor_svgene_list[[1]][[7]],data_dlpfc_samedonor_svgene_list[[2]][[7]]
                           ,data_dlpfc_samedonor_svgene_list[[3]][[7]],data_dlpfc_samedonor_svgene_list[[4]][[7]]
                           ,data_dlpfc_samedonor_svgene_list[[5]][[7]])))

genes_cauchy_in_three_or_more <- names(cauchy_list[cauchy_list >= 3])
length(intersect(genes_cauchy_in_three_or_more,unionproblem))#274

hmp_list=table(unlist(list(data_dlpfc_samedonor_svgene_list[[1]][[8]],data_dlpfc_samedonor_svgene_list[[2]][[8]]
                              ,data_dlpfc_samedonor_svgene_list[[3]][[8]],data_dlpfc_samedonor_svgene_list[[4]][[8]]
                              ,data_dlpfc_samedonor_svgene_list[[5]][[8]])))
genes_hmp_in_three_or_more <- names(hmp_list[hmp_list >= 3])
length(intersect(genes_hmp_in_three_or_more,unionproblem))#427


#4.weak genes can not be identified by other methods
all_union_union_list_nogaston=Reduce(union,list(data_dlpfc_samedonor_svgene_list[[1]][[5]],data_dlpfc_samedonor_svgene_list[[3]][[5]],
                          data_dlpfc_samedonor_svgene_list[[2]][[5]],data_dlpfc_samedonor_svgene_list[[4]][[5]],
                          data_dlpfc_samedonor_svgene_list[[5]][[5]],data_dlpfc_samedonor_svgene_list[[8]][[5]]))

weakproblem=setdiff(data_dlpfc_samedonor_svgene_list[[6]],all_union_union_list_nogaston)#76

#5.PASTE method's incompatibility
#86 genes identified by nnSVG-PASTE 
length(intersect(data_dlpfc_samedonor_svgene_list[[5]][[9]],all_union_inter_sv))#86
#286 genes identified by StarTrail-PASTE
length(intersect(data_dlpfc_samedonor_svgene_list[[8]][[7]],all_union_inter_sv))#263

#Finally, we choose the gene based on the results of "all_inter_sv", "interproblem", "unionproblem" and "weakproblem".

#######################################################################################################################################
                                #comparative analysis for dlpfc dataset with samedonor: Figure2

#######################################################################################################################################
library(ggplot2)
library(patchwork)

#function for ploting
meta_process =function(position1,count1){
  meta=cbind(position1,t(count1))
  colnames(meta)[c(1:2)] = c('x','y')
  meta=as.data.frame(meta)
  meta$x <- (meta$x - min(meta$x))/(max(meta$x) - min(meta$x))
  meta$y <- (meta$y - min(meta$y))/(max(meta$y) - min(meta$y))
  return(meta)
}

pattern_plot2 <- function(pltdat, igene, xy = TRUE, main = FALSE, titlesize = 2, 
                          pointsize = 3, xpand = 0, ypand = 1, title = NULL, legend = TRUE) {
  if (!xy) {
    xy <- matrix(as.numeric(do.call(rbind, strsplit(as.character(pltdat[, 1]), split = "x")), ncol = 2))
    rownames(xy) <- as.character(pltdat[, 1])
    colnames(xy) <- c("x", "y")
    pd <- cbind.data.frame(xy, pltdat[, 2:ncol(pltdat)])
  } else {
    pd <- pltdat
  }
  
  
  target_col <- pd[, igene + 2]
  min_val <- min(target_col, na.rm = TRUE)
  max_val <- max(target_col, na.rm = TRUE)
  pd$rel_value <- (target_col - min_val) / (max_val - min_val) * 100  
  
  
  pal <- colorRampPalette(c("mediumseagreen", "lightyellow2", "deeppink"))
  gpt <- ggplot(pd, aes(x = x, y = y, color = rel_value)) + 
    geom_point(size = pointsize) + 
    scale_color_gradientn(
      colours = pal(5),
      limits = c(0, 100), 
      labels = scales::percent_format(scale = 1)  
    ) +
    scale_x_discrete(expand = c(xpand, ypand)) + 
    scale_y_discrete(expand = c(xpand, ypand)) + 
    coord_equal() + 
    theme_bw() + 
    labs(color = "Relative Expression Level ")  
  
  if (legend) {
    gpt <- gpt + 
      theme(legend.position = "bottom",
            legend.justification = "center",
            legend.text = element_text(
              size = 8,
              margin = margin(t = 0, b = 0, unit = "cm"),  
              hjust = 0,                                           
              vjust = 0.5                                        
            )) +  
      guides(color = guide_colorbar(
        direction = "horizontal",
        barwidth = unit(5, "cm"),
        barheight = unit(0.3, "cm"),
        title.position = "left",
        label.position = "bottom"
      ))
  } else {
    gpt <- gpt + theme(legend.position = "none")
  }
  
  if (main) {
    out <- gpt + labs(title = title %||% colnames(pd)[igene + 2], x = NULL, y = NULL) + 
      theme(plot.title = element_text(hjust = 0.5, size = rel(titlesize)))
  } else {
    out <- gpt + labs(title = NULL, x = NULL, y = NULL)
  }
  
  return(out)
}
#load the data
dataset="samedonor"
position1=as.matrix(read.csv(here::here(paste0("data/Realdataset/dlpfc ",dataset),paste0("matrix1_position_",dataset,".csv")),row.names = 1))
count1=as.matrix(read.csv(here::here(paste0("data/Realdataset/dlpfc ",dataset),paste0("matrix1_count_",dataset,".csv")),row.names = 1))
position2=as.matrix(read.csv(here::here(paste0("data/Realdataset/dlpfc ",dataset),paste0("matrix2_position_",dataset,".csv")),row.names = 1))
count2=as.matrix(read.csv(here::here(paste0("data/Realdataset/dlpfc ",dataset),paste0("matrix2_count_",dataset,".csv")),row.names = 1))
position3=as.matrix(read.csv(here::here(paste0("data/Realdataset/dlpfc ",dataset),paste0("matrix3_position_",dataset,".csv")),row.names = 1))
count3=as.matrix(read.csv(here::here(paste0("data/Realdataset/dlpfc ",dataset),paste0("matrix3_count_",dataset,".csv")),row.names = 1))
position4=as.matrix(read.csv(here::here(paste0("data/Realdataset/dlpfc ",dataset),paste0("matrix4_position_",dataset,".csv")),row.names = 1))
count4=as.matrix(read.csv(here::here(paste0("data/Realdataset/dlpfc ",dataset),paste0("matrix4_count_",dataset,".csv")),row.names = 1))

#choose the problem gene
total_gene_choose=c("SCGB1D2","GFAP","KRT19","AKR1B1", "FBXO9","ROGDI","ACER3", "B3GAT2", "NR2F2","FAM81A","NAA20","PSMA4")
count1=count1[total_gene_choose,]
count2=count2[total_gene_choose,]
count3=count3[total_gene_choose,]
count4=count4[total_gene_choose,]
meta1=meta_process(position1,count1)
meta2=meta_process(position2,count2)
meta3=meta_process(position3,count3)
meta4=meta_process(position4,count4)

#plot and layout
combined_list=list()
for(genenum in c(1:12)){
  p1 <- pattern_plot2(meta1, genenum, xpand = 0.01, ypand = 0.01, main = FALSE, pointsize = 1,legend = FALSE)
  p2 <- pattern_plot2(meta2, genenum, xpand = 0.01, ypand = 0.01, main = FALSE, pointsize = 1,legend = FALSE)
  p3 <- pattern_plot2(meta3, genenum, xpand = 0.01, ypand = 0.01, main = FALSE, pointsize = 1,legend = FALSE)
  p4 <- pattern_plot2(meta4, genenum, xpand = 0.01, ypand = 0.01, main = FALSE, pointsize = 1,legend = FALSE)
  
  combined_list[[genenum]] <- wrap_elements((p1 | p2) / (p3 | p4)) + 
    labs(title = total_gene_choose[genenum]) + 
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0),
    )
}

legend_space <- plot_spacer()
combined_plot <- 
  (wrap_elements(combined_list[[1]]/combined_list[[2]]/combined_list[[3]]) + 
     labs(tag = "A") + 
     theme(plot.tag = element_text(size = 25, face = "bold", margin = margin(0)),
           plot.tag.position = c(0.01, 0.99))) |
  (wrap_elements(combined_list[[4]]/combined_list[[5]]/combined_list[[6]]) + 
     labs(tag = "B") + 
     theme(plot.tag = element_text(size = 25, face = "bold", margin = margin(0)),
           plot.tag.position = c(0.01, 0.99))) | 
  (wrap_elements(combined_list[[7]]/combined_list[[8]]/combined_list[[9]]) + 
     labs(tag = "C") + 
     theme(plot.tag = element_text(size = 25, face = "bold", margin = margin(0)),
           plot.tag.position = c(0.01,0.99))) |
  (wrap_elements(combined_list[[10]]/combined_list[[11]]/combined_list[[12]]) + 
     labs(tag = "D") + 
     theme(plot.tag = element_text(size = 25, face = "bold", margin = margin(0)),
           plot.tag.position = c(0.01, 0.99))) |
  plot_layout(heights = c(1.2, 1.2, 1.2,1.2))  

#print(combined_plot)
ggsave(paste0("combined_dlpfc_samedonor_verti.png"), plot = combined_plot, width = 24, height = 16, dpi = 300)


#######################################################################################################################################
                                     #comparative analysis for dlpfc dataset across donors: Figures12

#######################################################################################################################################

#load the data
dataset="acrossdonor"
position1=as.matrix(read.csv(here::here(paste0("data/Realdataset/dlpfc ",dataset),paste0("matrix1_position_",dataset,".csv")),row.names = 1))
count1=as.matrix(read.csv(here::here(paste0("data/Realdataset/dlpfc ",dataset),paste0("matrix1_count_",dataset,".csv")),row.names = 1))
position2=as.matrix(read.csv(here::here(paste0("data/Realdataset/dlpfc ",dataset),paste0("matrix2_position_",dataset,".csv")),row.names = 1))
count2=as.matrix(read.csv(here::here(paste0("data/Realdataset/dlpfc ",dataset),paste0("matrix2_count_",dataset,".csv")),row.names = 1))
position3=as.matrix(read.csv(here::here(paste0("data/Realdataset/dlpfc ",dataset),paste0("matrix3_position_",dataset,".csv")),row.names = 1))
count3=as.matrix(read.csv(here::here(paste0("data/Realdataset/dlpfc ",dataset),paste0("matrix3_count_",dataset,".csv")),row.names = 1))
position4=as.matrix(read.csv(here::here(paste0("data/Realdataset/dlpfc ",dataset),paste0("matrix4_position_",dataset,".csv")),row.names = 1))
count4=as.matrix(read.csv(here::here(paste0("data/Realdataset/dlpfc ",dataset),paste0("matrix4_count_",dataset,".csv")),row.names = 1))

#choose the problem gene
total_gene_choose=c("COX6C","MBP","SPARCL1","JAKMIP1", "GPRASP2","NFE2L1","VAV3", "TAGLN", "DNAJC2","AKT3","NOVA1","EIF2S2")

count1=count1[total_gene_choose,]
count2=count2[total_gene_choose,]
count3=count3[total_gene_choose,]
count4=count4[total_gene_choose,]
meta1=meta_process(position1,count1)
meta2=meta_process(position2,count2)
meta3=meta_process(position3,count3)
meta4=meta_process(position4,count4)

#plot and layout
for(genenum in c(1:12)){
  p1 <- pattern_plot2(meta1, genenum, xpand = 0.01, ypand = 0.01, main = FALSE, pointsize = 1,legend = FALSE)
  p2 <- pattern_plot2(meta2, genenum, xpand = 0.01, ypand = 0.01, main = FALSE, pointsize = 1,legend = FALSE)
  p3 <- pattern_plot2(meta3, genenum, xpand = 0.01, ypand = 0.01, main = FALSE, pointsize = 1,legend = FALSE)
  p4 <- pattern_plot2(meta4, genenum, xpand = 0.01, ypand = 0.01, main = FALSE, pointsize = 1,legend = FALSE)
  
  combined_list[[genenum]] <- wrap_elements((p1 | p2) / (p3 | p4)) + 
    labs(title = total_gene_choose[genenum]) + 
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0),
    )
}

combined_plot <- 
  (wrap_elements(combined_list[[1]]|combined_list[[2]]|combined_list[[3]]) + 
     labs(tag = "A") + 
     theme(plot.tag = element_text(size = 25, face = "bold", margin = margin(0)),
           plot.tag.position = c(0.01, 1))) /
  (wrap_elements(combined_list[[4]]|combined_list[[5]]|combined_list[[6]]) + 
     labs(tag = "B") + 
     theme(plot.tag = element_text(size = 25, face = "bold", margin = margin(0)),
           plot.tag.position = c(0.01, 1))) /
  (wrap_elements(combined_list[[7]]|combined_list[[8]]|combined_list[[9]]) + 
     labs(tag = "C") + 
     theme(plot.tag = element_text(size = 25, face = "bold", margin = margin(0)),
           plot.tag.position = c(0.01, 1))) /
  (wrap_elements(combined_list[[10]]|combined_list[[11]]|combined_list[[12]]) + 
     labs(tag = "D") + 
     theme(plot.tag = element_text(size = 25, face = "bold", margin = margin(0)),
           plot.tag.position = c(0.01, 1))) +
  plot_layout(heights = c(1.2, 1.2,1.2,1.2,0.2)) 

#print(combined_plot)

ggsave(paste0("combined_dlpfc_acrossdonor.png"), plot = combined_plot, width = 15, height = 24, dpi = 300)



#######################################################################################################################################
                                            #comparative analysis for scc dataset: Figures13

#######################################################################################################################################

#load the data
dataset="scc"
position1=as.matrix(read.csv(here::here(paste0("data/Realdataset/",dataset),paste0("matrix1_position_",dataset,".csv")),row.names = 1))
count1=as.matrix(read.csv(here::here(paste0("data/Realdataset/",dataset),paste0("matrix1_count_",dataset,".csv")),row.names = 1))
position2=as.matrix(read.csv(here::here(paste0("data/Realdataset/",dataset),paste0("matrix2_position_",dataset,".csv")),row.names = 1))
count2=as.matrix(read.csv(here::here(paste0("data/Realdataset/",dataset),paste0("matrix2_count_",dataset,".csv")),row.names = 1))
position3=as.matrix(read.csv(here::here(paste0("data/Realdataset/",dataset),paste0("matrix3_position_",dataset,".csv")),row.names = 1))
count3=as.matrix(read.csv(here::here(paste0("data/Realdataset/",dataset),paste0("matrix3_count_",dataset,".csv")),row.names = 1))

#choose the problem gene
total_gene_choose=c("LGALS7","SPRR1B","S100A2","AP2A2","CWC15","DNM2","PSCA","CSGALNACT1","LCN2","EGFR", "HIF1A","KLF5")
count1=count1[total_gene_choose,]
count2=count2[total_gene_choose,]
count3=count3[total_gene_choose,]
meta1=meta_process(position1,count1)
meta2=meta_process(position2,count2)
meta3=meta_process(position3,count3)

#plot and layout
combined_list=list()
for(genenum in c(1:12)){
  p1 <- pattern_plot2(meta1, genenum, xpand = 0.01, ypand = 0.01, main = FALSE, pointsize = 2,legend = FALSE)
  p2 <- pattern_plot2(meta2, genenum, xpand = 0.01, ypand = 0.01, main = FALSE, pointsize = 2,legend = FALSE)
  p3 <- pattern_plot2(meta3, genenum, xpand = 0.01, ypand = 0.01, main = FALSE, pointsize = 2,legend = FALSE)
  
  combined_list[[genenum]] <- wrap_elements(p1 | p2 | p3) + 
    labs(title = total_gene_choose[genenum]) + 
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, vjust = 0),
    )
}
combined_plot <- 
  ((wrap_elements(combined_list[[1]]/combined_list[[2]]/combined_list[[3]]) + 
      labs(tag = "A") + 
      theme(plot.tag = element_text(size = 25, face = "bold", margin = margin(0)),
            plot.tag.position = c(0.01, 0.99))) |
     (wrap_elements(combined_list[[4]]/combined_list[[5]]/combined_list[[6]]) + 
        labs(tag = "B") + 
        theme(plot.tag = element_text(size = 25, face = "bold", margin = margin(0)),
              plot.tag.position = c(0.01, 0.99)))) /
  ((wrap_elements(combined_list[[7]]/combined_list[[8]]/combined_list[[9]]) + 
      labs(tag = "C") + 
      theme(plot.tag = element_text(size = 25, face = "bold", margin = margin(0)),
            plot.tag.position = c(0.01, 1))) |
     (wrap_elements(combined_list[[10]]/combined_list[[11]]/combined_list[[12]]) + 
        labs(tag = "D") + 
        theme(plot.tag = element_text(size = 25, face = "bold", margin = margin(0)),
              plot.tag.position = c(0.01, 1)))) +
  plot_layout(heights = c(1.2, 1.2,0.1)) 

ggsave(paste0("combined_scc.png"), plot = combined_plot, width = 14, height = 18, dpi = 900)












