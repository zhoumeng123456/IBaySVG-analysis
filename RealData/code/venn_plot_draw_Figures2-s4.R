library(ggvenn)

#function for plot venn
make_venn_plot <- function(data_list, i, namelist,section_num=3, text_size = 5) {
  
  # 当前方法的多个 section 数据
  sections <- data_list[[i]][1:section_num]
  
  # 自动生成 Section 名称
  section_names <- paste0("Section ", seq_along(sections))
  names(sections) <- section_names
  
  # 根据 section 数量自动选择颜色
  base_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#f9d580")
  fill_colors <- base_colors[1:length(sections)]
  
  venn_plot <- ggvenn(
    sections,
    fill_color = fill_colors,
    stroke_size = 0,
    set_name_size = 0,
    show_percentage = FALSE,
    text_size = text_size
  ) +
    coord_fixed(ratio = 1) +
    ggtitle(namelist[i]) +
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5, vjust = -0.5),
      text = element_text(size = 4)
    )
  
  return(venn_plot)
}


#read the data for scc
dir="RealData/result_data/realdata svgene/"
load(here(dir,"data_scc_svgene_list.RData"))

#layout and save
venn_list=list()
namelist=names(data_scc_svgene_list)
for(i in c(1:5,8,9)){
  venn_list[[i]] <- make_venn_plot(data_scc_svgene_list,i,namelist,section_num=3,text_size = 5)
}

combined_plot <- (wrap_elements(venn_list[[5]] | venn_list[[3]] | venn_list[[4]]|venn_list[[1]])
                  /wrap_elements(venn_list[[2]] | venn_list[[9]] | venn_list[[8]])) 
print(combined_plot)

ggsave(paste0("combined_venn_scc.png"), plot = combined_plot, width = 14, height =8, dpi = 300)



#read the data for dlpfc_samedonor
load(here(dir,"data_dlpfc_samedonor_svgene_list.RData"))
venn_list=list()
namelist=names(data_dlpfc_samedonor_svgene_list)
for(i in c(1:5,8,9)){
  venn_list[[i]] <- make_venn_plot(data_dlpfc_samedonor_svgene_list,i,namelist,section_num=4,text_size = 4)
}

#layout and save
combined_plot <- (wrap_elements(venn_list[[5]] | venn_list[[3]] | venn_list[[4]]|venn_list[[1]])
                  /wrap_elements(venn_list[[2]] | venn_list[[9]] | venn_list[[8]])) 
print(combined_plot)

ggsave(paste0("combined_venn_dlpfc_samedonor.png"), plot = combined_plot, width = 12, height =6, dpi = 300)



#read the data for dlpfc acrossdonor
load(here(dir,"data_dlpfc_acrossdonor_svgene_list.RData"))
venn_list=list()
namelist=names(data_dlpfc_acrossdonor_svgene_list)
for(i in c(1:5,8,9)){
  venn_list[[i]] <- make_venn_plot(data_dlpfc_acrossdonor_svgene_list,i,namelist,section_num=4,text_size = 4)
}

#layout and save
combined_plot <- (wrap_elements(venn_list[[5]] | venn_list[[3]] | venn_list[[4]]|venn_list[[1]])
                  /wrap_elements(venn_list[[2]] | venn_list[[9]] | venn_list[[8]])) 
print(combined_plot)

ggsave(paste0("/Users/zhoum/Desktop/combined_venn_scc.png"), plot = combined_plot, width = 12, height =6, dpi = 300)






















