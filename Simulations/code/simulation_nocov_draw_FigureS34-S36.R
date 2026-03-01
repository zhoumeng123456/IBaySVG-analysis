

#load the data,include four index:recall(TPR),FPR,precision and F1 
total_result=c()
for(i in c("IBaySVG","DESpace","SPARK-AVE","SPARKX-AVE","HEARTSVG-AVE","nnSVG-AVE","spVC-AVE","GASTON-AVE",
           "SPARK-Union","SPARKX-Union","HEARTSVG-Union","nnSVG-Union","spVC-Union","GASTON-Union",
           "SPARK-Inter","SPARKX-Inter","HEARTSVG-Inter","nnSVG-Inter", "spVC-Inter","GASTON-Inter",
           "SPARK-Cauchy","SPARKX-Cauchy","HEARTSVG-Cauchy","nnSVG-Cauchy","spVC-Cauchy","SPARK-HMP","SPARKX-HMP",
           "HEARTSVG-HMP","nnSVG-HMP","spVC-HMP","SPARK-PASTE","SPARKX-PASTE","HEARTSVG-PASTE","nnSVG-PASTE","spVC-PASTE","GASTON-PASTE")){
  for(pattern_choose in c("linear","focal","period")){
    for(inf_choose in c("inf3")){
      for(inte_situ_choose in c("scenario 1","scenario 2","scenario 3","scenario 4")){
        if(inf_choose=="inf1"){
          infname="low_dropout"
        }else if(inf_choose=="inf3"){
          infname="middle_dropout"
        }else{
          infname="high_dropout"
        }
        if(inte_situ_choose =="scenario 1"){
          setname="setting1"
        }else if(inte_situ_choose =="scenario 2"){
          setname="setting2"
        }else if(inte_situ_choose =="scenario 3"){
          setname="setting3"
        }else{
          setname="setting4"
        }
        value_result=read.csv(here::here(paste0("Simulations/result_data/simulation result/basic simulation results without covariates/",pattern_choose,"/",setname,"_",infname),paste0(i,"_nocov_simu.csv")))
        TPR_value=value_result$TPR
        FPR_value=value_result$FPR
        F1_value=value_result$F1
        result_type_domain=c("TPR","FPR","F1")
        for(result_type in seq(result_type_domain)){
          temp_result=data.frame(
            method = i,
            set = result_type_domain[result_type],
            value = value_result[,result_type],
            inf = inf_choose,
            inte_situ = inte_situ_choose,
            pattern = pattern_choose
          )
          total_result=rbind(total_result,temp_result)
        }
      }
    }
  }
}



#choose the color
custom_colors <- c(
  "IBaySVG" = "red",             # Red
  "DESpace" = "#FFEB5E", 
  "spVC-AVE" = "black", 
  "spVC-Union" = "#2A2A2A",       
  "spVC-Inter" = "#555555",       
  "spVC-Cauchy" = "#808080",
  "spVC-HMP" = "#AAAAAA",  
  "spVC-PASTE" = "#D3D3D3",   
  "nnSVG-AVE" = "#F27400",            # Orange (Adjusted)
  "nnSVG-Union" = "#F68719",      # Slightly lighter orange
  "nnSVG-Inter" = "#F99B33",            # Orange (Adjusted)
  "nnSVG-Cauchy" = "#FCAE4C",      # Slightly lighter orange
  "nnSVG-HMP" = "#FDBF63",      # Slightly lighter orange
  "nnSVG-PASTE" = "#FFC775",      # Slightly lighter orange
  "SPARK-AVE" = "#0D47A1",            # Light blue
  "SPARK-Union" = "#2C64B3",      # Greyish blue
  "SPARK-Inter" = "#4C82C5",            # Light blue
  "SPARK-Cauchy" = "#6BAFD9",      # Greyish blue
  "SPARK-HMP" = "#8ACCF0",      # Greyish blue
  "SPARK-PASTE" = "#96CDF8",      # Greyish blue
  "SPARKX-AVE" = "#1B5E20",           # Light Green
  "SPARKX-Union" = "#3F8146",     # Slightly lighter green
  "SPARKX-Inter" = "#64A46D",           # Light Green
  "SPARKX-Cauchy" = "#88C793",     # Slightly lighter green
  "SPARKX-HMP" = "#AADDB0",     # Slightly lighter green
  "SPARKX-PASTE" = "#B3DDB6",     # Slightly lighter green
  "HEARTSVG-AVE" = "#4B0082",          # Purple
  "HEARTSVG-Union" = "#6A3296",    # Light purple
  "HEARTSVG-Inter" = "#8964AB",          # Purple
  "HEARTSVG-Cauchy" = "#A897BF",    # Light purple
  "HEARTSVG-HMP" = "#C0AFE1",    # Light purple
  "HEARTSVG-PASTE" = "#C4B6ED",    # Light purple
  "GASTON-AVE" = "#B22222",
  "GASTON-Union" = "#CC4C4C",
  "GASTON-Inter" = "#E06666",
  "GASTON-PASTE" = "#FF9999"
)


# function for plot
plot_uni_inter_single2 <- function(data, cluster_type,inf_type,result_type,legend_position="none",title="linear",x_title="",y_title="F1 Value") {
  ggplot(data %>% filter(inte_situ == cluster_type)%>%filter(inf== inf_type)%>%filter(set== result_type)%>%
           mutate(method = factor(method, levels = c("IBaySVG","DESpace",
                                                     "SPARK-AVE","SPARKX-AVE","HEARTSVG-AVE","nnSVG-AVE","spVC-AVE","GASTON-AVE",
                                                     "SPARK-Union","SPARKX-Union","HEARTSVG-Union","nnSVG-Union","spVC-Union","GASTON-Union",
                                                     "SPARK-Inter","SPARKX-Inter","HEARTSVG-Inter","nnSVG-Inter", "spVC-Inter","GASTON-Inter",
                                                     "SPARK-Cauchy","SPARKX-Cauchy",
                                                     "HEARTSVG-Cauchy"
                                                     ,"nnSVG-Cauchy"
                                                     ,"spVC-Cauchy","SPARK-HMP","SPARKX-HMP",
                                                     "HEARTSVG-HMP"
                                                     ,"nnSVG-HMP"
                                                     ,"spVC-HMP","SPARK-PASTE","SPARKX-PASTE",
                                                     "HEARTSVG-PASTE"
                                                     ,"nnSVG-PASTE"
                                                     ,"spVC-PASTE","GASTON-PASTE")),
                  method_group = case_when(
                    method %in% c("IBaySVG", "DESpace") ~ "Multi",
                    grepl("AVE", method) ~ "AVE",
                    grepl("Union", method) ~ "Union",
                    grepl("Inter", method) ~ "Inter",
                    grepl("Cauchy", method) ~ "Cauchy",
                    grepl("HMP", method) ~ "HMP",
                    grepl("PASTE", method) ~ "PASTE",
                    TRUE ~ "Other"
                  ),
                  method_group = factor(method_group, levels = c("Multi", "AVE", "Union", "Inter", "Cauchy","HMP","PASTE")),
                  pattern = factor(str_to_title(pattern), levels = c("Linear", "Focal", "Period"))),
         aes(y = value,x=method_group, fill = method, color = method)) +
    facet_grid(cols = vars(pattern),space = "free_x", scales = "free_x") +
    geom_boxplot(position = position_dodge2(preserve = "single", padding = 0.4), 
                 width = 0.9, outlier.size = 0.5, outlier.shape = NA) +
    labs(title =title,
         x = x_title,
         y = y_title) +
    theme_minimal() +
    theme(
      plot.margin = margin(0.2,0.2,0.2,0.2),
      plot.title = element_text(size = 12,face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.ticks.x = element_line(),  
      axis.text.x = element_text(size = 9),
      legend.position = legend_position,
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      strip.text = element_text(size = 12),
      panel.spacing.x = unit(0.5, "lines"),
      plot.tag = element_text(size = 15, face = "bold")
    ) +
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(values = custom_colors)
}


#three index you need to choose,respectively
result_type="F1"
y_title="F1 value"

result_type="TPR"
y_title="TPR"

result_type="FPR"
y_title="FPR"

plot_total1=plot_uni_inter_single2(total_result, "scenario 1","inf3",title="Setting 1",legend_position="right",result_type=result_type,y_title=y_title)
plot_total2=plot_uni_inter_single2(total_result, "scenario 2","inf3",title="Setting 2",result_type=result_type,y_title=y_title)
plot_total3=plot_uni_inter_single2(total_result, "scenario 3","inf3",title="Setting 3",result_type=result_type,y_title=y_title)
plot_total4=plot_uni_inter_single2(total_result, "scenario 4","inf3",title="Setting 4",result_type=result_type,y_title=y_title)

combined_plot <- ((plot_total1)/(plot_total2)/(plot_total3)/(plot_total4)) +
  plot_layout(guides = "collect", heights = c(1.2, 1.2, 1.2, 1.2)) + 
  plot_annotation(
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    tag_levels = 'A',  
    tag_prefix = "",  
    tag_sep = "",  
    theme = theme(plot.tag = element_text(size = 12, face = "bold"))  
  ) &
  theme(legend.position = "right") &  
  guides(fill = guide_legend(ncol = 1))

#print(combined_plot)

ggsave(paste0("simu_",result_type,"_no_cov_inf3.png"), plot = combined_plot, width = 12, height = 12, dpi = 300)






























