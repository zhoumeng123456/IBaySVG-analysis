library(ggplot2)

##function for handle the plot
sensitive_plot=function(sensitive_list){
  p_list=list()
  for(paras in c("abpi","abphi","etapsi_sigma","gamma12","cdq","cdp","akdomain")){
    xsize=8
    if(paras=="gamma12"){
      xsize=6
    }
    if(paras !="akdomain"){
      df <- melt(sensitive_list[[paras]])
      if(paras=="abpi"){
        xname=expression(a[pi])
        yname=expression(b[pi])
        param_x <- c("0.6", "0.8", "1", "1.2", "1.4")
        param_y <- c("0.6", "0.8", "1", "1.2", "1.4")
      }else if(paras=="abphi"){
        xname <- expression(a["\u03D5"])  # ϕ
        yname <- expression(b["\u03D5"])
        param_x <- c("0.0005"," 0.0008","0.001", "0.003", "0.005")
        param_y <- c("0.0005"," 0.0008","0.001", "0.003", "0.005")
      }else if(paras=="etapsi_sigma"){
        xname=expression(sigma[eta])
        yname=expression(sigma[psi])
        param_x <- c("0.6", "0.8", "1", "1.2", "1.4")
        param_y <- c("0.6", "0.8", "1", "1.2", "1.4")
      }else if(paras=="gamma12"){
        xname=expression(Gamma[1])
        yname=expression(Gamma[2])
        param_x <- c(expression(sqrt(0.0006)), expression(sqrt(0.0008)), expression(sqrt(0.001)),expression(sqrt(0.0012)), expression(sqrt(0.0014)))
        param_y <- c("0.005", "0.0075", "0.01", "0.0125", "0.015")
      }else if(paras=="cdp"){
        xname=expression(c[p])
        yname=expression(d[p])
        param_x <- c("0.1", "0.15", "0.2", "0.25", "0.3")
        param_y <- c("0.9", "1.35", "1.8", "2.25", "2.7")
      }else if(paras=="cdq"){
        xname=expression(c[q])
        yname=expression(d[q])
        param_x <- c("0.6", "0.8", "1", "1.2", "1.4")
        param_y <- c("0.6", "0.8", "1", "1.2", "1.4")
      }
      
      colnames(df) <- c( "X_param","Y_param", "Value")
      df$X_param <- factor(df$X_param, labels = param_x)
      df$Y_param <- factor(df$Y_param, labels = param_y)
      p1=ggplot(df, aes(x = X_param, y = Y_param)) +
        geom_point(aes(fill = Value), shape = 21, size = 12,color="transparent") + 
        scale_fill_gradient(
          low = "white", 
          high = "#e3716e",
          limits = c(0, 1),         
          name = "Jaccard"          
        ) +
        geom_text(aes(label = round(Value, 3)), color = "black", size = 3) + 
        theme_minimal() +
        labs(x =  xname, y =  yname, fill = "Jaccard")+
        theme(axis.text.x = element_text(size = xsize),legend.position = "none")
        p_list[[paras]]=p1
    }else{
      results_matrix <- sensitive_list[["akdomain"]]
      rownames(results_matrix) <- c("Degree1","Degree2","Degree3","Degree4")
      
      param_values <- list(
        Param1 = c(0.06, 0.07, 0.08, 0.09, 0.10),
        Param2 = c(0.04, 0.045, 0.05, 0.055, 0.06),
        Param3 = c(0.03, 0.035, 0.04, 0.045, 0.05),
        Param4 = c(0.02,0.025,0.03,0.035,0.04)
      )
      

      df <- data.frame(
        Parameter = rep(rownames(results_matrix), each = 5),
        Value = unlist(param_values),       
        Result = as.vector(t(results_matrix))
      )
      

      p1=ggplot(df, aes(x = Value, y = Parameter, fill = Result)) +
        geom_point(aes(fill = Result),shape = 21, size = 11, color = "transparent") +  
        scale_fill_gradient(    low = "white", 
                                high = "#e3716e",
                                limits = c(0, 1),          
                                name = "Jaccard",guide = guide_colorbar(
                                  barwidth  = 12,   
                                  barheight = 0.8,    title.position = "left",  
                                  title.vjust = 1            
                                )  ) +
        geom_text(aes(label = round(Result, 3)), color = "black", size = 3) +
        theme_minimal() +
        labs(x = "Ak",y="", fill = "Metric")+
        theme(legend.position = "bottom")

        p_list[[paras]]=p1
    }
  }
  return(p_list)
}

##choose the dataset
dataset="dlpfc_samedonor" # alternative choose: "dlpfc_acrossdonor", "scc"

##load the sensitive result and plot
sensitive_list=list()
for(paras in c("abpi","abphi","etapsi_sigma","gamma12","cdp","cdq","akdomain") ){
  sensitive_list[[paras]]=as.matrix(read.csv(here::here(paste0("RealData/result_data/sensitive analysis/",dataset),paste0(paras,"_",dataset,"_sensi.csv"))))
}
p_list=sensitive_plot(sensitive_list)#

##layout and save
combined_plot <- 
  (wrap_elements(p_list[[1]] | p_list[[2]] | p_list[[3]]) + 
     labs(tag = "A") + 
     theme(
       plot.tag = element_text(size = 14, face = "bold", margin = margin(0)),
       plot.tag.position = c(0.01, 1)
     )) /
  (wrap_elements(p_list[[4]] | p_list[[5]] | p_list[[6]]) + 
     labs(tag = "B") + 
     theme(
       plot.tag = element_text(size = 14, face = "bold", margin = margin(0)),
       plot.tag.position = c(0.01, 1)
     )) /
  (wrap_elements(p_list[[7]]) + 
     labs(tag = "C") + 
     theme(
       plot.tag = element_text(size = 14, face = "bold", margin = margin(0)),
       plot.tag.position = c(0.01, 1)
     )) +
  plot_layout(
    heights = c(1, 1, 1.2),
    guides = "collect"  
  ) &
  theme(
    legend.position = "bottom", 
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )

plot(combined_plot)  

ggsave(paste0("sensitive_",dataset,".png"), plot = combined_plot, width = 10, height = 10, dpi = 300)




