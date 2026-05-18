
##load the data
sensitive_list <- list()
pattern_domain=c("linear","focal","period")
for(paras in c("api_bpi","aphi_bphi","eta_sigma_psi_sigma","gamma1_gamma2","cp_dp","cq_dq","akdomain") ){
  sensitive_list[[paras]] <- list()  
  for(pattern_idx in seq_along(pattern_domain)){
    pattern_name <- pattern_domain[pattern_idx]
    seed_list <- list()
    for(seed1 in 1:10){
      file_path <- here::here(
        "Simulations/result_data/sensitivity analysis of hyperparameters",
        pattern_name,
        paras,
        paste0(paras, "_sensi_", pattern_name, "_seed", seed1, ".csv")
      )
      seed_list[[seed1]] <- read.csv(file_path, check.names = FALSE)[,-1]
    }
    seed_array <- array(unlist(seed_list), 
                        dim = c(nrow(seed_list[[1]]), ncol(seed_list[[1]]), 10))
    mean_matrix <- apply(seed_array, c(1,2), mean)
    sensitive_list[[paras]][[pattern_name]] <- mean_matrix
  }
}




##make the picture
parasdomain=c("api_bpi","aphi_bphi","eta_sigma_psi_sigma","gamma1_gamma2","cp_dp","cq_dq","akdomain")
pattern_domain=c("linear","focal","period")
p_total_list=list()
for(paras in parasdomain){
  xsize=8
  if(paras=="gamma1_gamma2"){
    xsize=5
  }
  p_list=replicate(3,list(), simplify = FALSE)
  if(paras !="akdomain"){
    for(pattern in seq(pattern_domain)){
      sensitive_matrix=as.matrix(sensitive_list[[paras]][[pattern]])
      df <- reshape2::melt(sensitive_matrix)
      
      if(paras=="api_bpi"){
        xname=expression(a[pi])
        yname=expression(b[pi])
        param_x <- c("0.6", "0.8", "1", "1.2", "1.4")
        param_y <- c("0.6", "0.8", "1", "1.2", "1.4")
      }else if(paras=="aphi_bphi"){
        xname <- expression(a["\u03D5"])  # ϕ
        yname <- expression(b["\u03D5"])
        param_x <- c("0.0005"," 0.0008","0.001", "0.003", "0.005")
        param_y <- c("0.0005"," 0.0008","0.001", "0.003", "0.005")
      }else if(paras=="eta_sigma_psi_sigma"){
        xname=expression(sigma[eta])
        yname=expression(sigma[psi])
        param_x <- c("0.6", "0.8", "1", "1.2", "1.4")
        param_y <- c("0.6", "0.8", "1", "1.2", "1.4")
      }else if(paras=="gamma1_gamma2"){
        xname=expression(Gamma[1])
        yname=expression(Gamma[2])
        param_x <- c(expression(sqrt(0.0006)), expression(sqrt(0.0008)), expression(sqrt(0.001)),expression(sqrt(0.0012)), expression(sqrt(0.0014)))
        param_y <- c("0.005", "0.0075", "0.01", "0.0125", "0.015")
      }else if(paras=="cq_dq"){
        xname=expression(c[q])
        yname=expression(d[q])
        param_x <- c("0.6", "0.8", "1", "1.2", "1.4")
        param_y <- c("0.6", "0.8", "1", "1.2", "1.4")
      }else if(paras=="cp_dp"){
        xname=expression(c[p])
        yname=expression(d[p])
        param_x <- c("0.1", "0.15", "0.2", "0.25", "0.3")
        param_y <- c("0.9", "1.35", "1.8", "2.25", "2.7")
      }
      
      colnames(df) <- c( "X_param","Y_param", "Value")
      
      df$X_param <- factor(df$X_param, labels = param_x)
      df$Y_param <- factor(df$Y_param, labels = param_y)
      
      p1=ggplot(df, aes(x = X_param, y = Y_param)) +
        geom_point(aes(fill = Value), shape = 21, size = 12,color="transparent") + 
        scale_fill_gradient(
          low = "white", 
          high = "#e3716e",
          limits = c(0.5, 1),        
          name = "F1 Value"         
        ) +
        geom_text(aes(label = round(Value, 3)), color = "black", size = 3) +
        theme_minimal() +
        labs(x =  xname, y =  yname, fill = "F1 Value")+
        theme(axis.text.x = element_text(size = xsize),legend.position = "none")
      
      p_list[[pattern]]=p1
    }
  }else{
    for(pattern in seq(pattern_domain)){
      results_matrix <- as.matrix(sensitive_list[[paras]][[pattern]])
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
                                name = "F1 Value",guide = guide_colorbar(
                                  barwidth  = 12,   
                                  barheight = 0.8,    title.position = "left",  
                                  title.vjust = 1            
                                )  ) +
        geom_text(aes(label = round(Result, 3)), color = "black", size = 3) +
        theme_minimal() +
        labs(x = "Ak", y ="", fill = "F1 Value")+theme(legend.position = "bottom")
      
      p_list[[pattern]]=p1
    }
  }
  p_total_list[[paras]]=p_list
}


##Layout and save 
for(i in seq(pattern_domain)){
  combined_plot <- 
    (wrap_elements(p_total_list[[1]][[i]] | p_total_list[[2]][[i]] | p_total_list[[3]][[i]]) + 
       labs(tag = "A") + 
       theme(
         plot.tag = element_text(size = 14, face = "bold", margin = margin(0)),
         plot.tag.position = c(0.01, 1)
       )) /
    (wrap_elements(p_total_list[[4]][[i]] | p_total_list[[6]][[i]] | p_total_list[[5]][[i]]) + 
       labs(tag = "B") + 
       theme(
         plot.tag = element_text(size = 14, face = "bold", margin = margin(0)),
         plot.tag.position = c(0.01, 1)
       )) /
    (wrap_elements(p_total_list[[7]][[i]]) + 
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
  ggsave(paste0("simulation_sensi_",pattern_domain[i],".png"), plot = combined_plot, width = 10, height = 10, dpi = 300)
  
}
























