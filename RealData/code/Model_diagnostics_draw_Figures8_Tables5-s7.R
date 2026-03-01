library(nonnest2)
library(future)
library(future.apply)
library(parallel)
library(splines)
library(pscl)
library(MASS)
library(ggplot2)

############################################################################################################################################
                                                        #1.produce the result

############################################################################################################################################

##compute function
#a modified version of CTIG() of IBaySVG to account for linear/focal/period spatial structure.
CTIG_modified<-function(spelist,pattern="linear"){
  #obtain the number of datasets
  m=length(spelist)
  Y.list<-list()
  coord.list<-list()
  samplenumber<-c()
  
  for(i in 1:m){
    Y.list[[i]]<-spelist[[i]][[1]]#spot*gene
    samplenumber[i]<-nrow(Y.list[[i]])## get sample number of each datasets
    coord<-spelist[[i]][[2]]
    coord <-coord-colMeans(coord)#Centercoordinates of spots
    coord <- coord / apply(coord,2,sd)# normalize coordinates of spots
    coord.list[[i]]<-coord
  }
  
  ## gene numbers
  G <- ncol(Y.list[[1]])
  
  ##spline transform to coord.list for 1:4
  co_spline_list_list <- list()
  for(j in c(1:4)){
    co_spline_list <- list()
    for(i in 1:m){
      coord_spline1 <- bs(x=coord.list[[i]][,1],df = j,degree=j)
      coord_spline2 <- bs(x = coord.list[[i]][,2],df = j,degree=j)
      co_spline_list[[i]] <- cbind(coord_spline1, coord_spline2)
    }
    co_spline_list_list[[j]] <- co_spline_list
  }
  
  co_spline_list <- list()
  
  for(i in 1:m){
    if(pattern=="linear"){
      co_spline_list[[i]]<-coord.list[[i]]
    }else if(pattern=="focal"){
      co_spline_list[[i]]<-exp(-coord.list[[i]]^2/2)
    }else if(pattern=="period"){
      co_spline_list[[i]]<-cos(2*pi*coord.list[[i]])
    }
  }
  
  co_spline_list_list[[5]] <- co_spline_list
  
  
  result <- list(Y.list,co_spline_list_list,samplenumber,G)
  return(result)
}

#function for vuong test
capture_vuong <- function(m1, m2) {
  out <- capture.output(pscl::vuong(m1, m2))
  
  target_lines <- out[ grep("Raw|AIC-corrected|BIC-corrected", out) ]
  
  parse_line <- function(line) {
    name <- sub("\\s+([A-Z].*)", "", line)
    name <- trimws(strsplit(name, " ")[[1]][1])
    
    # |Vuong z-statistic| H_A | p-value|
    z <- as.numeric(regmatches(line, regexpr("-?\\d+\\.\\d+", line)))
    
    # H_A: model1 > model2 OR model2 > model1
    ha <- regmatches(line, regexpr("model[12] > model[12]", line))
    
    # p-value
    p_pat <- gregexpr("[0-9]\\.?[0-9eE\\-]*$", line)
    p <- as.numeric(regmatches(line, p_pat)[[1]])
    
    list(z = z, HA = ha, p = p)
  }
  
  parsed <- lapply(target_lines, parse_line)
  names(parsed) <- c("raw", "aic", "bic")
  
  return(parsed)
}

#funtion to fit the best zinb/nb model across samples.
MODEL_CHECKING <- function(spelist, calpha.list=NULL, g = 1, result, model_type = "zinb") {
  ## test length
  n1 <- length(spelist)
  n2 <- length(calpha.list)
  if (n1 != n2) {
    print("warning! the length of expression matrix do not match the calpha's")
  }
  
  if (!model_type %in% c("zinb", "nb")) {
    stop("model_type must be either 'zinb' or 'nb'")
  }
  
  tunning_choose <- rep(1, n1)
  log_likelihood_list <- rep(NA, n1)
  aic_list <- rep(NA, n1)             
  
  for (iii in c(1:n1)) {
    
    y <- unname(spelist[[iii]][[1]][, g])
    zero_ratio <- sum(y == 0) / length(y)
    
    aic_value <- rep(Inf, n1)
    log_lik_values <- rep(NA, n1)     
    aic_values <- rep(NA, n1)        
    
    if (zero_ratio > 0.9) {
      min_position <- 1
      
    } else {
      
      
      for (kkk in 1:4) {
        splinelevel <- kkk
        x <- result[[2]][[kkk]][[iii]]
        
        celltypenum <- ncol(calpha.list[[iii]])
        calpha_names <- paste0("calpha", c(1:celltypenum))
        rows_with_sum_1 <- rowSums(calpha.list[[iii]]) == 1
        calpha.list[[iii]][rows_with_sum_1, 1] <- calpha.list[[iii]][rows_with_sum_1, 1] - 0.01
        calpha.list[[iii]][calpha.list[[iii]] < 0] <- 0
        data <- data.frame(x = x, calpha = calpha.list[[iii]], y = y)
        x_names = paste0("x", c(1:(2 * splinelevel)))
        data_name = c(x_names, calpha_names, "y")
        names(data) = data_name
        formula_str <- paste(data_name[1:(2*splinelevel+celltypenum)], collapse = " + ")
        
        if (model_type == "zinb") {
          formula_count <- paste("y ~", formula_str)
          formula_zero <- paste("|", formula_str)
          formula_full <- as.formula(paste(formula_count, formula_zero))
          model_fit <- zeroinfl(formula_full, data = data, dist = "negbin")
          
        } else {
          formula_full <- as.formula(paste("y ~", formula_str))
          model_fit <- glm.nb(formula_full, data = data)
        }
        
        aic_value[kkk] <- AIC(model_fit)
        log_lik_values[kkk] <- as.numeric(logLik(model_fit)) 
        aic_values[kkk] <- aic_value[kkk]                   
      }
      
      min_position <- which.min(aic_value)
      
      log_likelihood_list[iii] <- log_lik_values[min_position]
      aic_list[iii] <- aic_values[min_position]
    }
    
    tunning_choose[iii] <- min_position
    cat(sprintf("Dataset %d/%d: selected k = %d, log-likelihood = %.3f, AIC = %.3f, model: %s\n", 
                iii, n1, min_position, log_likelihood_list[iii], aic_list[iii], toupper(model_type)))
  }
  
  tunning_final <- max(tunning_choose)
  print(paste0("Choosing b-spline degree: ", tunning_final, ", using ", toupper(model_type), " model"))
  
  return(list(
    best_degree = tunning_final,
    individual_degrees = tunning_choose,
    log_likelihoods = log_likelihood_list,
    aic_values = aic_list
  ))
}

#function to compare the proposed model with NB or specified spatial model(linear/focal/period) using voung test or likelihood test
compare_models <- function(result_nb, result_zinb, spelist, calpha.list, g, result, gene_index = 1) {
  
  # Fit the best degree for the specified gene under two models
  y <- unname(spelist[[gene_index]][[1]][, g])
  best_k_nb <- result_nb$individual_degrees[gene_index]
  best_k_zinb <- result_zinb$individual_degrees[gene_index]
  
  ## ----- Fit NB model -----
  splinelevel <- best_k_nb
  x <- result[[2]][[splinelevel]][[gene_index]]
  
  celltypenum <- ncol(calpha.list[[gene_index]])
  calpha_names <- paste0("calpha", c(1:celltypenum))
  
  rows_with_sum_1 <- rowSums(calpha.list[[gene_index]]) == 1
  calpha.list[[gene_index]][rows_with_sum_1, 1] <- 
    calpha.list[[gene_index]][rows_with_sum_1, 1] - 0.01
  calpha.list[[gene_index]][calpha.list[[gene_index]] < 0] <- 0
  
  data <- data.frame(x = x, calpha = calpha.list[[gene_index]], y = y)
  
  x_names <- paste0("x", c(1:(2 * splinelevel)))
  data_name <- c(x_names, calpha_names, "y")
  names(data) <- data_name
  
  formula_str <- paste(data_name[1:(2*splinelevel + celltypenum)], collapse = " + ")
  formula_full <- as.formula(paste("y ~", formula_str))
  
  model_nb <- glm.nb(formula_full, data = data)
  
  ## ----- Fit ZINB model -----
  splinelevel <- best_k_zinb
  x <- result[[2]][[splinelevel]][[gene_index]]
  
  data <- data.frame(x = x, calpha = calpha.list[[gene_index]], y = y)
  x_names <- paste0("x", c(1:(2 * splinelevel)))
  data_name <- c(x_names, calpha_names, "y")
  names(data) <- data_name
  
  formula_str <- paste(data_name[1:(2*splinelevel + celltypenum)], collapse = " + ")
  formula_zinb <- as.formula(paste("y ~", formula_str, "|", formula_str))
  
  model_zinb <- zeroinfl(formula_zinb, data = data, dist = "negbin")
  
  ## ----- Fit Linear (fixed degree = 1 spline) ZINB model -----
  x <- result[[2]][[5]][[gene_index]]
  
  data <- data.frame(x = x, calpha = calpha.list[[gene_index]], y = y)
  x_names <- paste0("x", c(1:2))
  data_name <- c(x_names, calpha_names, "y")
  names(data) <- data_name
  
  formula_str <- paste(data_name[1:(2 + celltypenum)], collapse = " + ")
  formula_linear <- as.formula(paste("y ~", formula_str, "|", formula_str))
  
  model_linear <- zeroinfl(formula_linear, data = data, dist = "negbin")
  
  ## ----- Basic comparison -----
  cat("=== Model Comparison Results ===\n")
  cat(sprintf("Gene index: %d\n", gene_index))
  cat(sprintf("Sample size: %d\n", length(y)))
  cat(sprintf("Zero proportion: %.3f\n", sum(y == 0) / length(y)))
  cat(sprintf("Best NB degree: %d, Best ZINB degree: %d\n", best_k_nb, best_k_zinb))
  
  cat("\n--- Model Fit Metrics ---\n")
  cat(sprintf("NB log-likelihood: %.3f, AIC: %.3f\n", logLik(model_nb), AIC(model_nb)))
  cat(sprintf("ZINB log-likelihood: %.3f, AIC: %.3f\n", logLik(model_zinb), AIC(model_zinb)))
  cat(sprintf("ΔAIC (ZINB - NB): %.3f\n", AIC(model_zinb) - AIC(model_nb)))
  
  ## ----- Compare ZINB vs NB -----
  cat("\n--- Vuong Test ---\n")
  options(nonnest2.vuong.eps = 1e-6)
  vuong_test <- capture_vuong(model_zinb, model_nb)
  print(vuong_test)
  
  cat("\n--- Likelihood Ratio Test (Manual Calculation) ---\n")
  ll_nb <- as.numeric(logLik(model_nb))
  ll_zinb <- as.numeric(logLik(model_zinb))
  lr_stat <- 2 * (ll_zinb - ll_nb)
  df_diff <- length(coef(model_zinb)) - length(coef(model_nb))
  p_value <- 1 - pchisq(lr_stat, df_diff)
  
  cat(sprintf("LR statistic: %.3f\n", lr_stat))
  cat(sprintf("Degrees of freedom difference: %d\n", df_diff))
  cat(sprintf("p-value: %.4f\n", p_value))
  
  if (p_value < 0.05) {
    cat("Conclusion: ZINB model significantly outperforms NB model (p < 0.05)\n")
  } else {
    cat("Conclusion: No significant difference between ZINB and NB models\n")
  }
  
  delta_aic <- AIC(model_zinb) - AIC(model_nb)
  
  comparison_nb <- list(
    delta_aic = delta_aic,
    lr_statistic = lr_stat,
    lr_pvalue = p_value,
    vuong = vuong_test
  )
  
  ## ----- Compare ZINB vs Linear -----
  cat("\n--- Vuong Test ---\n")
  vuong_test_linear <- capture_vuong(model_zinb, model_linear)
  print(vuong_test_linear)
  
  cat("\n--- Likelihood Ratio Test (Manual Calculation) ---\n")
  ll_linear <- as.numeric(logLik(model_linear))
  lr_stat <- 2 * (ll_zinb - ll_linear)
  df_diff <- length(coef(model_zinb)) - length(coef(model_linear))
  p_value <- 1 - pchisq(lr_stat, df_diff)
  
  cat(sprintf("LR statistic: %.3f\n", lr_stat))
  cat(sprintf("Degrees of freedom difference: %d\n", df_diff))
  cat(sprintf("p-value: %.4f\n", p_value))
  
  if (p_value < 0.05) {
    cat("Conclusion: ZINB model significantly outperforms Linear model (p < 0.05)\n")
  } else {
    cat("Conclusion: No significant difference between ZINB and Linear models\n")
  }
  
  delta_aic <- AIC(model_zinb) - AIC(model_linear)
  
  comparison_linear <- list(
    delta_aic = delta_aic,
    lr_statistic = lr_stat,
    lr_pvalue = p_value,
    vuong = vuong_test_linear
  )
  
  return(
    list(
      best_k_zinb = best_k_zinb,
      best_k_nb = best_k_nb,
      comparison_nb = comparison_nb,
      comparison_linear = comparison_linear
    )
  )
}



##load the data for dlpfc
dir="RealData/result_data/realdata dataset/dlpfc samedonor" #You may switch to the across-donor dataset, as the processing steps are identical except for the number of genes.
subgene=c(1:2400)#We recommend selecting gene sets in batches, as the detection function has limitations on the size of the input dataset.
matrix4=as.matrix(read.csv(here(dir,"matrix4_count_samedonor.csv"),row.names = 1,check.names = FALSE))[subgene,]
matrix3=as.matrix(read.csv(here(dir,"matrix3_count_samedonor.csv"),row.names = 1,check.names = FALSE))[subgene,]
matrix2=as.matrix(read.csv(here(dir,"matrix2_count_samedonor.csv"),row.names = 1,check.names = FALSE))[subgene,]
matrix1=as.matrix(read.csv(here(dir,"matrix1_count_samedonor.csv"),row.names = 1,check.names = FALSE))[subgene,]

position4=as.matrix(read.csv(here(dir,"matrix4_position_samedonor.csv"),row.names = 1))
position3=as.matrix(read.csv(here(dir,"matrix3_position_samedonor.csv"),row.names = 1))
position2=as.matrix(read.csv(here(dir,"matrix2_position_samedonor.csv"),row.names = 1))
position1=as.matrix(read.csv(here(dir,"matrix1_position_samedonor.csv"),row.names = 1))

calpha4=as.matrix(read.csv(here(dir,"matrix4_celltype_samedonor.csv"),row.names = 1))
calpha3=as.matrix(read.csv(here(dir,"matrix3_celltype_samedonor.csv"),row.names = 1))
calpha2=as.matrix(read.csv(here(dir,"matrix2_celltype_samedonor.csv"),row.names = 1))
calpha1=as.matrix(read.csv(here(dir,"matrix1_celltype_samedonor.csv"),row.names = 1))

spelist<-list(list(t(matrix1),position1),list(t(matrix2),position2),
              list(t(matrix3),position3),list(t(matrix4),position4))
c_alpha<-list(calpha1,calpha2,calpha3,calpha4)#the default input for function CTIG_modified()

##or you can load data for scc
dir="RealData/result_data/realdata dataset/scc" 
subgene=c(1:2400)#We recommend selecting gene sets in batches, as the detection function has limitations on the size of the input dataset.

matrix3=as.matrix(read.csv(here(dir,"matrix3_count_scc.csv"),row.names = 1,check.names = FALSE))[subgene,]
matrix2=as.matrix(read.csv(here(dir,"matrix2_count_scc.csv"),row.names = 1,check.names = FALSE))[subgene,]
matrix1=as.matrix(read.csv(here(dir,"matrix1_count_scc.csv"),row.names = 1,check.names = FALSE))[subgene,]

position3=as.matrix(read.csv(here(dir,"matrix3_position_scc.csv"),row.names = 1))
position2=as.matrix(read.csv(here(dir,"matrix2_position_scc.csv"),row.names = 1))
position1=as.matrix(read.csv(here(dir,"matrix1_position_scc.csv"),row.names = 1))

calpha3=as.matrix(read.csv(here(dir,"matrix3_celltype_scc.csv"),row.names = 1))
calpha2=as.matrix(read.csv(here(dir,"matrix2_celltype_scc.csv"),row.names = 1))
calpha1=as.matrix(read.csv(here(dir,"matrix1_celltype_scc.csv"),row.names = 1))

spelist<-list(list(t(matrix1),position1),list(t(matrix2),position2),
              list(t(matrix3),position3))
c_alpha<-list(calpha1,calpha2,calpha3)#the default input for function CTIG_modified()


##produce the result after loading the data
pattern="focal"#alternative choose:linear/period
genedomain=c(1:24)#modify based on the genesize you choose

result_ctig <- CTIG_modified(spelist = spelist,pattern = pattern)
spelist_global <- spelist
c_alpha_global <- c_alpha
result_ctig_global <- result_ctig
plan(multisession, workers = 18)  
for(i in genedomain){
  if(i<=47){
    begin=100*i-99
    end=100*i
    chunk_data=c(begin:end)
  }else{
    begin=100*i-99
    end=100*i+8
    chunk_data=c(begin:end)
  }
  all_results <- future_lapply(chunk_data, function(g) {
    
    result_nb <- MODEL_CHECKING(
      spelist = spelist_global,
      calpha.list = c_alpha_global,
      g = g,
      result = result_ctig_global,
      model_type = "nb"
    )
    
    result_zinb <- MODEL_CHECKING(
      spelist = spelist_global,
      calpha.list = c_alpha_global,
      g = g,
      result = result_ctig_global,
      model_type = "zinb"
    )
    
    one_gene_result <- lapply(1:length(spelist), function(gene_index) {
      compare_models(
        result_nb = result_nb,
        result_zinb = result_zinb,
        spelist = spelist_global,
        calpha.list = c_alpha_global,
        g = g,
        result = result_ctig_global,
        gene_index = gene_index
      )
    })
    
    one_gene_result
  })
  save(all_results,file=paste0((chunk_data[1]),"_",(chunk_data[100]),"_model_checking_",pattern,".RData"))
}




############################################################################################################################################
                                               #2.consider the result: ZINB vs NB:Figures8

############################################################################################################################################
#load the result
dir="RealData/result_data/model check"
load(here(dir,"total_list_model_check_zinb_vs_nb.RData"))#This result contains the comparison between the ZINB and NB models across multiple slices from three datasets.

#handle the result 
datasetname=c("DLPFC-same donor","DLPFC-across donors","SCC")
df_long <- do.call(
  rbind,
  lapply(1:3, function(ds) {
    do.call(
      rbind,
      lapply(1:length(total_list_model_check[[ds]]), function(sl) {
        data.frame(
          dataset = datasetname[ds],
          slice   = sl,
          pval    = -log(total_list_model_check[[ds]][[sl]]$pval_bh,base=10),#choose the p-value after the BH procedure
          row.names = NULL
        )
      })
    )
  })
)
df_long$slice <- factor(df_long$slice)
df_long$pval[is.infinite(df_long$pval)] <- 10# truncated for visualization purposes.
df_long$pval[is.na(df_long$pval)] <- 10 #truncated for visualization purposes.
df_long$dataset <- factor(
  df_long$dataset,
  levels = c(
    "DLPFC-same donor",
    "DLPFC-across donors",
    "SCC"
  )
)

#plot,layout and save
p1=ggplot(df_long, aes(x = slice, y = pval, fill = slice)) +
  geom_boxplot(
    outlier.size = 0.6,
    alpha = 0.8,
    color = "grey30"
  ) +
  geom_hline(yintercept = 1.3, #set the predefined significance threshold
             linetype = "dashed",
             color = "red",
             linewidth = 0.6) +
  facet_wrap(~ dataset, scales = "free_x") +
  theme_bw() +
  scale_fill_manual(
    values = c("#4E79A7", "#59A14F", "#F28E2B", "#E15759")
  ) +
  labs(
    x = "Slice",
    y =  expression(-log[10](adj.~p))
  ) +
  theme(legend.position = "none")

ggsave("model_checking_noinf.png",plot=p1,width = 8, height = 4)

############################################################################################################################################
                                    #3.consider the result: proposed vs linear/focal/period:TableS5-S7

############################################################################################################################################

#function for comparision the result of linear/focal/period
count_model_comparison <- function(result, mode_names = c("linear","focal","period"),n_section = 4,alpha = 0.05){
  
  num_mode_result <- array(
    0,
    dim = c(n_section, 3, length(mode_names)),
    dimnames = list(
      paste0("sec", 1:n_section),
      c("Significant advantage",
        "No significant difference",
        "Opponent significant advantage"),
      mode_names
    )
  )
  
  # Initialize joint result matrix
  num_result <- matrix(
    0,
    nrow = n_section,
    ncol = 3,
    dimnames = list(
      paste0("sec", 1:n_section),
      c("Significant advantage",
        "No significant difference",
        "Opponent significant advantage")
    )
  )
  
  # Total number of genes per mode
  N_total <- length(result[[1]][[mode_names[1]]]) / n_section
  
  # Loop over sections
  for(k in 1:n_section){
    
    for(g in 1:N_total){
      
      idx <- (g-1)*n_section + k
      
      # Extract direction and BH-adjusted p-values
      models <- sapply(mode_names, function(mode) result[[1]][[mode]][idx])
      pvals  <- sapply(mode_names, function(mode) result[[3]][[mode]][idx])
      
      # ---- Single-mode statistics ----
      for(m in seq_along(mode_names)){
        
        if(models[m] == "model1 > model2" && pvals[m] < alpha){
          num_mode_result[k, 1, m] <- num_mode_result[k, 1, m] + 1
          
        } else if(models[m] == "model2 > model1" && pvals[m] < alpha){
          num_mode_result[k, 3, m] <- num_mode_result[k, 3, m] + 1
          
        } else {
          num_mode_result[k, 2, m] <- num_mode_result[k, 2, m] + 1
        }
      }
      
      # ---- Joint statistics across all modes ----
      if(all(models == "model1 > model2" & pvals < alpha)){
        
        num_result[k, 1] <- num_result[k, 1] + 1
        
      } else if(any(models == "model2 > model1" & pvals < alpha)){
        
        num_result[k, 3] <- num_result[k, 3] + 1
        
      } else {
        
        num_result[k, 2] <- num_result[k, 2] + 1
      }
    }
  }
  
  return(list(
    single_mode = num_mode_result,
    joint_mode  = num_result
  ))
}


#load the result
load(here(dir,"voung_spatial_dlpfc_samedonor.RData"))#name:voung_spatial_dlpfc_samedonor
#load(here(dir,"voung_spatial_dlpfc_acrossdonor.RData"))#name:voung_spatial_dlpfc_acrossdonor
#load(here(dir,"voung_spatial_dlpfc_scc.RData"))#name:voung_spatial_dlpfc_scc

res <- count_model_comparison(
  result = voung_spatial_dlpfc_samedonor,#alternative:voung_spatial_dlpfc_acrossdonor,voung_spatial_dlpfc_scc
  mode_names = c("linear","focal","period"),
  n_section = 4,#for dlpfc dataset, n_section=4 and for scc dataset, n_section=3.
  alpha = 0.05
)

res$single_mode   
res$joint_mode    



