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
dir="data/Realdataset/dlpfc samedonor" #You may switch to the across-donor dataset, as the processing steps are identical except for the number of genes.
subgene=c(1:2400)#We recommend selecting gene sets in batches, as the detection function has limitations on the size of the input dataset.
matrix4=as.matrix(read.csv(here::here(dir,"matrix4_count_samedonor.csv"),row.names = 1,check.names = FALSE))[subgene,]
matrix3=as.matrix(read.csv(here::here(dir,"matrix3_count_samedonor.csv"),row.names = 1,check.names = FALSE))[subgene,]
matrix2=as.matrix(read.csv(here::here(dir,"matrix2_count_samedonor.csv"),row.names = 1,check.names = FALSE))[subgene,]
matrix1=as.matrix(read.csv(here::here(dir,"matrix1_count_samedonor.csv"),row.names = 1,check.names = FALSE))[subgene,]

position4=as.matrix(read.csv(here::here(dir,"matrix4_position_samedonor.csv"),row.names = 1))
position3=as.matrix(read.csv(here::here(dir,"matrix3_position_samedonor.csv"),row.names = 1))
position2=as.matrix(read.csv(here::here(dir,"matrix2_position_samedonor.csv"),row.names = 1))
position1=as.matrix(read.csv(here::here(dir,"matrix1_position_samedonor.csv"),row.names = 1))

calpha4=as.matrix(read.csv(here::here(dir,"matrix4_celltype_samedonor.csv"),row.names = 1))
calpha3=as.matrix(read.csv(here::here(dir,"matrix3_celltype_samedonor.csv"),row.names = 1))
calpha2=as.matrix(read.csv(here::here(dir,"matrix2_celltype_samedonor.csv"),row.names = 1))
calpha1=as.matrix(read.csv(here::here(dir,"matrix1_celltype_samedonor.csv"),row.names = 1))

spelist<-list(list(t(matrix1),position1),list(t(matrix2),position2),
              list(t(matrix3),position3),list(t(matrix4),position4))
c_alpha<-list(calpha1,calpha2,calpha3,calpha4)#the default input for function CTIG_modified()

##or you can load data for scc
dir="data/Realdataset/scc" 
subgene=c(1:2400)#We recommend selecting gene sets in batches, as the detection function has limitations on the size of the input dataset.

matrix3=as.matrix(read.csv(here::here(dir,"matrix3_count_scc.csv"),row.names = 1,check.names = FALSE))[subgene,]
matrix2=as.matrix(read.csv(here::here(dir,"matrix2_count_scc.csv"),row.names = 1,check.names = FALSE))[subgene,]
matrix1=as.matrix(read.csv(here::here(dir,"matrix1_count_scc.csv"),row.names = 1,check.names = FALSE))[subgene,]

position3=as.matrix(read.csv(here::here(dir,"matrix3_position_scc.csv"),row.names = 1))
position2=as.matrix(read.csv(here::here(dir,"matrix2_position_scc.csv"),row.names = 1))
position1=as.matrix(read.csv(here::here(dir,"matrix1_position_scc.csv"),row.names = 1))

calpha3=as.matrix(read.csv(here::here(dir,"matrix3_celltype_scc.csv"),row.names = 1))
calpha2=as.matrix(read.csv(here::here(dir,"matrix2_celltype_scc.csv"),row.names = 1))
calpha1=as.matrix(read.csv(here::here(dir,"matrix1_celltype_scc.csv"),row.names = 1))

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
                                               #2.consider the result: ZINB vs NB:Figures34

############################################################################################################################################
#load the result
dir="RealData/result_data/model check"
load(here::here(dir,"total_list_model_check_zinb_vs_nb.RData"))#This result contains the comparison between the ZINB and NB models across multiple slices from three datasets.

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
                                    #3.consider the result: proposed vs linear/focal/period:Table S9-S11

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
load(here::here(dir,"voung_spatial_dlpfc_samedonor.RData"))#name:voung_spatial_dlpfc_samedonor
#load(here::here(dir,"voung_spatial_dlpfc_acrossdonor.RData"))#name:voung_spatial_dlpfc_acrossdonor
#load(here::here(dir,"voung_spatial_dlpfc_scc.RData"))#name:voung_spatial_dlpfc_scc

res <- count_model_comparison(
  result = voung_spatial_dlpfc_samedonor,#alternative:voung_spatial_dlpfc_acrossdonor,voung_spatial_dlpfc_scc
  mode_names = c("linear","focal","period"),
  n_section = 4,#for dlpfc dataset, n_section=4 and for scc dataset, n_section=3.
  alpha = 0.05
)

res$single_mode   
res$joint_mode    


#################################################################################################################
#bspline-check for Figures_S35-s39
#################################################################################################################
#loading packages
library(FNN)

######################################## Figure S35 : the certification of discontinuity###############################################################
######################################## Figure S35 : the certification of discontinuity###############################################################
######################################## Figure S35 : the certification of discontinuity###############################################################
#Note: Because the full computational procedure and resulting objects require substantial memory, 
#we recommend reproducing the analysis directly from Section 3 using the preprocessed results provided
#in this repository.

#function
#the procedure to obtain the file "paste0("_zinb_fit_", dataname, "\\.RData$")" in function "interface_meanbias_analysis"
#is illustrated in outlier_detection_Tables_S7-S8_Figure_S33.R
interface_meanbias_analysis <- function(dataname, inter_gene_num, spelist, layer_barcode,
                                        fitdir = "",k_nn = 6, min_diff_prop = 1/6,min_scale = 0.02) {
  M <- length(spelist)
  iface_list <- lapply(seq_len(M), function(m)
    make_interface_by_knn(spelist[[m]][[2]], layer_barcode[[m]], k_nn, min_diff_prop))
  nn_list <- lapply(seq_len(M), function(m)
    FNN::get.knn(as.matrix(spelist[[m]][[2]]), k = k_nn)$nn.index)   
  
  grp_idx_list <- lapply(iface_list, function(ifc_m) {
    bnd <- ifc_m$boundary; ifc <- ifc_m$interface; lay <- ifc_m$layer
    grp <- ifelse(is.na(bnd), NA_character_, ifelse(bnd, ifc, lay))
    keep <- !is.na(grp); split(which(keep), grp[keep])
  })
  
  files <- list.files(fitdir, pattern = paste0("_zinb_fit_", dataname, "\\.RData$"),
                      full.names = TRUE)
  stopifnot(length(files) > 0)
  
  per_gene <- list()
  for (f in files) {
    load(f)                                   # zinb_summary, lambda_list, pi_list
    gene_ids <- as.integer(colnames(lambda_list[[1]]))
    for (m in seq_len(M)) {
      Lm <- lambda_list[[m]]; Pm <- pi_list[[m]]; Ym <- spelist[[m]][[1]]
      grp_idx <- grp_idx_list[[m]]; n_grp <- lengths(grp_idx); nn <- nn_list[[m]]
      gnames <- names(grp_idx)
      for (j in seq_along(gene_ids)) {
        lam <- Lm[, j]; if (all(is.na(lam))) next
        g_abs <- gene_ids[j]; rel <- g_abs - inter_gene_num
        pp <- Pm[, j]; y <- Ym[, rel]
        fit <- (1 - pp) * lam                              # ZINB mean count
        scale_g <- mean(fit, na.rm = TRUE)                 # overall level 
        if (!is.finite(scale_g) || scale_g < min_scale) next
        obs_loc <- rowMeans(cbind(y, matrix(y[nn], nrow = length(y))))  
        adl <- abs(obs_loc - fit)                          # spot |real−fit| by gene after denoisy
        abs_diff=abs(y - fit)
        #residual
        mobs <- sapply(grp_idx, function(ix) mean(y[ix],   na.rm = TRUE))   
        mfit <- sapply(grp_idx, function(ix) mean(fit[ix], na.rm = TRUE))  
        madl <- sapply(grp_idx, function(ix) mean(adl[ix], na.rm = TRUE))   
        mabs_diff <- sapply(grp_idx, function(ix) mean(abs_diff[ix], na.rm = TRUE))   
        
        per_gene[[length(per_gene)+1]] <- data.frame(
          gene = g_abs, sample = m, group = gnames,
          type = ifelse(grepl("\\|", gnames), "interface", "layer"),
          mobs = mobs, mfit = mfit,
          bias_scaled = (mobs - mfit) / scale_g,           
          mad_scaled  = madl / scale_g,                    
          mabs_diff_scaled  = mabs_diff / scale_g,
          n = n_grp[gnames], row.names = NULL, stringsAsFactors = FALSE)
      }
    }
  }
  gene_tab <- do.call(rbind, per_gene)
  save(gene_tab, iface_list,
       file = file.path(fitdir, paste0("interface_meanbias_", dataname, ".RData")))
  list(gene_tab = gene_tab, iface_list = iface_list)
}

make_interface_by_knn <- function(position, layer, k = 6, min_diff_prop = 1/6) {
  pos <- as.matrix(position); layer <- as.character(layer); n <- nrow(pos)
  nn  <- FNN::get.knn(pos, k = k)$nn.index
  diff_prop <- numeric(n); interface <- rep(NA_character_, n)
  for (i in seq_len(n)) {
    own <- layer[i]; nb <- layer[nn[i, ]]
    if (is.na(own)) { diff_prop[i] <- NA; next }
    valid <- !is.na(nb)
    diff_prop[i] <- if (any(valid)) mean(nb[valid] != own) else NA
    if (!is.na(diff_prop[i]) && diff_prop[i] >= min_diff_prop) {
      diffn <- nb[valid][nb[valid] != own]
      if (length(diffn)) {
        adj <- names(sort(table(diffn), decreasing = TRUE))[1]   
        interface[i] <- paste(sort(c(own, adj)), collapse = "|") # L2|L3 == L3|L2
      }
    }
  }
  data.frame(spot = seq_len(n), layer = layer,
             boundary = diff_prop >= min_diff_prop,
             interface = interface, boundary_score = diff_prop,
             stringsAsFactors = FALSE)
}

jump_per_gene_3seg <- function(gt2, interface) {
  AB <- strsplit(interface, "\\|")[[1]]; A <- AB[1]; B <- AB[2]
  key <- paste(gt2$gene, gt2$sample)
  getv <- function(grp){ v<-gt2$mobs[gt2$group==grp]; names(v)<-key[gt2$group==grp]; v }
  mA <- getv(A); mAB <- getv(interface); mB <- getv(B)
  common <- Reduce(intersect, list(names(mA), names(mAB), names(mB)))
  if (length(common)==0) return(NULL)
  jA <- abs(mA[common]-mAB[common])/(0.5*(mA[common]+mAB[common])+1e-6)
  jB <- abs(mB[common]-mAB[common])/(0.5*(mB[common]+mAB[common])+1e-6)
  rbind(
    data.frame(interface=interface, side="left side",  rel_jump=as.numeric(jA)),
    data.frame(interface=interface, side="right side", rel_jump=as.numeric(jB))
  )
}

#1.load the data
#DLPFC-same donor
dir="data/Realdataset/dlpfc samedonor" #You may switch to the across-donor dataset
matrix4=as.matrix(read.csv(here::here(dir,"matrix4_count_samedonor.csv"),row.names = 1,check.names = FALSE))
matrix3=as.matrix(read.csv(here::here(dir,"matrix3_count_samedonor.csv"),row.names = 1,check.names = FALSE))
matrix2=as.matrix(read.csv(here::here(dir,"matrix2_count_samedonor.csv"),row.names = 1,check.names = FALSE))
matrix1=as.matrix(read.csv(here::here(dir,"matrix1_count_samedonor.csv"),row.names = 1,check.names = FALSE))
position4=as.matrix(read.csv(here::here(dir,"matrix4_position_samedonor.csv"),row.names = 1))
position3=as.matrix(read.csv(here::here(dir,"matrix3_position_samedonor.csv"),row.names = 1))
position2=as.matrix(read.csv(here::here(dir,"matrix2_position_samedonor.csv"),row.names = 1))
position1=as.matrix(read.csv(here::here(dir,"matrix1_position_samedonor.csv"),row.names = 1))
calpha4=as.matrix(read.csv(here::here(dir,"matrix4_celltype_samedonor.csv"),row.names = 1))
calpha3=as.matrix(read.csv(here::here(dir,"matrix3_celltype_samedonor.csv"),row.names = 1))
calpha2=as.matrix(read.csv(here::here(dir,"matrix2_celltype_samedonor.csv"),row.names = 1))
calpha1=as.matrix(read.csv(here::here(dir,"matrix1_celltype_samedonor.csv"),row.names = 1))
spelist_samedonor <- list(list(t(matrix1),position1),list(t(matrix2),position2),
                          list(t(matrix3),position3),list(t(matrix4),position4))
calpha_samedonor <- list(calpha1,calpha2,calpha3,calpha4)

#or DLPFC across donors
dir="data/Realdataset/dlpfc acrossdonor" 
matrix4=as.matrix(read.csv(here::here(dir,"matrix4_count_acrossdonor.csv"),row.names = 1,check.names = FALSE))
matrix3=as.matrix(read.csv(here::here(dir,"matrix3_count_acrossdonor.csv"),row.names = 1,check.names = FALSE))
matrix2=as.matrix(read.csv(here::here(dir,"matrix2_count_acrossdonor.csv"),row.names = 1,check.names = FALSE))
matrix1=as.matrix(read.csv(here::here(dir,"matrix1_count_acrossdonor.csv"),row.names = 1,check.names = FALSE))
position4=as.matrix(read.csv(here::here(dir,"matrix4_position_acrossdonor.csv"),row.names = 1))
position3=as.matrix(read.csv(here::here(dir,"matrix3_position_acrossdonor.csv"),row.names = 1))
position2=as.matrix(read.csv(here::here(dir,"matrix2_position_acrossdonor.csv"),row.names = 1))
position1=as.matrix(read.csv(here::here(dir,"matrix1_position_acrossdonor.csv"),row.names = 1))
calpha4=as.matrix(read.csv(here::here(dir,"matrix4_celltype_acrossdonor.csv"),row.names = 1))
calpha3=as.matrix(read.csv(here::here(dir,"matrix3_celltype_acrossdonor.csv"),row.names = 1))
calpha2=as.matrix(read.csv(here::here(dir,"matrix2_celltype_acrossdonor.csv"),row.names = 1))
calpha1=as.matrix(read.csv(here::here(dir,"matrix1_celltype_acrossdonor.csv"),row.names = 1))
spelist_acrossdonor <- list(list(t(matrix1),position1),list(t(matrix2),position2),
                            list(t(matrix3),position3),list(t(matrix4),position4))
calpha_acrossdonor <- list(calpha1,calpha2,calpha3,calpha4)

# or SCC
dir=paste0("data/Realdataset/scc")
matrix3=as.matrix(read.csv(here::here(dir,"matrix3_count_scc.csv"),row.names = 1,check.names = FALSE))
matrix2=as.matrix(read.csv(here::here(dir,"matrix2_count_scc.csv"),row.names = 1,check.names = FALSE))
matrix1=as.matrix(read.csv(here::here(dir,"matrix1_count_scc.csv"),row.names = 1,check.names = FALSE))
position3=as.matrix(read.csv(here::here(dir,"matrix3_position_scc.csv"),row.names = 1))
position2=as.matrix(read.csv(here::here(dir,"matrix2_position_scc.csv"),row.names = 1))
position1=as.matrix(read.csv(here::here(dir,"matrix1_position_scc.csv"),row.names = 1))
calpha3=as.matrix(read.csv(here::here(dir,"matrix3_celltype_scc.csv"),row.names = 1))
calpha2=as.matrix(read.csv(here::here(dir,"matrix2_celltype_scc.csv"),row.names = 1))
calpha1=as.matrix(read.csv(here::here(dir,"matrix1_celltype_scc.csv"),row.names = 1))
spelist_scc <- list(list(t(matrix1),position1),list(t(matrix2),position2),
                    list(t(matrix3),position3))
calpha_scc <- list(calpha1,calpha2,calpha3)

#2.compute the result
#DLPFC-same donor
load(here::here("RealData/result_data/spotcluster","spot_real_dlpfc_samedonor.RData"))
res2 <- interface_meanbias_analysis("dlpfcsamedonor", inter_gene_num = 0,spelist = spelist_samedonor,layer_barcode = layer_barcode)


#DLPFC-across donors
load(here::here("RealData/result_data/spotcluster","spot_real_dlpfc_acrossdonor.RData"))
res2 <- interface_meanbias_analysis("dlpfcacrossdonor", inter_gene_num = 0,spelist = spelist_acrossdonor,layer_barcode = layer_barcode)

#scc
load(here::here("RealData/result_data/spotcluster","spot_IBaySVG_scc.RData"))
res2 <- interface_meanbias_analysis("scc", inter_gene_num = 0,spelist = spelist_scc,layer_barcode = result_total, min_diff_prop = 1/6)

#3.choose the dataset for plot based on the result in section 1 and 2
##########################################DLPFC#######################################
dataname="dlpfcsamedonor" #or you can replace it with dlpfcacrossdonor
load(here::here("RealData/result_data/model check",paste0("interface_meanbias_", dataname, ".RData")))
res2=list("gene_tab"=gene_tab,"iface_list"=iface_list)
gt2 <- res2$gene_tab 
ifaces <- grep("\\|", unique(gt2$group), value=TRUE)
ifaces <- ifaces[!grepl("others", ifaces)]
jump_long <- do.call(rbind, lapply(ifaces, function(q) jump_per_gene_3seg(gt2, q)))
#dlpfc
ord_if <- c("L1|L2","L2|L3","L3|L4","L4|L5","L5|L6","L6|WM")
jump_long <- jump_long[jump_long$interface %in% ord_if, ]

##########################################SCC#########################################
dataname="scc" 
load(here::here("RealData/result_data/model check",paste0("interface_meanbias_", dataname, ".RData")))
res2=list("gene_tab"=gene_tab,"iface_list"=iface_list)
gt2 <- res2$gene_tab 
ifaces <- grep("\\|", unique(gt2$group), value=TRUE)
ifaces <- ifaces[!grepl("others", ifaces)]
jump_long <- do.call(rbind, lapply(ifaces, function(q) jump_per_gene_3seg(gt2, q)))
ord_if <- c("1|5","1|7","2|7","2|8","3|8","4|5","4|7","4|8","5|7","7|8")
jump_long <- jump_long[jump_long$interface %in% ord_if, ]

#4.plot by the ggplot function. you are advised to generate p1/p2/p3 plot derived from "dlpfcsamedonor","dlpfcacrossdonor" and "scc" datasets in section 3.
#Remember to replace the title name in "labs".
jump_long$interface <- factor(jump_long$interface, levels=ord_if)
jump_long$side <- factor(jump_long$side, levels=c("left side","right side"))
p1=ggplot(jump_long, aes(x=interface, y=rel_jump, fill=side)) +
  stat_boxplot(geom = "errorbar", width = 0.4,        
               position = position_dodge(0.8)) +
  geom_boxplot(outlier.size=0.4, outlier.alpha=0.25, position=position_dodge(0.8))+
  scale_fill_manual(values=c("left side"="#F4A582","right side"="#92C5DE")) +
  labs(x=NULL, y="relative expression jump", fill=NULL,
       title="SCC") +
  coord_cartesian(ylim=c(0, 2)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 0),
    plot.title = element_text(
      hjust = 0.5,      
      size = 12
    )
  )

#5.combined the plot and save
combined_plot <- (p1 /p2 / p3) +
  plot_layout(
    heights = c(1.5, 1.5, 1.5),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(size = 15, face = "bold"),
    plot.tag.position = c(0.01, 0.99)
  )

#plot(combined_plot)
ggsave("C:/Users/zhoum/Desktop/real_expression_jump.png",
       plot = combined_plot, width = 8, height = 10, dpi = 300)


#################################Figure S36: the overall fit performence based on the RQR #######################
#################################Figure S36: the overall fit performence based on the RQR #######################
#################################Figure S36: the overall fit performence based on the RQR #######################
#Note: Because the full computational procedure and resulting objects require substantial memory, 
#we recommend reproducing the analysis directly from Section 2 using the preprocessed results provided
#in this repository. The procedure to obtain the file "paste0("_zinb_fit_", dataname, "\\.RData$")" 
#in function "interface_rqr_analysis" is illustrated in outlier_detection_Tables_S7-S8_Figure_S33.R

#function
rqr_rand <- function(y, lambda, pi_i, theta, B) {
  n  <- length(y)
  Fy <- pi_i + (1 - pi_i) * pnbinom(y,     mu = lambda, size = theta)
  Fl <- ifelse(y >= 1, pi_i + (1 - pi_i) * pnbinom(y - 1, mu = lambda, size = theta), 0)
  lo <- pmin(Fl, Fy); hi <- pmax(Fl, Fy)
  u  <- matrix(runif(B * n, rep(lo, each = B), rep(hi, each = B)), nrow = B)
  qnorm(pmin(pmax(u, 1e-10), 1 - 1e-10))
}

interface_rqr_analysis <- function(dataname, inter_gene_num, spelist, layer_barcode,fitdir = "",
                                   k_nn = 6, min_diff_prop = 1/6, B = 20) {
  M <- length(spelist)
  iface_list <- lapply(seq_len(M), function(m)
    make_interface_by_knn(spelist[[m]][[2]], layer_barcode[[m]], k_nn, min_diff_prop))
  
  grp_idx_list <- lapply(iface_list, function(ifc_m) {
    bnd <- ifc_m$boundary; ifc <- ifc_m$interface; lay <- ifc_m$layer
    grp <- ifelse(is.na(bnd), NA_character_, ifelse(bnd, ifc, lay))
    keep <- !is.na(grp)
    split(which(keep), grp[keep])
  })
  
  files <- list.files(fitdir, pattern = paste0("_zinb_fit_", dataname, "\\.RData$"),
                      full.names = TRUE)
  stopifnot(length(files) > 0)
  
  per_gene <- list()
  for (f in files) {
    load(f)
    gene_ids <- as.integer(colnames(lambda_list[[1]]))
    for (m in seq_len(M)) {
      Lm <- lambda_list[[m]]; Pm <- pi_list[[m]]; Ym <- spelist[[m]][[1]]
      grp_idx <- grp_idx_list[[m]]; n_grp <- lengths(grp_idx)
      stopifnot(nrow(Lm) == nrow(Ym))
      rel_all <- gene_ids - inter_gene_num
      stopifnot(min(rel_all) >= 1, max(rel_all) <= ncol(Ym))
      
      for (j in seq_along(gene_ids)) {
        lam <- Lm[, j]; if (all(is.na(lam))) next
        g_abs <- gene_ids[j]; rel <- g_abs - inter_gene_num
        pp <- Pm[, j]; y <- Ym[, rel]
        th <- zinb_summary$theta[zinb_summary$gene == g_abs & zinb_summary$sample == m]
        if (length(th) != 1 || is.na(th)) next
        
        aR <- abs(rqr_rand(y, lam, pp, th, B))
        ar <- colMeans(aR); a2 <- colMeans(aR > 2); a3 <- colMeans(aR > 3); a4 <- colMeans(aR > 4)
        
        gnames <- names(grp_idx)
        per_gene[[length(per_gene)+1]] <- data.frame(
          gene = g_abs, sample = m, group = gnames,
          type = ifelse(grepl("\\|", gnames), "interface", "layer"),
          absrqr = sapply(grp_idx, function(ix) mean(ar[ix])),
          gt2 = sapply(grp_idx, function(ix) mean(a2[ix])),
          gt3 = sapply(grp_idx, function(ix) mean(a3[ix])),
          gt4 = sapply(grp_idx, function(ix) mean(a4[ix])),
          n = n_grp[gnames], row.names = NULL, stringsAsFactors = FALSE)
      }
    }
  }
  gene_tab <- do.call(rbind, per_gene)
  save(gene_tab, iface_list,
       file = file.path(fitdir, paste0("interface_rqr_", dataname, ".RData")))
  list(gene_tab = gene_tab, iface_list = iface_list)
}

summarise_interface <- function(gt, Delta = 0.04, min_gene = 20) {
  parts <- split(gt, gt$interface)
  out <- lapply(names(parts), function(q) {
    d <- na.omit(parts[[q]]$delta)
    if (length(d) < min_gene) return(NULL)
    p_up <- t.test(d, mu =  Delta, alternative = "less")$p.value
    p_lo <- t.test(d, mu = -Delta, alternative = "greater")$p.value
    data.frame(interface = q, n_gene = length(d),
               n_spot_med = median(parts[[q]]$n_I, na.rm = TRUE),
               delta_med  = median(d), frac_pos = mean(d > 0),
               absrqr_I = mean(parts[[q]]$absrqr_I, na.rm = TRUE),
               absrqr_N = mean(parts[[q]]$absrqr_N, na.rm = TRUE),
               gt4_I = mean(parts[[q]]$gt4_I, na.rm = TRUE),
               cohens_d = mean(d) / sd(d),
               tost_equiv = max(p_up, p_lo) < 0.05)      
  })
  res <- do.call(rbind, out)
  res[order(-abs(res$delta_med)), ]
}

gradient_summary <- function(gt, Delta = 0.04, min_gene = 20) {
  gt$key <- paste(gt$gene, gt$sample, sep = "_")
  val <- setNames(gt$absrqr, paste(gt$key, gt$group, sep = "@@"))   
  getv <- function(keys, g) val[paste(keys, g, sep = "@@")]
  ifaces <- grep("\\|", unique(gt$group), value = TRUE)
  out <- lapply(ifaces, function(q) {
    AB <- strsplit(q, "\\|")[[1]]; A <- AB[1]; B <- AB[2]
    keys <- gt$key[gt$group == q]
    vAB <- getv(keys, q); vA <- getv(keys, A); vB <- getv(keys, B)
    ok <- is.finite(vAB) & is.finite(vA) & is.finite(vB)
    if (sum(ok) < min_gene) return(NULL)
    bump <- vAB[ok] - 0.5*(vA[ok] + vB[ok])                
    p_up <- t.test(bump, mu= Delta, alternative="less")$p.value
    p_lo <- t.test(bump, mu=-Delta, alternative="greater")$p.value
    data.frame(interface=q, n=sum(ok),
               flankA=mean(vA[ok]), boundary=mean(vAB[ok]), flankB=mean(vB[ok]),  
               bump_025=unname(quantile(bump,0.25, na.rm = TRUE)),bump_med=median(bump),bump_075=unname(quantile(bump,0.75, na.rm = TRUE)),frac_pos=mean(bump>0),
               cohens_d=mean(bump)/sd(bump), tost_equiv=max(p_up,p_lo)<0.05)
  })
  res <- do.call(rbind, out)
  res[order(-res$bump_med), ]
}

wmean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

#1. compute the RQR result
#The results have been stored in “RealData/result_data/model check”
#DLPFC-same donor
load(here::here("RealData/result_data/spotcluster","spot_real_dlpfc_samedonor.RData"))
res <- interface_rqr_analysis("dlpfcsamedonor", inter_gene_num = 0,
                              spelist = spelist_samedonor,
                              layer_barcode = layer_barcode)
#save(res,file=here::here(“RealData/result_data/model check”,"res_mean_boundary_acrossdonor.RData"))

#DLPFC-across donor
load(here::here("RealData/result_data/spotcluster","spot_real_dlpfc_acrossdonor.RData"))
res <- interface_rqr_analysis("dlpfcacrossdonor", inter_gene_num = 0,
                              spelist = spelist_acrossdonor,
                              layer_barcode = layer_barcode)
#save(res,file=here::here(“RealData/result_data/model check”,"res_mean_boundary_acrossdonor.RData"))

#SCC
load(here::here("RealData/result_data/spotcluster","spot_IBaySVG_scc.RData"))
res <- interface_rqr_analysis("scc", inter_gene_num = 0,
                              spelist = spelist_scc,
                              layer_barcode = result_total)
#save(res,file=here::here(“RealData/result_data/model check”,"res_mean_boundary_scc.RData"))


#2.summary the result 
dataname="samedonor" #or you can replace it with "acrossdonor" and "scc".
load(here::here("RealData/result_data/model check",paste0("res_mean_boundary_",dataname,".RData")))
gradient_summary(res$gene_tab,Delta = 0.05)       
aggregate(absrqr ~ group, res$gene_tab, mean)   
# gt <- res$gene_tab #total gene
gt <- res$gene_tab
gname <- colnames(spelist_samedonor[[1]][[1]])
gt$gene_name <- gname[gt$gene]
sv_names <- data_dlpfc_samedonor_svgene_list[[6]]#the SV gene identified
gt <- gt[gt$gene_name %in% sv_names, ]    
gt2 <- do.call(rbind,lapply(split(gt, interaction(gt$gene, gt$sample)), function(d){
  data.frame(
    gene      = d$gene[1],sample    = d$sample[1],gene_name = d$gene_name[1],
    absrqr_B = weighted.mean(d$absrqr[d$type == "interface"],w = d$n[d$type == "interface"],na.rm = TRUE),
    absrqr_N = weighted.mean(d$absrqr[d$type == "layer"], w = d$n[d$type == "layer"],na.rm = TRUE))})
)
gt2$delta <- gt2$absrqr_B - gt2$absrqr_N
summary(gt2$absrqr_B)
summary(gt2$absrqr_N)
summary(gt2$delta)
mean(gt2$delta > 0, na.rm = TRUE)

#3.choose the dataset for plot based on the result in section 1 and 2
#dlpfc-samedonor
dataname="samedonor" 
load(here::here("RealData/result_data/model check",paste0("res_mean_boundary_",dataname,".RData")))
gt <- res$gene_tab
gname <- colnames(spelist_samedonor[[1]][[1]])
gt$gene_name <- gname[gt$gene]
load(here::here("RealData/result_data/realdata svgene","data_dlpfc_samedonor_svgene_list.RData"))
sv_names <- data_dlpfc_samedonor_svgene_list[[6]]
gt <- gt[gt$gene_name %in% sv_names, ]  
ord <- c("L1","L1|L2","L2","L2|L3","L3","L3|L4","L4","L4|L5","L5","L5|L6","L6","L6|WM","WM")

#dlpfc-acrossdonor
dataname="acrossdonor" 
load(here::here("RealData/result_data/model check",paste0("res_mean_boundary_",dataname,".RData")))
gt <- res$gene_tab
gname <- colnames(spelist_acrossdonor[[1]][[1]])
gt$gene_name <- gname[gt$gene]
load(here::here("RealData/result_data/realdata svgene","data_dlpfc_acrossdonor_svgene_list.RData"))
sv_names <- data_dlpfc_acrossdonor_svgene_list[[6]]
gt <- gt[gt$gene_name %in% sv_names, ]  
ord <- c("L1","L1|L2","L2","L2|L3","L3","L3|L4","L4","L4|L5","L5","L5|L6","L6","L6|WM","WM")

#scc
dataname="scc" 
load(here::here("RealData/result_data/model check",paste0("res_mean_boundary_",dataname,".RData")))
gt <- res$gene_tab
gname <- colnames(spelist_scc[[1]][[1]])
gt$gene_name <- gname[gt$gene]
load(here::here("RealData/result_data/realdata svgene","data_scc_svgene_list.RData"))
sv_names <- data_scc_svgene_list[[6]]
gt <- gt[gt$gene_name %in% sv_names, ]
ord <- c("1","1|5","1|7","2","2|7","2|8","3","3|8","4","4|5","4|7","4|8","5","5|7","7","7|8","8")
#NOTE: you are advised to modify the coord_cartesian(ylim=c(0.65, 1.4)) in ggplot function when you are using the SCC dataset.

#4.plot by the ggplot function. you are advised to generate p1/p2/p3 plot derived from "dlpfcsamedonor","dlpfcacrossdonor" and "scc" datasets in section 3.
ord <- ord[ord %in% unique(gt$group)]
gt <- gt[gt$group %in% ord, ]
gt$group <- factor(gt$group, levels = ord)
gt$type <- ifelse(grepl("\\|", gt$group), "Interface (boundary)", "Layer (interior)")
gt$type <- factor(gt$type, levels = c("Layer (interior)", "Interface (boundary)"))
p1 <- ggplot(gt, aes(x = group, y = absrqr, fill = type)) +
  stat_boxplot(geom = "errorbar", width = 0.4) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = sqrt(2/pi), linetype = 2, color = "grey30") +
  scale_fill_manual(values = c("Layer (interior)" = "#92C5DE",
                               "Interface (boundary)" = "#F4A582")) +
  labs(x = NULL, y = "mean |RQR| per gene", fill = NULL,
       title = dataname) +
  coord_cartesian(ylim=c(0.65, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 11))
p1

#5.combined the plots when you produce the p1/p2/p3 and save.
combined_plot <- (p1 /p2 / p3) +
  plot_layout(
    heights = c(1.5, 1.5, 1.5),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(size = 15, face = "bold"),
    plot.tag.position = c(0.01, 0.99)
  )

#plot(combined_plot)
ggsave("real_rqr_estimated.png",plot = combined_plot, width = 8, height = 10, dpi = 300)


#######################################Figures S37-S39: details of specific genes############################################################
#######################################Figures S37-S39: details of specific genes############################################################
#######################################Figures S37-S39: details of specific genes############################################################
#function
plot_gene_profile <- function(gene_name, gt2, spelist, show_ylab=TRUE, show_xlab=TRUE) {
  g_id <- which(colnames(spelist[[1]][[1]]) == gene_name)
  sub <- gt2[gt2$gene == g_id, ]
  ord <- c("L1","L1|L2","L2","L2|L3","L3","L3|L4","L4","L4|L5","L5","L5|L6","L6","L6|WM","WM")
  ord <- ord[ord %in% sub$group]
  sub$group <- factor(sub$group, levels=ord)
  agg <- aggregate(mobs ~ group, sub, mean)
  agg$group <- factor(agg$group, levels=ord); agg <- agg[order(agg$group),]
  plot(as.numeric(agg$group), agg$mobs, type="b", pch=16, lwd=2,
       xaxt="n", xlab="", ylab=if(show_ylab) "mean count" else "", main=gene_name)
  if (show_xlab) axis(1, at=seq_along(ord), labels=ord, las=2)
  else axis(1, at=seq_along(ord), labels=FALSE)   # 只画刻度不画标签
}
plot_gene_rqr_profile <- function(gene_name, gt, ord, show_ylab=TRUE, show_xlab=TRUE,ylim=c(0.5,1.2)) {
  sub <- gt[gt$gene_name == gene_name, ]
  sub <- sub[sub$group %in% ord, ]
  agg <- aggregate(absrqr ~ group, sub, mean)
  agg$group <- factor(agg$group, levels=ord); agg <- agg[order(agg$group),]
  plot(as.numeric(agg$group), agg$absrqr, type="b", pch=16, lwd=2,
       xaxt="n", xlab="", ylab=if(show_ylab) "mean |RQR|" else "",
       main=gene_name, ylim=ylim)
  abline(h=sqrt(2/pi), lty=2, col="red")            # N(0,1) 参考线
  if (show_xlab) axis(1, at=seq_along(ord), labels=ord, las=2)
  else axis(1, at=seq_along(ord), labels=FALSE)
}
plot_gene_profile_scc <- function(gene_name, gt2, spelist, show_ylab=TRUE, show_xlab=TRUE) {
  g_id <- which(colnames(spelist[[1]][[1]]) == gene_name)
  sub <- gt2[gt2$gene == g_id, ]
  ord <- c(c("1","1|5","1|7","2","2|7","2|8","3","3|8","4","4|5","4|7","4|8","5","5|7","7","7|8","8"))
  ord <- ord[ord %in% sub$group]
  sub$group <- factor(sub$group, levels=ord)
  agg <- aggregate(mobs ~ group, sub, mean)
  agg$group <- factor(agg$group, levels=ord); agg <- agg[order(agg$group),]
  plot(as.numeric(agg$group), agg$mobs, type="b", pch=16, lwd=2,
       xaxt="n", xlab="", ylab=if(show_ylab) "mean count" else "", main=gene_name)
  if (show_xlab) axis(1, at=seq_along(ord), labels=ord, las=2)
  else axis(1, at=seq_along(ord), labels=FALSE)   # 只画刻度不画标签
}

#1.DLPFC-samedonor
dataname="samedonor" 
load(here::here("RealData/result_data/model check",paste0("interface_meanbias_dlpfc", dataname, ".RData")))
res2=list("gene_tab"=gene_tab,"iface_list"=iface_list)
gt2 <- res2$gene_tab 
load(here::here("RealData/result_data/model check",paste0("res_mean_boundary_",dataname,".RData")))
gt <- res$gene_tab
gname <- colnames(spelist_samedonor[[1]][[1]])
gt$gene_name <- gname[gt$gene]
genes <- c("CXCL14","CUX2","CARTPT","HPCAL1","RORB","PCP4","MBP","PLP1")
ord <- c("L1","L1|L2","L2","L2|L3","L3","L3|L4","L4","L4|L5","L5","L5|L6","L6","L6|WM","WM")
png("profile_rqr_acrossdonor.png", width=1200, height=1200, res=130)
layout(matrix(c(
  0, 0, 0,0,
  1,  2,  3,  4,
  5,  6,  7,  8,
  0,  0,  0,  0,   
  9, 10, 11, 12,
  13, 14, 15, 16
), nrow = 6, byrow = TRUE),
heights = c(0.2,1, 1, 0.2, 1, 1))
par(mar = c(5,4,2,1), oma = c(0,3,0,0))
# row1-2：profile
for (i in seq_along(genes)) {
  show_ylab <- (i %% 4 == 1)
  show_xlab <- (i > 4)          
  plot_gene_profile(genes[i], gt2, spelist_acrossdonor, show_ylab, show_xlab)
}
# row3-4：RQR
for (i in seq_along(genes)) {
  show_ylab <- (i %% 4 == 1)
  show_xlab <- (i > 4)
  plot_gene_rqr_profile(genes[i], gt, ord, show_ylab, show_xlab)
}
mtext("A", side=2, outer=TRUE, at=0.98, line=0.5, cex=1.2, font=2, las=1, adj=1)
mtext("B", side=2, outer=TRUE, at=0.50, line=0.5, cex=1.2, font=2, las=1, adj=1)
dev.off()


#2.DLPFC-across donors
dataname="acrossdonor" 
load(here::here("RealData/result_data/model check",paste0("interface_meanbias_dlpfc", dataname, ".RData")))
res2=list("gene_tab"=gene_tab,"iface_list"=iface_list)
gt2 <- res2$gene_tab 
load(here::here("RealData/result_data/model check",paste0("res_mean_boundary_",dataname,".RData")))
gt <- res$gene_tab
gname <- colnames(spelist_acrossdonor[[1]][[1]])
gt$gene_name <- gname[gt$gene]
genes <- c("CXCL14","CUX2","CARTPT","HPCAL1","RORB","PCP4","MBP","PLP1")
ord <- c("L1","L1|L2","L2","L2|L3","L3","L3|L4","L4","L4|L5","L5","L5|L6","L6","L6|WM","WM")
png("profile_rqr_acrossdonor.png", width=1200, height=1200, res=130)
layout(matrix(c(
  0, 0, 0,0,
  1,  2,  3,  4,
  5,  6,  7,  8,
  0,  0,  0,  0,   
  9, 10, 11, 12,
  13, 14, 15, 16
), nrow = 6, byrow = TRUE),
heights = c(0.2,1, 1, 0.2, 1, 1))
par(mar = c(5,4,2,1), oma = c(0,3,0,0))
# row1-2：profile
for (i in seq_along(genes)) {
  show_ylab <- (i %% 4 == 1)
  show_xlab <- (i > 4)          
  plot_gene_profile(genes[i], gt2, spelist_acrossdonor, show_ylab, show_xlab)
}
# row3-4：RQR
for (i in seq_along(genes)) {
  show_ylab <- (i %% 4 == 1)
  show_xlab <- (i > 4)
  plot_gene_rqr_profile(genes[i], gt, ord, show_ylab, show_xlab)
}
mtext("A", side=2, outer=TRUE, at=0.98, line=0.5, cex=1.2, font=2, las=1, adj=1)
mtext("B", side=2, outer=TRUE, at=0.50, line=0.5, cex=1.2, font=2, las=1, adj=1)
dev.off()


#SCC
dataname="scc"
load(here::here("RealData/result_data/model check",paste0("interface_meanbias_", dataname, ".RData")))
res2=list("gene_tab"=gene_tab,"iface_list"=iface_list)
gt2 <- res2$gene_tab 
load(here::here("RealData/result_data/model check",paste0("res_mean_boundary_scc.RData")))
gt <- res$gene_tab
gname <- colnames(spelist_scc[[1]][[1]])
gt$gene_name <- gname[gt$gene]
genes <- c("TP63","DCN","PECAM1","LUM","LAMC2","MMP10","COL1A2","PTPRC")
ord <- c("1","1|5","1|7","2","2|7","2|8","3","3|8","4","4|5","4|7","4|8","5","5|7","7","7|8","8")
png("profile_rqr_scc.png", width=1400, height=1300, res=130)
layout(matrix(c(
  0, 0, 0,0,
  1,  2,  3,  4,
  5,  6,  7,  8,
  0,  0,  0,  0,   
  9, 10, 11, 12,
  13, 14, 15, 16
), nrow = 6, byrow = TRUE),
heights = c(0.2,1, 1, 0.2, 1, 1))
par(mar = c(5,4,2,1), oma = c(0,3,0,0))
for (i in seq_along(genes)) {
  show_ylab <- (i %% 4 == 1)
  show_xlab <- (i > 4)          
  plot_gene_profile_scc(genes[i], gt2, spelist_scc, show_ylab, show_xlab)
}
for (i in seq_along(genes)) {
  show_ylab <- (i %% 4 == 1)
  show_xlab <- (i > 4)
  plot_gene_rqr_profile(genes[i], gt, ord, show_ylab, show_xlab,ylim=c(0.4,1.4))
}
mtext("A", side=2, outer=TRUE, at=0.98, line=0.5, cex=1.2, font=2, las=1, adj=1)
mtext("B", side=2, outer=TRUE, at=0.50, line=0.5, cex=1.2, font=2, las=1, adj=1)
dev.off()








