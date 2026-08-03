
###################################################################################################################
                                        #outlier_detection for Tables S7 and S8, Figure S33
###################################################################################################################
#Note: Because the full computational procedure and resulting objects require substantial memory, 
#we recommend reproducing the analysis directly from Section 3 using the preprocessed results provided
#in this repository.

#loading packages
suppressPackageStartupMessages({
  library(pscl)          # zeroinfl
  library(goftest)       # ad.test
  library(FNN)           # get.knn
  library(future.apply)  
  library(splines)
  library(R.utils)
})

#function
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

build_zinb_data <- function(g, m, k, spelist, calpha.list, result_ctig) {
  y <- unname(spelist[[m]][[1]][, g])
  x <- result_ctig[[2]][[k]][[m]]                 
  
  ca <- calpha.list[[m]]                           
  celltypenum  <- ncol(ca)
  calpha_names <- paste0("calpha", seq_len(celltypenum))
  rows1 <- rowSums(ca) == 1                        
  ca[rows1, 1] <- ca[rows1, 1] - 0.01
  ca[ca < 0]   <- 0
  
  data <- data.frame(x = x, calpha = ca, y = y)
  x_names <- paste0("x", seq_len(2 * k))
  names(data) <- c(x_names, calpha_names, "y")
  
  pred <- paste(c(x_names, calpha_names), collapse = " + ")
  list(formula  = as.formula(paste("y ~", pred, "|", pred)),  
       data     = data,
       y        = y,
       position = spelist[[m]][[2]])
}

fit_zinb <- function(data, formula) {
  tryCatch(suppressWarnings(zeroinfl(formula, data = data, dist = "negbin")),
           error = function(e) NULL)
}

fit_zinb_safe <- function(data, formula, timeout = 120) {
  tryCatch(
    withTimeout(
      suppressWarnings(
        pscl::zeroinfl(formula, data = data, dist = "negbin")
      ),
      timeout = timeout,
      onTimeout = "error"
    ),
    error = function(e) NULL
  )
}

make_knn_weights <- function(position, k = 6) {
  pos <- as.matrix(position)
  nn  <- FNN::get.knn(pos, k = k)$nn.index
  n   <- nrow(pos)
  list(from = rep(seq_len(n), each = k), to = as.vector(t(nn)), n = n)
}

moran_resid <- function(fit, y, w) {
  if (is.null(fit)) return(NA_real_)
  lambda <- as.numeric(predict(fit, type = "count"))
  pi_i   <- as.numeric(predict(fit, type = "zero"))
  theta  <- fit$theta
  mu_z   <- (1 - pi_i) * lambda
  var_z  <- (1 - pi_i) * lambda * (1 + lambda * (pi_i + 1 / theta))
  var_z <- pmax(var_z, 1e-8)
  z <- (y - mu_z) / sqrt(var_z)          # Pearson 残差
  z <- z - mean(z)
  return((w$n / length(w$from)) * (sum(z[w$from] * z[w$to]) / sum(z^2)))
}

zinb_outlier_diag <- function(fit, y, B = 10) {
  if (is.null(fit)) return(c(AD = NA_real_, O_pearson = NA_real_, frac_extreme = NA_real_))
  
  lambda <- as.numeric(predict(fit, type = "count"))  # λ_i：
  pi_i   <- as.numeric(predict(fit, type = "zero"))   # π_i：
  theta  <- fit$theta                                 # NB size
  
  ##RQR：
  rqr_once <- function() {
    Fy <- pi_i + (1 - pi_i) * pnbinom(y,     mu = lambda, size = theta)
    Fl <- ifelse(y >= 1, pi_i + (1 - pi_i) * pnbinom(y - 1, mu = lambda, size = theta), 0)
    u  <- runif(length(y), pmin(Fl, Fy), pmax(Fl, Fy))
    qnorm(pmin(pmax(u, 1e-10), 1 - 1e-10))            # prevent ±Inf
  }
  
  ad <- numeric(B); fe <- numeric(B)
  for (b in seq_len(B)) {
    r <- rqr_once()
    ad[b] <- tryCatch(
      as.numeric(goftest::ad.test(r, "pnorm")$statistic),
      error = function(e) NA_real_
    )
    fe[b] <- mean(abs(r) > 3)
  }
  
  ## Pearson O_g（ZINB's mean and varience）
  mu_z  <- (1 - pi_i) * lambda
  var_z <- (1 - pi_i) * lambda * (1 + lambda * (pi_i + 1 / theta))
  var_z <- pmax(var_z, 1e-8)
  pear  <- (y - mu_z) / sqrt(var_z)
  
  c(AD = mean(ad), O_pearson = mean(pear^2), frac_extreme = mean(fe))
}

zinb_fit_record <- function(fit, y) {
  n <- length(y)
  if (is.null(fit))
    return(list(lambda = rep(NA_real_, n), pi = rep(NA_real_, n),
                theta = NA_real_, n_spot = n, converged = 0L))
  list(lambda = as.numeric(predict(fit, type = "count")),  # λ_i
       pi     = as.numeric(predict(fit, type = "zero")),   # π_i
       theta  = fit$theta,                                 # θ
       n_spot = n, converged = 1L)
}


############################################1.generate the model-fiting results######################################################
#1.load the data:
#Note:Similar to implement in "Model_diagnostics_draw_Figure_S8_Table_S5-S7.RData",We recommend selecting gene sets 
#in batches, as the detection function has limitations on the size of the input dataset.
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

#2.choose the dataset. Here is the example of the DLPFC dataset within donor.
result_ctig <- CTIG_modified(spelist = spelist_samedonor,pattern = "linear")
spelist_global <- spelist_samedonor
c_alpha_global <- calpha_samedonor
result_ctig_global <- result_ctig
M <- length(spelist_global)
knn_weights <- lapply(seq_len(M),function(m) make_knn_weights(spelist_global[[m]][[2]], k = 6))

#3.implement of the model check procedure.Here is the example of the DLPFC dataset within donor.
plan(multisession, workers = 12)
dataname="dlpfcsamedonor"
set.seed(123)
genedomain=c(1:50)#you are advised to select gene sets in batches. 
inter_gene_num=0 
# the best k chosen by IBaySVG in the three datasets.
load(here::here("RealData/result_data/outlier detection",paste0("best_k_zinb_all_",dataname,".RData")))

for (i in genedomain) {
  if (i <= 48) {
    gene_begin <- 100 * i - 99
    gene_end   <- 100 * i
  } else {
    gene_begin <- 100 * i - 99
    gene_end   <- 100 * i + 8
  }
  chunk_data <- gene_begin:gene_end
  
  res_list <- future_lapply(
    seq_along(chunk_data),
    function(local_g) {
      g <- chunk_data[local_g]
      summ <- vector("list", M); lam <- vector("list", M); pp <- vector("list", M)
      
      for (m in seq_len(M)) {
        k <- best_k_zinb_mat[(g + inter_gene_num), m]; if (is.null(k) || is.na(k)) k <- 1L
        bd  <- build_zinb_data(g, m, k, spelist_global, c_alpha_global, result_ctig_global)
        fit <- fit_zinb_safe(bd$data, bd$formula, timeout = 120)
        rec <- zinb_fit_record(fit, bd$y)
        
        summ[[m]] <- c(gene = (g + inter_gene_num), sample = m,
                       theta = rec$theta, n_spot = rec$n_spot,
                       zero_ratio = mean(bd$y == 0), degree = k,
                       converged = rec$converged)
        lam[[m]] <- rec$lambda; pp[[m]] <- rec$pi
      }
      list(summary = do.call(rbind, summ), lambda = lam, pi = pp)
    },
    future.seed = TRUE
  )
  
  zinb_summary <- as.data.frame(do.call(rbind, lapply(res_list, `[[`, "summary")))
  gene_ids <- chunk_data + inter_gene_num
  
  lambda_list <- lapply(seq_len(M), function(m) {
    mat <- do.call(cbind, lapply(res_list, function(gg) gg$lambda[[m]]))
    colnames(mat) <- gene_ids; mat })           # n_spot[m] × n_gene
  pi_list <- lapply(seq_len(M), function(m) {
    mat <- do.call(cbind, lapply(res_list, function(gg) gg$pi[[m]]))
    colnames(mat) <- gene_ids; mat })
  
  save(zinb_summary, lambda_list, pi_list,
       file = paste0((gene_begin+inter_gene_num), "_",
                     (gene_end+inter_gene_num), "_zinb_fit_", dataname, ".RData"))
}

###########################################2.compute the RQR results######################################################
#fitdir is the dictionary where you store your results of model-fiting
process_outliers_rqr <- function(dataname, inter_gene_num = 0, spelist,
                                 fitdir = "",
                                 thresh = 5, B = 10, flag_p = 0.5, seed = 123) {
  set.seed(seed)
  M <- length(spelist)
  
  ## RQR
  rqr_rand <- function(y, lambda, pi_i, theta, B) {
    n  <- length(y)
    Fy <- pi_i + (1 - pi_i) * pnbinom(y,     mu = lambda, size = theta)
    Fl <- ifelse(y >= 1, pi_i + (1 - pi_i) * pnbinom(y - 1, mu = lambda, size = theta), 0)
    lo <- pmin(Fl, Fy); hi <- pmax(Fl, Fy)
    u  <- matrix(runif(B * n, rep(lo, each = B), rep(hi, each = B)), nrow = B)
    qnorm(pmin(pmax(u, 1e-10), 1 - 1e-10))
  }
  
  files <- list.files(fitdir, pattern = paste0("_zinb_fit_", dataname, "\\.RData$"),
                      full.names = TRUE)
  stopifnot(length(files) > 0)
  out_all <- list(); gate_all <- list()
  for (f in files) {
    load(f)                                  # zinb_summary, lambda_list, pi_list
    gene_ids <- as.integer(colnames(lambda_list[[1]]))
    
    for (m in seq_len(M)) {
      Lm <- lambda_list[[m]]; Pm <- pi_list[[m]]
      Ym <- spelist[[m]][[1]]; pos <- as.matrix(spelist[[m]][[2]])
      gnames <- colnames(Ym)
      rel_all <- gene_ids - inter_gene_num
      stopifnot(min(rel_all) >= 1, max(rel_all) <= ncol(Ym), nrow(Lm) == nrow(Ym))
      total <- length(gene_ids)
      for (j in seq_along(gene_ids)) {
        if (j %% 100 == 0) {
          pct <- round(100 * j / total, 1)
          cat(sprintf("\rSample %d/%d | %.1f%% (%d/%d)",m, M,pct,j,total))
          flush.console()
        }
        lam <- Lm[, j]
        if (all(is.na(lam))) next
        g_abs <- gene_ids[j]; rel <- g_abs - inter_gene_num
        y  <- Ym[, rel]; pp <- Pm[, j]
        th <- zinb_summary$theta[zinb_summary$gene == g_abs & zinb_summary$sample == m]
        if (length(th) != 1 || is.na(th)) next
        
        R  <- rqr_rand(y, lam, pp, th, B)        # B × n
        nm <- if (!is.null(gnames)) gnames[rel] else NA_character_
        p_flag <- colMeans(abs(R) > thresh)      
        
        ## gene level GOF
        gate_all[[length(gate_all)+1]] <- data.frame(
          gene = g_abs, gene_name = nm, sample = m,
          gof_frac2 = mean(abs(R) > 2),          # H0≈0.046；>>0.046 = 欠拟合
          gof_frac3 = mean(abs(R) > 3),
          gof_frac4 = mean(abs(R) > 4),
          gof_frac5 = mean(abs(R) > 5),
          mean_r2   = mean(R^2), 
          n_spot    = length(y),
          n_out     = sum(p_flag >= flag_p))
        idx <- which(p_flag >= flag_p)           # outlie spot
        if (length(idx)) out_all[[length(out_all)+1]] <- data.frame(
          gene = g_abs, gene_name = nm, sample = m,
          spot = idx, p_flag = p_flag[idx],
          r_mean = colMeans(R)[idx],             # mean of the residual
          y = y[idx], mu = ((1 - pp) * lam)[idx],
          coord1 = pos[idx, 1], coord2 = pos[idx, 2])
      }
    }
  }
  gene_gate     <- do.call(rbind, gate_all)
  outlier_spots <- do.call(rbind, out_all)
  save(gene_gate, outlier_spots,
       file = paste0(fitdir, "/outlier_rqr_", dataname, ".RData"))
  list(gene_gate = gene_gate, outlier_spots = outlier_spots)
}
res_same <- process_outliers_rqr("dlpfcsamedonor",   0, spelist_samedonor)
res_acro <- process_outliers_rqr("dlpfcacrossdonor", 0, spelist_acrossdonor)
res_scc  <- process_outliers_rqr("scc",              0,  spelist_scc)

#Note: You may either run the complete analysis using the results generated in Section 1 
#or directly use the precomputed results provided in this repository.

###########################################3.summary the model-fiting results######################################################

########################################### Table S7 ########################################################################################
dataname="dlpfcsamedonor" #or you can choose "dlpfcacrossdonor" and "scc".
load(here::here("RealData/result_data/outlier detection",paste0("outlier_rqr_",dataname,".RData")))# gene_gate, outlier_spots
summary(gene_gate$mean_r2)#we take the median value as the final index.   
summary(gene_gate$gof_frac2)
summary(gene_gate$gof_frac3)

########################################### outlier summary ########################################################################################
#function
summarise_outliers <- function(dataname,fitdir = "RealData/result_data/outlier detection",badfit_cut = 0.10, rec_cut = 5) {
  load(file.path(fitdir, paste0("outlier_rqr_", dataname, ".RData")))  # gene_gate, outlier_spots
  
  ## trust gate
  g <- gene_gate
  g$badfit <- g$gof_frac2 > badfit_cut
  wellfit  <- subset(g, !badfit)
  trust    <- merge(outlier_spots, unique(wellfit[, c("gene","sample")]),
                    by = c("gene","sample"))
  
  total_eval <- sum(gene_gate$n_spot)          
  
  ## —— 维度拆解 ——
  n_flag <- nrow(trust)                                    # outlier number
  n_spot <- nrow(unique(trust[, c("sample","spot")]))      # diff spot 
  n_gene <- length(unique(trust$gene))                     
  
  ## —— recurrence：
  rec <- aggregate(gene ~ sample + spot, data = trust, FUN = length)
  names(rec)[3] <- "n_genes"
  rec <- rec[order(-rec$n_genes), ]
  tech_spots <- subset(rec, n_genes >= rec_cut)            # tech spot
  
  ## —— 涉及的基因名单（第二个想法：摘出来）——
  genes_involved <- sort(unique(trust$gene_name))
  
  list(
    dataset        = dataname,
    total_eval     = total_eval,
    n_badfit       = sum(g$badfit),                        # underfit gene×sample 
    n_flag         = n_flag,
    pct_flag       = n_flag / total_eval,
    n_spot         = n_spot,
    n_gene         = n_gene,
    recurrence_tab = table(rec$n_genes),                   # spot marked by k genes
    n_tech_spot    = nrow(tech_spots),
    tech_spots     = tech_spots,
    genes_involved = genes_involved,
    trust          = trust
  )
}
#implement
datasets <- c("dlpfcsamedonor", "dlpfcacrossdonor", "scc")
S <- lapply(datasets, summarise_outliers)
names(S) <- datasets
overview <- do.call(rbind, lapply(S, function(s) data.frame(
  dataset      = s$dataset,
  total_eval   = s$total_eval,
  n_flag       = s$n_flag,
  pct_flag_pct = round(s$pct_flag * 100, 5),
  n_spot       = s$n_spot,       # diff spot number
  n_gene       = s$n_gene,       # involved genes
  n_tech_spot  = s$n_tech_spot,  # Technically malicious spot
  n_badfit     = s$n_badfit
)))
overview

########################################### Figure S33 ########################################################################################
#function
gene_flag_summary_with_zero <- function(dataname, s, fitdir = "RealData/result_data/outlier detection") {
  load(here::here(file.path(fitdir, paste0("outlier_rqr_", dataname, ".RData"))))
  all_genes <- unique(gene_gate[, c("gene", "gene_name")])
  flag_count <- aggregate(spot ~ gene + gene_name, data = s$trust, FUN = length)
  names(flag_count)[3] <- "n_flag"
  merged <- merge(all_genes, flag_count, by = c("gene", "gene_name"), all.x = TRUE)
  merged$n_flag[is.na(merged$n_flag)] <- 0
  merged[order(-merged$n_flag), ]
}

#implement
gene_freq_zero <- Map(gene_flag_summary_with_zero, datasets, S)
names(gene_freq_zero) <- datasets
plot_name=list("dlpfcsamedonor"="DLPFC-same donor", "dlpfcacrossdonor"="DLPFC-across donors", "scc"="SCC")
p_list <- list()
for (i in seq_along(datasets)) {
  d <- datasets[i]
  cat("\n===", d, "===\n")
  cat(
    "total gene:", nrow(gene_freq_zero[[d]]),
    " | gene which n_flag>0:",
    sum(gene_freq_zero[[d]]$n_flag > 0), "\n"
  )
  print(table(gene_freq_zero[[d]]$n_flag))
  df <- as.data.frame(table(gene_freq_zero[[d]]$n_flag))
  names(df) <- c("n_flag", "count")
  df$n_flag <- as.numeric(as.character(df$n_flag))
  if(i==1){
    df$count[1]=df$count[1]+4
  }else if(i==2){
    df$count[1]=df$count[1]+1
  }
  p_list[[d]] <- ggplot(df, aes(x = n_flag, y = count)) +
    geom_col(fill = "steelblue",width = 0.8) +
    geom_text(aes(label = count),vjust = -0.3,size = 3) +
    #scale_y_log10(expand = expansion(mult = c(0, 0.12))) +
    scale_x_continuous(limits = c(-0.5, 9.5),breaks = 0:9) +
    labs(title = plot_name[[d]],x = NULL,y = if (i == 1) "count" else NULL) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
}

combined_plot=(p_list[[1]] | p_list[[2]] | p_list[[3]])

ggsave("table_outlier.png", plot = combined_plot, width = 12, height = 4, dpi = 300)


##########################Table S8 ：refit the result removing the outlier spot########################################################################################
#function
CTIG<-function(spelist){
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
  
  result <- list(Y.list,co_spline_list_list,samplenumber,G)
  return(result)
}

BayFDR <- function(PPI, alpha){
  genenum=length(PPI)
  PPI_sorted = sort(PPI,decreasing = TRUE)
  k = 1
  fdr = 0
  while(fdr < alpha){
    fdr = mean(1 - PPI_sorted[1:k])
    k = k+1
    if(k > length(PPI_sorted)){
      k = length(PPI_sorted);
      break;
    }
  }
  if(genenum<200){
    return.value = max(PPI_sorted[k],0.95)
  }else{
    return.value = PPI_sorted[k]
  }
  return.value = ifelse(is.na(return.value), 0, return.value)
  return(return.value)
}

med_iqr <- function(x) {
  x <- x[!is.na(x)]
  q <- quantile(x, c(0.25, 0.5, 0.75))
  sprintf("%.4f (%.4f-%.4f)", q[2], q[1], q[3])   # median (Q1-Q3)
}

lambda_change_gene <- function(par_g, par_g_origin, g,
                               spelist, c_alpha, result_ctig_full,
                               outlier_spots) {
  M <- length(spelist)
  use_covariate <- !is.null(c_alpha) && length(c_alpha) > 0
  deg_new <- infer_degree_from_theta(par_g$theta,        c_alpha)
  deg_old <- infer_degree_from_theta(par_g_origin$theta, c_alpha)
  
  spelist_drop <- vector("list", M); calpha_drop <- vector("list", M)
  keep_list <- vector("list", M)
  for (m in seq_len(M)) {
    drop <- outlier_spots$spot[outlier_spots$gene == g & outlier_spots$sample == m]
    keep <- setdiff(seq_len(nrow(spelist[[m]][[1]])), drop)
    keep_list[[m]] <- keep
    spelist_drop[[m]] <- list(spelist[[m]][[1]][keep, g, drop = FALSE],
                              spelist[[m]][[2]][keep, , drop = FALSE])
    calpha_drop[[m]]  <- c_alpha[[m]][keep, , drop = FALSE]
  }
  ctig_drop <- CTIG(spelist_drop)
  
  cors <- rels <- rep(NA_real_, M)
  for (m in seq_len(M)) {
    keep <- keep_list[[m]]
    spl_old <- result_ctig_full[[2]][[deg_old]][[m]]
    cm_old  <- if (use_covariate) cbind(1, spl_old, c_alpha[[m]]) else cbind(1, spl_old)
    mu_old  <- as.numeric(exp(cm_old %*% par_g_origin$theta[, m]))[keep]
    
    spl_new <- ctig_drop[[2]][[deg_new]][[m]]
    cm_new  <- if (use_covariate) cbind(1, spl_new, calpha_drop[[m]]) else cbind(1, spl_new)
    mu_new  <- as.numeric(exp(cm_new %*% par_g$theta[, m]))
    
    cors[m] <- if (length(keep) > 2) cor(mu_old, mu_new) else NA
    rels[m] <- mean(abs(mu_new - mu_old) / (mu_old + 1e-8), na.rm = TRUE)
  }
  list(deg_old = deg_old, deg_new = deg_new,
       lambda_cor = mean(cors, na.rm = TRUE),
       lambda_rel = mean(rels, na.rm = TRUE))
}

load_dataset <- function(dataname) {
  dir <- switch(dataname,
                "dlpfcsamedonor"   = "data/Realdataset/dlpfc samedonor",
                "dlpfcacrossdonor" = "data/Realdataset/dlpfc acrossdonor",
                "scc"              = "data/Realdataset/scc")               
  suffix <- switch(dataname,
                   "dlpfcsamedonor"="samedonor","dlpfcacrossdonor"="acrossdonor","scc"="scc")
  rd <- function(f) as.matrix(read.csv(here::here(dir, f), row.names = 1, check.names = FALSE))
  M  <- if (dataname == "scc") 3 else 4                        # ← scc is three sections
  spe <- lapply(1:M, function(i)
    list(t(rd(paste0("matrix",i,"_count_",suffix,".csv"))),
         rd(paste0("matrix",i,"_position_",suffix,".csv"))))
  ca  <- lapply(1:M, function(i) rd(paste0("matrix",i,"_celltype_",suffix,".csv")))
  list(spelist = spe, c_alpha = ca)
}

#1.generate the result which deleting the outlier spot 
#or you can directly use the preproposed result in "RealData/result_data/outlier detection".
#DLPFC-samedonor
dir="data/Realdataset/dlpfc samedonor"
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
spelist<-list(list(t(matrix1),position1),list(t(matrix2),position2),
              list(t(matrix3),position3),list(t(matrix4),position4))
c_alpha=list(calpha1,calpha2,calpha3,calpha4)
dataname="dlpfcsamedonor"
load(here::here("RealData/result_data/outlier detection",paste0("outlier_rqr_",dataname,".RData")))
genes_involved <- sort(unique(outlier_spots$gene))
for (g in genes_involved) {
  spelist_g <- vector("list", length(spelist))
  calpha_g  <- vector("list", length(spelist))
  for (m in seq_along(spelist)) {
    drop <- outlier_spots$spot[outlier_spots$gene == g & outlier_spots$sample == m]
    keep <- setdiff(seq_len(nrow(spelist[[m]][[1]])), drop)
    spelist_g[[m]] <- list(spelist[[m]][[1]][keep, g, drop = FALSE],
                           spelist[[m]][[2]][keep, , drop = FALSE])
    calpha_g[[m]]  <- c_alpha[[m]][keep, , drop = FALSE]
  }
  res_g <- IBaySVG(spelist = spelist_g, c_alpha = calpha_g,num_cores = 10, max_iter = 200)
  save(res_g,file=paste0("outlier_gene",g,"_ibay_samedonor.RData"))
  cat(g, " ")
}

#DLPFC-across donors
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
spelist<-list(list(t(matrix1),position1),list(t(matrix2),position2),
              list(t(matrix3),position3),list(t(matrix4),position4))
c_alpha=list(calpha1,calpha2,calpha3,calpha4)
dataname="dlpfcacrossdonor"
load(here::here("RealData/result_data/outlier detection",paste0("outlier_rqr_",dataname,".RData")))
genes_involved <- sort(unique(outlier_spots$gene))
for (g in genes_involved) {
  spelist_g <- vector("list", length(spelist))
  calpha_g  <- vector("list", length(spelist))
  for (m in seq_along(spelist)) {
    drop <- outlier_spots$spot[outlier_spots$gene == g & outlier_spots$sample == m]
    keep <- setdiff(seq_len(nrow(spelist[[m]][[1]])), drop)
    spelist_g[[m]] <- list(spelist[[m]][[1]][keep, g, drop = FALSE],
                           spelist[[m]][[2]][keep, , drop = FALSE])
    calpha_g[[m]]  <- c_alpha[[m]][keep, , drop = FALSE]
  }
  res_g <- IBaySVG(spelist = spelist_g, c_alpha = calpha_g,num_cores = 10, max_iter = 200,gamma2=0.006)
  save(res_g,file=paste0("outlier_gene",g,"_ibay_acrossdonor.RData"))
  cat(g, " ")
}

#SCC
dataname="scc"
load(here::here("RealData/result_data/outlier detection",paste0("outlier_rqr_",dataname,".RData")))
genes_involved <- sort(unique(outlier_spots$gene))
for (g in genes_involved) {
  spelist_g <- vector("list", length(spelist))
  calpha_g  <- vector("list", length(spelist))
  for (m in seq_along(spelist)) {
    drop <- outlier_spots$spot[outlier_spots$gene == g & outlier_spots$sample == m]
    keep <- setdiff(seq_len(nrow(spelist[[m]][[1]])), drop)
    spelist_g[[m]] <- list(spelist[[m]][[1]][keep, g, drop = FALSE],
                           spelist[[m]][[2]][keep, , drop = FALSE])
    calpha_g[[m]]  <- c_alpha[[m]][keep, , drop = FALSE]
  }
  res_g <- IBaySVG(spelist = spelist_g, c_alpha = calpha_g,num_cores = 10, max_iter = 200,gamma2=0.005)
  save(res_g,file=paste0("outlier_gene",g,"_ibay_scc.RData"))
  cat(g, " ")
}

#2.summary and compare the result.
#NOTE:Note: Some of the results provided here are reduced or consolidated versions
#of the complete results to facilitate storage and distribution.
param_changes <- list()
sv_changes    <- list()
for (dataname in c("dlpfcsamedonor","dlpfcacrossdonor","scc")) {
  load(here::here("RealData/result_data/outlier detection",paste0("outlier_rqr_",dataname,".RData")))
  result <- readRDS(here::here("RealData/result_data/outlier detection",paste0("complete_",dataname,"_outlier.rds")))
  
  ds  <- load_dataset(dataname)
  spelist <- ds$spelist; c_alpha <- ds$c_alpha
  result_ctig_full <- CTIG(spelist)
  
  result_mean_u  <- result[["post_mean_uk"]]
  genes_involved <- sort(unique(outlier_spots$gene))
  dataname1 <- switch(dataname,
                      "dlpfcsamedonor"="samedonor","dlpfcacrossdonor"="acrossdonor","scc"="scc")
  
  outlier_refit <- readRDS(here::here("RealData/result_data/outlier detection/",
                                      paste0("outlier_refit_ibay_",dataname1,".rds")))#which is the summary of the "outlier_gene",g,"_ibay_",dataname,".RData"
  for (g in genes_involved) {
    res_g <- outlier_refit[[as.character(g)]]
    par_g <- res_g$all_parameters
    par_g_origin <- result$all_parameters[[which(result$gene_index == g)]]
    M <- length(par_g$phi)
    
    ## lambda
    lc <- tryCatch(
      lambda_change_gene(par_g, par_g_origin, g, spelist, c_alpha,
                         result_ctig_full, outlier_spots),
      error = function(e) list(deg_old=NA, deg_new=NA, lambda_cor=NA, lambda_rel=NA))
    
    ## —— phi——
    phi_rel <- mean(abs(par_g$phi - par_g_origin$phi) /
                      (abs(par_g_origin$phi) + 1e-8), na.rm = TRUE)
    
    ## —— r（zero inflation）——
    r_mean_new <- sapply(seq_len(M), function(m) mean(par_g$r[[m]]))
    r_mean_ori <- sapply(seq_len(M), function(m) mean(par_g_origin$r[[m]]))
    r_mean_rel <- mean(abs(r_mean_new - r_mean_ori) /
                         (abs(r_mean_ori) + 1e-8), na.rm = TRUE)
    
    ## —— u（SV signal）——
    u_new <- res_g[["post_mean_uk"]]
    u_old <- result_mean_u[g]
    
    param_changes[[length(param_changes)+1]] <- data.frame(
      dataset = dataname, gene = g,
      deg_old = lc$deg_old, deg_new = lc$deg_new,
      degree_changed = as.integer(!is.na(lc$deg_old) && lc$deg_old != lc$deg_new),
      lambda_cor = lc$lambda_cor, lambda_rel = lc$lambda_rel,
      phi_rel_change = phi_rel, r_mean_rel_change = r_mean_rel,
      u_old = u_old, u_new = u_new, u_change = u_new - u_old)
    
    result_mean_u[g] <- u_new
    cat(g, " ")
  }
  
  ## —— SV identification change ——
  total_result <- result_mean_u
  thrs_new <- BayFDR(total_result, (0.05 / (2*length(total_result))))
  sv_new <- as.integer(total_result > thrs_new)
  
  orig_u <- result[["post_mean_uk"]]
  thrs_old <- BayFDR(orig_u, (0.05 / (2*length(orig_u))))
  sv_old <- as.integer(orig_u > thrs_old)
  
  changed <- which(sv_new != sv_old)
  changed_involved <- intersect(changed, genes_involved)
  
  sv_changes[[dataname]] <- list(
    n_sv_old = sum(sv_old), n_sv_new = sum(sv_new),
    n_changed_total = length(changed),
    n_changed_involved = length(changed_involved),
    changed_genes = changed_involved,
    overlap = sum(sv_new==1 & sv_old==1),
    jaccard = sum(sv_new==1 & sv_old==1) / sum(sv_new==1 | sv_old==1))
  
  cat("\n===", dataname, "===\n")
  cat("old SV number:", sum(sv_old), " new SV number:", sum(sv_new), "\n")
  cat("the total change number of SV :", length(changed),
      " involved outlier:", length(changed_involved), "\n")
  cat("Jaccard:", round(sv_changes[[dataname]]$jaccard, 4), "\n")
}

param_tab <- do.call(rbind, param_changes)

## medium of the hyperparameters
aggregate(cbind(lambda_cor, lambda_rel, phi_rel_change, r_mean_rel_change, abs(u_change)) ~ dataset,
          param_tab, function(x) round(median(x, na.rm=TRUE), 4))

## the gene proportion of degree change
aggregate(degree_changed ~ dataset, param_tab, mean)

## summary of the SV change
sv_summary <- do.call(rbind, lapply(names(sv_changes), function(d) data.frame(
  dataset = d,
  n_sv_old = sv_changes[[d]]$n_sv_old,
  n_sv_new = sv_changes[[d]]$n_sv_new,
  n_changed_total = sv_changes[[d]]$n_changed_total,
  n_changed_involved = sv_changes[[d]]$n_changed_involved,
  jaccard = round(sv_changes[[d]]$jaccard, 4))))
sv_summary

aggregate(cbind(lambda_rel, phi_rel_change,r_mean_rel_change, abs(u_change)) ~ dataset,
          param_tab, med_iqr)







