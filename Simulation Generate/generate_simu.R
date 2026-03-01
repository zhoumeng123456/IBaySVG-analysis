library(MCMCpack)
library(parallel)

#' Simulation function for synthetic spatial transcriptomic data under ZINB framework
#'
#' @param gene_size Integer. Total number of genes generated in the simulation.
#' @param svgene_size Numeric (0–1). Proportion of spatially variable (SV) genes 
#'        among all generated genes.
#' @param sv_mark Numeric vector of length 2. Signal strength parameters for SV genes 
#'        along the x- and y-coordinate spatial effect functions.
#' @param no_sv_mark Numeric vector of length 2. Signal strength parameters for 
#'        non-SV genes (baseline spatial effect; typically set to c(0,0)).
#' @param inf_size Numeric (0–1). Zero-inflation probability in the ZINB model, 
#'        controlling the dropout rate.
#' @param phi Numeric. Dispersion (shape) parameter of the zero-inflated negative 
#' @param etamean Numeric. Baseline mean expression level in the log-scale mean 
#'        component of the ZINB model. Default is 2.
#' @param xspace Character. Spatial effect function applied to the x-coordinate. 
#'        Options include three canonical spatial patterns (e.g., linear, focal, 
#'        periodic) and additional heterogeneous spatial structures.
#' @param yspace Character. Spatial effect function applied to the y-coordinate, 
#'        defined similarly to `xspace`.
#' @param seed Integer. Random seed for reproducibility.
#' @param domainnum Integer. Specifies the covariate distribution scenario 
#'        (e.g., different Dirichlet configurations). Default is 5.
#'
#' @param use_covariate Logical. Whether to include covariate effects in the 
#'        mean expression model.
sim_create <- function(gene_size =100,svgene_size=0.1,sv_mark=c(1,1),no_sv_mark = c(0,0),inf_size=0.1,phi=15,etamean=2,
                       xspace="linear",yspace="linear",seed=1,domainnum=5,use_covariate=TRUE){
  ###
  ###
  set.seed(seed)
  x_coords <- rep(0:31, each = 32)
  y_coords <- rep(0:31, times = 32)
  if(domainnum==1){
    z_type <- c(rep(rep(1:2,each=16),9),rep(1,96),rep(c(rep(3,12),rep(4,20)),10),rep(3,128),rep(rep(3:4,each=16),6))#横平竖直的分成了四个区域
  }else if(domainnum==2){
    z_type <- rep(1,1024)
  }else if(domainnum==3){
    z_type <- c(rep(1,512),rep(2,512))
  }else if(domainnum==4){
    z_type <- c(rep(rep(1:2,each=16),16),rep(3,512))
  }else if(domainnum==5){
    z_type <- c(rep(c(rep(1,20),rep(2,12)),16),rep(c(rep(3,12),rep(1,8),rep(4,12)),4),rep(rep(3:4,each=16),12))
  }else if(domainnum==6){
    z_type <- c(rep(c(rep(2,12),rep(1,20)),16),rep(c(rep(2,12),rep(3,20)),6),rep(3,320))
  }else if(domainnum==7){
    z_type <- c(rep(c(rep(2,12),rep(1,20)),16),rep(c(rep(2,12),rep(3,20)),6),rep(5,320))
  }else if(domainnum==8){
    z_type <- c(rep(c(rep(2,12),rep(1,20)),16),rep(c(rep(2,12),rep(3,20)),6),rep(4,320))
  }else if(domainnum==9){
    z_type <- c(rep(c(rep(1,20),rep(2,12)),16),rep(c(rep(3,12),rep(1,8),rep(4,12)),8),rep(rep(3:4,each=16),8))
  }else if(domainnum==10){
    z_type <- c(rep(c(rep(1,20),rep(2,12)),16),rep(c(rep(3,4),rep(1,16),rep(4,12)),4),rep(rep(3:4,each=16),12))
  }else if(domainnum==11){
    #
    z_type <- c(rep(rep(2:3,each=16),16),rep(1,512))
  }
  
  
  location <- as.matrix(data.frame(x = x_coords, y = y_coords,z=z_type))
  
  x <- location[,1]
  y <- location[,2]
  
  x <- x-mean(x)
  x <- x/sd(x)
  y <- y-mean(y)
  y <- y/sd(y)
  
  x <- switch(xspace, "focal" = exp(-x^2/2),
              "period" = cos(2*pi*x),
              "linear" = x,
              "sigmoid"=1/(1+exp(-x)),
              "polynomial1"=0.5 * (x + 1) * (x - 0.8) * (x - 1.6),
              "polynomial2"=-0.5*(x^3)+0.3*x,
              "polynomial3"=0.15*(x^4)-0.1*x^2+0.7,
              "polynomial4"=0.25 * x^3 + 0.1 * x^2 - 0.15 * x + 0.3,
              stop("Invalid xspace!"))
  
  y <- switch(yspace,
              "focal" = exp(-y^2/2),
              "period" = cos(2*pi*y),
              "linear" = y,
              "sigmoid"=1/(1+exp(-y)),
              "polynomial1"=0.5 * (y + 1) * (y - 0.8) * (y - 1.6),
              "polynomial2"=-0.5*(y^3)+0.3*y,
              "polynomial3"=0.15*(y^4)-0.1*y^2+0.7,
              "polynomial4"=0.25 * y^3 + 0.1 * y^2 - 0.15 * y + 0.3,
              stop("Invalid yspace!"))
  
  kern_coord<-cbind(x,y)
  
  npoints <- nrow(location)
  rownames(location) = paste('spot', 1:npoints, sep = '')
  expres_marx = as.data.frame(matrix(NA, nrow = npoints, ncol = gene_size))
  rownames(expres_marx) = paste('spot', 1:npoints, sep = '')
  colnames(expres_marx) = paste('gene', 1:gene_size, sep = '')
  
  sv_points=svgene_size*gene_size
  sv_gene <- c(1:sv_points)
  no_sv_gene <- setdiff(1:gene_size, sv_gene)
  
  eta <- rnorm(gene_size,mean = 2,sd = 0.5)
  
  cell_matrix <- matrix(NA,nrow = npoints,ncol = 6)
  
  for(i in 1:npoints){
    if(z_type[i]==1){
      cell_matrix[i,] <- MCMCpack::rdirichlet(1,alpha = c(1,1,1,1,1,1))
    }else if(z_type[i]==2){
      cell_matrix[i,] <- MCMCpack::rdirichlet(1,alpha = c(1,3,5,7,9,11))
    }else if(z_type[i]==3){
      cell_matrix[i,] <- MCMCpack::rdirichlet(1,alpha = c(14,12,10,8,6,4))
    }else if(z_type[i]==4){
      cell_matrix[i,] <- MCMCpack::rdirichlet(1,alpha = c(1,4,4,4,4,1))
    }else{
      cell_matrix[i,] <- MCMCpack::rdirichlet(1,alpha = c(18,16,14,12,10,8))
    }
  }
  
  cell_mark <- rnorm(6,mean=0,sd=1)
  
  for(i in sv_gene){
    for(t in 1:npoints){
      if(use_covariate){
        mu = exp(eta[i]+sum(kern_coord[t,]*sv_mark)+sum(cell_matrix[t,]*cell_mark))
      }else{
        mu = exp(eta[i]+sum(kern_coord[t,]*sv_mark))
      }
      expres_marx[t,i] <- rnbinom(1,mu=mu,size=phi)
    }
  }
  for(i in no_sv_gene){
    for(t in 1:npoints){
      if(use_covariate){
        mu= exp(eta[i]+sum(kern_coord[t,]*no_sv_mark)+sum(cell_matrix[t,]*cell_mark))
      }else{
        mu= exp(eta[i]+sum(kern_coord[t,]*no_sv_mark))
      }
      expres_marx[t,i] <- rnbinom(1,mu=mu,size=phi)
    }
  }
  
  
  expres_marx <- as.matrix(expres_marx)
  total.size <- npoints*gene_size
  zero_size<-floor(npoints*gene_size*inf_size)
  zeroin<-sample(c(1:total.size),zero_size)
  expres_marx[zeroin]<-0
  location <- as.matrix(location[,1:2])
  cell_matrix=as.matrix(cell_matrix)
  rownames(cell_matrix)=paste("spot",c(1:npoints))
  rownames(location)=paste("spot",c(1:npoints))
  rownames(expres_marx)=paste("spot",c(1:npoints))
  colnames(expres_marx)=paste("gene",c(1:gene_size))
  colnames(location)=c("x","y")
  colnames(cell_matrix)=paste0("cell",c(1:6))
  if(!use_covariate){
    cell_matrix=NULL
  }
  spe <- list(expres_marx,location,cell_matrix)
  
  return(spe)
}



#' Simulation function for synthetic spatial transcriptomic data under Gaussian process framework
#' @param gene_size Integer. Total number of genes generated in the simulation.
#' @param svgene_size Numeric (0–1). Proportion of spatially variable (SV) genes 
#'        among all generated genes.
#' @param inf_size Numeric (0–1). Zero-inflation probability in the ZINB model, 
#'        controlling the dropout rate.
#' @param seed Integer. Random seed for reproducibility.
#' @param domainnum Integer. Specifies the covariate distribution scenario 
#'        (e.g., different Dirichlet configurations). Default is 1.
#' @param use_covariate Logical. Whether to include covariate effects in the 
#'        mean expression model.
#' @param num_cores Number of CPU cores used for parallel simulation.
#' @param kernel Type of spatial covariance kernel for SV genes.
#'                    Options:
#'                      - "exp" : Exponential kernel (smooth decay)
#'                      - "cos" : Cosine kernel (periodic pattern)
#'
#' @param tau1 Strength of spatial differential expression effect where 9 (high signal), 4 (moderate signal), 
#'   1 (low signal), 0.25 (very weak signal).
#' @param tau2 Baseline variance (nugget effect). 
#' 
normal_create <- function(gene_size = 100,svgene_size = 0.1,inf_size = 0,kernel = "exp",seed = 1,domainnum = 1,
                          tau1 = 1,tau2 = 1,num_cores = 8,use_covariate = TRUE){
  
  set.seed(seed)
  
  # ----------------------------------------------------------
  # Step 1: Generate spatial grid (32 × 32 lattice)
  # ----------------------------------------------------------
  x_coords <- rep(0:31, each = 32)
  y_coords <- rep(0:31, times = 32)
  
  # ----------------------------------------------------------
  # Step 2: Define spatial domain structures
  # Different domainnum values correspond to different
  # spatial covariate distribution scenarios
  # ----------------------------------------------------------
  
  if(domainnum==1){
    z_type <- c(rep(rep(1:2,each=16),9),rep(1,96),
                rep(c(rep(3,12),rep(4,20)),10),
                rep(3,128),rep(rep(3:4,each=16),6))
  }else if(domainnum==2){
    z_type <- rep(1,1024)
  }else if(domainnum==3){
    z_type <- c(rep(1,512),rep(2,512))
  }else if(domainnum==4){
    z_type <- c(rep(rep(1:2,each=16),16),rep(3,512))
  }else if(domainnum==5){
    z_type <- c(rep(c(rep(1,20),rep(2,12)),16),
                rep(c(rep(3,12),rep(1,8),rep(4,12)),4),
                rep(rep(3:4,each=16),12))
  }else if(domainnum==6){
    z_type <- c(rep(c(rep(2,12),rep(1,20)),16),
                rep(c(rep(2,12),rep(3,20)),6),
                rep(3,320))
  }else if(domainnum==7){
    z_type <- c(rep(c(rep(2,12),rep(1,20)),16),
                rep(c(rep(2,12),rep(3,20)),6),
                rep(5,320))
  }else if(domainnum==8){
    z_type <- c(rep(c(rep(2,12),rep(1,20)),16),
                rep(c(rep(2,12),rep(3,20)),6),
                rep(4,320))
  }else if(domainnum==9){
    z_type <- c(rep(c(rep(1,20),rep(2,12)),16),
                rep(c(rep(3,12),rep(1,8),rep(4,12)),8),
                rep(rep(3:4,each=16),8))
  }else if(domainnum==10){
    z_type <- c(rep(c(rep(1,20),rep(2,12)),16),
                rep(c(rep(3,4),rep(1,16),rep(4,12)),4),
                rep(rep(3:4,each=16),12))
  }else if(domainnum==11){
    z_type <- c(rep(rep(2:3,each=16),16),rep(1,512))
  }
  
  location <- as.matrix(data.frame(x = x_coords, y = y_coords, z = z_type))
  
  # Standardize spatial coordinates
  x <- scale(location[,1])
  y <- scale(location[,2])
  
  # ----------------------------------------------------------
  # Step 3: Construct spatial covariance matrix
  # ----------------------------------------------------------
  
  kernel_matrix <- matrix(NA, 1024, 1024)
  
  for (i in 1:1024) {
    for (j in 1:1024) {
      # Euclidean distance
      kernel_matrix[i, j] <- sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2)
    }
  }
  
  base_corr <- tau2 * diag(rep(1,1024))
  
  if(kernel=="exp"){
    kernel_matrix <- exp(-kernel_matrix^2/2)*tau1 + base_corr
  }else if(kernel=="cos"){
    kernel_matrix <- cos(2*pi*kernel_matrix)*tau1 + base_corr
  }
  
  kernel_matrix <- round(kernel_matrix,3)
  kernel_matrix <- make_positive_definite(kernel_matrix)
  
  # ----------------------------------------------------------
  # Step 4: Generate covariate matrix (cell composition)
  # ----------------------------------------------------------
  
  npoints <- nrow(location)
  
  cell_matrix <- matrix(NA, nrow = npoints, ncol = 6)
  
  for(i in 1:npoints){
    if(z_type[i]==1){
      cell_matrix[i,] <- rdirichlet(1,c(1,1,1,1,1,1))
    }else if(z_type[i]==2){
      cell_matrix[i,] <- rdirichlet(1,c(1,3,5,7,9,11))
    }else if(z_type[i]==3){
      cell_matrix[i,] <- rdirichlet(1,c(14,12,10,8,6,4))
    }else if(z_type[i]==4){
      cell_matrix[i,] <- rdirichlet(1,c(1,4,4,4,4,1))
    }else{
      cell_matrix[i,] <- rdirichlet(1,c(18,16,14,12,10,8))
    }
  }
  
  # ----------------------------------------------------------
  # Step 5: Covariate effect (optional)
  # ----------------------------------------------------------
  
  if(use_covariate){
    cell_mark <- rnorm(6, mean=0, sd=1)
  }else{
    cell_mark <- rep(0, 6)
  }
  
  nor_mean <- cell_matrix %*% cell_mark
  nor_mean <- as.vector(nor_mean)
  
  # ----------------------------------------------------------
  # Step 6: Generate gene expression
  # SV genes follow spatial GP
  # Non-SV genes follow independent noise
  # ----------------------------------------------------------
  
  expres_marx <- matrix(NA, nrow = npoints, ncol = gene_size)
  
  sv_points <- floor(svgene_size * gene_size)
  sv_gene <- 1:sv_points
  no_sv_gene <- setdiff(1:gene_size, sv_gene)
  
  num_cores <- min(num_cores, parallel::detectCores() - 2)
  cl <- makeCluster(num_cores)
  clusterEvalQ(cl, library(MASS))
  
  sv_results <- parLapply(cl, sv_gene, function(i) {
    mvrnorm(n = 1, mu = nor_mean, Sigma = kernel_matrix)
  })
  
  no_sv_results <- parLapply(cl, no_sv_gene, function(i) {
    mvrnorm(n = 1, mu = nor_mean, Sigma = base_corr)
  })
  
  stopCluster(cl)
  
  for (i in seq_along(sv_gene)) {
    expres_marx[, sv_gene[i]] <- sv_results[[i]]
  }
  
  for (i in seq_along(no_sv_gene)) {
    expres_marx[, no_sv_gene[i]] <- no_sv_results[[i]]
  }
  
  # Truncate negative values to zero
  expres_marx[expres_marx < 0] <- 0
  expres_marx <- round(expres_marx)
  
  # ----------------------------------------------------------
  # Step 7: Add zero-inflation
  # ----------------------------------------------------------
  
  total.size <- npoints * gene_size
  zero_size <- floor(total.size * inf_size)
  zeroin <- sample(1:total.size, zero_size)
  expres_marx[zeroin] <- 0
  
  # ----------------------------------------------------------
  # Output
  # ----------------------------------------------------------
  
  if(!use_covariate){
    cell_matrix=NULL
  }
  
  return(list(
    count_matrix = expres_marx,
    location = location[,1:2],
    covariate_matrix = cell_matrix
  ))
}



#' Parallel simulation data generation wrapper under ZINB framework
#'
#' This function serves as a parallelized wrapper of `sim_create()`, 
#' enabling batch generation of synthetic spatial transcriptomic datasets 
#' under multiple signal situations.
#'
#' @param k1 Character. Spatial effect function for the x-coordinate.
#' @param k2 Character. Spatial effect function for the y-coordinate.
#'
#' @param seeddomain Integer vector. Random seeds used for generating 
#'        multiple simulation domains.
#'
#' @param random1 Numeric vector. Parameter configurations controlling 
#'        domain-level random settings (passed to `sim_create()`).
#'
#' @param svmark Numeric vector. Signal strength levels for spatially 
#'        variable (SV) genes.
#'
#' @param infdomain Numeric vector. Zero-inflation probabilities for 
#'        different simulation scenarios.
#'
#' @param genesize Integer. Total number of genes generated per simulation.
#'
#' @param use_covariate Logical. Whether to include covariate effects 
#'        in the simulation model.
#'
#' @param curdir Character. Directory path for saving generated datasets.
#'
intedata <- function(k1 = "linear",k2="linear", seeddomain = c(1:5), random1 = rep(5, 8), svmark = c(0.8, 0.5, 0.2, 0.05), infdomain=c(0.1, 0.3, 0.5),
                     genesize = 5000,use_covariate=TRUE,curdir="/gemini/code/simulation/") {
  for(seed1 in seeddomain){
    sv1 <- rep(svmark[1], 2)
    sv2 <- rep(svmark[2], 2)
    sv3 <- rep(svmark[3], 2)
    sv4 <- rep(svmark[4], 2) 
    svgenesize <- genesize * 0.1
    
    seeds <- seed1*8 + 0:7
    svmarks <- list(sv1, sv2, sv2, sv2, sv3, sv3, sv4, sv4)
    
    # Define a function for generating a single dataset
    generate_data <- function(i1, i) {
      sim_create(
        svgene_size = 0.1,
        sv_mark = svmarks[[i]],
        no_sv_mark = rep(0, 2),
        inf_size = i1,
        phi = 15,
        gene_size = genesize,
        xspace = k1,
        yspace = k2,
        seed = seeds[i],
        domainnum = random1[i],
        use_covariate=use_covariate
      )
    }
    
    if(use_covariate){
      covar_choose="cov"
    }else{
      covar_choose="nocov"
    }
    all_results_list <- list()
    
    for (i1 in infdomain) {
      # Run the simulations in parallel using mclapply
      results <- mclapply(1:8, function(i) generate_data(i1, i), mc.cores = 8)
      
      # Extract data and calpha lists
      data_list <- lapply(results, function(spe) as.data.frame(spe[[1]]))
      if(use_covariate){
        calpha_list <- lapply(results, function(spe) as.data.frame(spe[[3]]))
      }else{
        calpha_list=NULL
      }
      # Create a list to store both data_list and calpha_list for this i1 value
      dataset_info_list <- list(
        data_list = data_list,
        calpha_list = calpha_list
      )
      all_results_list[[paste0("inf", (i1*10))]] <- dataset_info_list
      
      for (i in 1:8) {
        print(dim(data_list[[i]]))
      }
      
      print(paste0("Successfully generated datasets for i1 = ", i1))
    }
    
    
    save(all_results_list, file = paste0(curdir, "data_",k1,k2,"_",covar_choose,"_seed", seed1,".RData"))
  }
}



#Internal helper function used inside normal_create().
make_positive_definite <- function(kernel_matrix) {
  epsilon <- 1e-6  # Small diagonal increment
  while (TRUE) {
    # Check whether the covariance matrix is positive definite
    eigen_values <- eigen(kernel_matrix)$values
    if (all(eigen_values > 0)) {
      # If positive definite, return the matrix
      return(kernel_matrix)
    } else {
      # Otherwise, increase diagonal entries slightly
      kernel_matrix <- kernel_matrix + diag(epsilon, nrow(kernel_matrix))
    }
  }
}



##1.Example under ZINB framework
#Generation of Linear Spatial Pattern with Low Zero-Inflation and High Signal Strength
result=sim_create(gene_size =100,svgene_size=0.1,sv_mark=c(1,1),no_sv_mark = c(0,0),inf_size=0.1,phi=15,etamean=2,
            xspace="linear",yspace="linear",seed=1,domainnum=5,use_covariate=TRUE)

##2.Example under Gaussian process framework
result=normal_create(gene_size = 100,svgene_size = 0.1,inf_size = 0.1,kernel = "exp",seed = 1,domainnum = 1,
                     tau1 = 1,tau2 = 1,num_cores = 8,use_covariate = TRUE)

head(result[[1]])  
# Count matrix (N × G): spot-by-gene expression counts generated from the ZINB model/Gaussian process framework.

head(result[[2]])  
# Spatial location matrix (N × 2): two-dimensional coordinates for each spot.

head(result[[3]])  
# Covariate matrix (N × J, if use_covariate = TRUE): spot-level covariates 

##3.Example for reproducing simulation settings in paper: linear pattern/middle drop out/high signal
#load the packages according to the github library
library(truncnorm)
library(nlme)
library(spatstat.data)
library(spatstat.geom)
library(spatstat.random)
library(spatstat.explore)
library(spatstat.model)
library(spatstat.linnet)
library(spatstat)
library(splines)
library(mvtnorm)
library(Rcpp)
library(statmod)
library(R.utils)
library(MCMCpack)
library(Matrix)
library(MASS)
library(parallel)
library(pscl)
source(here("Simulations/code","IBaySVG_main.R"))#load the function

#2.generate the data
result1=sim_create(gene_size =5000,svgene_size=0.1,sv_mark=c(0.8,0.8),no_sv_mark = c(0,0),inf_size=0.3,phi=15,etamean=2,
                  xspace="linear",yspace="linear",seed=1,domainnum=5,use_covariate=TRUE)
result2=sim_create(gene_size =5000,svgene_size=0.1,sv_mark=c(0.5,0.5),no_sv_mark = c(0,0),inf_size=0.3,phi=15,etamean=2,
                  xspace="linear",yspace="linear",seed=2,domainnum=5,use_covariate=TRUE)
result3=sim_create(gene_size =5000,svgene_size=0.1,sv_mark=c(0.5,0.5),no_sv_mark = c(0,0),inf_size=0.3,phi=15,etamean=2,
                  xspace="linear",yspace="linear",seed=3,domainnum=5,use_covariate=TRUE)
result4=sim_create(gene_size =5000,svgene_size=0.1,sv_mark=c(0.5,0.5),no_sv_mark = c(0,0),inf_size=0.3,phi=15,etamean=2,
                  xspace="linear",yspace="linear",seed=4,domainnum=5,use_covariate=TRUE)
spelist <- list(list(result1[[1]], result1[[2]]),list(result2[[1]], result2[[2]]),
                list(result3[[1]], result3[[2]]),list(result4[[1]], result4[[2]])) 
c_alpha <- list(result1[[3]],  result2[[3]],result3[[3]],result4[[3]])# covariates
    
##3.process            
result <- IBaySVG(spelist = spelist,c_alpha = c_alpha,num_cores = 10)

#4.show the result  
#the max value of u_k
print(result[[2]])

#compute tp/fp/F1
tp=mean(result[[3]][1:500],na.rm=TRUE)
fp=mean(result[[3]][501:1000],na.rm=TRUE)
f1=2 * tp / (1 + 9 * fp + tp)
print(c(tp,fp,f1))










