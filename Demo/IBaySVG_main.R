#' Main Spatial Variable Gene Identification Function
#'
#' @param spelist A list of length m. Each element is a sublist containing:
#'   \itemize{
#'     \item{expression_matrix: Matrix (spots x genes)}
#'     \item{location_matrix: Matrix (spots x coordinates)}
#'   }
#' @param c_alpha A list of length m, each element being a covariates matrix. Default is NULL for the case of lack of cell-type covariates.
#' @param num_cores Number of CPU cores for parallelization (default: 1).
#' @param alpha : Significance level, default is 0.05.
#' @param acrate :the minimum threshold for EBLO convergence.
#' @param gamma1 :hyperparameter of the slab and spike prior of beta. Default is sqrt(0.001).
#' @param gamma2 :hyperparameter of the bi-level structure of alpha. We set gamma2 to 0.01, 0.005, and 0.0001 for more than four, three, and two or fewer slices, respectively.
#' @param min_iter : the minimum number of iteration. Default to 30.
#' @param max_iter : the maximum number of iteration. Default to 100.
#' @param speci_iter : Number of times the early convergence criterion is met.
#' @param ak_domain :Hyperparameter of varience of beta to control the slab and spike prior for different degree. Default are 0.08,0.05,0.04 and 0.03 for degree 1-4, respectively.
#' @param dp :hyperparameter of alpha's beta prior
#' @param cp :hyperparameter of alpha's beta prior
#' @param cq :hyperparameter of u's beta prior
#' @param dq :hyperparameter of u's beta prior
#' @param a_phi_single :Hyperparameter of the shape parameter $\phi$. Default is 0.001.
#' @param b_phi_single :Hyperparameter of the shape parameter $\phi$. Default is 0.001.
#' @param a_pi,b_pi Numeric. Hyperparameters for the Beta prior of the zero-inflation probability \eqn{\pi}. Defaults are 1 and 1.
#' @param eta_sigma_single :Hyperparameter of variance of $\eta$, the baseline gene expression. Default is 1.
#' @param psi_sigma_single :Hyperparameter of variance of $\psi$, controlling the contribution of covariates. Default is 1.
#'
#' @return A list with:
#'   \itemize{
#'     \item{initial_result: Posterior probability u in 2D coordinates}
#'     \item{post_mean_uk: Posterior probability u after maximization}
#'     \item{identify_number: Binary vector (1=SV gene, 0=non-SV gene)}
#'     \item{svgenename: Names of identified SV genes}
#'     \item{Posterior mean of remained parameters.}
#'   }
#' @export
#' @importFrom splines bs
#' @importFrom stats AIC as.formula na.omit rnbinom rnorm sd
#' @importFrom parallel makeCluster clusterEvalQ clusterExport parLapply stopCluster detectCores
#' @importFrom R.utils withTimeout
#' @importFrom MASS mvrnorm
#' @importFrom statmod gauss.quad.prob
#' @importFrom pscl hurdle zeroinfl
#'
#' @examples
#' # Simulate 4 datasets with different spatial patterns
#' \donttest{
#' seed <- 123
#' result1 <- sim_create(gene_size = 20, svgene_size = 0.5, sv_mark = c(0.8, 0.8), inf_size = 0.5, seed = seed)
#' result2 <- sim_create(gene_size = 20, svgene_size = 0.5, sv_mark = c(0.5, 0.5), inf_size = 0.5, seed = seed + 1)
#' result3 <- sim_create(gene_size = 20, svgene_size = 0.5, sv_mark = c(0.5, 0.5), inf_size = 0.5, seed = seed + 2)
#' result4 <- sim_create(gene_size = 20, svgene_size = 0.5, sv_mark = c(0.5, 0.5), inf_size = 0.5, seed = seed + 3)
#'
#' # Format input for IBaySVG()
#' spelist <- list(list(result1[[1]], result1[[2]]),
#'                 list(result2[[1]], result2[[2]]),
#'                 list(result3[[1]], result3[[2]]),
#'                 list(result4[[1]], result4[[2]])) #expression matrix and location
#' c_alpha <- list(result1[[3]],  result2[[3]],result3[[3]],result4[[3]])# covariates
#'# you can choose c_alpha <- NULL for the lack of cell-type covariates.
#' # Run analysis (parallel with 4 cores)
#' result=IBaySVG(spelist=spelist,c_alpha=c_alpha,num_cores=4)
#' # Print key outputs:
#' # 1. Posterior mean probabilities after maximization
#' print("Posterior mean probabilities:")
#' print(result[[2]])
#'
#' # 2. the identified SV genes
#' print("the identified SV genes:")
#' print(result[[4]])
#'
#' # 3. the posterior mean of parameters
#' names(result[[5]][[1]])
#' print("the posterior mean of parameters theta:")
#' result[[5]][[1]][["theta"]]
#' }
IBaySVG <- function(spelist,c_alpha=NULL,num_cores=1,acrate=0.01,alpha=0.05,gamma1=sqrt(0.001),gamma2=0.01,min_iter=30,max_iter=200,speci_iter=12,
                    cp=0.2,dp=1.8,cq=1,dq=1,ak_domain=c(0.08,0.05,0.04,0.03),a_phi_single=0.001,b_phi_single=0.001,
                    a_pi=1,b_pi=1,eta_sigma_single=1,psi_sigma_single=1){

  validate_data_consistency_strict(spelist, c_alpha)
  genenum <- length(colnames(spelist[[1]][[1]]))
  genename <- colnames(spelist[[1]][[1]])
  use_covariate <- !is.null(c_alpha) && length(c_alpha) > 0# logical variable
  # initial
  result_ctig <- CTIG(spelist = spelist)

  available_cores <- detectCores()-1
  cl <- makeCluster(min(num_cores,available_cores))
  clusterEvalQ(cl, {
    library(truncnorm)
    library(nlme)
    library(spatstat.data)
    library(spatstat.geom)
    library(spatstat.random)
    library(spatstat.explore)
    library(spatstat.model)
    library(spatstat.linnet)
    library(splines)
    library(mvtnorm)
    library(statmod)
    library(R.utils)
    library(MCMCpack)
    library(Matrix)
    library(MASS)
    library(pscl)
  })

  # Import necessary functions and objects
  clusterExport(cl, c("sim_create", "CTIG", "f", "gausslq", "compu_phi",
                      "safe_compu_phi", "log_compu_phi", "safe_log_compu_phi",
                      "BayFDR", "compute_diag_ABC_corrected",
                      "VIZINB", "citgtest", "c_alpha", "spelist","tun_spl"))


  # Parallel processing
  input_data <- 1:genenum
  chunk_size <- 100
  total_result <- matrix(NA, nrow = genenum, ncol = 2)
  all_parameters <- vector("list", genenum)
  progress_intervals <- seq(0.1, 1, by = 0.1)
  current_progress <- 0

  params <- list(acrate_val = acrate, gamma1_val = gamma1,gamma2_val = gamma2,
                 min_iter_val=min_iter,max_iter_val=max_iter,speci_iter_val=speci_iter,
                 a_pi_val =a_pi,b_pi_val =b_pi,cp_val=cp,dp_val=dp,cq_val=cq,dq_val=dq,ak_domain_val=ak_domain,
                 a_phi_single_val=a_phi_single,b_phi_single_val=b_phi_single,
                 eta_sigma_single_val=eta_sigma_single,psi_sigma_single_val=psi_sigma_single)


  for (i in seq(1, length(input_data), by = chunk_size)) {
    chunk_data <- input_data[i:min(i + chunk_size - 1, length(input_data))]
    result_chunk <- parLapply(cl, chunk_data, function(x, params,use_covariate) {
      acrate_val <- params$acrate_val
      gamma1_val <- params$gamma1_val
      gamma2_val<-  params$gamma2_val
      min_iter_val<-  params$min_iter_val
      max_iter_val<-  params$max_iter_val
      speci_iter_val<-  params$speci_iter_val
      a_pi_val<-params$a_pi_val
      b_pi_val<-params$b_pi_val
      cp_val <- params$cp_val
      dp_val <- params$dp_val
      cq_val <- params$cq_val
      dq_val <- params$dq_val
      ak_domain_val <- params$ak_domain_val
      a_phi_single_val <- params$a_phi_single_val
      b_phi_single_val <- params$b_phi_single_val
      eta_sigma_single_val <- params$eta_sigma_single_val
      psi_sigma_single_val <- params$psi_sigma_single_val

      result <- tryCatch({
        withTimeout({
          citgtest(x, result = result_ctig, c_alpha = c_alpha, spelist = spelist,
                   acrate = acrate_val, gamma1 = gamma1_val,gamma2 = gamma2_val,
                   min_iter=min_iter_val,max_iter=max_iter_val,speci_iter=speci_iter_val,use_covariate = use_covariate,
                   a_pi=a_pi_val,b_pi=b_pi_val,cp=cp_val,dp=dp_val,cq=cq_val,dq=dq_val,ak_domain=ak_domain_val,
                   a_phi_single=a_phi_single_val,b_phi_single=b_phi_single_val,
                   eta_sigma_single=eta_sigma_single_val,psi_sigma_single=psi_sigma_single_val)
        }, timeout = 600, onTimeout = "error")
      }, error = function(e) {
        return(list(
          "u" = c(NA, NA),
          "g" = NULL, "r" = NULL, "phi" = NULL, "sigma_beta" = NULL,
          "a_beta" = NULL, "alpha" = NULL, "q" = NULL, "p" = NULL,
          "theta" = NULL, "sigma_theta" = NULL
        ))
      })
      return(result)
    },params,use_covariate)

    for (j in seq_along(chunk_data)) {
      idx <- chunk_data[j]
      total_result[idx, ] <- result_chunk[[j]]$u   # only u (2 numbers)
      all_parameters[[idx]] <- result_chunk[[j]]   # full parameter list
    }

    completed_fraction <- i / length(input_data)
    if (completed_fraction >= progress_intervals[current_progress + 1]) {
      current_progress <- current_progress + 1
      print(paste0(current_progress * 10, "% completed"))
    }
  }

  out1=total_result
  total_result <- apply(total_result, 1, max)
  out2<-total_result
  thrs <- BayFDR(total_result, (alpha / (2*length(total_result))))
  total_result[total_result > thrs] <- 1
  total_result[total_result <= thrs] <- 0
  out3<-total_result
  svgenename <- as.vector(na.omit(genename[total_result == 1]))
  out4=svgenename
  result_NBIMSVG <- list("initial_result"=out1, "post_mean_uk"=out2,"the identify number"=out3, "svgenename"=out4,"all_parameters"  = all_parameters)
  stopCluster(cl)
  return(result_NBIMSVG)
}

#' function to produce simulation
#'
#' @param gene_size :the number of gene numbers you need
#' @param svgene_size :the proportion of SV gene among the gene_size you need
#' @param sv_mark :the coefficients of spatial effect for SV genes in 2-dimention direction you need
#' @param no_sv_mark :the coefficients of spatial effect for non-SV genes in 2-dimention direction you need.Default is 0.
#' @param inf_size :the probability of dropout
#' @param phi :the shape parameter of NB distribution
#' @param etamean :the mean value of the baseline eta.Default is 2.
#' @param xspace :the spatial pattern in x coordinate.Default situations are "linear","focal" and "period".
#' @param yspace :the spatial pattern in x coordinate.Default situations are "linear","focal" and "period".
#' @param seed :seed number
#' @param cell_dist :the alpha of the Dirichlet distribution of the coefficients of covariates. We default the covariate's number =6.
#' @param domainnum :different covariate domain structures. The value domain is from 1 to 11 and we default it to 5.
#' @param use_covariate :the logistic variable of using covariate. Default is TRUE which use covariates.
#'
#' @return a list consist of expression matrix, spot location and covariate matrix.
#' @export
#'
sim_create <- function(gene_size =100,svgene_size=0.1,sv_mark=c(1,1),no_sv_mark = c(0,0),inf_size=0.1,phi=15,etamean=2,
                       xspace="linear",yspace="linear",seed=1,cell_dist=rep(1,6),domainnum=5,use_covariate=TRUE){
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


#' implement function
#'
#' @param y ：list consist of a gene expression in m datasets
#' @param splinevalue :list consist of basic functions in m datasets with predetermined degree
#' @param samplenumber :vector consist of gene numbers in m datasets
#' @param splinelevel :predetermined degree
#' @param acrate :threshold value for determining ELBO convergence
#' @param ak :hyperparameter of varience of beta to control the slab and spike prior
#' @param dp :hyperparameter of alpha's beta prior
#' @param cp :hyperparameter of alpha's beta prior
#' @param c_alpha :list consist of coveriates in m datasets
#' @param knum :the dimension of coveriates
#' @param cq :hyperparameter of u's beta prior
#' @param dq :hyperparameter of u's beta prior
#' @param gamma1 :hyperparameter of the slab and spike prior of beta. Default is sqrt(0.001).
#' @param gamma2 :hyperparameter of the bi-level structure of alpha. We set gamma2 to 0.01, 0.005, and 0.001 to represent scenarios integrating more than four, three, and two or fewer slices, respectively.
#' @param min_iter : the minimum number of iteration.
#' @param max_iter : the maximum number of iteration.
#' @param speci_iter : Number of times the early convergence criterion is met
#' @param use_covariate :the logistic variable of using covariate. Default is TRUE which use covariates.
#' @param a_phi_single :Hyperparameter of the shape parameter $\phi$. Default is 0.001.
#' @param b_phi_single :Hyperparameter of the shape parameter $\phi$. Default is 0.001.
#' @param a_pi,b_pi Numeric. Hyperparameters for the Beta prior of the zero-inflation probability \eqn{\pi}. Defaults are 1 and 1.
#' @param eta_sigma_single :Hyperparameter of variance of $\eta$, the baseline gene expression. Default is 1.
#' @param psi_sigma_single :Hyperparameter of variance of $\psi$, controlling the contribution of covariates. Default is 1.
#'
#' @return the posterior mean of all the parameters
#' @export
#'
#'
VIZINB <- function(y,splinevalue,samplenumber,splinelevel,acrate=0.01,ak=1,dp=1.8,cp=0.2,c_alpha,knum=6,cq=1,dq=1,
                   a_phi_single=0.001,b_phi_single=0.001,a_pi=1,b_pi=1,eta_sigma_single=1,psi_sigma_single=1,
                   gamma1=sqrt(0.001),gamma2=0.01,min_iter=30,max_iter=200,speci_iter=12,use_covariate=TRUE){
  ##handle the input and initial the parameters
  # 添加参数检查
  stopifnot(length(y) == length(splinevalue))
  stopifnot(length(y) == length(samplenumber))

  M=length(samplenumber)
  G=splinelevel

  #the shape parameter
  u_phi <- rep(100,M)
  phiold <- rep(15,M)
  logphiold <- rep(1,M)

  #hyperparameter of shape parameter
  a_phi=rep(a_phi_single,M)
  b_phi=rep(b_phi_single,M)

  #dropout probability
  u_pi <- rep(list(),M)
  for(i in 1:M){
    u_pi[[i]] <- rep(0.5,samplenumber[i])
  }
  for(i in 1:M){
    for(j in 1:samplenumber[i]){
      if(y[[i]][j]==0){
        u_pi[[i]][j]=0.5
      }else{
        u_pi[[i]][j]=0
      }
    }
  }


  #the poission parameter gi
  u_g <- rep(list(),M)
  for(i in 1:M){
    u_g[[i]] <- rep(3,samplenumber[i])
  }

  #the log of poission parameter gi
  u_log_g <- rep(list(),M)
  for(i in 1:M){
    u_log_g[[i]] <- rep(2,samplenumber[i])
  }

  #the exp of poission parameter gi
  u_exp_g <- rep(list(),M)
  for(i in 1:M){
    u_exp_g[[i]] <- rep(5,samplenumber[i])
  }

  #the coefficient of the log-mean term n^m*(2L+1+J)
  c_m <- rep(list(),M)
  if(use_covariate){
    for(i in 1:M){
      c_m[[i]] <- cbind(rep(1,samplenumber[i]),splinevalue[[i]],c_alpha[[i]])
    }
  }else{
    for(i in 1:M){
      c_m[[i]] <- cbind(rep(1,samplenumber[i]),splinevalue[[i]])
    }
  }

  #the mean of integrated variables
  u_theta <- matrix(rep(0.2,(2*G+1+knum)*M),nrow=(2*G+1+knum),ncol=M)
  u_theta[1,] <- 2
  #the varience of integrated variables
  sigma_theta <- rep(list(),M)
  for(i in 1:M){
    sigma_theta[[i]] <- diag(rep(1,(2*G+1+knum)))
  }

  #the varience of beta
  u_sigmaT_beta=matrix(rep(1,2*M),nrow = M,ncol = 2)

  #hyperparameter of integrated variables
  eta_sigma = rep(eta_sigma_single,M)
  alpha_sigma=rep(psi_sigma_single,M)
  u_abetaT=matrix(rep(1,2*M),nrow = M,ncol = 2)
  A_k=rep(ak,2)

  #the indicator of each slice
  u_alpha=matrix(rep(0,2*M),nrow = M,ncol = 2)

  Gamma1=gamma1
  Gamma2=gamma2

  #the parameter of the total indicator
  u_qk <- rep(0.5,2)
  u_pk <- rep(0.5,2)
  u_uk <- rep(0.5,2)
  u_uk_old <- u_uk
  #the hyperparameter of the bi-level structure
  c_q=cq
  d_q=dq
  c_p=cp
  d_p=dp

  #the ELBO
  ELBO=0
  ELBO_OLD=1000
  signo=100
  num=0

  ##count the number when u_uk is stability
  stability_counter <- 0

  while((signo > acrate) || (num < min_iter)){

    #constant
    lg_gamma2=log(Gamma2)
    lg_1_gamma2=log(1-Gamma2)
    lg_gamma1=log(Gamma1)
    beta1 <- beta(a_pi+1, b_pi)
    beta2 <- beta(a_pi, b_pi+1)

    ##iterate u_uk
    for(k in 1:2){
      #iterate u_qk
      q1=u_uk[k]*sum(u_alpha[,k])+c_q
      q2=u_uk[k]*M+d_q-u_uk[k]*sum(u_alpha[,k])
      u_qk[k]=q1/(q1+q2)
      #iterate u_pk
      p1=c_p+u_uk[k]
      p2=1+d_p-u_uk[k]
      u_pk[k]=p1/(p1+p2)
      #iterate u_uk
      dg_q1 <- digamma(q1)
      dg_q2 <- digamma(q2)
      dg_qsum <- digamma(q1 + q2)
      dg_p1 <- digamma(p1)
      dg_p2 <- digamma(p2)
      dg_psum <- digamma(p1 + p2)
      upart1=1
      upart2=1
      for(i in 1:M){
        upart1=upart1*exp(u_alpha[i,k]*(dg_q1-dg_qsum)+(1-u_alpha[i,k])*(dg_q2-dg_qsum))
        upart2=upart2*exp(u_alpha[i,k]*lg_gamma2+(1-u_alpha[i,k])*lg_1_gamma2)
      }
      upart1=upart1*exp(dg_p1-dg_psum)
      upart2=upart2*exp(dg_p2-dg_psum)

      u_uk[k]=upart1/(upart2+upart1)
    }


    ##iterate u_theta and sigma_theta（eta,beta1,beta2) and u_alpha
    for(i in 1:M){
      #compute M_q_sigma
      one_minus_pi <- (1 - u_pi[[i]])
      M_q_1=u_alpha[i,1]*u_sigmaT_beta[i,1]+(1-u_alpha[i,1])/(Gamma1^2)
      M_q_2=u_alpha[i,2]*u_sigmaT_beta[i,2]+(1-u_alpha[i,2])/(Gamma1^2)
      if(use_covariate){
        M_q = diag(c(1/((eta_sigma[i])^2),rep(M_q_1,G),rep(M_q_2,G),rep(1/(alpha_sigma[i]^2),knum)))
      }else{
        M_q = diag(c(1/((eta_sigma[i])^2),rep(M_q_1,G),rep(M_q_2,G)))
      }
      #compute the mean term
      theta1=as.vector(exp(-c_m[[i]]%*%u_theta[,i]+unname(compute_diag_ABC_corrected(c_m[[i]],sigma_theta[[i]]))/2))
      temp1=one_minus_pi*u_g[[i]]*theta1
      D_u_theta= u_phi[i]*crossprod(c_m[[i]], temp1 - (one_minus_pi))-(M_q%*%u_theta[,i])

      #compute the varience
      vec_D_sigma_theta=u_phi[i]*crossprod(c_m[[i]]*temp1, c_m[[i]])+M_q


      #keep regularity
      rigd=1e-8
      while(abs(det(vec_D_sigma_theta))<.Machine$double.eps){
        vec_D_sigma_theta=vec_D_sigma_theta+rigd*diag(1+2*G+knum)
        rigd <- rigd*10
      }

      #inverse matrix
      sigma_theta[[i]] = solve(vec_D_sigma_theta)
      u_theta[,i] = u_theta[,i]+sigma_theta[[i]]%*%D_u_theta


      c_m_theta_i=c_m[[i]]%*%u_theta[,i]
      if(num<=2){
        theta1=as.vector(exp(-c_m_theta_i+unname(compute_diag_ABC_corrected(c_m[[i]],sigma_theta[[i]]))/2))
      }

      #iterate u_g
      g1=one_minus_pi+u_phi[i]*one_minus_pi*theta1
      g2=(y[[i]]+u_phi[i]-1)*one_minus_pi+1
      u_g[[i]] <- g2/g1
      u_log_g[[i]] <- digamma(g2)-log(g1)
      u_exp_g[[i]] <- (g1/(g1+1))^g2


      #iterate ri(m)/u_pi,
      idx <- which(y[[i]] == 0)
      if(length(idx) > 0){
        u_pi[[i]][idx] <- beta1 / (beta1 + beta2 * u_exp_g[[i]][idx])
      }

      ##iterate u_phi
      #compute c1
      c1=as.numeric(t(one_minus_pi)%*%(c_m_theta_i-u_log_g[[i]])+b_phi[i]+
                      t(one_minus_pi*u_g[[i]])%*%theta1)
      #compute N_pi
      N_pi = sum(one_minus_pi)

      #compute the posterior using numerical integration
      u_phi[i]=safe_compu_phi(n=200,s=N_pi,t=c1,ap=a_phi[i],old = u_phi[i])
      phiold[i] <- u_phi[i]

      #iterate u_sigmaT_beta M*2 matrix
      digg_sigma=unname(diag(sigma_theta[[i]]))
      theta_sample1=sum(u_theta[2:(G+1),i]^2)+sum(digg_sigma[2:(G+1)])
      sigma1=2*u_abetaT[i,1]+u_alpha[i,1]*theta_sample1
      u_sigmaT_beta[i,1]= (G*u_alpha[i,1]+1)/sigma1
      theta_sample2=sum(u_theta[(G+2):(2*G+1),i]^2)+sum(digg_sigma[(G+2):(2*G+1)])
      sigma2=2*u_abetaT[i,2]+u_alpha[i,2]*theta_sample2
      u_sigmaT_beta[i,2]= (G*u_alpha[i,2]+1)/sigma2
      sigma12=list(sigma1,sigma2)
      theta_samples <- c(theta_sample1, theta_sample2)
      u_abetaT[i,1]=1/(1/(A_k[1]^2)+u_sigmaT_beta[i,1])
      u_abetaT[i,2]=1/(1/(A_k[2]^2)+u_sigmaT_beta[i,2])

      #iterate u_alpha M*2 matrix
      for(k in 1:2){
        q1=u_uk[k]*sum(u_alpha[,k])+c_q
        q2=u_uk[k]*M+d_q-u_uk[k]*sum(u_alpha[,k])
        dg_q1 <- digamma(q1)
        dg_q2 <- digamma(q2)
        dg_qsum <- digamma(q1 + q2)
        theta_sample <- theta_samples[k]
        log_part1=-0.5*theta_sample*u_sigmaT_beta[i,k]+u_uk[k]*(dg_q1-dg_qsum)+0.5*G*(digamma(0.5*G*u_alpha[i,k]+0.5)-log(0.5*sigma12[[k]]))+(1-u_uk[k])*lg_gamma2
        log_part2=-0.5*theta_sample/(Gamma1^2)+u_uk[k]*(dg_q2-dg_qsum)+(1-u_uk[k])*lg_1_gamma2-G*lg_gamma1
        u_alpha[i,k]=1 / (1 + exp(log_part2 - log_part1))
      }

      #compute ELBO
      ##first part
      ELBO=t(one_minus_pi)%*%(-lgamma(y[[i]]+1)-u_phi[i]*c_m_theta_i)
      -0.5*sum(diag(M_q%*%(u_theta[,i]%*%t(u_theta[,i])+sigma_theta[[i]])))
      -(1-u_alpha[i,1])*G*lg_gamma1-(1-u_alpha[i,2])*G*lg_gamma1
      +(1-u_uk[1])*(u_alpha[i,1]*lg_gamma2+(1-u_alpha[i,1])*lg_1_gamma2)
      +(1-u_uk[2])*(u_alpha[i,2]*lg_gamma2+(1-u_alpha[i,2])*lg_1_gamma2)
      -b_phi[i]*u_phi[i]
      if (samplenumber[i] > 0) {
        u_pi_i <- u_pi[[i]]
        pi_terms <- sum(u_pi_i * log(beta1) - log(u_pi_i^u_pi_i) +
                          (1 - u_pi_i) * log(beta2) -log((1 - u_pi_i)^(1 - u_pi_i)))
        ELBO <- ELBO + pi_terms
      }

      ##second part
      ELBO=ELBO+sum(-u_log_g[[i]]+lgamma(g2)
                    -g2*log(one_minus_pi*(u_phi[i]*theta1+1)))
      +c1*u_phi[i]+safe_log_compu_phi(s=N_pi,t=c1,ap=a_phi[i])+0.5*log(det(sigma_theta[[i]]))

      for(k in 1:2){
        e_betak=theta_samples[k]
        u_alpha_ik=u_alpha[i,k]
        ELBO=ELBO-(u_alpha_ik*G+1)/2*log(u_abetaT[i,k]+u_alpha_ik*e_betak/2)+lgamma((u_alpha_ik*G+1)/2)
        +(u_alpha_ik*e_betak/2+u_abetaT[i,k])*u_sigmaT_beta[i,k]
        -log(u_sigmaT_beta[i,k]+1/A_k[k]^2)
        -log(u_alpha_ik^u_alpha_ik)-log((1-u_alpha_ik)^(1-u_alpha_ik))
      }

    }
    ELBO=ELBO-log(u_uk[1]^u_uk[1])-log((1-u_uk[1])^((1-u_uk[1])))-log(u_uk[2]^u_uk[2])-log((1-u_uk[2])^((1-u_uk[2])))+
      log(beta(u_uk[1]*sum(u_alpha[,1])+c_q,u_uk[1]*(M-sum(u_alpha[,1]))+d_q))+log(beta(u_uk[1]+c_p,d_p-u_uk[1]+1))+
      log(beta(u_uk[2]*sum(u_alpha[,2])+c_q,u_uk[2]*(M-sum(u_alpha[,2]))+d_q))+log(beta(u_uk[2]+c_p,d_p-u_uk[2]+1))
    signo=abs(ELBO-ELBO_OLD)
    ELBO_OLD=ELBO

    if (all(abs(u_uk - u_uk_old) < 1e-9)) {
      stability_counter <- stability_counter + 1
    } else {
      stability_counter <- 0
    }

    # updata u_uk
    u_uk_old <- u_uk

    if (stability_counter > speci_iter && signo < 1) {
      break
    }

    num=num+1
    if(num>max_iter){
      break
    }
  }
  final_result=list("u"=u_uk,
                    "g"=u_g,
                    "r"=u_pi,
                    "phi"=u_phi,
                    "sigma_beta"=u_sigmaT_beta,
                    "a_beta"=u_abetaT,
                    "alpha"=u_alpha,
                    "q"=u_qk,
                    "p"=u_pk,
                    "theta"=u_theta,
                    "sigma_theta"=sigma_theta
                    )
  return(final_result)
}




#' choose the basic function degree
#'
#' @param spelist A list of length m. Each element is a sublist containing:
#'   \itemize{
#'     \item{expression_matrix: Matrix (spots x genes)}
#'     \item{location_matrix: Matrix (spots x  coordinates)}
#'   }
#' @param calpha.list :same to c_alpha
#' @param g :gene number
#' @param use_covariate :the logistic variable of using covariate. Default is TRUE which use covariates.
#'
#' @return the number of basic function degree
#' @export
tun_spl <- function(spelist, calpha.list=NULL, g = 1, result, use_covariate=TRUE) {
  ## test length
  n1 <- length(spelist)
  n2 <- if(use_covariate) length(calpha.list) else n1
  if (n1 != n2) {
    print("warning! the length of expression matrix do not match the calpha's")
  }

  tunning_choose <- rep(1, n1)
  current_max_k <- 1  # 当前最大k值

  for (iii in c(1:n1)) {
    # 如果已经找到k=4，直接跳出循环
    if (current_max_k == 4) {
      cat(sprintf("Early termination: k=4 already found, skipping remaining datasets\n"))
      tunning_choose[iii:n1] <- 4  # 剩余数据集都设为4
      break
    }

    y <- unname(spelist[[iii]][[1]][, g])
    zero_ratio <- sum(y == 0) / length(y)
    samplenum <- length(y)
    if(use_covariate && !is.null(calpha.list)){
      celltypenum <- ncol(calpha.list[[iii]])
      calpha_names <- paste0("calpha", c(1:celltypenum))
    }else{
      celltypenum <- 0
    }

    if (zero_ratio > 0.9) {
      min_position <- 1
    } else {
      # 只计算从当前最大k到4的范围
      k_range <- c(current_max_k:4)
      aic_value <- rep(Inf, 4)  # 初始化为无穷大

      for (kkk in k_range) {
        splinelevel <- kkk
        x <- result[[2]][[kkk]][[iii]]

        if (use_covariate && !is.null(calpha.list)) {
          rows_with_sum_1 <- rowSums(calpha.list[[iii]]) == 1
          calpha.list[[iii]][rows_with_sum_1, 1] <- calpha.list[[iii]][rows_with_sum_1, 1] - 0.01
          calpha.list[[iii]][calpha.list[[iii]] < 0] <- 0
          data <- data.frame(x = x, calpha = calpha.list[[iii]], y = y)
          x_names = paste0("x", c(1:(2 * splinelevel)))
          data_name = c(x_names, calpha_names, "y")
          names(data) = data_name
          formula_str <- paste(data_name[1:(2*splinelevel+celltypenum)], collapse = " + ")
        } else {
          # 不使用协变量
          data <- data.frame(x = x, y = y)
          x_names <- paste0("x", 1:(2 * splinelevel))
          data_name <- c(x_names, "y")
          names(data) <- data_name
          formula_str <- paste(x_names, collapse = " + ")
        }

        formula_count <- paste("y ~", formula_str)
        formula_zero <- paste("|", formula_str)
        formula_full <- as.formula(paste(formula_count, formula_zero))

        tryCatch({
          if (zero_ratio < 0.7 && samplenum > 500) {
            # model with zeroinfl
            model_zinb <- zeroinfl(formula_full, data = data, dist = "negbin")
          } else {
            # model with hurdle
            model_zinb <- hurdle(formula_full, data = data, dist = "negbin")
          }
          aic_value[kkk] <- AIC(model_zinb)
        }, error = function(e) {
          aic_value[kkk] <- Inf
          cat(sprintf("Error fitting model for dataset %d, k=%d: %s\n", iii, kkk, e$message))
        })
      }

      min_position <- which.min(aic_value)
      # iterate the max
      current_max_k <- max(current_max_k, min_position)
    }

    tunning_choose[iii] <- min_position
    cat(sprintf("Dataset %d/%d: selected k = %d, current max k = %d\n", iii, n1, min_position, current_max_k))

    # 如果当前数据集选择了4，直接跳出循环
    if (min_position == 4) {
      if (iii < n1) {
        tunning_choose[(iii+1):n1] <- 4  # 剩余数据集都设为4
        cat(sprintf("Early termination: k=4 found at dataset %d, skipping remaining %d datasets\n",
                    iii, n1 - iii))
      }
      break
    }
  }

  tunning_final <- max(tunning_choose)
  print(paste0("Choosing b-spline degree: ", tunning_final))
  return(tunning_final)
}



#' data preprocess
#'
#' @param spelist A list of length m. Each element is a sublist containing:
#'   \itemize{
#'     \item{expression_matrix: Matrix (spots x genes)}
#'     \item{location_matrix: Matrix (spots x coordinates)}
#'   }
#'
#' @return: preprocessed data
#' @export

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

#' single gene process
#'
#' @param g :gene number
#' @param result : preprocessed data produced by CTIG
#' @param spelist A list of length m. Each element is a sublist containing:
#'   \itemize{
#'     \item{expression_matrix: Matrix (spots x genes)}
#'     \item{location_matrix: Matrix (spots x coordinates)}
#'   }
#' @param c_alpha A list of length m, each element being a covariates matrix.
#' @param acrate :the minimum threshold for EBLO convergence.
#' @param gamma1 :hyperparameter of the slab and spike prior of beta. Default is sqrt(0.001).
#' @param gamma2 :hyperparameter of the bi-level structure of alpha. We set gamma2 to 0.01, 0.005, and 0.001 for more than four, three, and two or fewer slices, respectively.
#' @param min_iter : the minimum number of iteration.
#' @param max_iter : the maximum number of iteration.
#' @param speci_iter : Number of times the early convergence criterion is met.
#' @param use_covariate :the logistic variable of using covariate. Default is TRUE which use covariates.
#' @param ak_domain :Hyperparameter of varience of beta to control the slab and spike prior for different degree. Default are 0.08,0.05,0.04 and 0.03 for degree 1-4, respectively.
#' @param dp :hyperparameter of alpha's beta prior
#' @param cp :hyperparameter of alpha's beta prior
#' @param cq :hyperparameter of u's beta prior
#' @param dq :hyperparameter of u's beta prior
#' @param use_covariate :the logistic variable of using covariate. Default is TRUE which use covariates.
#' @param a_phi_single :Hyperparameter of the shape parameter $\phi$. Default is 0.001.
#' @param b_phi_single :Hyperparameter of the shape parameter $\phi$. Default is 0.001.
#' @param a_pi,b_pi Numeric. Hyperparameters for the Beta prior of the zero-inflation probability \eqn{\pi}. Defaults are 1 and 1.
#' @param eta_sigma_single :Hyperparameter of variance of $\eta$, the baseline gene expression. Default is 1.
#' @param psi_sigma_single :Hyperparameter of variance of $\psi$, controlling the contribution of covariates. Default is 1.
#' @return the posterior mean of all parameters for gene g.
#' @export
#'
#'
#'
#'
#'
#'
#'
citgtest <- function(g,result,c_alpha,spelist,acrate,gamma1,gamma2,min_iter,max_iter,speci_iter,use_covariate=TRUE,
                     a_pi,b_pi,cp,dp,cq,dq,ak_domain,a_phi_single,b_phi_single,eta_sigma_single,psi_sigma_single){
  Ysingle.list<-list()
  for(i in 1:length(result[[3]])){
    Ysingle.list[[i]]<-as.vector(result[[1]][[i]][,g])
  }
  begin <- Sys.time()
  #obtain the b-spline level
  tun_choose <- tun_spl(spelist = spelist,calpha.list = c_alpha,g=g,use_covariate=use_covariate,result=result)
  if(tun_choose==1){
    ak=ak_domain[tun_choose]
  }else if(tun_choose==2){
    ak=ak_domain[tun_choose]
  }else if(tun_choose==3){
    ak=ak_domain[tun_choose]
  }else{
    ak=ak_domain[tun_choose]
  }
  if (use_covariate) {
    knum1 <- dim(c_alpha[[1]])[2]
  } else {
    knum1 <- 0
  }
  bbb <- VIZINB(Ysingle.list,result[[2]][[tun_choose]],result[[3]],tun_choose,acrate=acrate,ak=ak,dp=dp,cp=cp,cq=cq,dq=dq,c_alpha=c_alpha,knum=knum1,
                a_phi_single=a_phi_single,b_phi_single=b_phi_single,a_pi=a_pi,b_pi=b_pi,eta_sigma_single=eta_sigma_single,psi_sigma_single=psi_sigma_single,
                gamma1=gamma1,gamma2=gamma2,min_iter=min_iter,max_iter=max_iter,speci_iter=speci_iter,use_covariate=use_covariate)
  end <- Sys.time()
  return(bbb)
}


#'data Detection
#'
#' @param spelist A list of length m. Each element is a sublist containing:
#'   \itemize{
#'     \item{expression_matrix: Matrix (spots x genes)}
#'     \item{location_matrix: Matrix (spots x coordinates)}
#'   }
#' @param c_alpha A list of length m, each element being a covariates matrix.
#'
#' @return TRUE or warnings
#' @export
validate_data_consistency_strict <- function(spelist, c_alpha=NULL) {
  # 1. Gene name Detection
  ref_genes <- colnames(spelist[[1]][[1]])  # Taking the first dataset as reference

  for (m in seq_along(spelist)) {
    current_genes <- colnames(spelist[[m]][[1]])
    if (!identical(current_genes, ref_genes)) {
      stop(sprintf(
        "[Gene name Detection] dataset ",m," : number or order of genes does not match.Stop!"
      ))
    }
  }

  #2. Coveriate name Detection
  use_covariate <- !is.null(c_alpha) && length(c_alpha) > 0
  if (use_covariate) {
    ref_calpha_cols <- colnames(c_alpha[[1]])#Taking the first dataset as reference
    for (m in seq_along(c_alpha)) {
      current_cols <- colnames(c_alpha[[m]])
      if (!identical(current_cols, ref_calpha_cols)) {
        stop(sprintf(
          "[Coveriate name Detection] dataset ",m," : number or order of coveriates does not match.Stop!"
        ))
      }
    }
  }

  #3. spot name consistency Detection
  for (m in seq_along(spelist)) {
    spots_expr <- rownames(spelist[[m]][[1]])
    spots_coord <- rownames(spelist[[m]][[2]])

    if (use_covariate) {
      spots_calpha <- rownames(c_alpha[[m]])
      if (any(sapply(list(spots_expr, spots_coord, spots_calpha), is.null))) {
        stop(sprintf(
          "[Spot name consistency Detection] One dataset is missing spot name. Stop!"
        ))
      }
      if (!identical(spots_expr, spots_coord) || !identical(spots_expr, spots_calpha)) {
        stop(sprintf(
          "[Spot name consistency Detection] Number or order of spots does not match. Stop!"
        ))
      }
    } else {
      # 不使用协变量的情况，只检查表达矩阵与坐标是否匹配
      if (any(sapply(list(spots_expr, spots_coord), is.null))) {
        stop(sprintf(
          "[Spot name consistency Detection] One dataset (expr/coord) is missing spot name. Stop!"
        ))
      }
      if (!identical(spots_expr, spots_coord)) {
        stop(sprintf(
          "[Spot name consistency Detection] Number or order of spots between expr and coord does not match. Stop!"
        ))
      }
    }
  }

  return(TRUE)
}


#' function to compute complex function
#'
#' @param x :value
#' @param p :parameters of function
#' @param q :parameters of function
#' @param r :parameters of function
#' @param s :parameters of function
#' @param t :parameters of function
#' @param h :parameters of function
#'
#' @return value of the function f(x)
#' @export
f <- function(x, p, q, r, s, t,h=0) {
  result <- ((log(1+r*x))^q)*exp(p*log(x)+s*x*log(x)-s*lgamma(x)-t*x)
  return(result)
}


#' function of Numerical integration
#'
#' @param n :parameters of algorithm
#' @param p :parameters of function
#' @param q :parameters of function
#' @param r :parameters of function
#' @param s :parameters of function
#' @param t :parameters of function
#' @param h :parameters of function
#'
#' @return result of the Numerical integration
#' @export
gausslq <- function(n=100,p,q,r,s,t,h=0){
  out <- gauss.quad.prob(n,"gamma",alpha =1,beta = 1 )
  return(sum(f(out[[1]],p,q,r,s,t,h)*out[[2]]))
}

#' function of Numerical integration2
#'
#' @param n :parameters of algorithm
#' @param s :parameters of function
#' @param t :parameters of function
#' @param ap :parameters of function
#' @param max_attempts :max attempts
#'
#' @return result of the Numerical integration
#' @export
compu_phi <- function(n = 200, s, t, ap, max_attempts = 5) {
  h <- 0
  result <- NA
  attempt <- 1

  while(attempt <= max_attempts &&
        (is.na(result) || !is.finite(result) || result == 0)) {
    result <- tryCatch({
      int1 <- gausslq(n, p = ap, q = 0, r = 1, s = s, t = t - 1, h = h)
      int2 <- gausslq(n, p = ap - 1, q = 0, r = 1, s = s, t = t - 1, h = h)

      if(is.nan(int2) || int2 == 0) {
        NA
      } else {
        int1 / int2
      }
    }, error = function(e) NA)

    h <- h + 100
    attempt <- attempt + 1
  }

  if(is.na(result) || !is.finite(result)) {
    return(1)  # 返回默认值
  } else {
    return(result)
  }
}

#'function of Numerical integration3
#'
#' @param n :parameters of algorithm
#' @param s :parameters of function
#' @param t :parameters of function
#' @param ap :parameters of function
#' @param timeout :the maximum process time
#' @param old : post value of parameter
#'
#' @return result of the Numerical integration
#' @export
safe_compu_phi <- function(n=200,s,t,ap,timeout=0.05,old=1){
  result <- old
  tryCatch(
    {result <- withTimeout(compu_phi(n=200,s=s,t=t,ap=ap),timeout = timeout)},
    TimeoutException = function(ex){result <- old}
  )
  return(result)
}

#' function of Numerical integration4
#'
#' @param s :parameters of function
#' @param t :parameters of function
#' @param ap :parameters of function
#' @param max_attempts :max attempts
#'
#' @return result of the Numerical integration
#' @export
log_compu_phi <- function(s, t, ap, max_attempts = 5) {
  for (j in 0:(max_attempts - 1)) {
    h <- j * 100  # 大幅减少h增量

    result <- tryCatch({
      log_val <- log(gausslq(n = 100, p = ap - 1, q = 0, r = 1, s = s, t = t - 1, h = h))
      if (is.finite(log_val) && log_val != 0) {
        return(log_val - 100 * j)  # 对应调整修正值
      }
      NA
    }, error = function(e) NA)
  }
  return(0)
}


#' function of Numerical integration5
#'
#' @param s :parameters of function
#' @param t :parameters of function
#' @param ap :parameters of function
#' @param timeout :the maximum process time
#' @param old : post value of parameter
#'
#' @return result of the Numerical integration
#' @export
safe_log_compu_phi <- function(s,t,ap,timeout=0.05,old=0){
  result <- old
  tryCatch(
    {result <- withTimeout(log_compu_phi(s=s,t=t,ap=ap),timeout = timeout)},
    TimeoutException = function(ex){result <- old}
  )
  return(result)
}



#'function of matrix compution
#'
#' @param A :Matrix1
#' @param B :Matrix2
#'
#' @return result of multi-matrix compution
#' @export
compute_diag_ABC_corrected <- function(A, B) {
  n <- nrow(A)
  diag_elements <- numeric(n)
  for (i in 1:n) {
    diag_elements[i] <- A[i, ] %*% B %*% t(A[i, , drop=FALSE])
  }
  return(diag_elements)
}


#' function of BFDR
#'
#' @param PPI :the sequence of the result
#' @param alpha :the significance level
#'
#' @return ：the threshold determining the SV and non-SV genes
#' @export
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









