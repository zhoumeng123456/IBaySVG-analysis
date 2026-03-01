############################################################################################################################################
                                                  #time cost for tables4

############################################################################################################################################
#function for produce the experiment of time cost
time_compute_create <- function(spot_num=32,gene_size =1000,svgene_size=0.1,sv_mark=c(1,1),no_sv_mark = c(0,0),inf_size=0.1,phi=15,etamean=2,
                                xspace="linear",yspace="linear",seed=1,cell_dist=rep(1,6),domainnum=5,use_covariate=TRUE){
  ###
  ###
  set.seed(seed)
  x_coords <- rep(0:(spot_num-1), each = spot_num)
  y_coords <- rep(0:(spot_num-1), times = spot_num)
  z_type=rep(1,(spot_num*spot_num))
  
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

#example for produce 
seed1=1
spot_num=32
inf_size=0.5
result1=time_compute_create(spot_num=spot_num,gene_size = 1000,sv_mark = c(0.8,0.8),inf_size = inf_size,seed = (seed1*4+1),use_covariate=FALSE)
result2=time_compute_create(spot_num=spot_num,gene_size = 1000,sv_mark = c(0.5,0.5),inf_size = inf_size,seed = (seed1*4+1),use_covariate=FALSE)
result3=time_compute_create(spot_num=spot_num,gene_size = 1000,sv_mark = c(0.5,0.5),inf_size = inf_size,seed = (seed1*4+2),use_covariate=FALSE)
result4=time_compute_create(spot_num=spot_num,gene_size = 1000,sv_mark = c(0.5,0.5),inf_size = inf_size,seed = (seed1*4+3),use_covariate=FALSE)

matrix1=as.matrix(result1[[1]])
matrix2=as.matrix(result2[[1]])
matrix3=as.matrix(result3[[1]])
matrix4=as.matrix(result4[[1]])

position1=as.matrix(result1[[2]])
position2=as.matrix(result2[[2]])
position3=as.matrix(result3[[2]])
position4=as.matrix(result4[[2]])

##运行函数
spelist <- list(list(matrix1,position1),list(matrix2,position2),list(matrix3,position3),list(matrix4,position4))
c_alpha=NULL
begin=Sys.time()
result <- NBIMSVG(spelist = spelist,c_alpha = c_alpha,num_cores =11,max_iter=200)
end=Sys.time()
time_value=as.numeric((end-begin), units = "mins")
result=list(result[[1]],result[[2]],result[[3]],result[[4]],time_value)
print(time_value)


##load the result
# This matrix records the computational time of IBaySVG for analyzing 1,000 genes.
# Columns correspond to different zero-inflation rates (0.1, 0.3, 0.5, 0.7, and 0.9).
# Rows represent increasing numbers of spatial sample points (from low to high).
TIME_MATRIX=read.csv(here("RealData/result_data","compute time.csv"))#

print(TIME_MATRIX/1000)
