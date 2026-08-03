#####################################################################################################################
                                      #1. zeroinflation_motivation for Figure S4
#####################################################################################################################
library(ggplot2)
library(patchwork)
library(splines)

#1.function for plot
zinb_promotion <- function(counts,zeroylim=0.4) {
  n <- ncol(counts)
  gene_mean <- rowMeans(counts)
  gene_var  <- (rowSums(counts^2) - n * gene_mean^2) / (n - 1)
  gene_zero <- rowMeans(counts == 0)
  df_zero <- data.frame(gene_zero = gene_zero)
  df_mv <- data.frame(
    gene_mean = gene_mean,
    gene_var  = gene_var
  )
  p_zero <- ggplot(df_zero, aes(x = gene_zero)) +
    geom_histogram(
      breaks = seq(0, 1, by = 0.05),
      aes(y = after_stat(count / sum(count))),
      fill = "grey80",
      color = "grey35",
      linewidth = 0.35
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, zeroylim)) +
    labs(
      x = "Zero proportion",
      y = "Proportion"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 10)
    )
  p_mv <- ggplot(df_mv, aes(x = log1p(gene_mean), y = log1p(gene_var))) +
    geom_point(color = "steelblue",size = 0.6,alpha = 1) +
    geom_abline(
      intercept = 0,
      slope = 1,
      color = "red",
      linewidth = 0.6
    ) +
    coord_cartesian(xlim = c(0, 4), ylim = c(0, 4)) +
    labs(
      x = "Mean ",
      y = "Variance"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 10)
    )
  list(p_zero = p_zero, p_mv = p_mv)
}

#2.load the data
#DLPFC within donor
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
counts_samedonor <- cbind(matrix1, matrix2, matrix3, matrix4)

#DLPFC across donors
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
counts_acrossdonor<- cbind(matrix1, matrix2, matrix3, matrix4)

#SCC
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
counts_scc <- cbind(matrix1, matrix2, matrix3)

#3. Generate plots
same  <- zinb_promotion(counts_samedonor,zeroylim=0.35)
across <- zinb_promotion(counts_acrossdonor,zeroylim=0.35)
scc   <- zinb_promotion(counts_scc,zeroylim=0.5)

p1 <- same$p_zero
p2 <- same$p_mv
p3 <- across$p_zero
p4 <- across$p_mv
p5 <- scc$p_zero
p6 <- scc$p_mv

row_A <- wrap_elements(p1 | p2) +
  labs(tag = "A", title = "DLPFC-same donor") +
  theme(plot.tag = element_text(size=12, face="bold"),
        plot.tag.position = c(0.01, 0.99),
        plot.title = element_text(size=12, hjust=0.5, margin=margin(b=2)))

row_B <- wrap_elements(p3 | p4) +
  labs(tag = "B", title = "DLPFC-across donors") +
  theme(plot.tag = element_text(size=12, face="bold"),
        plot.tag.position = c(0.01, 0.99),
        plot.title = element_text(size=12, hjust=0.5, margin=margin(b=2)))

row_C <- wrap_elements(p5 | p6) +
  labs(tag = "C", title = "SCC") +
  theme(plot.tag = element_text(size=12, face="bold"),
        plot.tag.position = c(0.01, 0.99),
        plot.title = element_text(size=12, hjust=0.5, margin=margin(b=2)))

combined_plot <- row_A / row_B / row_C +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(
    title = NULL,
    theme = theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
    )
  )
combined_plot

ggsave("zeroinf_promotion.png",plot = combined_plot,width = 7,height = 9,dpi = 300)


#####################################################################################################################
                                    #2. bilevel_motivation for Figures S4-S8
#####################################################################################################################
#function for plot
pattern_plot3 <- function(pltdat, igene, xy = T, main = F, titlesize = 2,
                          pointsize = 3, xpand = 0, ypand = 1, title = NULL,
                          color_limits = NULL, show_legend = FALSE) {   # 新增两个参数
  if (!xy) {
    xy <- matrix(as.numeric(do.call(rbind, strsplit(as.character(pltdat[,1]), split="x"))), ncol=2)
    rownames(xy) <- as.character(pltdat[,1]); colnames(xy) <- c("x","y")
    pd <- cbind.data.frame(xy, pltdat[, 2:ncol(pltdat)])
  } else pd <- pltdat
  
  pal <- colorRampPalette(c("mediumseagreen", "lightyellow2", "deeppink"))
  gpt <- ggplot(pd, aes(x=x, y=y, color=pd[,igene+2])) + geom_point(size=pointsize) +
    scale_color_gradientn(colours=pal(5), limits=color_limits) +   # ← 加 limits
    scale_x_discrete(expand=c(xpand,ypand)) + scale_y_discrete(expand=c(xpand,ypand)) +
    coord_equal() + theme_void()
  
  if (main) {
    if (is.null(title)) title = colnames(pd)[igene+2]
    out = gpt + labs(title=title, x=NULL, y=NULL, color=NULL) +
      theme(plot.title=element_text(hjust=0.5, size=rel(titlesize)),
            legend.position = if(show_legend) "right" else "none")
  } else {
    out = gpt + labs(title=NULL, x=NULL, y=NULL, color=NULL) +
      theme(legend.position = if(show_legend) "right" else "none")
  }
  return(out)
}

meta_process =function(position1,count1){
  meta=cbind(position1,t(count1))
  colnames(meta)[c(1:2)] = c('x','y')
  meta=as.data.frame(meta)
  meta$x <- (meta$x - min(meta$x))/(max(meta$x) - min(meta$x))
  meta$y <- (meta$y - min(meta$y))/(max(meta$y) - min(meta$y))
  return(meta)
}

simulate_pattern_expr <- function(pattern, m, spelist, coef = 0.5,intercept = 0, phi = 1, pi0 = 0.3, seed = 1) {
  set.seed(seed)
  # 生成该模式的基（调 CTIG_modified 对应 pattern）
  ctig_p <- CTIG_modified1(spelist, pattern = pattern)
  basis <- ctig_p[[2]][[5]][[m]]           # n×2，该模式的基（两个方向）
  
  # 空间效应 = 基 × 系数（两方向都用 coef）
  spatial <- as.numeric(basis %*% rep(coef, ncol(basis)))
  mu <- exp(intercept + spatial)           # NB 均值
  
  n <- length(mu)
  is0 <- rbinom(n, 1, pi0)                  # 零膨胀
  ysim <- ifelse(is0 == 1, 0, rnbinom(n, mu = mu, size = phi))
  ysim
}

CTIG_modified1<-function(spelist,pattern="linear"){
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
      co_spline_list[[i]]<-cbind(coord.list[[i]][,1],(-coord.list[[i]][,2]))
    }else if(pattern=="focal"){
      co_spline_list[[i]]<-exp(-coord.list[[i]]^2/2)
    }else if(pattern=="period"){
      co_spline_list[[i]]<-cos(3*pi*coord.list[[i]]/2)
    }
  }
  
  co_spline_list_list[[5]] <- co_spline_list
  
  
  result <- list(Y.list,co_spline_list_list,samplenumber,G)
  return(result)
}


###########################################Figure 4###################################################
#1.load the data
dir="data/Realdataset/dlpfc samedonor" #You may switch to the across-donor dataset
matrix4=as.matrix(read.csv(here::here(dir,"matrix4_count_samedonor.csv"),row.names = 1,check.names = FALSE))
matrix3=as.matrix(read.csv(here::here(dir,"matrix3_count_samedonor.csv"),row.names = 1,check.names = FALSE))
matrix2=as.matrix(read.csv(here::here(dir,"matrix2_count_samedonor.csv"),row.names = 1,check.names = FALSE))
matrix1=as.matrix(read.csv(here::here(dir,"matrix1_count_samedonor.csv"),row.names = 1,check.names = FALSE))
position4=as.matrix(read.csv(here::here(dir,"matrix4_position_samedonor.csv"),row.names = 1))
position3=as.matrix(read.csv(here::here(dir,"matrix3_position_samedonor.csv"),row.names = 1))
position2=as.matrix(read.csv(here::here(dir,"matrix2_position_samedonor.csv"),row.names = 1))
position1=as.matrix(read.csv(here::here(dir,"matrix1_position_samedonor.csv"),row.names = 1))
spelist_samedonor <- list(list(t(matrix1),position1),list(t(matrix2),position2),
                          list(t(matrix3),position3),list(t(matrix4),position4))

#2.choose the determined gene for "linear","focal","period". 
rep_genes <- c(linear="CNP", focal="COX6C", period="PRNP") 
patterns <- c("linear","focal","period")
mlist <- c(4,3,3)
plots_row <- list()

# 3.plot the component plot for each row 
# row1: the true expression of the gene
for (i in seq_along(patterns)) {
  gname <- rep_genes[patterns[i]]
  g <- which(colnames(spelist_samedonor[[1]][[1]]) == gname)
  ytrue <- as.numeric(spelist_samedonor[[mlist[i]]][[1]][, g])
  lim <- range(ytrue)
  count_true <- matrix(ytrue, nrow=1, dimnames=list(gname, NULL))
  meta_true <- meta_process(spelist_samedonor[[mlist[i]]][[2]], count_true)
  p <- pattern_plot3(meta_true, 1, xpand=-2, ypand=2,
                     main = TRUE, title = gname, titlesize = 1,color_limits = lim)
  p <- p + theme(plot.title = element_text(hjust=0.5, size=rel(1),
                                           margin = margin(b = 8),   
                                           vjust = 1))
  plots_row[[i+6]] <- p        
}

# row2：sampled expression from the simulated model
coef_list <- c(linear = 1, focal = 6, period = 3)   
intercept_list <- c(linear = 0, focal = 2, period = 0)
for (i in seq_along(patterns)) {
  ysim <- simulate_pattern_expr(patterns[i], mlist[i], spelist_samedonor,
                                coef = coef_list[i],
                                intercept = intercept_list[i],
                                phi = 30, pi0 = 0.01)
  count_sim <- matrix(ysim, nrow=1, dimnames=list(patterns[i], NULL))
  meta_sim <- meta_process(spelist_samedonor[[mlist[i]]][[2]], count_sim)
  plots_row[[i+3]] <- pattern_plot3(meta_sim, 1, xpand=-2, ypand=2, main=FALSE)
}

# row3:the curve for the expression
s <- seq(-2.5, 2.5, length = 300)
func_val <- list(linear = s, focal = exp(-s^2/2), period = cos(2*pi*s))
titles <- list(
  linear = expression(paste("Linear:  ", f(s) == s)),
  focal  = expression(paste("Focal:  ", f(s) == e^{-s^2/2})),
  period = expression(paste("Period:  ", f(s) == cos(2*pi*s))))

for (i in seq_along(patterns)) {
  df <- data.frame(s = s, value = func_val[[patterns[i]]])
  plots_row[[i]] <- ggplot(df, aes(s, value)) +
    geom_line(linewidth=1, color="steelblue") +
    geom_hline(yintercept=0, linetype=3, color="grey60") +
    labs(title = titles[[patterns[i]]],
         y = if(i==1) "Function Value" else NULL) +
    theme_minimal() +
    theme(plot.title=element_text(hjust=0.5, size=11),
          axis.text.y=element_blank(), axis.ticks.y=element_blank())
}

plots_row[[3]] <- plots_row[[3]] + theme(
  legend.position = c(0.98, 0.02),      # 右上角（x=0.98, y=0.98）
  legend.justification = c(1, 1),        # 图例右上角对齐到该点
  legend.background = element_rect(fill=alpha("white",0.7), color=NA),
  legend.title = element_blank())
plots_row[[1]] <- plots_row[[1]] + theme(legend.position = "none")
plots_row[[2]] <- plots_row[[2]] + theme(legend.position = "none")


row_A <- (wrap_elements(plots_row[[7]]) | wrap_elements(plots_row[[8]]) | wrap_elements(plots_row[[9]]))
row_A <- wrap_elements(row_A) + 
  labs(tag = "A", title = "Observed expression") +
  theme(plot.tag = element_text(size=12, face="bold"),
        plot.tag.position = c(0.01, 0.99),
        plot.title = element_text(size=12, hjust=0.5, margin=margin(b=2)))

row_B <- (wrap_elements(plots_row[[4]]) | wrap_elements(plots_row[[5]]) | wrap_elements(plots_row[[6]]))
row_B <- wrap_elements(row_B) + 
  labs(tag = "B", title = "Model-simulated expression") +
  theme(plot.tag = element_text(size=12, face="bold"),
        plot.tag.position = c(0.01, 0.99),
        plot.title = element_text(size=12, hjust=0.5, margin=margin(b=2)))

row_C <- (wrap_elements(plots_row[[1]]) | wrap_elements(plots_row[[2]]) | wrap_elements(plots_row[[3]]))
row_C <- wrap_elements(row_C) + 
  labs(tag = "C", title = "Fitted spatial-effect curves") +
  theme(plot.tag = element_text(size=12, face="bold"),
        plot.tag.position = c(0.01, 0.99),
        plot.title = element_text(size=12, hjust=0.5, margin=margin(b=2)))


combined_plot <- row_A / row_B / row_C + plot_layout(heights = c(1.6, 1.4, 1.55))

ggsave("C:/Users/zhoum/Desktop/spatial_baisc.png",plot = combined_plot, width = 10, height = 11, dpi = 300)

###################################Figures S5, S6 and S7##########################################################
#function
plot_nonpara_function=function(spelist_samedonor,calpha_samedonor,choose_genes,seedlist,mlist,dataname,result){
  result_ctig <- CTIG_modified(spelist = spelist_samedonor,pattern = "linear")
  gene_names_all <- colnames(spelist_samedonor[[1]][[1]])
  gene_info <- lapply(choose_genes, function(gn) {
    g <- match(gn, gene_names_all)
    
    if (is.na(g)) {
      stop("Gene not found in expression matrix: ", gn)
    }
    
    if (is.null(result$all_parameters[[gn]])) {
      stop("Parameters not found in reduced result object: ", gn)
    }
    
    list(
      name  = gn,
      g     = g,
      par_g = result$all_parameters[[gn]]
    )
  })
  
  plots <- list()
  for (i in seq_along(gene_info)) {
    gi <- gene_info[[i]]
    p <- spatial_curve_gg(gi$par_g, gi$g, mlist[i], spelist_samedonor, calpha_samedonor,
                          result_ctig, title = NULL,titlesize = 12)
    p <- p + theme(axis.text = element_blank(), axis.ticks = element_blank())
    if (i != 1) {
      p <- p + labs(y = NULL)          # 中间、右边去掉 y 轴标签 "spatial effect"
    }
    # x 轴标签 "spatial coordinate" 三个都留也行，或只留中间：
    if (i != 2) {
      p <- p + labs(x = NULL)     
    }
    plots[[i]] <- p
  }
  
  ## —— 行B：抽样表达；行C：真实表达；同基因统一色标 ——
  sim_list <- true_list <- vector("list", length(gene_info))
  for (i in seq_along(gene_info)) {
    gi <- gene_info[[i]]
    sim_list[[i]]  <- simulate_count_row(gi$par_g, gi$g, mlist[i], spelist_samedonor,
                                         calpha_samedonor, result_ctig,seed=seedlist[i])
    true_list[[i]] <- as.numeric(spelist_samedonor[[mlist[i]]][[1]][, gi$g])
  }
  
  for (i in seq_along(gene_info)) {                       # 行B：抽样
    gi <- gene_info[[i]]
    lim <- range(c(sim_list[[i]], true_list[[i]]))        # 该基因抽样+真实的公共色标
    count_sim <- matrix(sim_list[[i]], nrow=1, dimnames=list(gi$name, NULL))
    meta_sim  <- meta_process(spelist_samedonor[[mlist[i]]][[2]], count_sim)
    plots[[3 + i]] <- pattern_plot3(meta_sim, 1, xpand=-2, ypand=2, main=FALSE,
                                    color_limits = lim)
  }
  for (i in seq_along(gene_info)) {                       # 行C：真实
    gi <- gene_info[[i]]
    lim <- range(c(sim_list[[i]], true_list[[i]]))
    count_true <- matrix(true_list[[i]], nrow=1, dimnames=list(gi$name, NULL))
    meta_true  <- meta_process(spelist_samedonor[[mlist[i]]][[2]], count_true)
    plots[[6 + i]] <- pattern_plot3(meta_true, 1, xpand=-2, ypand=2, main=TRUE,title=gi$name,titlesize = 1,
                                    color_limits = NULL)+theme(plot.title = element_text(hjust=0.5, size=rel(1),
                                                                                         margin = margin(b = 8),   # 标题和图之间加间距
                                                                                         vjust = 1))
  }
  
  plots[[3]] <- plots[[3]] + theme(
    legend.position = c(0.98, 0.02),      # 右上角（x=0.98, y=0.98）
    legend.justification = c(1, 1),        # 图例右上角对齐到该点
    legend.background = element_rect(fill=alpha("white",0.7), color=NA),
    legend.title = element_blank())
  plots[[1]] <- plots[[1]] + theme(legend.position = "none")
  plots[[2]] <- plots[[2]] + theme(legend.position = "none")
  
  row_C <- (plots[[1]] | plots[[2]] | plots[[3]])
  row_C <- wrap_elements(row_C) + 
    labs(tag = "C", title = "Fitted spatial-effect curves") +
    theme(plot.tag = element_text(size=12, face="bold"),
          plot.tag.position = c(0.01, 0.99),
          plot.title = element_text(size=12, hjust=0.5, margin=margin(b=2)))
  
  row_B <- (wrap_elements(plots[[4]]) | wrap_elements(plots[[5]]) | wrap_elements(plots[[6]]))
  row_B <- wrap_elements(row_B) + 
    labs(tag = "B", title = "Model-simulated expression") +
    theme(plot.tag = element_text(size=12, face="bold"),
          plot.tag.position = c(0.01, 0.99),
          plot.title = element_text(size=12, hjust=0.5, margin=margin(b=2)))
  
  row_A <- (wrap_elements(plots[[7]]) | wrap_elements(plots[[8]]) | wrap_elements(plots[[9]]))
  row_A <- wrap_elements(row_A) + 
    labs(tag = "A", title = "Observed expression") +
    theme(plot.tag = element_text(size=12, face="bold"),
          plot.tag.position = c(0.01, 0.99),
          plot.title = element_text(size=12, hjust=0.5, margin=margin(b=2)))
  
  combined_plot <- row_A / row_B / row_C + plot_layout(heights = c(1.6, 1.4, 1.4))
  ggsave(paste0("C:/Users/zhoum/Desktop/nonparameter_",dataname,".png"),
         plot = combined_plot, width = 10, height = 11, dpi = 300)
  return(combined_plot)
}

spatial_curve_gg <- function(par_g, g, m, spelist, c_alpha, result_ctig,mode = "loess", span = 0.1,
                             infer_degree_fn = infer_degree_from_theta,title = NULL,titlesize=12) {
  J <- infer_degree_fn(par_g$theta, c_alpha)
  spline_m <- result_ctig[[2]][[J]][[m]]
  theta_m  <- par_g$theta[, m]
  beta  <- theta_m[2:(1 + 2*J)]
  beta1 <- beta[1:J]; beta2 <- beta[(J+1):(2*J)]
  X1 <- spline_m[, 1:J, drop=FALSE]; X2 <- spline_m[, (J+1):(2*J), drop=FALSE]
  s1 <- spelist[[m]][[2]][,1]; s2 <- spelist[[m]][[2]][,2]
  f1 <- as.numeric(X1 %*% matrix(beta1, ncol=1))
  f2 <- as.numeric(X2 %*% matrix(beta2, ncol=1))
  
  if (mode == "loess") {
    x1g <- seq(min(s1),max(s1),length=200); x2g <- seq(min(s2),max(s2),length=200)
    df <- rbind(
      data.frame(x=x1g, y=predict(loess(f1~s1, span=span), x1g), dir="b1(s1)"),
      data.frame(x=x2g, y=predict(loess(f2~s2, span=span), x2g), dir="b2(s2)"))
  } else {
    df <- rbind(
      data.frame(x=s1, y=f1, dir="b1(s1)"),
      data.frame(x=s2, y=f2, dir="b2(s2)"))
  }
  
  ggplot(df, aes(x, y, color=dir)) +
    geom_line(linewidth=1) +
    geom_hline(yintercept=0, linetype=3, color="grey60") +
    scale_color_manual(values=c("b1(s1)"="steelblue","b2(s2)"="darkorange")) +
    labs(title = title,
         x="spatial coordinate", y="spatial effect", color=NULL) +
    theme_minimal() +
    theme(legend.position="top", plot.title=element_text(hjust=0.5, size=titlesize))
}

infer_degree_from_theta <- function(theta_mat, c_alpha = NULL) {
  p_theta <- nrow(theta_mat)
  
  q <- if (!is.null(c_alpha) && length(c_alpha) > 0) {
    ncol(c_alpha[[1]])
  } else {
    0
  }
  
  k <- (p_theta - 1 - q) / 2
  
  if (abs(k - round(k)) > 1e-8) {
    stop("Cannot infer spline degree from theta dimension.")
  }
  
  as.integer(round(k))
}

simulate_count_row <- function(par_g, g, m, spelist, c_alpha, result_ctig,
                               infer_degree_fn = infer_degree_from_theta, seed=1) {
  set.seed(seed)
  J <- infer_degree_fn(par_g$theta, c_alpha)
  spline_m <- result_ctig[[2]][[J]][[m]]
  c_m <- cbind(1, spline_m, c_alpha[[m]])
  theta_m <- par_g$theta[, m]
  #theta_m[c(2:(1+2*J))]=4*theta_m[c(2:(1+2*J))]
  theta_m=2*theta_m
  mu  <- as.numeric(exp(c_m %*% theta_m))
  pi  <- as.numeric(par_g$r[[m]])
  phi <- par_g$phi[m]
  n <- length(mu)
  is0 <- rbinom(n, 1, pi)
  ysim <- ifelse(is0==1, 0, rnbinom(n, mu=mu, size=phi))
  ysim   # 长度 n 的抽样计数
}

#1.load the data
#DLPFC within donor
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

#2.plot for each datasets based on the chosen genes

load(here::here("RealData/result_data/exploratory_analysis","part_dlpfcsamedonor_result.RData"))
plot_nonpara_function(spelist_samedonor,calpha_samedonor,choose_genes= c("YWHAH","SYT1","SPARCL1"),
                      seedlist=c(13,5,5),mlist=c(3,2,4),dataname="samedonor",result=result_small)

load(here::here("RealData/result_data/exploratory_analysis","part_dlpfcacrossdonor_result.RData"))
plot_nonpara_function(spelist_acrossdonor,calpha_acrossdonor,choose_genes= c("HOPX","THY1","PLP1"),
                      seedlist=c(4,5,3),mlist=c(4,3,4),dataname="acrossdonor",result=result_small)

load(here::here("RealData/result_data/exploratory_analysis","part_scc_result.RData"))
plot_nonpara_function(spelist_scc,calpha_scc,choose_genes= c("KRT10","SPRR1B","IFI27"),
                      seedlist=c(2,2,1),mlist=c(2,1,2),dataname="scc",result=result_small)



###############################################################################################################################
                                          #3.bilevel structure for Figure S9
##############################################################################################################################
#function
plot_gene_slices <- function(gene, spelist, high_q = 0.85, gene_title = TRUE) {
  M <- length(spelist)
  gene_names <- colnames(spelist[[1]][[1]])
  
  if (is.character(gene)) {
    g <- which(gene_names == gene)
    if (length(g) == 0) stop("基因 '", gene, "' 不在该数据集中。相似的: ",
                             paste(grep(substr(gene,1,3), gene_names, value=TRUE)[1:5], collapse=", "))
  } else g <- gene
  
  cells <- lapply(1:M, function(m) {
    y <- as.numeric(spelist[[m]][[1]][, g])
    pos <- spelist[[m]][[2]]
    df <- data.frame(x=pos[,1], y=pos[,2], expr=y)
    thr <- quantile(y[y>0], high_q); df$high <- df$expr >= thr
    ggplot(df, aes(x,y)) +
      geom_point(aes(color=expr), size=0.8) +
      scale_color_gradientn(colours=colorRampPalette(
        c("mediumseagreen","lightyellow2","deeppink"))(5)) +
      geom_density_2d(data=subset(df,high), aes(x,y),
                      color="red", linewidth=0.4, bins=3) +
      coord_equal() + theme_void() + theme(legend.position="none")
  })
  
  p <- wrap_plots(cells, ncol = 2)      # 固定 2 列：M=4→2×2，M=3→2×2(第4格空)
  if (gene_title)
    p <- p + plot_annotation(title=gene,
                             theme=theme(plot.title=element_text(hjust=0.5, size=12, face="bold")))
  p
}

#1.plot based on the datasets in section  “2.bilevel_motivation for Figures S4-S8”
plot1=plot_gene_slices("CNP", spelist_samedonor)
plot2=plot_gene_slices("MAPK10", spelist_samedonor)
plot3=plot_gene_slices("MGP", spelist_samedonor)
plot4=plot_gene_slices("PCP4", spelist_acrossdonor)
plot5=plot_gene_slices("RNASE1", spelist_acrossdonor)
plot6=plot_gene_slices("MGP", spelist_acrossdonor)
plot7=plot_gene_slices("LGALS7", spelist_scc)
plot8=plot_gene_slices("CSTA", spelist_scc)
plot9=plot_gene_slices("CALML5", spelist_scc)

row_A <- (wrap_elements(plot1) | wrap_elements(plot2) | wrap_elements(plot3))
row_A <- wrap_elements(row_A) + 
  labs(tag = "A", title = NULL) +
  theme(plot.tag = element_text(size=15, face="bold"),
        plot.tag.position = c(0.01, 0.99),
        plot.title = element_text(size=12, hjust=0.5, margin=margin(b=2)))

row_B <- (wrap_elements(plot4) | wrap_elements(plot5) | wrap_elements(plot6))
row_B <- wrap_elements(row_B) + 
  labs(tag = "B", title = NULL) +
  theme(plot.tag = element_text(size=15, face="bold"),
        plot.tag.position = c(0.01, 0.99),
        plot.title = element_text(size=12, hjust=0.5, margin=margin(b=2)))

row_C <- (wrap_elements(plot7) | wrap_elements(plot8) | wrap_elements(plot9))
row_C <- wrap_elements(row_C) + 
  labs(tag = "C", title = NULL) +
  theme(plot.tag = element_text(size=15, face="bold"),
        plot.tag.position = c(0.01, 0.99),
        plot.title = element_text(size=12, hjust=0.5, margin=margin(b=2)))

combined_plot <- row_A / row_B / row_C + plot_layout(heights = c(1.5, 1.5, 1.5))

ggsave("bilevel_promotion.png",plot = combined_plot, width = 12, height = 13, dpi = 400)






