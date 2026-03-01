library(dynamicTreeCut)
library(mclust)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

##1.funtion for gene cluster
svg.clust <- function(expression, coor, svg, deepsplit=1,hcl_med="ward.D",dis_method="manhattan",method = 'h', num_cluster = NULL, plot = FALSE) {
  # Check input validity
  if (length(svg) == 0) {
    stop("svg must contain at least one gene")
  }
  # Extract coordinates and scale the data
  ct <- as.data.frame(scale(expression, center = TRUE, scale = TRUE))
  counts <- subset(ct, select = svg)
  
  if (method == 'h') {
    # Calculate distance matrix
    co <- dist(t(counts),method = dis_method)
    
    # Perform hierarchical clustering
    hh <- hclust(co, method = hcl_med)
    
    # Determine clusters using Dynamic Tree Cut if num_cluster is not provided
    if (is.null(num_cluster)) {
      c6 <- cutreeDynamic(dendro = hh, distM = as.matrix(co), method = "hybrid", deepSplit = deepsplit)
    } else {  # Use specified number of clusters
      c6 <- cutree(hh, k = num_cluster)
    }
    
    # Optionally plot the dendrogram
    if (plot) {
      plot(hh, main = "Dendrogram of Hierarchical Clustering", xlab = "Genes", ylab = "Height")
      rect.hclust(hh, k = length(unique(c6)), border = "red")
    }
    
    svg.cl <- data.frame(gene = hh$labels, cl = c6)
  }
  
  if (method == 'k') {
    if (is.null(num_cluster)) {  # Number of clusters must be provided
      stop('Missing the parameter num_cluster!')
    } else {
      km <- kmeans(t(counts), centers = num_cluster)  # kmeans defaults to clustering rows, so transpose
      svg.cl <- data.frame(gene = names(km$cluster), cl = km$cluster)
    }
  }
  
  colnames(svg.cl) <- c('gene', 'cluster')
  return(svg.cl)
}

#'
#' @param svgenename A character vector of spatially variable gene names 
#'        used for clustering.
#'
#' @param data_dir A character string specifying the directory 
#'        containing count and position matrices.
#'
#' @param dis_method Character string specifying the distance metric used 
#'        in hierarchical clustering. Common options include:
#'        "euclidean" (default), "manhattan", "maximum", 
#'        "canberra", "binary", and "minkowski".
#'
#' @param deepsplit Integer controlling the sensitivity of dynamic tree 
#'        cutting.Default is 2.
#'
#' @param hcl_med Character string specifying the agglomeration (linkage) 
#'        method used in hierarchical clustering. Common options include:
#'        "average" (default), "complete", "single", 
#'        "ward.D", "ward.D2", "mcquitty", "median", and "centroid".
#'
#' @return A list of gene cluster assignments for each section.
run_svg_clustering <- function(svgenename,data_dir,sectnumber=4,dis_method = "euclidean",deepsplit_domain = 2,hcl_med = "average",dataset="samedonor"){
  result_gene_cluster_list <- list()
  for(sectnum in 1:sectnumber){
    count1 <- as.matrix(read.csv(here(data_dir, paste0("matrix", sectnum, "_count_",dataset,".csv")),row.names = 1))
    count1 <- t(count1)
    
    # ---- Read position matrix ----
    position1 <- as.matrix(read.csv(here(data_dir, paste0("matrix", sectnum, "_position_",dataset,".csv")),row.names = 1))
    
    # ---- SVG clustering ----
    result_svg_clust <- svg.clust(count1,position1,svgenename,dis_method = dis_method,deepsplit = deepsplit_domain[sectnum],hcl_med = hcl_med)
    expression_df <- as.data.frame(count1)
    gene_clusters <- result_svg_clust
    clusters <- sort(unique(gene_clusters$cluster))
    
    # ---- Compute cluster mean expression ----
    cluster_means <- data.frame(rowname = rownames(expression_df))
    for (cl in clusters) {
      genes_in_cluster <- gene_clusters$gene[gene_clusters$cluster == cl]
      mean_expression <- rowMeans(expression_df[, genes_in_cluster, drop = FALSE],na.rm = TRUE)
      cluster_means[[paste0("Cluster", cl)]] <- mean_expression
    }
    meta <- meta_process(position1, t(cluster_means[, -1]))
    result_gene_cluster_list[[sectnum]] <- gene_clusters
  }
  
  return(result_gene_cluster_list)
}


##2.Run spatial gene clustering for multiple tissue sections
#Example:for dlpfc dataset with same donor.
load(here("RealData/result_data/realdata svgene","data_dlpfc_samedonor_svgene_list.RData"))
svgene_list=data_dlpfc_samedonor_svgene_list
result_gene_cluster_list <- run_svg_clustering(svgenename = svgene_list[[6]],
                                               data_dir = "RealData/result_data/realdata dataset/dlpfc samedonor",
                                               dataset="samedonor",
                                               sectnumber=4,
                                               deepsplit_domain =c(2,2,2,2),
                                               dis_method = "euclidean",
                                               hcl_med = "average"
)#The same results can be obtained with minor modifications.

###########################################################################################################################################
                                            #3.plot for gene cluster expression :(S14,S16,S18)

############################################################################################################################################
#function for plot
f_genecluster=function(result_gene_cluster_list,count_list,position_list,secnum=1,cluster_num=1,main=F,title = NULL,titlesize = 2,hjust=0.5,legend.position = "none"){
  expre_mat=t(count_list[[secnum]])
  position_mat=position_list[[secnum]]
  gene_clusters=result_gene_cluster_list[[secnum]]
  genes_in_cluster <- gene_clusters$gene[gene_clusters$cluster == cluster_num]
  
  mean_expression <- rowMeans(expre_mat[, genes_in_cluster, drop = FALSE], na.rm = TRUE)
  meta=meta_process(position_mat,t(mean_expression))
  plot=pattern_plot2(meta,1,xpand=0.01,ypand = 0.01,pointsize = 2,main=main,title = title,titlesize = titlesize,hjust=hjust,legend.position = legend.position)
  return(plot)
}

meta_process =function(position1,count1){
  meta=cbind(position1,t(count1))
  colnames(meta)[c(1:2)] = c('x','y')
  meta=as.data.frame(meta)
  meta$x <- (meta$x - min(meta$x))/(max(meta$x) - min(meta$x))
  meta$y <- (meta$y - min(meta$y))/(max(meta$y) - min(meta$y))
  return(meta)
}

pattern_plot2 <- function(pltdat, igene, xy = T, main = F, titlesize = 2, 
                          pointsize = 3, xpand = 0, ypand = 1, title = NULL,hjust=0.5,legend.position = "none") {
  if (!xy) {
    xy <- matrix(as.numeric(do.call(rbind, strsplit(as.character(pltdat[, 
                                                                        1]), split = "x"))), ncol = 2)
    rownames(xy) <- as.character(pltdat[, 1])
    colnames(xy) <- c("x", "y")
    pd <- cbind.data.frame(xy, pltdat[, 2:ncol(pltdat)])
  } else {
    pd <- pltdat
  }
  
  # pal <- colorRampPalette(c('seagreen1','orangered')) pal <-
  # colorRampPalette(c('#00274c','#ffcb05')) pal <-
  # colorRampPalette(c('deepskyblue','goldenrod1')) pal <-
  # colorRampPalette(c('deepskyblue','deeppink'))
  pal <- colorRampPalette(c("mediumseagreen", "lightyellow2", "deeppink"))
  gpt <- ggplot(pd, aes(x = x, y = y, color = pd[, igene + 2])) + geom_point(size = pointsize) + 
    # scale_color_gradientn(colours=pal(5))+
    scale_color_gradientn(colours = pal(5)) + scale_x_discrete(expand = c(xpand, 
                                                                          ypand)) + scale_y_discrete(expand = c(xpand, ypand)) + coord_equal() + 
    # labs(title = colnames(pd)[igene+2], x = NULL, y = NULL)+
    theme_bw()
  if (main) {
    if (is.null(title)) {
      title = colnames(pd)[igene + 2]
    }
    out = gpt + labs(title = title, x = NULL, y = NULL) + theme(legend.position = legend.position, 
                                                                plot.title = element_text(hjust = hjust,face="bold",size = rel(titlesize)))
  } else {
    out = gpt + labs(title = NULL, x = NULL, y = NULL) + theme(legend.position = legend.position)
  }
  return(out)
}


#(1)for the dlpfc dataset with same donor
#choose the svgene list and  corresponse gene cluster result
dataset="samedonor" #alternative choose:acrossdonor,scc
load(here("RealData/result_data/genecluster",paste0("gene_clusters_",dataset,".RData")))
load(here("RealData/result_data/realdata svgene",paste0("data_dlpfc_",dataset,"_svgene_list.RData")))
svgene_list=data_dlpfc_samedonor_svgene_list

position1=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix1_position_",dataset,".csv")),row.names = 1))
count1=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix1_count_",dataset,".csv")),row.names = 1))
position2=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix2_position_",dataset,".csv")),row.names = 1))
count2=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix2_count_",dataset,".csv")),row.names = 1))
position3=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix3_position_",dataset,".csv")),row.names = 1))
count3=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix3_count_",dataset,".csv")),row.names = 1))
position4=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix4_position_",dataset,".csv")),row.names = 1))
count4=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix4_count_",dataset,".csv")),row.names = 1))

count1=count1[svgene_list[[6]],]
count2=count2[svgene_list[[6]],]
count3=count3[svgene_list[[6]],]
count4=count4[svgene_list[[6]],]
count_list=list(count1,count2,count3,count4)
position_list=list(position1,position2,position3,position4)

p1=f_genecluster(secnum=1,cluster=1,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 1",titlesize = 0.8,hjust=0.5)
p2=f_genecluster(secnum=2,cluster=1,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 1",titlesize = 0.8,hjust=0.5)
p3=f_genecluster(secnum=3,cluster=1,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 1",titlesize = 0.8,hjust=0.5)
p4=f_genecluster(secnum=4,cluster=1,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 1",titlesize = 0.8,hjust=0.5)
p5=f_genecluster(secnum=1,cluster=2,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 2",titlesize = 0.8,hjust=0.5)
p6=f_genecluster(secnum=2,cluster=2,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 2",titlesize = 0.8,hjust=0.5)
p7=f_genecluster(secnum=3,cluster=2,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 2",titlesize = 0.8,hjust=0.5)
p8=f_genecluster(secnum=4,cluster=2,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 2",titlesize = 0.8,hjust=0.5)
p9=f_genecluster(secnum=1,cluster=4,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 3",titlesize = 0.8,hjust=0.5)
p10=f_genecluster(secnum=2,cluster=4,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 3",titlesize = 0.8,hjust=0.5)
p11=f_genecluster(secnum=3,cluster=4,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 3",titlesize = 0.8,hjust=0.5)
p12=f_genecluster(secnum=4,cluster=4,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 3",titlesize = 0.8,hjust=0.5)
p13=f_genecluster(secnum=1,cluster=3,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 4",titlesize = 0.8,hjust=0.5)
p14=f_genecluster(secnum=2,cluster=3,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 4",titlesize = 0.8,hjust=0.5)
p15=f_genecluster(secnum=3,cluster=3,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 4",titlesize = 0.8,hjust=0.5)
p16=f_genecluster(secnum=4,cluster=3,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 4",titlesize = 0.8,hjust=0.5)

combined_plot <- 
  ((wrap_elements((p1 | p5) /( p9 | p13)) + 
      labs(tag = "A") + 
      theme(plot.tag = element_text(size = 14, face = "bold", margin = margin(0)),
            plot.tag.position = c(0.01, 1))) |
     (wrap_elements((p2 | p6) / ( p10 | p14)) + 
        labs(tag = "B") + 
        theme(plot.tag = element_text(size = 14, face = "bold", margin = margin(0)),
              plot.tag.position = c(0.01, 1)))) /
  ((wrap_elements((p3 | p7) /( p11 | p15)) + 
      labs(tag = "C") + 
      theme(plot.tag = element_text(size = 14, face = "bold", margin = margin(0)),
            plot.tag.position = c(0.01, 1))) |
     (wrap_elements((p4 | p8) /( p12 | p16)) + 
        labs(tag = "D") + 
        theme(plot.tag = element_text(size = 14, face = "bold", margin = margin(0)),
              plot.tag.position = c(0.01, 1)))) +
  plot_layout(heights = c(1, 1,0.2)) 
print(combined_plot)
ggsave(paste0("genecluster_samedonor.png"), plot = combined_plot, width = 9, height = 10, dpi = 300)


#(2)for the dlpfc dataset with across donor
dataset="acrossdonor" 
load(here("RealData/result_data/genecluster",paste0("gene_clusters_",dataset,".RData")))
load(here("RealData/result_data/realdata svgene",paste0("data_dlpfc_",dataset,"_svgene_list.RData")))
svgene_list=data_dlpfc_acrossdonor_svgene_list

position1=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix1_position_",dataset,".csv")),row.names = 1))
count1=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix1_count_",dataset,".csv")),row.names = 1))
position2=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix2_position_",dataset,".csv")),row.names = 1))
count2=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix2_count_",dataset,".csv")),row.names = 1))
position3=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix3_position_",dataset,".csv")),row.names = 1))
count3=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix3_count_",dataset,".csv")),row.names = 1))
position4=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix4_position_",dataset,".csv")),row.names = 1))
count4=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix4_count_",dataset,".csv")),row.names = 1))

count1=count1[svgene_list[[6]],]
count2=count2[svgene_list[[6]],]
count3=count3[svgene_list[[6]],]
count4=count4[svgene_list[[6]],]
count_list=list(count1,count2,count3,count4)
position_list=list(position1,position2,position3,position4)

p1=f_genecluster(secnum=1,cluster=1,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 1",titlesize = 0.8,hjust=0.5)
p2=f_genecluster(secnum=2,cluster=1,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 1",titlesize = 0.8,hjust=0.5)
p3=f_genecluster(secnum=3,cluster=1,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 1",titlesize = 0.8,hjust=0.5)
p4=f_genecluster(secnum=4,cluster=1,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 1",titlesize = 0.8,hjust=0.5)
p5=f_genecluster(secnum=1,cluster=2,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 2",titlesize = 0.8,hjust=0.5)
p6=f_genecluster(secnum=2,cluster=2,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 2",titlesize = 0.8,hjust=0.5)
p7=f_genecluster(secnum=3,cluster=2,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 2",titlesize = 0.8,hjust=0.5)
p8=f_genecluster(secnum=4,cluster=2,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 2",titlesize = 0.8,hjust=0.5)
p9=f_genecluster(secnum=1,cluster=3,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 3",titlesize = 0.8,hjust=0.5)
p10=f_genecluster(secnum=2,cluster=3,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 3",titlesize = 0.8,hjust=0.5)
p11=f_genecluster(secnum=3,cluster=3,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 3",titlesize = 0.8,hjust=0.5)
p12=f_genecluster(secnum=4,cluster=3,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 3",titlesize = 0.8,hjust=0.5)


combined_plot <- 
  ((wrap_elements((p1 | p5| p9 )) + 
      labs(tag = "A") + 
      theme(plot.tag = element_text(size = 14, face = "bold", margin = margin(0)),
            plot.tag.position = c(0.01, 1))) /
     (wrap_elements((p2 | p6| p10)) + 
        labs(tag = "B") + 
        theme(plot.tag = element_text(size = 14, face = "bold", margin = margin(0)),
              plot.tag.position = c(0.01, 1)))) /
  (wrap_elements((p3 | p7|p11)) + 
     labs(tag = "C") + 
     theme(plot.tag = element_text(size = 14, face = "bold", margin = margin(0)),
           plot.tag.position = c(0.01, 1))) /
  (wrap_elements((p4 | p8| p12 )) + 
     labs(tag = "D") + 
     theme(plot.tag = element_text(size = 14, face = "bold", margin = margin(0)),
           plot.tag.position = c(0.01, 1))) +
  plot_layout(heights = c(1,1,1,1,0.2)) 

print(combined_plot)

ggsave(paste0("genecluster_dlpfc_34711.png"), plot = combined_plot, width = 8, height = 12, dpi = 300)


#(3)for the scc dataset
dataset="scc" 
load(here("RealData/result_data/genecluster",paste0("gene_clusters_",dataset,".RData")))
load(here("RealData/result_data/realdata svgene",paste0("data_",dataset,"_svgene_list.RData")))
svgene_list=data_scc_svgene_list

position1=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/",dataset),paste0("matrix1_position_",dataset,".csv")),row.names = 1))
count1=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/",dataset),paste0("matrix1_count_",dataset,".csv")),row.names = 1))
position2=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/",dataset),paste0("matrix2_position_",dataset,".csv")),row.names = 1))
count2=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/",dataset),paste0("matrix2_count_",dataset,".csv")),row.names = 1))
position3=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/",dataset),paste0("matrix3_position_",dataset,".csv")),row.names = 1))
count3=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/",dataset),paste0("matrix3_count_",dataset,".csv")),row.names = 1))

count1=count1[svgene_list[[6]],]
count2=count2[svgene_list[[6]],]
count3=count3[svgene_list[[6]],]
count_list=list(count1,count2,count3)
position_list=list(position1,position2,position3)

p1=f_genecluster(secnum=1,cluster=1,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 1",titlesize = 1,hjust=0.5)
p2=f_genecluster(secnum=1,cluster=2,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 2",titlesize = 1,hjust=0.5)
p3=f_genecluster(secnum=1,cluster=3,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 3",titlesize = 1,hjust=0.5)
p4=f_genecluster(secnum=2,cluster=1,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 1",titlesize = 1,hjust=0.5)
p5=f_genecluster(secnum=2,cluster=2,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 2",titlesize = 1,hjust=0.5)
p6=f_genecluster(secnum=2,cluster=3,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 3",titlesize = 1,hjust=0.5)
p7=f_genecluster(secnum=3,cluster=1,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 1",titlesize = 1,hjust=0.5)
p8=f_genecluster(secnum=3,cluster=2,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 2",titlesize = 1,hjust=0.5)
p9=f_genecluster(secnum=3,cluster=3,result_gene_cluster_list=result_gene_cluster_list,count_list=count_list,position_list=position_list,main=TRUE,title ="Cluster 3",titlesize = 1,hjust=0.5)

combined_plot <- 
  (wrap_elements(p1 | p2|p3 ) + 
     labs(tag = "A") + 
     theme(plot.tag = element_text(size = 15, face = "bold", margin = margin(0)),
           plot.tag.position = c(0.01, 1))) /
  (wrap_elements(p4|p5 |p6) + 
     labs(tag = "B") + 
     theme(plot.tag = element_text(size = 15, face = "bold", margin = margin(0)),
           plot.tag.position = c(0.01, 1))) /
  (wrap_elements(p7 | p8 | p9) + 
     labs(tag = "C") + 
     theme(plot.tag = element_text(size = 15, face = "bold", margin = margin(0)),
           plot.tag.position = c(0.01, 1))) +
  plot_layout(heights = c(1, 1, 1,0.2)) 
print(combined_plot)
ggsave(paste0("genecluster_scc.png"), plot = combined_plot, width = 9, height = 10, dpi = 300)


###########################################################################################################################################
                                                  #4.heatmap plot for gene cluster:s15,s17,s19

############################################################################################################################################
#function for heatmap
heatmapplot <- function(result_gene_cluster_list, i, j, show_legend = TRUE) {
  # 计算重叠比例矩阵
  table_mat <- table(result_gene_cluster_list[[i]][[2]], result_gene_cluster_list[[j]][[2]])
  overlap_mat <- sweep(table_mat, 1, rowSums(table_mat), FUN = "/")
  
  # 设置切片标签
  x_slice <- ifelse(j == 1, "A", ifelse(j == 2, "B",ifelse(j==3,"C","D")))
  y_slice <- ifelse(i == 1, "A", ifelse(i == 2, "B",ifelse(i==3,"C","D")))
  
  # 转换数据为长格式
  df_long <- reshape2::melt(overlap_mat)
  colnames(df_long) <- c("Cluster1", "Cluster2", "Overlap")
  
  # 设置因子水平（保证顺序一致）
  df_long <- df_long[
    as.numeric(as.character(df_long$Cluster1)) %in% 1:4 &
      as.numeric(as.character(df_long$Cluster2)) %in% 1:4,
  ]
  df_long$Cluster2 <- factor(
    df_long$Cluster2,
    levels = colnames(overlap_mat),
    labels = paste0("cluster", colnames(overlap_mat))
  )
  df_long$Cluster1 <- factor(
    df_long$Cluster1,
    levels = rownames(overlap_mat),
    labels = paste0("cluster", rownames(overlap_mat))
  )
  
  # 绘图
  p1 <- ggplot(df_long, aes(x = Cluster2, y = Cluster1, fill = Overlap)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "#deebf7", high = "#3182bd", 
                        limits = c(0, 1),  # 固定颜色范围
                        name = "Overlap Ratio") +
    geom_text(aes(label = round(Overlap, 2)), color = "black", size = 5, fontface = "bold") +
    labs(
      x = paste0("slice ", x_slice),
      y = paste0("slice ", y_slice)
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_text(size = 10)
    )
  
  # 根据参数控制图例显示
  if (!show_legend) {
    p1 <- p1 + theme(legend.position = "none")
  } else {
    p1 <- p1 + guides(fill = guide_colorbar(barwidth = 1, barheight = 10))
  }
  
  return(p1)
}

#(1)for dlpfc with same donor
dataset="samedonor" 
load(here("RealData/result_data/genecluster",paste0("gene_clusters_",dataset,".RData")))
p1=heatmapplot(result_gene_cluster_list,1,2,show_legend = FALSE)
p2=heatmapplot(result_gene_cluster_list,1,3,show_legend = FALSE)
p3=heatmapplot(result_gene_cluster_list,1,4,show_legend = FALSE)
p4=heatmapplot(result_gene_cluster_list,2,3,show_legend = FALSE)
p5=heatmapplot(result_gene_cluster_list,2,4,show_legend = FALSE)
p6=heatmapplot(result_gene_cluster_list,3,4,show_legend = FALSE)
combined_heatmap=((p1 + p2 + p3)/(p4+p5+p6)) + plot_layout(guides = "collect") & theme(legend.position = "right")
plot(combined_heatmap)
ggsave(paste0("dlpfc_cluster_heatmap_samedonor.png"), plot = combined_heatmap, width = 13, height = 8, dpi = 300)

#(2)for dlpfc with across donor
dataset="acrossdonor" 
load(here("RealData/result_data/genecluster",paste0("gene_clusters_",dataset,".RData")))
p1=heatmapplot(result_gene_cluster_list,2,1,show_legend = FALSE)
p2=heatmapplot(result_gene_cluster_list,3,1,show_legend = FALSE)
p3=heatmapplot(result_gene_cluster_list,1,4,show_legend = FALSE)
p4=heatmapplot(result_gene_cluster_list,2,3,show_legend = FALSE)
p5=heatmapplot(result_gene_cluster_list,2,4,show_legend = FALSE)
p6=heatmapplot(result_gene_cluster_list,3,4,show_legend = FALSE)
combined_heatmap=((p1 + p2 + p3)/(p4+p5+p6)) + plot_layout(guides = "collect") & theme(legend.position = "right")
plot(combined_heatmap)
ggsave(paste0("dlpfc_cluster_heatmap_acrossdonor.png"), plot = combined_heatmap, width = 13, height = 8, dpi = 300)

#(3)for scc dataset
dataset="scc" 
load(here("RealData/result_data/genecluster",paste0("gene_clusters_",dataset,".RData")))
p1=heatmapplot(result_gene_cluster_list,2,1,show_legend = FALSE)
p2=heatmapplot(result_gene_cluster_list,3,1,show_legend = FALSE)
p3=heatmapplot(result_gene_cluster_list,2,3,show_legend = FALSE)
combined_heatmap=(p1 + p2 + p3) + plot_layout(guides = "collect") & theme(legend.position = "right")
plot(combined_heatmap)
ggsave(paste0("/Users/zhoum/Desktop/scc_cluster_heatmap.png"), plot = combined_heatmap, width = 10, height = 3, dpi = 300)


###########################################################################################################################################
                                              #5.enrichment for gene cluster:Table S17-S19

############################################################################################################################################
#function for gene cluster
enrich_function=function(result_gene_cluster_list){
  result_enrich_list=replicate(length(result_gene_cluster_list), replicate(4, list(), simplify = FALSE), simplify = FALSE)
  for(sectnum in c(1:length(result_gene_cluster_list))){
    clusters <- sort(unique(result_gene_cluster_list[[sectnum]]$cluster))
    for (cl in seq(clusters)){
      genes_in_cluster <- result_gene_cluster_list[[sectnum]]$gene[result_gene_cluster_list[[sectnum]]$cluster == clusters[cl]]
      gene_list <- bitr(genes_in_cluster, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      go_enrich_bp <- enrichGO(gene = gene_list$ENTREZID,
                               OrgDb = org.Hs.eg.db,
                               keyType = "ENTREZID",
                               ont = "BP", 
                               pAdjustMethod = "BH", 
                               qvalueCutoff = 0.05,
                               readable = TRUE)
      result_enrich_list[[sectnum]][[cl]]=go_enrich_bp
    }
  }
}

#produce the enrichment result and obtain the first five cluster of each dataset
for(dataset in c("samedonor","acrossdonor","scc")){
  load(here("RealData/result_data/genecluster",paste0("gene_clusters_",dataset,".RData")))
  result_enrich_list=enrich_function(result_gene_cluster_list)
  save(result_enrich_list,file=here("RealData/result_data/genecluster",paste0("enrich_result_",dataset,".RData")))#the result have been saved in correspond dictionary.
  cluster1_results <- list()
  for(clusternum in c(1:length(unique(result_gene_cluster_list[[1]]$cluster)))){
    for(sect in c(1:length(result_enrich_list))) {
      enrich_obj <- result_enrich_list[[sect]][[clusternum]]  
      if(!is.null(enrich_obj) && nrow(enrich_obj) > 0) {
        go_df <- as.data.frame(enrich_obj)[, c("ID", "Description", "p.adjust")]
        go_df$Section <- sect
        cluster1_results[[sect]] <- go_df
      }
    }
    combined_cluster1 <- do.call(rbind, cluster1_results)
    max_padjust_by_go <- aggregate(p.adjust ~ ID + Description, 
                                   data = combined_cluster1, 
                                   FUN = max, 
                                   na.rm = TRUE)
    max_padjust_by_go$Rank <- rank(max_padjust_by_go$p.adjust, ties.method = "min")
    final_result <- max_padjust_by_go[order(max_padjust_by_go$Rank), ]
    print(final_result[c(1:5),])
  }
}


















