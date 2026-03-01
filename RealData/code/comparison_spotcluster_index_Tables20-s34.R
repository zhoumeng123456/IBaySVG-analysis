##################################################################################################################################################
                                                        #1.generate the spot cluster results


##################################################################################################################################################
##example to produce the spot cluster results based on the identified SV genes
library(BASS)#a efficient algorithm to spatial clustering. Code and detailed implement can be avaiable at https://github.com/xzhoulab/BASS
library(fpc)
library(mclust)
library(aricode)


#load the data
dataset="samedonor"
matrix1=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix1_count_",dataset,".csv")),row.names = 1,check.names = FALSE))
matrix2=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix2_count_",dataset,".csv")),row.names = 1,check.names = FALSE))
matrix3=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix3_count_",dataset,".csv")),row.names = 1,check.names = FALSE))
matrix4=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix4_count_",dataset,".csv")),row.names = 1,check.names = FALSE))
position1=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix1_position_",dataset,".csv")),row.names = 1))
position2=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix2_position_",dataset,".csv")),row.names = 1))
position3=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix3_position_",dataset,".csv")),row.names = 1))
position4=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix4_position_",dataset,".csv")),row.names = 1))
calpha1=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix1_celltype_",dataset,".csv")),row.names = 1))
calpha2=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix2_celltype_",dataset,".csv")),row.names = 1))
calpha3=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix3_celltype_",dataset,".csv")),row.names = 1))
calpha4=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix4_celltype_",dataset,".csv")),row.names = 1))

colnames(position1)=c("row","col")
colnames(position2)=c("row","col")
colnames(position3)=c("row","col")
colnames(position4)=c("row","col")

rowname <- c("nnSVG","SPARKX","SPARK","spVC", "HEARTSVG","IBaySVG","DESpace","StarTrail","GASTON")
C <- 20
R <- 7
set.seed(0)
load(here("RealData/result_data/realdata svgene",paste0("data_dlpfc_",dataset,"_svgene_list.RData")))
svgene_list=data_dlpfc_samedonor_svgene_list
  
#implement
for(i in c(1:9)){
  if(i <=5){
    for(j in c(1:9)){
      select_gene=svgene_list[[i]][[j]]
      print(length(select_gene))
      choose_matrix1=matrix1[select_gene,]
      choose_matrix2=matrix2[select_gene,]
      choose_matrix3=matrix3[select_gene,]
      choose_matrix4=matrix4[select_gene,]
      
      sparse_matrix1 <- as(choose_matrix1, "dgCMatrix")
      sparse_matrix2 <- as(choose_matrix2, "dgCMatrix")
      sparse_matrix3 <- as(choose_matrix3, "dgCMatrix")
      sparse_matrix4 <- as(choose_matrix4, "dgCMatrix")
      
      cntm=list(sparse_matrix1,sparse_matrix2,sparse_matrix3,sparse_matrix4)
      xym=list(position1,position2,position3,position4)
      BASS <- createBASSObject(cntm, xym, C = C, R = R,
                               beta_method = "SW", init_method = "mclust", 
                               nsample = 10000)
      BASS <- BASS.preprocess(BASS, doLogNormalize = TRUE,
                              geneSelect = "sparkx", nSE = 3000, doPCA = TRUE, 
                              scaleFeature = FALSE, nPC = 20)
      BASS <- BASS.run(BASS)
      BASS <- BASS.postprocess(BASS)
      zlabels <- BASS@results$z
      result_total=zlabels
      save(result_total,file=paste0(rowname[i],j,"_bass_cluster.RData"))
    }
  }else if(i==6|i==7){
    select_gene=svgene_list[[i]]
    choose_matrix1=matrix1[select_gene,]
    choose_matrix2=matrix2[select_gene,]
    choose_matrix3=matrix3[select_gene,]
    choose_matrix4=matrix4[select_gene,]
    
    sparse_matrix1 <- as(choose_matrix1, "dgCMatrix")
    sparse_matrix2 <- as(choose_matrix2, "dgCMatrix")
    sparse_matrix3 <- as(choose_matrix3, "dgCMatrix")
    sparse_matrix4 <- as(choose_matrix4, "dgCMatrix")
    
    cntm=list(sparse_matrix1,sparse_matrix2,sparse_matrix3,sparse_matrix4)
    xym=list(position1,position2,position3,position4)
    BASS <- createBASSObject(cntm, xym, C = C, R = R,
                             beta_method = "SW", init_method = "mclust", 
                             nsample = 10000)
    BASS <- BASS.preprocess(BASS, doLogNormalize = TRUE,
                            geneSelect = "sparkx", nSE = 3000, doPCA = TRUE, 
                            scaleFeature = FALSE, nPC = 20)
    BASS <- BASS.run(BASS)
    BASS <- BASS.postprocess(BASS)
    zlabels <- BASS@results$z
    result_total=zlabels
    save(result_total,file=paste0(rowname[i],"_bass_cluster.RData"))
  }else{
    for(j in c(1:7)){
      select_gene=svgene_list[[i]][[j]]
      print(length(select_gene))
      choose_matrix1=matrix1[select_gene,]
      choose_matrix2=matrix2[select_gene,]
      choose_matrix3=matrix3[select_gene,]
      choose_matrix4=matrix4[select_gene,]
      
      sparse_matrix1 <- as(choose_matrix1, "dgCMatrix")
      sparse_matrix2 <- as(choose_matrix2, "dgCMatrix")
      sparse_matrix3 <- as(choose_matrix3, "dgCMatrix")
      sparse_matrix4 <- as(choose_matrix4, "dgCMatrix")
      
      cntm=list(sparse_matrix1,sparse_matrix2,sparse_matrix3,sparse_matrix4)
      xym=list(position1,position2,position3,position4)
      BASS <- createBASSObject(cntm, xym, C = C, R = R,
                               beta_method = "SW", init_method = "mclust", 
                               nsample = 10000)
      BASS <- BASS.preprocess(BASS, doLogNormalize = TRUE,
                              geneSelect = "sparkx", nSE = 3000, doPCA = TRUE, 
                              scaleFeature = FALSE, nPC = 20)
      BASS <- BASS.run(BASS)
      BASS <- BASS.postprocess(BASS)
      zlabels <- BASS@results$z
      result_total=zlabels
      save(result_total,file=paste0(rowname[i],j,"_bass_cluster.RData"))
    }
  }
}


##################################################################################################################################################
                                                      #2.ARI,NMI results for Table s20,S21,s26,S27


##################################################################################################################################################
#plot of real spot cluster in dlpfc dataset with same donor
dataset="dlpfc_samedonor"# Alternative: "dlpfc_acrossdonor"
load(here("RealData/result_data/spotcluster",paste0("spot_real_",dataset,".RData")))
load(here("RealData/result_data/spotcluster",paste0("total_bass_cluster_",dataset,".RData")))

for(j in c(1:4)){
  ari_matrix=matrix(NA,nrow=10,ncol=7)
  rownames(ari_matrix)=c("nnSVG","SPARKX","SPARK","spVC", "HEARTSVG","IBaySVG","DESpace","StarTrail","GASTON","Full trans")
  colnames(ari_matrix)=c("Multi-layer","Single-layer","Union","Intersect","Cauchy","HMP", "PASTE")
  for(i in c(1:10)){
    if(i <=5){
      ari_matrix[i,2]=adjustedRandIndex(result_cluster[[i]][[j]][[j]],layer_barcode[[j]])
      for(k in c(5:9)){
        
        ari_matrix[i,k-2]=adjustedRandIndex(result_cluster[[i]][[k]][[j]],layer_barcode[[j]])
      }
    }else if(i==6|i==7|i==10){
      ari_matrix[i,1]=adjustedRandIndex(result_cluster[[i]][[1]][[j]],layer_barcode[[j]])
    }else{
      ari_matrix[i,2]=adjustedRandIndex(result_cluster[[i]][[j]][[j]],layer_barcode[[j]])
      for(k in c(5:7)){
        ari_matrix[i,k-2]=adjustedRandIndex(result_cluster[[i]][[k]][[j]],layer_barcode[[j]])
      }
    }
  }
  print(ari_matrix)
}

#compute NMI
for(j in c(1:4)){
  nmi_matrix=matrix(NA,nrow=10,ncol=7)
  rownames(nmi_matrix)=c("nnSVG","SPARKX","SPARK","spVC", "HEARTSVG","IBaySVG","DESpace","StarTrail","GASTON","Full trans")
  colnames(nmi_matrix)=c("Multi-layer","Single-layer","Union","Intersect","Cauchy","HMP", "PASTE")
  for(i in c(1:10)){
    if(i <=5){
      nmi_matrix[i,2]=NMI(layer_barcode[[j]],result_cluster[[i]][[j]][[j]])
      for(k in c(5:9)){
        nmi_matrix[i,k-2]=NMI(layer_barcode[[j]],result_cluster[[i]][[k]][[j]])
      }
    }else if(i==6|i==7|i==10){
      nmi_matrix[i,1]=NMI(layer_barcode[[j]],result_cluster[[i]][[1]][[j]])
    }else{
      nmi_matrix[i,2]=NMI(layer_barcode[[j]],result_cluster[[i]][[j]][[j]])
      for(k in c(5:7)){
        nmi_matrix[i,k-2]=NMI(layer_barcode[[j]],result_cluster[[i]][[k]][[j]])
      }
    }
  }
  print(nmi_matrix)
}



##################################################################################################################################################
                                                     #3.normal index results for Table s22-S25 and s28-s34

##################################################################################################################################################
library(parallelDist)
library(cluster)
library(clusterCrit)
library(patchwork)
##1.produce the index function
calculate_dbi <- function(distance_matrix, cluster_labels) {
  clusters <- unique(cluster_labels)
  n_clusters <- length(clusters)
  
  scatter <- sapply(clusters, function(cluster) {
    points_in_cluster <- which(cluster_labels == cluster)
    cluster_points <- distance_matrix[points_in_cluster, points_in_cluster, drop = FALSE]
    mean(apply(cluster_points, 1, function(row) mean(row, na.rm = TRUE)), na.rm = TRUE)
  })
  
  dbi <- 0
  for (i in 1:n_clusters) {
    max_ratio <- -Inf
    for (j in 1:n_clusters) {
      if (i != j) {
        distance_ij <- mean(distance_matrix[which(cluster_labels == clusters[i]), 
                                            which(cluster_labels == clusters[j])])
        ratio <- (scatter[i] + scatter[j]) / distance_ij
        max_ratio <- max(max_ratio, ratio, na.rm = TRUE)
      }
    }
    dbi <- dbi + max_ratio
  }
  
  dbi / n_clusters  
}
symbolgenerate<-function(distance_matrix,cluster_result,sectnum=2,count1,idomain=c(1:8),sectnumber=4){
  ##compute silhouette score
  result_lunkuo_matrix=matrix(NA,nrow = 9,ncol = 9)
  rownames(result_lunkuo_matrix)=c("nnSVG","SPARKX","SPARK","spVC","HEARTSVG","Proposed","DEspace","startrail","gaston")
  for(i in idomain){
    if(i<=5){
      for(j in c(1:(sectnumber+5))){
        silhouette_score <- silhouette(as.numeric(cluster_result[[i]][[j]]), distance_matrix)  
        avg_silhouette_score <- mean(silhouette_score[, 3])  
        result_lunkuo_matrix[i,j]=round(avg_silhouette_score,3)
      }
    }else if(i==6|i==7){
        silhouette_score <- silhouette(as.numeric(cluster_result[[i]][[1]]), distance_matrix)  
        avg_silhouette_score <- mean(silhouette_score[, 3])  
        result_lunkuo_matrix[i,1]=round(avg_silhouette_score,3)
    }else{
      for(j in c(1:(sectnumber+3))){
        silhouette_score <- silhouette(as.numeric(cluster_result[[i]][[j]]), distance_matrix)  
        avg_silhouette_score <- mean(silhouette_score[, 3]) 
        result_lunkuo_matrix[i,j]=round(avg_silhouette_score,3)
      }
    }
  }
  print(result_lunkuo_matrix)

  ##compute dbi
  result_dbi_matrix=matrix(NA,nrow = 9,ncol = 9)
  rownames(result_dbi_matrix)=c("nnSVG","SPARKX","SPARK","spVC","HEARTSVG","Proposed","DEspace","startrail","gaston")
  for(i in idomain){
    if(i <=5){
      for(j in c(1:(sectnumber+5))){
        dbi_value <- calculate_dbi(as.matrix(distance_matrix), as.numeric(cluster_result[[i]][[j]]))
        result_dbi_matrix[i,j]=round(dbi_value,3) 
      }
    }else if(i==6|i==7){
      dbi_value <- calculate_dbi(as.matrix(distance_matrix), as.numeric(cluster_result[[i]][[1]]))
      result_dbi_matrix[i,1]=round(dbi_value,3) 
    }else{
      for(j in c(1:(sectnumber+3))){
        dbi_value <- calculate_dbi(as.matrix(distance_matrix), as.numeric(cluster_result[[i]][[j]]))
        result_dbi_matrix[i,j]=round(dbi_value,3) 
      }
    }
  }
  print(result_dbi_matrix)
  
  ##compute ch
  result_ch_matrix=matrix(NA,nrow = 9,ncol = 9)
  rownames(result_ch_matrix)=c("nnSVG","SPARKX","SPARK","spVC","HEARTSVG","Proposed","DEspace","startrail","gaston")
  for(i in idomain){
    if(i <=5){
      for(j in c(1:(sectnumber+5))){
        ch_index <- cluster.stats(as.matrix(distance_matrix), as.numeric(cluster_result[[i]][[j]]))$ch
        result_ch_matrix[i,j]=round(ch_index,3)
      }
    }else if(i==6|i==7){
      ch_index <- cluster.stats(as.matrix(distance_matrix), as.numeric(cluster_result[[i]][[1]]))$ch
      result_ch_matrix[i,1]=round(ch_index,3)
    }else{
      for(j in c(1:(sectnumber+3))){
        ch_index <- cluster.stats(as.matrix(distance_matrix), as.numeric(cluster_result[[i]][[j]]))$ch
        result_ch_matrix[i,j]=round(ch_index,3)
      }
    }
  }
  print(result_ch_matrix)
  
  #compute avona result 
  mean_count1=apply(count1,1,mean)
  result_avo_matrix=matrix(NA,nrow = 9,ncol = 9)
  rownames(result_avo_matrix)=c("nnSVG","SPARKX","SPARK","spVC","HEARTSVG","Proposed","DEspace","startrail","gaston")
  for(i in idomain){
    if(i <=5){
      for(j in c(1:(sectnumber+5))){
        unique_clusters <- unique(as.numeric(cluster_result[[i]][[j]]))
        regions <- list()
        for (cluster in unique_clusters) {
          region_data <- mean_count1[as.numeric(cluster_result[[i]][[j]]) == cluster]
          regions[[paste0("Region", cluster)]] <- region_data
        }
        
        combined_data <- data.frame(
          expression = unlist(regions), 
          region = factor(rep(names(regions), times = sapply(regions, length)))  
        )
        
        anova_result <- aov(expression ~ region, data = combined_data)
        result_avo_matrix[i,j]=round(summary(anova_result)[[1]]$`F value`[1],3)
      }
    }else if(i==6|i==7){
      unique_clusters <- unique(as.numeric(cluster_result[[i]][[1]]))
      regions <- list()
      for (cluster in unique_clusters) {
        region_data <- mean_count1[as.numeric(cluster_result[[i]][[1]]) == cluster]
        regions[[paste0("Region", cluster)]] <- region_data
      }
      
      combined_data <- data.frame(
        expression = unlist(regions),  
        region = factor(rep(names(regions), times = sapply(regions, length)))  
      )
      
      anova_result <- aov(expression ~ region, data = combined_data)
      result_avo_matrix[i,1]=round(summary(anova_result)[[1]]$`F value`[1],3)
    }else{
      for(j in c(1:(sectnumber+3))){
        unique_clusters <- unique(as.numeric(cluster_result[[i]][[j]]))
        regions <- list()
        for (cluster in unique_clusters) {
          region_data <- mean_count1[as.numeric(cluster_result[[i]][[j]]) == cluster]
          regions[[paste0("Region", cluster)]] <- region_data
        }
        
        combined_data <- data.frame(
          expression = unlist(regions),  
          region = factor(rep(names(regions), times = sapply(regions, length)))  
        )
        
        anova_result <- aov(expression ~ region, data = combined_data)
        result_avo_matrix[i,j]=round(summary(anova_result)[[1]]$`F value`[1],3)
      }
    }
  }
  print(result_avo_matrix)
  
}

##2.for dlpfc dataset
dataset="acrossdonor"#Alternative: "acrossdonor"
load(here("RealData/result_data/spotcluster",paste0("total_bass_cluster_dlpfc_",dataset,".RData")))

for(datanum in c(1:4)){
  count1=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/dlpfc ",dataset),paste0("matrix1_count_",dataset,".csv")),row.names = 1,check.names = FALSE))
  count1=t(count1)
  #distance_matrix <- parallelDist::parDist(as.matrix(count1), method = "euclidean", threads = 8)
  #save(distance_matrix,file=paste0("/Users/zhoum/Desktop/dlpfc_data/distance",datanum,"_matrix.RData"))
  load(paste0("C:/Users/zhoum/OneDrive/distance_matrix_0714_",datanum,".RData"))
  cluster_result=cluster_dlpfc[[datanum]]
  index_result=symbolgenerate(distance_matrix=distance_matrix,cluster_result=cluster_result,sectnum=datanum,count1=count1,idomain=c(1:9))
  print(index_result)
}

##for scc dataset
dataset="scc"
load(here("RealData/result_data/spotcluster",paste0("total_bass_cluster_",dataset,".RData")))
for(datanum in c(1:3)){
  count1=as.matrix(read.csv(here(paste0("RealData/result_data/realdata dataset/scc"),paste0("matrix",datanum,"_count_",dataset,".csv")),row.names = 1,check.names = FALSE))
  count1=t(count1)
  distance_matrix <- parallelDist::parDist(as.matrix(count1), method = "euclidean", threads = 8)
  cluster_result=cluster_scc[[datanum]]
  index_result=symbolgenerate(distance_matrix=distance_matrix,cluster_result=cluster_result,sectnum=datanum,count1=count1,idomain=c(1:9),sectnumber=3)
  print(index_result)
}






































