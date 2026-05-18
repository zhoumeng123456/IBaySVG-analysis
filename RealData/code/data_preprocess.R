#Here are the examples for proprocessing the raw data, while the processed datasets are available in "data/Realdataset/dlpfc acrossdonor",
#"data/Realdataset/dlpfc samedonor" and "data/Realdataset/scc".


##########################################################################################################
                                  #1.preprocess for the DLPFC data
        
##########################################################################################################

library(rhdf5)
library(Matrix)
library(Seurat)
library(Redeconve) #remotes::install_github("ZxZhou4150/Redeconve", build_vignettes = F)
library(HDF5Array)

##1.read and handle the single-cell data 
# Note: Due to the large size of the `assays.h5` file, this step may be memory-intensive and difficult to 
# execute on systems with limited RAM. Users are advised to ensure adequate local memory resources before 
# running this step.
# or you can directly use the pre-filtered scRNA data stored in "result_matrix.RData" and turn to step 2
load(here::here("data/Realdataset/raw_data/dlpfc_dataset","result_matrix.RData"))

#load the snRNA-seq cell-type and gene annotations
single1 <- readRDS(here::here("data/Realdataset/raw_data/dlpfc_dataset","se.rds"))

#load the snRNA-seq
#Note: Due to GitHub's file size limitations, ultra-large raw baseline references (`assays.h5`) exceed 
#the upload threshold and cannot be directly archived in this repository. Users **must** manually download
#the complete reference package via [Dropbox Link](https://www.dropbox.com/s/5919zt00vm1ht8e/sce_DLPFC_annotated.zip?dl=1)
#unzip it, and place the extracted files directly into the `data/Realdataset/raw_data/dlpfc_dataset/` directory 
#before executing the preprocessing pipeline.
file_path<-here::here("data/Realdataset/raw_data/dlpfc_dataset","assays.h5")
h5ls(file_path)
data <- h5read(file_path, "/assay001")

genelist=single1@rowRanges@ranges@NAMES
rownames(data)=genelist
barcode=single1@colData@listData[["Barcode"]]
colnames(data)=barcode

#choose the cell-type annotations
celllist=as.character(single1@colData@listData[["cellType_broad_k"]])
unique_types <- unique(celllist)
n=nrow(data)
k <- length(unique_types)
result_matrix <- matrix(NA, nrow = n, ncol = k)
for (i in seq_along(unique_types)) {
  cols <- which(celllist == unique_types[i])
  result_matrix[, i] <- rowMeans(data[, cols, drop = FALSE])
}
colnames(result_matrix) <- unique_types
rownames(result_matrix) <- genelist

# obtain the simplified scRNA-seq reference matrix
result_matrix <- Matrix(result_matrix, sparse = TRUE)


## 2.load the ST matrix and apply the redeconve
#The position{X}.txt is copied by the tissue_positions_list.txt file available at https://github.com/LieberInstitute/HumanPilot.
#The `.h5` files correspond to the 12 DLPFC samples, identified by their sample IDs (e.g., 151507–151676). 
#The files `position1`–`position12.txt` contain the corresponding spatial coordinate information for each sample.
#For example, `position1` corresponds to sample `151507_filtered_feature_bc_matrix.h5`, 
#`position12` corresponds to `151676_filtered_feature_bc_matrix.h5`, and so forth.

sample_ids <- c(
  "151507", "151508", "151509", "151510",
  "151669", "151670", "151671", "151672",
  "151673", "151674", "151675", "151676"
)

make_unique_gene_names <- function(gene_names) {
  make.unique(as.character(gene_names), sep = "_")
}

## function to generate Redeconve deconvolution results
load_dlpfc_sample <- function(sample_index,result_matrix,data_dir = here::here("data/Realdataset/raw_data/dlpfc_dataset"),ncores = 8) {
  
  sample_id <- sample_ids[sample_index]
  
  position_file <- file.path(data_dir, paste0("position", sample_index, ".txt"))
  h5_file <- file.path(data_dir, paste0(sample_id, "_filtered_feature_bc_matrix.h5"))
  
  message("Processing sample ", sample_index, ": ", sample_id)
  
  # read spatial coordinates
  position <- read.csv(position_file, header = FALSE)
  position <- position[position$V2 == 1, ]
  rownames(position) <- position$V1
  position <- position[, c("V5", "V6")]
  colnames(position) <- c("x", "y")
  
  # read 10X h5 matrix
  stdata <- h5read(h5_file, "/matrix/data")
  indices <- h5read(h5_file, "/matrix/indices")
  indptr <- h5read(h5_file, "/matrix/indptr")
  shape <- h5read(h5_file, "/matrix/shape")
  barcodes <- h5read(h5_file, "/matrix/barcodes")
  gene_names <- h5read(h5_file, "/matrix/features/name")
  
  gene_names <- make_unique_gene_names(gene_names)
  
  # transform CSR into sparse matrix
  matrix1 <- sparseMatrix(
    i = rep(seq_len(length(indptr) - 1), diff(indptr)),
    j = indices + 1,
    x = stdata,
    dims = c(shape[2], shape[1])
  )
  
  rownames(matrix1) <- barcodes
  colnames(matrix1) <- gene_names
  matrix1 <- t(matrix1)
  
  # deconvolution
  res <- deconvoluting(
    result_matrix,
    matrix1,
    genemode = "def",
    hpmode = "def",
    dopar = TRUE,
    ncores = ncores,
    normalize = TRUE
  )
  
  res <- res[rownames(res) != "drop", , drop = FALSE]
  
  # obtain spot × cell-type proportion matrix
  res.prop <- to.proportion(res)
  res.prop <- t(res.prop)
  
  return(list(
    expression = matrix1,
    position = position,
    deconvolution = res.prop
  ))
}

#for DLPFC dataset with same donor, we integrate the sample with index 1,2,3,4
#for DLPFC dataset with different donor, we integrate the sample with index 3,4,7,11
inte_sample_index=c(1,2,3,4)
dlpfc1 <- load_dlpfc_sample(inte_sample_index[1], result_matrix,data_dir = here::here("data/Realdataset/raw_data/dlpfc_dataset"))
dlpfc2 <- load_dlpfc_sample(inte_sample_index[2], result_matrix,data_dir = here::here("data/Realdataset/raw_data/dlpfc_dataset"))
dlpfc3 <- load_dlpfc_sample(inte_sample_index[3], result_matrix,data_dir = here::here("data/Realdataset/raw_data/dlpfc_dataset"))
dlpfc4 <- load_dlpfc_sample(inte_sample_index[4], result_matrix,data_dir = here::here("data/Realdataset/raw_data/dlpfc_dataset"))

                            
##3. variable selection and data cleaning

#choose the sample ID of integrated four samples and filter the common gene
file_path=here::here("data/Realdataset/raw_data/dlpfc_dataset",paste0(sample_ids[inte_sample_index[1]],"_filtered_feature_bc_matrix.h5"))
gene_names1 <- h5read(file_path, "/matrix/features/name")  
file_path=here::here("data/Realdataset/raw_data/dlpfc_dataset",paste0(sample_ids[inte_sample_index[2]],"_filtered_feature_bc_matrix.h5"))
gene_names2 <- h5read(file_path, "/matrix/features/name")  
file_path=here::here("data/Realdataset/raw_data/dlpfc_dataset",paste0(sample_ids[inte_sample_index[3]],"_filtered_feature_bc_matrix.h5"))
gene_names3 <- h5read(file_path, "/matrix/features/name")  
file_path=here::here("data/Realdataset/raw_data/dlpfc_dataset",paste0(sample_ids[inte_sample_index[4]],"_filtered_feature_bc_matrix.h5"))
gene_names4 <- h5read(file_path, "/matrix/features/name")  
inte_gene_name = Reduce(intersect, list(gene_names1,gene_names2, gene_names3,gene_names4))

#prefilter the count matrix
dlpfc1[[1]] <- dlpfc1[[1]][rownames(dlpfc1[[1]]) %in% inte_gene_name, , drop = FALSE]
dlpfc2[[1]] <- dlpfc2[[1]][rownames(dlpfc2[[1]]) %in% inte_gene_name, , drop = FALSE]
dlpfc3[[1]] <- dlpfc3[[1]][rownames(dlpfc3[[1]]) %in% inte_gene_name, , drop = FALSE]
dlpfc4[[1]] <- dlpfc4[[1]][rownames(dlpfc4[[1]]) %in% inte_gene_name, , drop = FALSE]

#further filtering based on the gene expression matrix
genefilter <- function(mat1,position1,res1,num_features=100,num_cells=100,nfeature=5000){
  gene_names=rownames(mat1)
  barcodes=colnames(mat1)
  
  rownames(mat1) <- paste("gene",c(1:length(rownames(mat1))))
  colnames(mat1) <- paste("cell",c(1:length(colnames(mat1))))
  
  seurat_object <- CreateSeuratObject(counts = mat1, project = "scRNAseq", min.cells = num_cells, min.features = num_features)#获取高可辨基因
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = nfeature)
  variable_genes <- Seurat::VariableFeatures(seurat_object)
  seurat_object_hvg <- subset(seurat_object, features = variable_genes)
  counts_matrix <- Seurat::GetAssayData(seurat_object_hvg, layer = "counts")
  
  st_gene_indices <- as.numeric(gsub("gene", "", rownames(counts_matrix)))
  rownames(counts_matrix) <- gene_names[st_gene_indices]
  st_cell_indices <- as.numeric(gsub("cell", "", colnames(counts_matrix)))
  colnames(counts_matrix) <- barcodes[st_cell_indices]
  
  position1 <- position1[rownames(position1) %in% colnames(counts_matrix), , drop = FALSE]
  position1 <- position1[colnames(counts_matrix), , drop = FALSE]
  
  res1 <- res1[rownames(res1) %in% colnames(counts_matrix), , drop = FALSE]
  res1 <- res1[colnames(counts_matrix), , drop = FALSE]
  
  counts_matrix=as.matrix(counts_matrix)
  return(list(counts_matrix,position1,res1))
}

result1=genefilter(dlpfc1[[1]],dlpfc1[[2]],dlpfc1[[3]],num_features = 500,num_cells=100,nfeature = 8000)
result2=genefilter(dlpfc2[[1]],dlpfc2[[2]],dlpfc2[[3]],num_features = 500,num_cells=100,nfeature = 8000)
result3=genefilter(dlpfc3[[1]],dlpfc3[[2]],dlpfc3[[3]],num_features = 500,num_cells=100,nfeature = 8000)
result4=genefilter(dlpfc4[[1]],dlpfc4[[2]],dlpfc4[[3]],num_features = 500,num_cells=100,nfeature = 8000)

dim(result1[[1]])
dim(result2[[1]])
dim(result3[[1]])
dim(result4[[1]])

#make sure the shared genes
totalgene=Reduce(intersect,list(rownames(result1[[1]]),rownames(result2[[1]]),rownames(result3[[1]]),rownames(result4[[1]])))
length(totalgene)


dlpfc1[[1]] <- dlpfc1[[1]][rownames(dlpfc1[[1]]) %in% totalgene, , drop = FALSE]
dlpfc2[[1]] <- dlpfc2[[1]][rownames(dlpfc2[[1]]) %in% totalgene, , drop = FALSE]
dlpfc3[[1]] <- dlpfc3[[1]][rownames(dlpfc3[[1]]) %in% totalgene, , drop = FALSE]
dlpfc4[[1]] <- dlpfc4[[1]][rownames(dlpfc4[[1]]) %in% totalgene, , drop = FALSE]

dlpfc1[[1]] <- dlpfc1[[1]][,colnames(dlpfc1[[1]]) %in% colnames(result1[[1]]), drop = FALSE]
dlpfc2[[1]] <- dlpfc2[[1]][,colnames(dlpfc2[[1]]) %in% colnames(result2[[1]]), drop = FALSE]
dlpfc3[[1]] <- dlpfc3[[1]][,colnames(dlpfc3[[1]]) %in% colnames(result3[[1]]), drop = FALSE]
dlpfc4[[1]] <- dlpfc4[[1]][,colnames(dlpfc4[[1]]) %in% colnames(result4[[1]]), drop = FALSE]


##4.output the result
write.csv(as.matrix(dlpfc1[[1]]), file = here::here("data/Realdataset/raw_data/dlpfc_dataset","matrix1_count_dlpfc.csv"), row.names = TRUE)
write.csv(result1[[2]], file = here::here("data/Realdataset/raw_data/dlpfc_dataset","matrix1_position_dlpfc.csv"), row.names = TRUE)
write.csv(result1[[3]], file = here::here("data/Realdataset/raw_data/dlpfc_dataset","matrix1_celltype_dlpfc.csv"), row.names = TRUE)
write.csv(as.matrix(dlpfc2[[1]]), file = here::here("data/Realdataset/raw_data/dlpfc_dataset","matrix2_count_dlpfc.csv"), row.names = TRUE)
write.csv(result2[[2]], file = here::here("data/Realdataset/raw_data/dlpfc_dataset","matrix2_position_dlpfc.csv"), row.names = TRUE)
write.csv(result2[[3]], file = here::here("data/Realdataset/raw_data/dlpfc_dataset","matrix2_celltype_dlpfc.csv"), row.names = TRUE)
write.csv(as.matrix(dlpfc3[[1]]), file = here::here("data/Realdataset/raw_data/dlpfc_dataset","matrix3_count_dlpfc.csv"), row.names = TRUE)
write.csv(result3[[2]], file = here::here("data/Realdataset/raw_data/dlpfc_dataset","matrix3_position_dlpfc.csv"), row.names = TRUE)
write.csv(result3[[3]], file = here::here("data/Realdataset/raw_data/dlpfc_dataset","matrix3_celltype_dlpfc.csv"), row.names = TRUE)
write.csv(as.matrix(dlpfc4[[1]]), file = here::here("data/Realdataset/raw_data/dlpfc_dataset","matrix4_count_dlpfc.csv"), row.names = TRUE)
write.csv(result4[[2]], file = here::here("data/Realdataset/raw_data/dlpfc_dataset","matrix4_position_dlpfc.csv"), row.names = TRUE)
write.csv(result4[[3]], file = here::here("data/Realdataset/raw_data/dlpfc_dataset","matrix4_celltype_dlpfc.csv"), row.names = TRUE)



##########################################################################################################
                                          #2.preprocess for the SCC data

##########################################################################################################


##1.load the scRNA data and apply the redeconve
scc_sc_barcode <- read.table(here::here("data/Realdataset/raw_data/scc_dataset","GSE144236_patient_metadata_new.txt"),header = TRUE, sep = "\t", stringsAsFactors = FALSE)
scc_sc_barcode<-scc_sc_barcode[scc_sc_barcode$tum.norm=="Tumor"&scc_sc_barcode$patient=="P9",]#26299
scc_sc_barcode=scc_sc_barcode[,"level2_celltype",drop = FALSE]
# This "merge10pts_counts.txt" file is available at the GSE144236_cSCC_counts.txt.gz 
# Note: Due to the large size of the `merge10pts_counts.txt` file, this step may be memory-intensive and difficult to 
# execute on systems with limited RAM. Users are advised to ensure adequate local memory resources before 
# running this step.
scc_count <- read.table(here::here("data/Realdataset/raw_data/scc_dataset","merge10pts_counts.txt"),header = TRUE, sep = "\t", stringsAsFactors = FALSE)
scc_sc_count=scc_count[,rownames(scc_sc_barcode)]
colnames(scc_sc_count)=scc_sc_barcode$level2_celltype

#or you can directly use the pre-filtered scRNA data stored in "scc_sc_count.RData"
load(here::here("data/Realdataset/raw_data/scc_dataset","scc_sc_count.RData"))

# remove non-gene annotation rows
scc_sc_count=scc_sc_count[-c(1,2),]
colnames(scc_sc_count)=scc_sc_barcode[,1]
cleaned_names=scc_sc_barcode[,1]

# average expression within each cell type
scc_sc_avg <- t(sapply(unique(cleaned_names), function(ct) {
  cells_in_ct <- which(cleaned_names == ct)
  rowMeans(scc_sc_count[, cells_in_ct, drop = FALSE])
}))
scc_sc_avg <- t(scc_sc_avg)# gene × celltype 


##2.load the ST data and apply the Redeconve
#function for handling scc dataset
scc_handle <-function(mat){
  spotname=colnames(mat)
  x_coord <- as.integer(sub("(\\d+)x(\\d+)", "\\1", spotname))
  y_coord <- as.integer(sub("(\\d+)x(\\d+)", "\\2", spotname))
  position=matrix(NA,nrow = length(spotname),ncol = 2)
  position[,1]=x_coord
  position[,2]=y_coord
  spotnames=paste("spot",c(1:length(spotname)),sep = "")
  colnames(mat)=spotnames
  colnames(position)=c("x","y")
  rownames(position)=spotnames
  return(list(mat,position))
}
scc_handle2<-function(repnum,dir="data/Realdataset/raw_data/scc_dataset"){
  scc1 <- read.table(here::here(dir,paste0("scc_P9_rep",repnum,".tsv")), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  rownames(scc1)=make.unique(as.character(scc1[,1]))
  scc1=as.matrix(scc1[,-1])
  location1=read.table(here::here(dir,paste0("spot_P9_rep",repnum,".tsv")), header = TRUE, sep = "\t", stringsAsFactors = FALSE)[,c(1,2)]
  
  scc_mat1=scc_handle(t(scc1))[[1]]
  scc_loc1=scc_handle(t(scc1))[[2]]
  
  scc_key <- paste(scc_loc1[,1], scc_loc1[,2], sep = "_")
  loc_key <- paste(location1[,1], location1[,2], sep = "_")
  
  matrix1 <- scc_mat1[,scc_key %in% loc_key ]
  position1 <- scc_loc1[scc_key %in% loc_key, ]
  
  spotnum=as.numeric(dim(position1)[1])
  colnames(matrix1)=paste("spot",c(1:spotnum))
  rownames(position1)=paste("spot",c(1:spotnum))
  return(list(matrix1,position1))
}

matrix1=scc_handle2(1)[[1]]#gene × spot
matrix2=scc_handle2(2)[[1]]
matrix3=scc_handle2(3)[[1]]
position1=scc_handle2(1)[[2]]#spot × c(x,y)
position2=scc_handle2(2)[[2]]
position3=scc_handle2(3)[[2]]

#function for applying the Redeconve
cell_compute=function(sc,count){
  res <- deconvoluting(sc, count, genemode = "def", hpmode = "def", dopar = T, ncores = 8,normalize = T)
  rowname=rownames(res)
  colname=colnames(res)
  covar=as.matrix(res)
  rownames(covar)=rowname
  colnames(covar)=colname
  
  return(covar)
}
scc1_covar=cell_compute(scc_sc_avg, matrix1)
scc2_covar=cell_compute(scc_sc_avg, matrix2)
scc3_covar=cell_compute(scc_sc_avg, matrix3)

#delete the noisy celltype
scc1_covar=scc1_covar[!(rownames(scc1_covar) %in% "Multiplet"),]
scc2_covar=scc2_covar[!(rownames(scc2_covar) %in% "Multiplet"),]
scc3_covar=scc3_covar[!(rownames(scc3_covar) %in% "Multiplet"),]

# retain cell types with sufficiently large proportions across spots
select_cell=function(covar,probs1=0.9){
  covar <- to.proportion(covar)
  quantile_90 <- apply(covar, 1, function(x) quantile(x, probs = probs1))
  selected_cell <- rownames(covar[quantile_90 > 0.1, ])
  return(selected_cell)
}
union_cell=Reduce(union,list(select_cell(scc1_covar,probs1 = 0.9),select_cell(scc2_covar,probs1 = 0.9),select_cell(scc3_covar,probs1 = 0.9)))

scc1_covar=round(t(to.proportion(scc1_covar[union_cell,])),4)
scc2_covar=round(t(to.proportion(scc2_covar[union_cell,])),4)
scc3_covar=round(t(to.proportion(scc3_covar[union_cell,])),4)


##3. variable selection and data cleaning
scc1_filtered_list=genefilter(matrix1,position1,scc1_covar,num_features=100,num_cells=10,nfeature=8000)
scc2_filtered_list=genefilter(matrix2,position2,scc2_covar,num_features=100,num_cells=10,nfeature=8000)
scc3_filtered_list=genefilter(matrix3,position3,scc3_covar,num_features=100,num_cells=10,nfeature=8000)

dim(scc1_filtered_list[[1]])
dim(scc2_filtered_list[[1]])
dim(scc3_filtered_list[[1]])

totalgene=Reduce(intersect,list(rownames(scc1_filtered_list[[1]]),rownames(scc2_filtered_list[[1]]),rownames(scc3_filtered_list[[1]])))
length(totalgene)

scc1_filtered_list[[1]]=scc1_filtered_list[[1]][totalgene, , drop = FALSE]
scc2_filtered_list[[1]]=scc2_filtered_list[[1]][totalgene, , drop = FALSE]
scc3_filtered_list[[1]]=scc3_filtered_list[[1]][totalgene, , drop = FALSE]
sample_total=list(scc1_filtered_list,scc2_filtered_list,scc3_filtered_list)

#4.output the result
save(sample_total,file=here::here("data/Realdataset/raw_data/scc_dataset","scc_P9_new.RData"))

write.csv(scc1_filtered_list[[1]], file = here::here("data/Realdataset/raw_data/scc_dataset","matrix1_count_scc.csv"), row.names = TRUE) 
write.csv(scc1_filtered_list[[2]], file = here::here("data/Realdataset/raw_data/scc_dataset","matrix1_position_scc.csv"), row.names = TRUE) 
write.csv(scc1_filtered_list[[3]], file = here::here("data/Realdataset/raw_data/scc_dataset","matrix1_celltype_scc.csv"), row.names = TRUE) 
write.csv(scc2_filtered_list[[1]], file = here::here("data/Realdataset/raw_data/scc_dataset","matrix2_count_scc.csv"), row.names = TRUE) 
write.csv(scc2_filtered_list[[2]], file = here::here("data/Realdataset/raw_data/scc_dataset","matrix2_position_scc.csv"), row.names = TRUE) 
write.csv(scc2_filtered_list[[3]], file = here::here("data/Realdataset/raw_data/scc_dataset","matrix2_celltype_scc.csv"), row.names = TRUE) 
write.csv(scc3_filtered_list[[1]], file = here::here("data/Realdataset/raw_data/scc_dataset","matrix3_count_scc.csv"), row.names = TRUE) 
write.csv(scc3_filtered_list[[2]], file = here::here("data/Realdataset/raw_data/scc_dataset","matrix3_position_scc.csv"), row.names = TRUE) 
write.csv(scc3_filtered_list[[3]], file = here::here("data/Realdataset/raw_data/scc_dataset","matrix3_celltype_scc.csv"), row.names = TRUE) 















