####################################################################################################################################
                                                  #comparison_methods in python

####################################################################################################################################
#you should define the root library before processing the code
from pathlib import Path
import scanpy as sc

PROJECT_ROOT = Path(__file__).resolve().parent
data_dir = PROJECT_ROOT / "Demo" / "example"
####################################################################################################################################
                                                          #StarTrail

####################################################################################################################################
##1.StarTrail: the implement of curve fitting. Other details can be seen in file "Comparison_methods_R".
from sklearn.datasets import make_blobs
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import sys
import numpy as np
import pandas as pd
import torch
import os
from scipy.io import loadmat
from math import floor
import numpy as np
from sklearn import metrics
from sklearn.cluster import DBSCAN
import elpigraph

for datanum in range(1,5):
  for cellnumber in range(1,7):
    gradient_mesh = pd.read_csv(f"{datanum}_gradient_gene{cellnumber}.txt",sep='\t')
    print(f"Processing gene {cellnumber}", flush=True)
    coord_mesh = gradient_mesh[['s1', 's2']].to_numpy()
    pred = np.array(gradient_mesh['pred'])
    L2 = np.array((gradient_mesh.g1**2+gradient_mesh.g2**2)**0.5)
    index_keep = (L2 > np.quantile(L2, 0.85)) # use 0.9 quantile
    coord_mesh_grad = coord_mesh[index_keep, :]
    db = DBSCAN(eps=0.03, min_samples=1).fit(coord_mesh_grad)
    labels = db.labels_
    labels_unique, counts = np.unique(labels, return_counts=True)
    print("Cluster sizes:", counts)
    print("Max cluster size:", counts.max())
    min_num_spot =200
    labels_to_keep = labels_unique[counts > min_num_spot]
    gene_boundary = []
    for i in range(labels_to_keep.shape[0]):
      print(f"  Cluster {i}, n = {counts[counts > min_num_spot][i]}", flush=True)
      coord_final = coord_mesh_grad[labels==labels_to_keep[i],]
      coord_final[:,1] = -coord_final[:,1]
      PG = elpigraph.computeElasticPrincipalCurve(coord_final,NumNodes=int(counts[counts > min_num_spot][i]/50))[0]
      nodesp = PG["NodePositions"]
      Edges = PG["Edges"][0].T
      points_all = []
      for j in range(Edges.shape[1]):
        x_coo = np.concatenate((nodesp[Edges[0, j], [0]], nodesp[Edges[1, j], [0]]))
        y_coo = np.concatenate((nodesp[Edges[0, j], [1]], nodesp[Edges[1, j], [1]]))
        points_all.append(np.concatenate((x_coo,y_coo),axis=0))
        print(f"Cell {cellnumber} finished\n", flush=True)
      points_all = np.array(points_all)
      gene_boundary.append(gene_boundary)
      np.savetxt(f'{datanum}_cell{cellnumber}.txt',points_all)
        
        
        
####################################################################################################################################
                                                     #GASTON

####################################################################################################################################        
##2.GASTON: Avaiable on https://github.com/raphael-group/GASTON.
#(1)import the lib
import os
import numpy as np
import matplotlib.pyplot as plt
from glmpca import glmpca
from importlib import reload
import gaston
from gaston import neural_net,cluster_plotting, dp_related, segmented_fit, model_selection
from gaston import binning_and_plotting, isodepth_scaling, run_slurm_scripts
from gaston import spatial_gene_classification, plot_cell_types, filter_genes, process_NN_output
import pandas as pd
from torch.serialization import add_safe_globals
import torch
from gaston.neural_net import GASTON
from torch.nn.modules.container import Sequential
from torch.nn.modules.linear import Linear
from torch.nn import ReLU
from gaston import process_NN_output
add_safe_globals([GASTON])
add_safe_globals([Sequential])
add_safe_globals([Linear])
add_safe_globals([ReLU])


#(2)train the model
for t in ["count", "position"]:
    csv_file = data_dir / f"matrix1_{t}_example.csv"
    npy_file = data_dir / f"matrix1_{t}_example.npy"
    df = pd.read_csv(csv_file, index_col=0)
    np.save(npy_file, df.values)
    print(f"Saved {npy_file}")

csv_file = data_dir / "matrix1_count_example.csv"
npy_file = data_dir / "genename_example.npy"
matrix_df = pd.read_csv(csv_file, index_col=0)
colnames = matrix_df.index
np.save(npy_file, colnames)
gene_labels = np.load(npy_file, allow_pickle=True)  # array of names of G genes

path_to_glmpca = data_dir / "matrix1_glmpca_example.npy"
path_to_coords = data_dir / "matrix1_position_example.npy"
path_to_count = data_dir / "matrix1_count_example.npy"

coords_mat = np.load(path_to_coords)  # N x 2 spatial coordinate matrix
counts_mat = np.load(path_to_count)   # N x G UMI count array
gene_labels = np.load(npy_file, allow_pickle=True)

num_dims = 8  # 2 * number of clusters
penalty = 10  # may need to increase if too small
glmpca_res = glmpca.glmpca(counts_mat.T, num_dims, fam="poi", penalty=penalty, verbose=True)
A = glmpca_res['factors']  # N x num_dims
np.save(path_to_glmpca, A)
A=np.load(path_to_glmpca) 
S=np.load(path_to_coords)
S_torch, A_torch = neural_net.load_rescale_input_data(S,A)
isodepth_arch=[20,20] # architecture for isodepth neural network d(x,y) : R^2 -> R 
expression_fn_arch=[20,20] # architecture for 1-D expression function h(w) : R -> R^G
num_epochs = 2500 # number of epochs to train NN (NOTE: it is sometimes beneficial to train longer)
checkpoint = 500 # save model after number of epochs = multiple of checkpoint
out_dir= data_dir /"example_output" # folder to save model runs
optimizer = "adam"
num_restarts=5
device='cpu' # change to 'cpu' if you don't have a GPU
seed_list=range(num_restarts)
for seed in seed_list:
  print(f'training neural network for seed {seed}')
  out_dir_seed=f"{out_dir}/rep{seed}"
  os.makedirs(out_dir_seed, exist_ok=True)
  mod, loss_list = neural_net.train(S_torch, A_torch,
                                    S_hidden_list=isodepth_arch, A_hidden_list=expression_fn_arch, 
                                    epochs=num_epochs, checkpoint=checkpoint, device=device,
                                    save_dir=out_dir_seed, optim=optimizer, seed=seed, save_final=True)


#(3)process the result
gaston_model, A, S = process_NN_output.process_files(
  data_dir /"example_output""
)
  
coords_mat = np.load(data_dir / "matrix1_position_example.npy")  # N x 2 spatial coordinate matrix 
counts_mat = np.load(data_dir / "matrix1_count_example.npy").T   # N x G UMI count array
gene_labels = np.load(data_dir / "genename_example.npy", allow_pickle=True)  # array of names of G genes

num_layers=4
gaston_isodepth, gaston_labels=dp_related.get_isodepth_labels(gaston_model,A,S,num_layers)
gaston_isodepth=isodepth_scaling.adjust_isodepth(gaston_isodepth, gaston_labels, coords_mat, 
                                                 q_vals=[0.05, 0.05, 0.05, 0.05], scale_factor=0.05)

np.savetxt(data_dir / "gaston_labels.txt", gaston_labels)

umi_threshold = 0
pw_fit_dict=segmented_fit.pw_linear_fit(counts_mat, gaston_labels, gaston_isodepth,
                                        None, [], umi_threshold = umi_threshold)
binning_output=binning_and_plotting.bin_data(counts_mat, gaston_labels, gaston_isodepth, 
                                             None, gene_labels, num_bins_per_domain=[4,4,4,4], umi_threshold=umi_threshold)

q=0.9 # use 0.9 quantile for slopes, discontinuities
discont_genes_layer=spatial_gene_classification.get_discont_genes(pw_fit_dict, binning_output,q=q)
cont_genes_layer=spatial_gene_classification.get_cont_genes(pw_fit_dict, binning_output,q=q)

cont_svg = list(cont_genes_layer.keys())
discont_svg = list(discont_genes_layer.keys())
all_svg = list(set(cont_svg + discont_svg))

all_svg_sorted = sorted(all_svg)#the identified sv genes.
with open(data_dir / "GASTON_SVG_list.txt", "w", encoding="utf-8") as f:
  for gene in all_svg_sorted:
    f.write(gene + "\n")


####################################################################################################################################
                                                           #PASTE

####################################################################################################################################
##3. PASTE：integrating multi-sample datasets for single-sample method. Avaiable on https://github.com/raphael-group/paste.
import pandas as pd
import numpy as np
import scanpy as sc
import paste as pst
import rpy2.robjects as ro
from rpy2.robjects import r
from rpy2.robjects import pandas2ri




slice_1 = sc.read_csv(data_dir / "matrix1_count_example.csv").T
slice_2 = sc.read_csv(data_dir / "matrix2_count_example.csv").T
slice_3 = sc.read_csv(data_dir / "matrix3_count_example.csv").T
slice_4 = sc.read_csv(data_dir / "matrix4_count_example.csv").T

slice_1_coor = pd.read_csv(data_dir / "matrix1_position_example.csv")
slice_2_coor = pd.read_csv(data_dir / "matrix2_position_example.csv")
slice_3_coor = pd.read_csv(data_dir / "matrix3_position_example.csv")
slice_4_coor = pd.read_csv(data_dir / "matrix4_position_example.csv")

slice_1_coor=slice_1_coor[['x', 'y']].to_numpy()
slice_2_coor=slice_2_coor[['x', 'y']].to_numpy()
slice_3_coor=slice_3_coor[['x', 'y']].to_numpy()
slice_4_coor=slice_4_coor[['x', 'y']].to_numpy()

slice_1.obsm['spatial'] = slice_1_coor
slice_2.obsm['spatial'] = slice_2_coor
slice_3.obsm['spatial'] = slice_3_coor
slice_4.obsm['spatial'] = slice_4_coor

sc.pp.filter_genes(slice_1, min_counts=1)
sc.pp.filter_cells(slice_1, min_counts=1)
sc.pp.filter_genes(slice_2, min_counts=1)
sc.pp.filter_cells(slice_2, min_counts=1)
sc.pp.filter_genes(slice_3, min_counts=1)
sc.pp.filter_cells(slice_3, min_counts=1)
sc.pp.filter_genes(slice_4, min_counts=1)
sc.pp.filter_cells(slice_4, min_counts=1)

initial_slice = slice_1.copy()
slices = [slice_1, slice_2, slice_3, slice_4]

lmbda = len(slices)*[1/len(slices)]
center_slice, pis = pst.center_align(initial_slice, slices, lmbda)#choose the first slice as objective slice
W = center_slice.uns['paste_W']
H = center_slice.uns['paste_H']
exp_matrix = np.dot(W,H)
df = pd.DataFrame(np.around(exp_matrix, decimals=1))
output_path = data_dir /"example_paste.csv"
df.to_csv(output_path, index=False)













        
        
        
        
