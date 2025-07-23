
```
docker run -it -v /media/old_home/home/rlhoove1/U5:/files cplaisier/scrna_seq_velocity

pip3 install --upgrade pip
pip3 install --upgrade scvelo
pip install webcolors
pip install tqdm
pip install ipywidgets
python3
```
import scanpy as sc
import scvelo as scv
import pandas as pd
import numpy as np
import webcolors
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt

# Load Seurat object
adata_Seurat = scv.read('/files/U5_normalized_ensembl.h5ad')

print("Seurat:", adata_Seurat.obs_names[:5].tolist())

#set to 14!!!!!!!!!!!
scv.utils.clean_obs_names(adata_Seurat, id_length= 14)
# n_vars = 2962 \xd7 3000
print("Seurat:", adata_Seurat.obs_names[:5].tolist())


# Load Velocyto loom
adata_vc = scv.read_loom('/files/U5_velocyto.loom')
print("VC:", adata_vc.obs_names[:5].tolist())

adata_vc.obs_names = adata_vc.obs_names.astype(str)
scv.utils.clean_obs_names(adata_vc,id_length= 14)
print("VC:", adata_vc.obs_names[:5].tolist())

#get ensemble
adata_vc.var_names = adata_vc.var['Accession']
adata_vc.var = adata_vc.var.drop(columns=['Accession'])

print("VC genes:", adata_vc.var_names[:5].tolist())

adata_vc
#n_obs \xd7 n_vars = 3049 \xd7 32738

print("Seurat object genes:", adata_Seurat.var_names.shape[0])
print("Velocyto object genes:", adata_vc.var_names.shape[0])
shared_genes = adata_Seurat.var_names.intersection(adata_vc.var_names)
print("Number of shared genes:", len(shared_genes))

# Merge
adata = scv.utils.merge(adata_vc, adata_Seurat)

# Format and re-order cell cycle factor
adata.obs['ccAF'] = adata.obs['ccAF'].astype(str)
adata.obs['ccAF'] = adata.obs['ccAF'].astype('category')
adata.obs['ccAF'] = adata.obs['ccAF'].cat.reorder_categories(['Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1'], ordered=True)

#
scv.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
#

#scv.pp.neighbors(adata, n_pcs=19, n_neighbors=15)
#scv.pp.moments(adata, n_pcs=19, n_neighbors=15)


scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical', groupby='ccAF', groups=['Neural G0', 'G1','Late G1','S','S/G2','G2/M','M/Early G1'], groups_for_fit=['Neural G0', 'G1','Late G1','S','S/G2','G2/M','M/Early G1'])
scv.tl.velocity_graph(adata)

adata
# n_obs \xd7 n_vars = 2962 \xd7 3000

print(adata.obsm.keys())

# Get same colors for each cell cycle phase
from scipy.spatial import KDTree
hexnames = webcolors.CSS3_HEX_TO_NAMES
names = []
positions = []
for hex, name in hexnames.items():
    names.append(name)
    positions.append(webcolors.hex_to_rgb(hex))

spacedb = KDTree(positions)
result = []
for querycolor in [(123,175,65),(201,149,43),(243,118,110),(31,189,194),(166,129,186),(225,109,170),(43,181,103),(74,161,217)]:
    dist, index = spacedb.query(querycolor)
    result.append(names[index])


cmap1 = dict(zip(['G1/other', 'Neural G0', 'G1', 'Late G1', 'S', 'S/G2', 'G2/M', 'M/Early G1'],result))

# Plot stream lines
scv.pl.velocity_embedding_stream(adata, basis='X_umap', color=['ccAF'], save='/files/output/adata_velocity_stream_umap_9_7.png', dpi=300, palette = cmap1, title='', legend_loc='right margin', figsize=(3,3), size=70)


#get velocity_pca for cellpath to use later 
scv.pl.velocity_embedding_stream(adata, basis='X_pca', color=['ccAF'], save='/files/output/adata_velocity_stream_pca_9_7.png', dpi=300, palette = cmap1, title='', legend_loc='right margin', figsize=(3,3), size=70)

keepers = adata.var.index[pd.DataFrame(adata.layers['velocity']).T.dropna().index]
adata2 = adata[:,keepers]

adata2

adata

adata2.write('/files/U5_hNSC_keepers_velocity_july_23_2025.h5ad', compression='gzip')
adata.write('/files/U5_hNSC_velocity_july_23_2025.h5ad', compression='gzip')
