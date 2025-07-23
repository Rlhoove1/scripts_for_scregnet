```
docker run -it -v /home/rlhoove1/endocrin:/files cplaisier/scrna_seq_velocity

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



scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.settings.set_figure_params('scvelo')


adata = scv.datasets.pancreas()

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(adata)

scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)

adata.write('/files/endo_velocity_july_14_2025.h5ad', compression='gzip')

scv.pl.velocity_embedding_stream(adata, basis='umap',save='/files/endo_velocity_stream_umap_july_22.png')

scv.pl.velocity_embedding_stream(adata, basis='pca',save='/files/endo_velocity_stream_umap_july_22.png')

keepers = adata.var.index[pd.DataFrame(adata.layers['velocity']).T.dropna().index]
adata2 = adata[:,keepers]

adata2.write('/files/endo_keepers_velocity_july_22_2025.h5ad', compression='gzip')
