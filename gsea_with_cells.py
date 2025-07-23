import os
import pandas as pd
!pip install gseapy
import gseapy as gp
import matplotlib.pyplot as plt
!pip install scanpy
import scanpy as sc
!pip install scvelo
import scvelo as scv

#set up
seurat1 = sc.read('/content/U5_hNSC_velocity_july_11_2025.h5ad')
seurat1

trajectory_dir = "/content/traj"

filename = "cellpath_trajectories_kn_30_t_0.4_lb_0.3.csv"
traj_path = os.path.join(trajectory_dir, filename)

# Load file
traj_df = pd.read_csv(traj_path, index_col=0)
print(f"Loaded {filename}, shape = {traj_df.shape}")

# make ranked cell list 
trajectory_name = traj_df.columns[0]
rnk_series = traj_df[trajectory_name].dropna().sort_values(ascending=False)
rnk = pd.DataFrame(rnk_series)
rnk.columns = ['traj_value']


print(rnk.head())
# clip the "-1" of the end of the rownames 
rnk.index = rnk.index.str.rstrip('-1')
print(rnk.head())

# make "gene-set" in this case it actualty cells

cell_sets = seurat1.obs.groupby('ccAF').groups
gene_sets = {k: list(v) for k, v in cell_sets.items()}
#its a dictionary with cell states as keys 
print(gene_sets.keys())
print(gene_sets['Neural G0'])


# run gsea
pre_res = gp.prerank(
    rnk=rnk,
    gene_sets=gene_sets,
    min_size=1,
    max_size=3000,
    processes=4,
    permutation_num=1000,
    outdir=None,
    seed=42,
    verbose=True
)

gsea_results = pre_res.res2d
print(gsea_results.head())

# plot enrichment 

terms = pre_res.res2d.Term
axs = pre_res.plot(terms=terms[1])


terms = pre_res.res2d.Term
axs = pre_res.plot(terms=terms[3])

axs = pre_res.plot(terms=terms[1:8],

show_ranking=True,
figsize=(3,4)
)
# https://app.readthedocs.org/projects/gseapy/downloads/pdf/latest/
