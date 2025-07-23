!pip install scanpy
import scanpy as sc
import pandas as pd

adata = sc.read('/content/adata.h5ad')
adata

adata.var_names
adata.obs_names

# Export expression matrix (X)
expr_df = pd.DataFrame(
    seurat1.X.toarray() if hasattr(seurat1.X, "toarray") else seurat1.X,
    index=seurat1.obs_names,
    columns=seurat1.var_names
)

expr_df.to_csv('/content/expression_matrix.csv', compression='gzip')

# Export cell metadata (obs)
seurat1.obs.to_csv('/content/cell_metadata.csv')

gene_names = seurat1.var_names
print(gene_names)

gene_df = pd.DataFrame(gene_names)
print (gene_df)

#gene_df.to_csv("/content/gene_names_U5_hNSC_velocity_9_7_2022.csv", index=False)

# basic scanpy pre-processing

# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

sc.pl.umap(
    adata,
    color="cell_cycle_position",
    size=20
)


