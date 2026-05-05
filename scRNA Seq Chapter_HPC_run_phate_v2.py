import scanpy as sc
import pandas as pd
import phate
import matplotlib.pyplot as plt
import os
from scipy.io import mmread

# Setup folders immediately
bridge_dir = "celltypist_input/bridge_files/"
out_dir = "celltypist_input/final_outputs/"
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

print("--- PYTHON PHASE START ---")
print("Rebuilding AnnData...")

X = mmread(f"{bridge_dir}matrix.mtx").T.tocsr()
genes = pd.read_table(f"{bridge_dir}genes.tsv", header=None)[0].values
metadata = pd.read_csv(f"{bridge_dir}metadata.csv", index_col='barcode')

adata = sc.AnnData(X=X, obs=metadata)
adata.var_names = genes

# SMART COORDINATE PICKER
# This looks for UMAP_1, umap_1, or UMAP1
def get_coords(df, prefix):
    cols = [c for c in df.columns if c.lower().replace("_", "") ==
f"{prefix.lower()}1" or
                                     c.lower().replace("_", "") ==
f"{prefix.lower()}2"]
    return df[sorted(cols)].values

try:
    adata.obsm['X_umap'] = get_coords(metadata, "UMAP")
    print("Found UMAP coordinates.")
except:
    print("Warning: Could not find UMAP coordinates. Check metadata.csv columns.")

# Extract PCA (usually PC_1, PC_2...)
pca_cols = [c for c in metadata.columns if c.startswith('PC_') or c.startswith('PC')]
adata.obsm['X_pca'] = metadata[pca_cols].values

# Run PHATE
print("Running PHATE (this is the heavy lifting)...")
phate_op = phate.PHATE(n_jobs=-1, verbose=True)
adata.obsm['X_phate'] = phate_op.fit_transform(adata.obsm['X_pca'])

# Save and Plot
print(f"Saving plots to {out_dir}...")
sc.settings.figdir = out_dir
genes_to_plot = ["CDX2", "MUC2", "SOX2", "NKX6-3", "MUC6"]
sc.pl.embedding(adata, basis="phate", color="final_annotation", show=False,
frameon=False)
plt.savefig(f"{out_dir}PHATE_Identity.png", dpi=300, bbox_inches='tight')

sc.pl.embedding(adata, basis="phate", color=genes_to_plot, show=False,
color_map="Reds", frameon=False)
plt.savefig(f"{out_dir}PHATE_Markers.png", dpi=300, bbox_inches='tight')

adata.write(f"{out_dir}master_final_with_phate.h5ad")
print("--- PYTHON PHASE COMPLETE ---")
# 1. Define your custom muted palette for tissue groups
# Gastric: Muted Red, IM: Purple, Duodenal: Blue
tissue_palette = {
    "GAS": "#d65f5f",      # Muted Red
    "IM": "#9575cd",       # Muted Purple
    "DUO": "#4fc3f7"       # Muted Blue
}

# 2. Plotting PHATE colored by CytoTRACE 2 Score
print("Plotting CytoTRACE results on PHATE...")
sc.pl.embedding(
    adata,
    basis="phate",
    color="CytoTRACE2_Score",
    color_map="magma", # 'magma' or 'viridis' work well for potency gradients
    frameon=False,
    show=False)

plt.savefig(f"{out_dir}PHATE_CytoTRACE_Score.png", dpi=300, bbox_inches='tight')

# 3. Plotting PHATE colored by Tissue Group with your specific palette
# Assuming 'tissue_group' contains GAS, IM, and DUO
sc.pl.embedding(
    adata,
    basis="phate",
    color="tissue_group",
    palette=tissue_palette,
    frameon=False,
    show=False
)
plt.savefig(f"{out_dir}PHATE_Tissue_Group.png", dpi=300,
bbox_inches='tight')

# 4. Plotting Potency Categories (Categorical)
sc.pl.embedding(
    adata,
    basis="phate",
    color="CytoTRACE2_Potency",
    frameon=False,
    show=False
)
plt.savefig(f"{out_dir}PHATE_Potency_Class.png", dpi=300,
bbox_inches='tight')



