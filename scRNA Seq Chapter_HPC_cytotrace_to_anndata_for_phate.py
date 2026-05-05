import scanpy as sc
import anndata
import phate
import matplotlib.pyplot as plt
import os

# Fix for nullable strings error
anndata.settings.allow_write_nullable_strings = True

# 1. Load the converted AnnData
adata = sc.read_h5ad("seu_integrated_final.h5ad")

print(f"Available embeddings: {adata.obsm.keys()}")

# Check for harmony key (your R script maps it to X_harmony)
if 'X_harmony' in adata.obsm.keys():
    harmony_key = 'X_harmony'
elif 'harmony' in adata.obsm.keys():
    harmony_key = 'harmony'
else:
    raise KeyError(f"No harmony embedding found. Available: {adata.obsm.keys()}")

print(f"Running PHATE on {harmony_key} coordinates...")
phate_op = phate.PHATE(n_jobs=-1, verbose=True)
adata.obsm['X_phate'] = phate_op.fit_transform(adata.obsm[harmony_key])

custom_palette = {
    "GAS": "#d57a7a",
    "IM": "#9575cd",
    "DUO": "#5c92e8"
}

out_dir = "phate_final_results/"
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
sc.settings.figdir = out_dir

sc.pl.embedding(adata, basis="phate", color="tissue_group",
                palette=custom_palette, frameon=False, show=False)
plt.savefig(f"{out_dir}PHATE_Tissue_Integrated.png", dpi=300, bbox_inches='tight')

sc.pl.embedding(adata, basis="phate", color="CytoTRACE2_Score",
                color_map="magma", frameon=False, show=False)
plt.savefig(f"{out_dir}PHATE_CytoTRACE_Score.png", dpi=300, bbox_inches='tight')

sc.pl.embedding(adata, basis="phate", color="final_annotation",
                color_map="magma", frameon=False, show=False)
plt.savefig(f"{out_dir}PHATE_final_annotations.png", dpi=300, bbox_inches='tight')

# Fix for nullable strings — already set at top
adata.write(f"{out_dir}final_integrated_adata.h5ad")
print("--- PHATE ANALYSIS COMPLETE ---")
