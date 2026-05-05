import scanpy as sc
import scvelo as scv
import anndata
import matplotlib.pyplot as plt
import os
import pandas as pd
import numpy as np

# —— THE FIX FOR PANDAS 3.0 + SCVELO ———————————————————————————
# This forces Pandas to accept standard lists from scVelo without crashing
orig_unique = pd.unique
def patched_unique(values):
    if isinstance(values, list):
        values = np.array(values)
    return orig_unique(values)
pd.unique = patched_unique
# ——————————————————————————————————————————————————————————————————————————————

# THE MULTIPROCESSING GUARD
if __name__ == '__main__':

    anndata.settings.allow_write_nullable_strings = True

    scv.settings.verbosity = 3
    scv.settings.set_figure_params('scvelo')

    out_dir = "scvelo_results/"
    os.makedirs(out_dir, exist_ok=True)
    scv.settings.figdir = out_dir

    # —— 1. Load both objects ——————————————————————————————————————————————————————
    print("Loading AnnData objects...")
    adata_final  = sc.read_h5ad("phate_final_results/final_integrated_adata.h5ad")
    adata_layers = sc.read_h5ad("scvelo_raw_integrated.h5ad")

    print(f"adata_final:  {adata_final.shape}")
    print(f"adata_layers: {adata_layers.shape}")

    # —— 2. Find common cells and genes ———————————————————————————————————————————
    common_cells = adata_final.obs_names.intersection(adata_layers.obs_names)
    common_genes = adata_final.var_names.intersection(adata_layers.var_names)

    print(f"Common cells: {len(common_cells)}")
    print(f"Common genes: {len(common_genes)}")

    if len(common_cells) == 0:
        raise ValueError("No common cells found — check barcode formatting!")
    if len(common_genes) == 0:
        raise ValueError("No common genes found — check var_names!")

    adata_final  = adata_final[common_cells, common_genes].copy()
    adata_layers = adata_layers[common_cells, common_genes].copy()

    # —— 3. Transfer spliced/unspliced layers ——————————————————————————————————————
    print("Transferring spliced/unspliced layers...")
    adata_final.layers['spliced']   = adata_layers.layers['spliced']
    adata_final.layers['unspliced'] = adata_layers.layers['unspliced']

    # —— 3b. Proportions plot ——————————————————————————————————————————————————————
    print("Plotting spliced/unspliced proportions...")
    scv.pl.proportions(adata_final, groupby='tissue_group', show=False)
    plt.savefig(f"{out_dir}spliced_unspliced_proportions_by_tissue.png", dpi=300, bbox_inches='tigh$
    plt.close()

    if 'final_annotations' in adata_final.obs.columns:
        scv.pl.proportions(adata_final, groupby='final_annotations', show=False)
        plt.savefig(f"{out_dir}spliced_unspliced_proportions_by_annotation.png", dpi=300, bbox_inch$
        plt.close()

    scv.pl.proportions(adata_final, show=False)
    plt.savefig(f"{out_dir}spliced_unspliced_proportions_overall.png", dpi=300, bbox_inches='tight')
    plt.close()

    # —— 4. Set tissue colours —————————————————————————————————————————————————————
    tissue_colors = {'GAS': '#d65f5f', 'IM': '#9575cd', 'DUO': '#4e79a7'}
    if 'tissue_group' in adata_final.obs.columns:
        adata_final.uns['tissue_group_colors'] = [
            tissue_colors[cat] for cat in adata_final.obs['tissue_group'].cat.categories
        ]

    # —— 5. Preprocessing ——————————————————————————————————————————————————————————
    print("Preprocessing...")
    scv.pp.filter_genes(adata_final, min_shared_counts=20)
    scv.pp.normalize_per_cell(adata_final)

    # Ensure var_names are unique and strings
    adata_final.var_names_make_unique()
    adata_final.var_names = adata_final.var_names.astype(str)

    sc.pp.highly_variable_genes(adata_final, n_top_genes=2000, subset=True)
    sc.pp.log1p(adata_final)

    # Explicitly compute neighbors before moments to resolve DeprecationWarning
    print("Computing PCA and neighbors...")
    sc.pp.pca(adata_final)
    sc.pp.neighbors(adata_final, n_pcs=50, n_neighbors=30)

    print("Computing moments...")
    scv.pp.moments(adata_final)

    # —— 6. Dynamical model ————————————————————————————————————————————————————————
    print("Running recover_dynamics (this will take a while)...")
    scv.tl.recover_dynamics(adata_final, n_jobs=1, show_progress_bar=False)

    print("Computing velocity...")
    scv.tl.velocity(adata_final, mode='dynamical')
    scv.tl.velocity_graph(adata_final)

    # —— 7. Latent time ————————————————————————————————————————————————————————————
    print("Computing latent time...")
    scv.tl.latent_time(adata_final)
    scv.tl.terminal_states(adata_final)

    # —— 8. Embed velocities onto PHATE ————————————————————————————————————————————
    print("Embedding velocities onto PHATE...")
    scv.tl.velocity_embedding(adata_final, basis='phate')
    # —— 9. Plots ——————————————————————————————————————————————————————————————————
    print("Generating plots...")

    # Stream plot coloured by tissue
    scv.pl.velocity_embedding_stream(
        adata_final, basis='phate',
        color='tissue_group',
        frameon=False, show=False,
        title='RNA Velocity (PHATE) — Tissue'
    )
    plt.savefig(f"{out_dir}velocity_stream_tissue.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Stream plot coloured by CytoTRACE2
    if 'CytoTRACE2_Score' in adata_final.obs.columns:
        scv.pl.velocity_embedding_stream(
            adata_final, basis='phate',
            color='CytoTRACE2_Score',
            cmap='magma',
            frameon=False, show=False,
            title='RNA Velocity (PHATE) — CytoTRACE2'
        )
        plt.savefig(f"{out_dir}velocity_stream_cytotrace.png", dpi=300, bbox_inches='tight')
        plt.close()

    # Latent time
    scv.pl.scatter(
        adata_final, basis='phate',
        color='latent_time',
        color_map='gnuplot',
        frameon=False, show=False,
        title='Latent Time (PHATE)'
    )
    plt.savefig(f"{out_dir}latent_time_phate.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Root and end point cells
    scv.pl.scatter(
        adata_final, basis='phate',
        color=['root_cells', 'end_points'],
        frameon=False, show=False,
    )
    plt.savefig(f"{out_dir}root_end_points_phate.png", dpi=300, bbox_inches='tight')
    plt.close()

    # —— 10. Save ——————————————————————————————————————————————————————————————————
    print("Saving final object...")
    adata_final.write(f"{out_dir}final_integrated_with_velocity.h5ad")
    print("--- scVELO COMPLETE ---")



