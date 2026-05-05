import scanpy as sc
import scvelo as scv
import anndata as ad
import matplotlib.pyplot as plt
import os
import pandas as pd
import numpy as np

# --- GLOBAL SETTINGS ---
# This fixes the "RuntimeError: allow_write_nullable_strings"
ad.settings.allow_write_nullable_strings = True
scv.settings.verbosity = 3

if __name__ == '__main__':
    # --- SETUP ---
    out_dir = "scvelo_results/final_plots/"
    os.makedirs(out_dir, exist_ok=True)

    print("1. Loading both objects...")
    adata_final = sc.read_h5ad("phate_final_results/final_integrated_adata.h5ad")
    adata_velo = sc.read_h5ad("scvelo_results/RAW_velocity_computed.h5ad")

    print(f"Original final object shape: {adata_final.shape}")
    print(f"Velocity object shape: {adata_velo.shape}")

    # --- 2. ALIGNING DATA ---
    print("Aligning genes and cells...")
    common_genes = adata_final.var_names.intersection(adata_velo.var_names)
    common_cells = adata_final.obs_names.intersection(adata_velo.obs_names)

    adata_final = adata_final[common_cells, common_genes].copy()
    adata_velo = adata_velo[common_cells, common_genes].copy()

    # To be absolutely safe with writing, we convert the index to standard strings
    adata_final.obs_names = adata_final.obs_names.astype(str)

    print(f"New aligned shape: {adata_final.shape}")
    # --- 3. TRANSFERRING MATH & CLEANING ---
    print("Transferring velocity results...")
    adata_final.layers['velocity'] = adata_velo.layers['velocity']
    adata_final.layers['spliced'] = adata_velo.layers['spliced']
    adata_final.layers['unspliced'] = adata_velo.layers['unspliced']
    adata_final.layers['Ms'] = adata_velo.layers['Ms']
    adata_final.layers['Mu'] = adata_velo.layers['Mu']

    # Handle latent_time and reindexing
    adata_final.obs['latent_time'] = adata_velo.obs['latent_time'].reindex(adata_final.obs_names).f$
    adata_final.obs['root_cells'] = adata_velo.obs['root_cells'].reindex(adata_final.obs_names).fil$
    adata_final.obs['end_points'] = adata_velo.obs['end_points'].reindex(adata_final.obs_names).fil$

    adata_final.uns['velocity_graph'] = adata_velo.uns['velocity_graph']
    adata_final.uns['velocity_params'] = adata_velo.uns['velocity_params']

    # --- 4. VELOCITY PROJECTION ---
    print("Projecting velocity onto PHATE embeddings...")
    scv.tl.velocity_embedding(adata_final, basis='phate')

    # --- 5. COLOR SCHEME ---
    # Gastric = Muted Red, IM = Purple, Duodenal = Blue
    tissue_palette = {'GAS': '#d65f5f', 'IM': '#9575cd', 'DUO': '#4e79a7'}

    if 'tissue_group' in adata_final.obs.columns:
        adata_final.obs['tissue_group'] = adata_final.obs['tissue_group'].astype('category')
        categories = adata_final.obs['tissue_group'].cat.categories
        adata_final.uns['tissue_group_colors'] = [tissue_palette.get(c, '#808080') for c in categor$

    # --- 6. PLOTTING ---
    print("Generating final PNGs...")

    # Velocity Stream
    scv.pl.velocity_embedding_stream(
        adata_final, basis='phate', color='tissue_group',
        title='RNA Velocity (PHATE) - Tissue Lineages',
        frameon=False, show=False
    )
    plt.savefig(f"{out_dir}final_velocity_stream_phate.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Latent Time
    scv.pl.scatter(
        adata_final, basis='phate',
        color=adata_final.obs['latent_time'].values,
        color_map='gnuplot', title='Calculated Latent Time',
        frameon=False, show=False
    )
    plt.savefig(f"{out_dir}final_latent_time_phate.png", dpi=300, bbox_inches='tight')
    plt.close()

    # --- 7. SAVE ---
    print("Saving the integrated result...")
    # Double check that we are saving as standard types to avoid further HDF5 errors
    adata_final.write("scvelo_results/FULLY_INTEGRATED_velocity.h5ad")
    print("--- INTEGRATION COMPLETE ---")



