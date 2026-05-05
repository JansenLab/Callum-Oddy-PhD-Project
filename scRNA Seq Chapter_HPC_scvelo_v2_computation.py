import scanpy as sc
import scvelo as scv
import os
import gc

if __name__ == '__main__':

    scv.settings.verbosity = 3
    out_dir = "scvelo_results/"
    os.makedirs(out_dir, exist_ok=True)

    print("1. Loading raw layers for computation...")
    # This file should load fine as it was likely created with standard settings
    adata_velo = sc.read_h5ad("scvelo_raw_integrated.h5ad")

    print(f"Original shape: {adata_velo.shape}")

    print("2. Forcing memory copies to prevent read-only crashes...")
    # This ensures the math libraries have full write-access to the data
    adata_velo.layers['spliced'] = adata_velo.layers['spliced'].copy()
    adata_velo.layers['unspliced'] = adata_velo.layers['unspliced'].copy()
    gc.collect()

    print("3. Preprocessing velocity object...")
    # Standard scVelo pipeline: Filter to genes with enough spliced/unspliced signal
    scv.pp.filter_genes(adata_velo, min_shared_counts=20)
    scv.pp.normalize_per_cell(adata_velo)

    # Ensure gene names are unique and strings
    adata_velo.var_names_make_unique()
    adata_velo.var_names = adata_velo.var_names.astype(str)

    # Focus on the top 2000 genes to make the dynamical model feasible
    sc.pp.highly_variable_genes(adata_velo, n_top_genes=2000, subset=True)
    sc.pp.log1p(adata_velo)

    print("4. Computing PCA, neighbors, and moments...")
    sc.pp.pca(adata_velo)
    sc.pp.neighbors(adata_velo, n_pcs=50, n_neighbors=30)
    scv.pp.moments(adata_velo)

    print("5. Running recover_dynamics (The Heavy Lifting)...")
    # Using 4 cores since we are in the stable 'scvelo' environment
    scv.tl.recover_dynamics(adata_velo, n_jobs=4)

    print("6. Computing velocity and latent time...")
    scv.tl.velocity(adata_velo, mode='dynamical')
    scv.tl.velocity_graph(adata_velo)
    scv.tl.latent_time(adata_velo)

    print("7. Saving raw velocity results...")
    # We save this as a standalone file.
    # We will merge it with your PHATE coordinates in the next step.
    adata_velo.write(f"{out_dir}RAW_velocity_computed.h5ad")
    print("--- COMPUTATION COMPLETE ---")


