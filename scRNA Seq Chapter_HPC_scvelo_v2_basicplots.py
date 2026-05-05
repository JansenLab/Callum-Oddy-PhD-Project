import scvelo as scv
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

# -----------------------------------------------------------------------------
# 1. Setup and Data Loading
# -----------------------------------------------------------------------------
out_dir = "scvelo_results/final_plots/"
os.makedirs(out_dir, exist_ok=True)
scv.settings.figdir = out_dir
scv.settings.set_figure_params('scvelo', dpi=300, format='png')

print("Loading integrated data...")
adata = sc.read_h5ad("scvelo_results/FULLY_INTEGRATED_velocity.h5ad")

# -----------------------------------------------------------------------------
# 2. Robust Category and Color Handling
# -----------------------------------------------------------------------------
# Clean up labels and force specific order/colors
adata.obs['tissue_group'] = adata.obs['tissue_group'].astype(str).str.strip()
target_order = ['Gastric', 'IM', 'Duodenal']

adata.obs['tissue_group'] = pd.Categorical(
    adata.obs['tissue_group'],
    categories=target_order,
    ordered=True
)

# Muted Red, Muted Purple, Muted Blue
adata.uns['tissue_group_colors'] = ['#d65f5f', '#9575cd', '#4e79a7']
print(f"Categories successfully set: {adata.obs['tissue_group'].cat.categories.tolist()}")

# -----------------------------------------------------------------------------
# 3. Recompute Neighbors and Moments
# -----------------------------------------------------------------------------
print("Recomputing neighbours and moments in scvelo format...")
scv.pp.neighbors(adata)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)


# -----------------------------------------------------------------------------
# PART A: Interpreting Velocities (Phase Portraits)
# -----------------------------------------------------------------------------
print("Plotting sample phase portraits...")
sample_genes = ['MUC5AC', 'MUC2', 'CDX2', 'TFF3']
valid_genes = [g for g in sample_genes if g in adata.var_names]

if len(valid_genes) > 0:
    print(f"Found genes: {valid_genes}")
    # Phase Portraits
    scv.pl.velocity(adata, var_names=valid_genes, color='tissue_group',
                    ncols=2, basis='phate', save='phase_portraits.png')

    # Scatter Plots - Plotted individually to avoid the Pandas 'unique' list bug
    for gene in valid_genes:
        scv.pl.scatter(adata, basis=gene, color='tissue_group',
                       save=f'scatter_{gene}_tissue.png')
        scv.pl.scatter(adata, basis=gene, color='velocity',
                       save=f'scatter_{gene}_velocity.png')
else:
    print("Warning: None of the sample genes were found in the dataset.")

# -----------------------------------------------------------------------------
# PART B: Rank Velocity Genes (Find Lineage Drivers)
# -----------------------------------------------------------------------------
print("Ranking velocity genes per lineage...")
scv.tl.rank_velocity_genes(adata, groupby='tissue_group', min_corr=.3)

df_ranked = pd.DataFrame(adata.uns['rank_velocity_genes']['names'])
df_ranked.to_csv(f"{out_dir}Top_Velocity_Driving_Genes.csv", index=False)

# -----------------------------------------------------------------------------
# PART C: Speed and Coherence
# -----------------------------------------------------------------------------
print("Calculating velocity speed and confidence...")
scv.tl.velocity_confidence(adata)

# Plotting metrics individually to bypass the Pandas TypeError in scvelo
scv.pl.scatter(adata, c='velocity_length', cmap='coolwarm',
               perc=[5, 95], basis='phate', save='velocity_speed.png')

scv.pl.scatter(adata, c='velocity_confidence', cmap='coolwarm',
               perc=[5, 95], basis='phate', save='velocity_confidence.png')

# -----------------------------------------------------------------------------
# PART D: Velocity Graph & Pseudotime
# -----------------------------------------------------------------------------
print("Calculating velocity pseudotime...")
scv.tl.velocity_pseudotime(adata)

scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot',
               basis='phate', save='velocity_pseudotime_phate.png')

# -----------------------------------------------------------------------------
# PART E: PAGA Velocity Graph (Trajectory Abstraction)
# -----------------------------------------------------------------------------
print("Running directed PAGA...")
# Sync neighbors structure for PAGA bugfix
if 'neighbors' in adata.uns:
    adata.uns['neighbors']['distances'] = adata.obsp['distances']
    adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='tissue_group')

# Export transition matrix
df_paga = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
df_paga.to_csv(f"{out_dir}PAGA_Transition_Confidence.csv")

# Superimpose PAGA on PHATE
scv.pl.paga(adata, basis='phate', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5,
            save='directed_PAGA_phate.png')

# -----------------------------------------------------------------------------
# Save final results
# -----------------------------------------------------------------------------
print("Saving final annotated object...")
adata.write("scvelo_results/SCVELO_ADVANCED_METRICS.h5ad")
print("--- ALL ADVANCED PLOTS COMPLETE ---")





