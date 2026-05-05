mport scanpy as sc
import numpy as np
import pandas as pd
import graphtools
from sklearn.preprocessing import MinMaxScaler
from scipy.interpolate import Rbf
import plotly.graph_objects as go

# ==========================================
# 1. CALCULATION FUNCTION
# ==========================================
def calculate_vr_score(adata, cluster_key='tissue_group', sample_key='patient'):
    scaler = MinMaxScaler()

    # Global Inverse Velocity (Stability Metric)
    v_mag = np.linalg.norm(adata.obsm['velocity_phate'], axis=1)
    inv_v = 1.0 / (v_mag + 1e-6)
    adata.obs['inv_velo_scaled'] = scaler.fit_transform(inv_v.reshape(-1, 1)).flatten()

    # Centrality Distance (Topographic Metric)
    adata.obs['dist_scaled'] = np.nan
    clusters = adata.obs[cluster_key].unique()

    for cluster in clusters:
        print(f"Computing Manhattan centrality for phenotype: {cluster}")
        mask = adata.obs[cluster_key] == cluster
        cell_index = adata.obs_names[mask]

        if len(cell_index) < 5:
            continue

        X_phate = adata[cell_index].obsm['X_phate']

        # Master Script Alignment: knn=50, distance='manhattan'
        G = graphtools.Graph(X_phate,
                             knn=min(50, len(cell_index)-1),
                             decay=None,
                             distance='manhattan',
                             verbose=False)

        dist_matrix = G.shortest_path(distance='data')
        median_dists = np.median(dist_matrix, axis=1)

        # Clip Q99 outliers
        q99 = np.percentile(median_dists, 99)
        cluster_median = np.median(median_dists)
        median_dists[median_dists > q99] = cluster_median

        dist_scaled = scaler.fit_transform(median_dists.reshape(-1, 1)).flatten()
        adata.obs.loc[cell_index, 'dist_scaled'] = dist_scaled
    # Create Integrated Dataframe
    df = pd.DataFrame({
        'PHATE1': adata.obsm['X_phate'][:, 0],
        'PHATE2': adata.obsm['X_phate'][:, 1],
        'Sample': adata.obs[sample_key].values,
        'Phenotype': adata.obs[cluster_key].values,
        'Annotation': adata.obs['final_annotation'].values,
        'CCAT': adata.obs['CytoTRACE2_Score'].values,
        'inv_velo_scaled': adata.obs['inv_velo_scaled'].values,
        'dist_scaled': adata.obs['dist_scaled'].values,
    }, index=adata.obs_names)

    # Master Script Median Logic
    df['invVELOmed'] = df.groupby(['Sample', 'Phenotype'])['inv_velo_scaled'].transform('median')

    # Final VR Score Formula
    df['VR_score'] = 0.9 * df['CCATmed'] + 0.1 * (df['invVELOmed'] * df['dist_scaled'])

    return df
# ==========================================
# 2. VISUALIZATION FUNCTIONS
# ==========================================
def plot_3d_landscape(df, target_phenotype, color_hex, filename):
    print(f"Plotting landscape for {target_phenotype}...")
    subset = df[df['Phenotype'] == target_phenotype].dropna(subset=['PHATE1', 'PHATE2', 'VR_score'])

    if len(subset) == 0:
        print(f"  -> Skipping {target_phenotype}: No cells found. Check spelling/case sensitivity.")
        return

    grid_x, grid_y = np.mgrid[subset['PHATE1'].min():subset['PHATE1'].max():100j,
                              subset['PHATE2'].min():subset['PHATE2'].max():100j]

    # RBF Surface Interpolation (Subsampled for Memory Safety)
    if len(subset) > 10000:
        surface_calc_subset = subset.sample(n=10000, random_state=42)
    else:
        surface_calc_subset = subset

    rbfi = Rbf(surface_calc_subset["PHATE1"], surface_calc_subset["PHATE2"], surface_calc_subset["V$
                smooth=1, kernel="thin_plate_spline")
    z_grid = rbfi(grid_x, grid_y)

    fig = go.Figure()
    fig.add_trace(go.Surface(
        x=grid_x, y=grid_y, z=z_grid,
        colorscale=[[0, 'blue'], [0.5, 'green'], [1, 'white']],
        opacity=0.8, showscale=False
    ))
    
    fig.add_trace(go.Scatter3d(
        x=subset["PHATE1"], y=subset["PHATE2"], z=subset["VR_score"] + 0.005,
        mode='markers', marker=dict(size=2, opacity=0.7, color=color_hex),
        name=f'{target_phenotype} Cells'
    ))
    fig.update_layout(title=f"VR Landscape: {target_phenotype}", scene=dict(xaxis_title="PHATE 1", $
    fig.write_html(filename)
    print(f"Saved: {filename}")

def plot_all_phenotypes_landscape(df, colors_dict, filename):
    print("Plotting Master Landscape for ALL cells...")
    subset = df.dropna(subset=['PHATE1', 'PHATE2', 'VR_score'])

    grid_x, grid_y = np.mgrid[subset['PHATE1'].min():subset['PHATE1'].max():100j,
                              subset['PHATE2'].min():subset['PHATE2'].max():100j]

    # Memory Safety Net: Subsample surface building to 10k cells
    if len(subset) > 10000:
        surface_calc_subset = subset.sample(n=10000, random_state=42)
    else:
        surface_calc_subset = subset

    rbfi = Rbf(surface_calc_subset["PHATE1"], surface_calc_subset["PHATE2"], surface_calc_subset["V$
                smooth=1, kernel="thin_plate_spline")
    z_grid = rbfi(grid_x, grid_y)

    fig = go.Figure()
    fig.add_trace(go.Surface(
        x=grid_x, y=grid_y, z=z_grid,
        colorscale=[[0, 'blue'], [0.5, 'green'], [1, 'white']],
        opacity=0.5, showscale=False
    ))
    # Plot cells for each specific color group
    for pheno, color in colors_dict.items():
        pheno_df = subset[subset['Phenotype'] == pheno]
        if not pheno_df.empty:
            fig.add_trace(go.Scatter3d(
                x=pheno_df["PHATE1"], y=pheno_df["PHATE2"], z=pheno_df["VR_score"] + 0.005,
                mode='markers', marker=dict(size=2, opacity=0.7, color=color),
                name=f'{pheno}'
            ))

    fig.update_layout(title="VR Landscape: All Phenotypes Integrated", scene=dict(xaxis_title="PHAT$
    fig.write_html(filename)
    print(f"Saved: {filename}")



# ==========================================
# 3. MAIN EXECUTION SCRIPT
# ==========================================
# Load Data
path_to_file = "scvelo_results/FULLY_INTEGRATED_velocity.h5ad"
print(f"Loading data from {path_to_file}...")
adata = sc.read_h5ad(path_to_file)

# Run Calculation
print("Analyzing cells...")
df_vr = calculate_vr_score(adata)

# Print unique phenotypes to debug any spelling/case issues
unique_phenotypes = df_vr['Phenotype'].dropna().unique()
print(f"Phenotypes found in your data: {unique_phenotypes}")

# Define Colors (Ensure the keys below EXACTLY match the spelling/casing in the print statement abo$
phenotype_colors = {
    "GAS": "#e5989b",   # Muted Red
    "IM": "#b79ced",        # Muted Purple
    "DUO": "#a2d2ff"   # Muted Blue
}

# Generate Individual Plots
print("\n--- Generating Plots ---")
for pheno, color in phenotype_colors.items():
    output_name = f"scvelo_results/landscape_v2_{pheno}.html"
    plot_3d_landscape(df_vr, target_phenotype=pheno, color_hex=color, filename=output_name)

# Generate Master Plot
plot_all_phenotypes_landscape(df_vr, phenotype_colors, "scvelo_results/landscape_v2_ALL_CELLS.html")

# Save Data
df_vr.to_csv("scvelo_results/VR_score_aligned_with_paper.csv")
print("All tasks complete.")




