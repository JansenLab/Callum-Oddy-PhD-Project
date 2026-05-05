import scanpy as sc
import numpy as np
import pandas as pd
import graphtools
import plotly.graph_objects as go
import plotly.io as pio
from sklearn.preprocessing import MinMaxScaler
from scipy.interpolate import Rbf

pio.renderers.default = "plotly_mimetype+notebook"

# ==========================================
# STEP 1 — Load Data
# ==========================================
adata = sc.read_h5ad("scvelo_results/FULLY_INTEGRATED_velocity.h5ad")

# ==========================================
# STEP 2 — Compute dist_scaled (Paper exact: knn=50, manhattan)
# ==========================================
def compute_distdeg(adata, cluster_key='tissue_group',
                    knn=50, distance_metric='manhattan',
                    scale=True, rmv_outliers=True):

    scaler = MinMaxScaler()
    adata.obs['dist_scaled'] = np.nan

    for cluster in adata.obs[cluster_key].unique():
        print(f"Computing centrality for: {cluster}")
        mask = adata.obs[cluster_key] == cluster
        cell_index = adata.obs_names[mask]
        X_phate = adata[cell_index].obsm['X_phate']

        G = graphtools.Graph(X_phate, knn=knn, decay=None,
                             distance=distance_metric, verbose=False)
        dist_matrix = G.shortest_path(distance='data')
        median_dists = np.median(dist_matrix, axis=1)

        if rmv_outliers:
            q99 = np.percentile(median_dists, 99)
            cluster_median = np.median(median_dists)
            median_dists[median_dists > q99] = cluster_median

        if scale:
            dist_scaled = scaler.fit_transform(
                median_dists.reshape(-1, 1)).flatten()
            adata.obs.loc[cell_index, 'dist_scaled'] = dist_scaled

    return adata

adata = compute_distdeg(adata)


# ==========================================
# STEP 3 — Compute inverse velocity (scaled globally)
# ==========================================
scaler = MinMaxScaler()
v_mag = np.linalg.norm(adata.obsm['velocity_phate'], axis=1)
inv_v = 1.0 / (v_mag + 1e-6)
adata.obs['inv_velo_scaled'] = scaler.fit_transform(inv_v.reshape(-1, 1)).flatten()

# ==========================================
# STEP 4 — Build dataframe & Export CSV
# ==========================================
df = pd.DataFrame({
    'PHATE1':            adata.obsm['X_phate'][:, 0],
    'PHATE2':            adata.obsm['X_phate'][:, 1],
    'Sample':            adata.obs['patient'].values,
    'Cluster':           adata.obs['tissue_group'].values,
    'CCAT':              adata.obs['CytoTRACE2_Score'].values,
    'inv_velo_scaled':   adata.obs['inv_velo_scaled'].values,
    'dist_scaled':       adata.obs['dist_scaled'].values,
    'final_annotation': adata.obs['final_annotation'].values
}, index=adata.obs_names)

# Cluster-level medians
df['CCATmed']    = df.groupby('Cluster')['CCAT'].transform('median')
df['invVELOmed'] = df.groupby('Cluster')['inv_velo_scaled'].transform('median')

# VR score — paper formula
df['VR_score'] = 0.9 * df['CCATmed'] + 0.1 * (df['invVELOmed'] * df['dist_scaled'])

print(f"VR_score range: {df['VR_score'].min():.6f} to {df['VR_score'].max():.6f}")

# Format and save the requested CSV for your collaborator
df_export = df[['PHATE1', 'PHATE2', 'VR_score', 'Cluster','final_annotation']].copy()
df_export.rename(columns={
    'PHATE1': 'Phate1',
    'PHATE2': 'Phate2',
    'VR_score': 'VR',
    'Cluster': 'organ of origin',
    'final_annotation': 'annotation of cell types (final_annotation)'
}, inplace=True)

df_export.to_csv("scvelo_results/VR_score_final.csv", index_label="Cell_ID")
print("Saved VR_score_final.csv")

# ==========================================
# STEP 5 — Plot (one landscape per organoid type)
# ==========================================
grid_x, grid_y = np.mgrid[-0.031:0.024:0.0024, -0.028:0.019:0.0024]

colours = {
    'GAS': '#D65F5F',
    'IM':  '#9575CD',
    'DUO': '#4E79A7'
}

df_dict = {
    'GAS': df[df['Cluster'] == 'GAS'].copy(),
    'IM':  df[df['Cluster'] == 'IM'].copy(),
    'DUO': df[df['Cluster'] == 'DUO'].copy()
}

for key, df_sub in df_dict.items():
    print(f"Plotting {key}...")
    df_sub = df_sub.dropna(subset=['PHATE1', 'PHATE2', 'VR_score'])

    rbf_sub = df_sub.sample(n=min(5000, len(df_sub)), random_state=42)

    rbfi = Rbf(rbf_sub['PHATE1'], rbf_sub['PHATE2'], rbf_sub['VR_score'],
               smooth=1, function='thin_plate')
    di = rbfi(grid_x, grid_y)

    di[di > (round(df_sub['VR_score'].max(), 3) + 0.002)] = None

    fig = go.Figure()
   fig.add_trace(go.Surface(
        x=grid_x, y=grid_y, z=di,
        colorscale=[
            [0,    'blue'],
            [0.25, 'green'],
            [0.5,  'yellow'],
            [0.75, '#5f3316'],
            [0.85, '#513b1b'],
            [1.0,  'white']
        ],
        cmin=0.02, cmax=0.11,
        opacity=0.9,
        showscale=True
    ))

    scatter_z = 0.1 * df_sub['CCAT'] + 0.9 * df_sub['VR_score'] + 0.012

    if len(df_sub) > 20000:
        idx = np.random.choice(len(df_sub), size=20000, replace=False)
        df_plot = df_sub.iloc[idx]
        sz = scatter_z.iloc[idx]
    else:
        df_plot = df_sub
        sz = scatter_z

    fig.add_trace(go.Scatter3d(
        x=df_plot['PHATE1'],
        y=df_plot['PHATE2'],
        z=sz,
        mode='markers',
        marker=dict(
            size=4,
            color=colours[key],
            opacity=0.42
        ),
 name=f'{key} cells'
    ))

    fig.update_layout(
        title=None,
        template='simple_white',
        scene=dict(
            xaxis=dict(nticks=5, title=''),
            yaxis=dict(nticks=5, title=''),
            zaxis=dict(nticks=5, range=[0.02, 0.13], title='')
        ),
        autosize=True,
        width=1000, height=1000,
        margin=dict(l=0, r=0, b=0, t=0),
        scene_aspectmode='cube',
        scene_camera=dict(
            up=dict(x=0, y=0, z=1),
            center=dict(x=0, y=0, z=0),
            eye=dict(x=1.5, y=1, z=2)
        )
    )

    fig.write_html(f"scvelo_results/landscape_{key}.html")
    print(f"Saved landscape_{key}.html")

# ==========================================
# MASTER PLOT — All three organoid types together
# ==========================================
print("Plotting master landscape...")

df_all = df.dropna(subset=['PHATE1', 'PHATE2', 'VR_score'])

rbf_sub = df_all.sample(n=5000, random_state=42)
rbfi = Rbf(rbf_sub['PHATE1'], rbf_sub['PHATE2'], rbf_sub['VR_score'],
           smooth=1, function='thin_plate')
di = rbfi(grid_x, grid_y)
di[di > (round(df_all['VR_score'].max(), 3) + 0.002)] = None

fig = go.Figure()

fig.add_trace(go.Surface(
    x=grid_x, y=grid_y, z=di,
    colorscale=[
        [0,    'blue'],
        [0.25, 'green'],
        [0.5,  'yellow'],
        [0.75, '#5f3316'],
        [0.85, '#513b1b'],
        [1.0,  'white']
    ],
    cmin=0.02, cmax=0.11,
    opacity=0.9,
    showscale=True
))
for key, colour in colours.items():
    df_sub = df_all[df_all['Cluster'] == key]
    scatter_z = 0.1 * df_sub['CCAT'] + 0.9 * df_sub['VR_score'] + 0.012

    if len(df_sub) > 20000:
        idx = np.random.choice(len(df_sub), size=20000, replace=False)
        df_plot = df_sub.iloc[idx]
        sz = scatter_z.iloc[idx]
    else:
        df_plot = df_sub
        sz = scatter_z

    fig.add_trace(go.Scatter3d(
        x=df_plot['PHATE1'],
        y=df_plot['PHATE2'],
        z=sz,
        mode='markers',
        marker=dict(size=4, color=colour, opacity=0.42),
        name=key
    ))
  
fig.update_layout(
    title=None,
    template='simple_white',
    scene=dict(
        xaxis=dict(nticks=5, title=''),
        yaxis=dict(nticks=5, title=''),
        zaxis=dict(nticks=5, range=[0.02, 0.13], title='')
    ),
    autosize=True,
    width=1000, height=1000,
    margin=dict(l=0, r=0, b=0, t=0),
    scene_aspectmode='cube',
    scene_camera=dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=1.5, y=1, z=2)
    )
)

fig.write_html("scvelo_results/landscape_ALL.html")
print("Saved landscape_ALL.html")









