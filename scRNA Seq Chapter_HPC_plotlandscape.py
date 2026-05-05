import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
from scipy.interpolate import Rbf

# --- 1. Load Data ---
df = pd.read_csv("scvelo_results/VR_score_per_cell.csv", index_col=0)

# --- 2. Grid (34x30 mesh matching paper, bounds from your PHATE range) ---
grid_x, grid_y = np.mgrid[-0.031:0.024:34j, -0.028:0.019:30j]

# --- 3. RBF interpolation of VR score onto grid (Memory Fix) ---
# Subsample data to max 5000 points for the surface calculation
interp_sub = df.sample(n=min(5000, len(df)), random_state=42)

# 'thin_plate_spline' was renamed to 'thin_plate' in newer scipy
rbfi = Rbf(interp_sub['PHATE1'],
           interp_sub['PHATE2'],
           interp_sub['VR_score'],
           smooth=1,
           function='thin_plate')

di = rbfi(grid_x, grid_y)

# Clip interpolated values that exceed data max (avoids artefacts at edges)
vr_max = round(df['VR_score'].max(), 3) + 0.002
di[di > vr_max] = None

# --- 4. Cluster colours ---
colours = {
    'GAS': '#D65F5F',  # Muted Red
    'IM':  '#9575CD',  # Muted Purple
    'DUO': '#4E79A7'   # Muted Blue
}

# --- 5. Scatter elevation: 0.9*VR + 0.1*CCAT + 0.012 ---
scatter_z = 0.9 * df['VR_score'] + 0.1 * df['CCAT'] + 0.012
# --- 6. Plotting ---
fig = go.Figure()

# Surface Trace
fig.add_trace(go.Surface(
    x=grid_x, y=grid_y, z=di,
    colorscale=[
        [0,    'blue'],
        [0.25, 'green'],
        [0.5,  'yellow'],
        [0.75, '#8B4513'],
        [1.0,  'white']
    ],
    cmin=df['VR_score'].min(),
    cmax=df['VR_score'].max(),
    opacity=0.9,
    showscale=True
))

# Scatter overlay Trace (downsample to 20k for UI speed)
sample_idx = np.random.choice(len(df), size=min(20000, len(df)), replace=False)
df_s = df.iloc[sample_idx]
sz = scatter_z.iloc[sample_idx]

fig.add_trace(go.Scatter3d(
    x=df_s['PHATE1'],
    y=df_s['PHATE2'],
    z=sz,
    mode='markers',
    marker=dict(
        size=2,
        color=df_s['Cluster'].map(colours),
        opacity=0.4
    ),
    name='Cells'
))


# Layout Updates
fig.update_layout(
    title='Waddington Landscape',
    template='simple_white',
    scene=dict(
        xaxis=dict(title='PHATE1', nticks=5),
        yaxis=dict(title='PHATE2', nticks=5),
        zaxis=dict(title='VR Score', nticks=5,
                   range=[df['VR_score'].min() - 0.005,
                          df['VR_score'].max() + 0.02])
    ),
    width=1000, height=1000,
    margin=dict(l=0, r=0, b=0, t=40),
    scene_aspectmode='cube',
    scene_camera=dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=1.5, y=1, z=2)
    )
)

# --- 7. Save Output ---
fig.write_html("scvelo_results/landscape.html")
print("Saved to scvelo_results/landscape.html")


