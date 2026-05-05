import scanpy as sc
import celltypist
import pandas as pd
import os
import matplotlib.pyplot as plt

# Path setup
model_dir = "Models/"
out_dir = "celltypist_input/lineage_results/CellTypist_Comparisons/"
os.makedirs(out_dir, exist_ok=True)

print("Loading data from seu_int.h5ad...")
# Updated to match your 'ls' output
adata = sc.read_h5ad("celltypist_input/seu_int.h5ad")

# Load and merge the Yue labels from R
print("Merging Yue et al. metadata labels...")
yue_meta = pd.read_csv("celltypist_input/yue_bridge_metadata.csv", index_col='barcode')

# Attach the labels. reindex ensures barcodes match perfectly.
adata.obs['Yue_Lineage'] = yue_meta.reindex(adata.obs_names)['Yue_Lineage']

# Define the models based on your Models/ folder
models = {
    "Stomach": "7_healthy_reference_AP_stomach_finalmodel.pkl",
    "Small_Intestine": "8_healthy_reference_AP_small_intestine_finalmodel.pkl",
    "All_Organ": "2_full_healthy_reference_AP_all_organs_finalmodel.pkl"
}

for label, model_file in models.items():
    model_path = os.path.join(model_dir, model_file)
    print(f"\n>>> Validating against {label} model...")

    # Run prediction
    predictions = celltypist.annotate(adata, model=model_path, majority_voting=True)

    # Generate the Comparison DotPlot
    # This maps your 'Yue_Lineage' (columns) vs CellTypist 'majority_voting' (rows)
    celltypist.dotplot(
        predictions,
        use_as_reference='Yue_Lineage',
        use_as_prediction='majority_voting',
        title=f"Yue et al. vs CellTypist {label} Model",
        show=False
    )

    # Save the plot
    save_path = os.path.join(out_dir, f"validation_dotplot_{label}.pdf")
    plt.savefig(save_path, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")

print(f"\nValidation complete! All plots are in {out_dir}")

