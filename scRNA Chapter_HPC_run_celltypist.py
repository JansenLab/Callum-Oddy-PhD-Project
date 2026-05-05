import celltypist
from celltypist import models
import anndata as ad
import pandas as pd
import os

# 1. Define Paths
# Path to your custom models
model_path = "Models/2_full_healthy_reference_AP_all_organs_finalmodel.pkl"
input_h5ad = "celltypist_input/seu_int.h5ad"

# 2. Load the data
print(f"Loading AnnData from {input_h5ad}...")
adata = ad.read_h5ad(input_h5ad)

# 3. Load the local model
print(f"Loading custom model from {model_path}...")
model = models.Model.load(model_path)

# 4. Run Prediction
# majority_voting=True uses your Seurat clusters to provide cleaner results
print("Running CellTypist annotation...")
predictions = celltypist.annotate(adata, model = model, majority_voting = True)

# 5. Save results
print("Saving results...")

# Export an updated h5ad that contains the new labels in adata.obs
adata_labeled = predictions.to_adata()
adata_labeled.write_h5ad("celltypist_input/seu_labeled.h5ad")

# Also save a CSV of just the barcodes and labels for easy R import later
# We grab 'predicted_labels' and 'majority_voting'
annotation_df = adata_labeled.obs[['predicted_labels', 'majority_voting']]
annotation_df.to_csv("celltypist_input/cell_labels.csv")

print("--- Annotation Complete ---")
print("Files created:")
print("- celltypist_input/seu_labeled.h5ad")
print("- celltypist_input/cell_labels.csv")
