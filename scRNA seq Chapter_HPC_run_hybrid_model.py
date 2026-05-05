import celltypist
from celltypist import models
import anndata as ad
import pandas as pd
import numpy as np

# 1. Load Data
print("Loading AnnData...")
adata = ad.read_h5ad("celltypist_input/seu_int.h5ad")

# 2. Run Stomach Model
print("Calculating Stomach probabilities...")
model_st = models.Model.load("Models/7_healthy_reference_AP_stomach_finalmodel.pkl")
pred_st = celltypist.annotate(adata, model=model_st, majority_voting=True)

st_scores = pred_st.probability_matrix.max(axis=1)
st_labels = pred_st.predicted_labels['majority_voting']

# 3. Run Intestine Model
print("Calculating Intestine probabilities...")
model_in = models.Model.load("Models/8_healthy_reference_AP_small_intestine_finalmodel.pkl")
pred_in = celltypist.annotate(adata, model=model_in, majority_voting=True)

in_scores = pred_in.probability_matrix.max(axis=1)
in_labels = pred_in.predicted_labels['majority_voting']

# 4. Create Hybrid Logic
print("Generating Hybrid Labels...")
# Compare scores: if Intestine is more confident, use Intestine label; otherwise use Stomach
hybrid_labels = np.where(in_scores > st_scores, in_labels, st_labels)

# 5. Save everything
results = pd.DataFrame({
    "stomach_majority": st_labels,
    "stomach_conf": st_scores,
    "intestine_majority": in_labels,
    "intestine_conf": in_scores,
    "hybrid_label": hybrid_labels
}, index=adata.obs_names)

results.to_csv("celltypist_input/hybrid_labels.csv")
print("Hybrid labels saved to celltypist_input/hybrid_labels.csv")
