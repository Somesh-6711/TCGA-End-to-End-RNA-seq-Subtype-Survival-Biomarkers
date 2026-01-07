# TCGA-End-to-End-RNA-seq-Subtype-Survival-Biomarkers
pipeline that pulls TCGA data and predicts cancer subtype and survival risk, then explains which genes/pathways drive risk.


## Results: Key Figures (Paper-ready)

> Tip: Click any figure to open it full-size.
---
## Results: Key Figures (Paper-ready)

## Results: Key Figures (Paper-ready)
Click any thumbnail to open the full-resolution figure.

| Figure | Summary |
|---|---|
| <a href="reports/figures/fig1_embedding_subtypes.png"><img src="reports/figures/fig1_embedding_subtypes.png" width="220"></a><br><b>Fig 1 — Embedding</b> | - Low-dimensional embedding of RNA-seq profiles.<br>- Checks clustering/separation by label.<br>- Useful sanity-check for outliers/batch effects. |
| <a href="reports/figures/fig2_confusion_matrix.png"><img src="reports/figures/fig2_confusion_matrix.png" width="220"></a><br><b>Fig 2 — Confusion Matrix</b> | - Class-wise performance for subtype prediction.<br>- Highlights dominant misclassification patterns. |
| <a href="reports/figures/fig3_km_risk_groups.png"><img src="reports/figures/fig3_km_risk_groups.png" width="220"></a><br><b>Fig 3 — Kaplan–Meier (Risk)</b> | - High vs Low risk survival curves from Cox risk score.<br>- Log-rank test quantifies separation significance. |
| <a href="reports/figures/fig4_forest_cox_biomarkers.png"><img src="reports/figures/fig4_forest_cox_biomarkers.png" width="220"></a><br><b>Fig 4 — Cox Biomarkers (Forest)</b> | - Top survival-associated genes ranked by effect size.<br>- Helps interpret model-driving biomarkers. |
| <a href="reports/figures/fig5_biomarker_heatmap_risk.png"><img src="reports/figures/fig5_biomarker_heatmap_risk.png" width="220"></a><br><b>Fig 5 — Biomarker Heatmap</b> | - Expression patterns across patients ordered by risk group.<br>- Links model output to biological signal. |
| <a href="reports/figures/fig6_volcano_oligo_vs_astro.png"><img src="reports/figures/fig6_volcano_oligo_vs_astro.png" width="220"></a><br><b>Fig 6 — Volcano (DE)</b> | - Differential expression (oligo vs astro).<br>- Highlights significant genes by effect size and FDR. |


### Results Summary + Key Figures

| Figure | Summary |
|---|---|
| **Overall Results** | Subtype classification achieved **Acc=0.585**, **Bal Acc=0.558**, **Macro F1=0.543**. Survival modeling achieved **C-index=0.902 (train)** and **0.863 (test)**. KM risk groups show significant separation (**log-rank p=3.52e-4**). |
| **Fig 1 — Embedding (PCA/UMAP)** | Visual separation of subtype structure in reduced dimensions using top-variance genes. |
| **Fig 2 — Confusion Matrix** | Shows subtype-specific errors and which classes are most frequently confused. |
| **Fig 3 — Kaplan–Meier Risk Groups** | High-risk vs low-risk survival curves with significant separation (log-rank p-value). |
| **Fig 4 — Forest Plot (Cox Biomarkers)** | Top survival-associated genes ranked by Cox effect size (hazard ratio / coefficient). |
| **Fig 5 — Biomarker Heatmap (Risk)** | Heatmap of top biomarkers across samples grouped by predicted risk for pattern discovery. |
| **Fig 6 — Volcano (DE Analysis)** | Differential expression (e.g., oligo vs astro): highlights significant genes by effect size and FDR. |
