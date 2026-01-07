# TCGA-End-to-End-RNA-seq-Subtype-Survival-Biomarkers
pipeline that pulls TCGA data and predicts cancer subtype and survival risk, then explains which genes/pathways drive risk.


## Results: Key Figures (Paper-ready)

> Tip: Click any figure to open it full-size.
## Results Summary + Key Figures
Click any thumbnail to open the full-resolution figure.

| Item | Summary |
|---|---|
| **Overall results (LGG)** | **Subtype classification (LogReg, 5-fold CV):** Acc **0.585**, Bal Acc **0.558**, Macro F1 **0.543**. <br><br> **Survival (CoxPH, patient-aggregated, 500 genes, penalizer=0.5):** C-index **0.902 (train)**, **0.863 (test)**. <br><br> **Risk stratification (KM, median split):** Log-rank p-value **3.52e-4** (significant separation). |
| <a href="reports/figures/fig1_embedding_subtypes.png"><img src="reports/figures/fig1_embedding_subtypes.png" width="240"></a><br><b>Fig 1 — Embedding (PCA/UMAP)</b> | Low-dimensional embedding of RNA-seq profiles using top-variance genes to visually assess class separation and potential outliers/batch effects. |
| <a href="reports/figures/fig2_confusion_matrix.png"><img src="reports/figures/fig2_confusion_matrix.png" width="240"></a><br><b>Fig 2 — Confusion Matrix</b> | Test-set confusion matrix showing class-wise errors and which histological groups are most frequently confused. |
| <a href="reports/figures/fig3_km_risk_groups.png"><img src="reports/figures/fig3_km_risk_groups.png" width="240"></a><br><b>Fig 3 — Kaplan–Meier (Risk Groups)</b> | Kaplan–Meier curves for predicted High vs Low risk patients from Cox risk scores; includes log-rank test for survival separation. |
| <a href="reports/figures/fig4_forest_cox_biomarkers.png"><img src="reports/figures/fig4_forest_cox_biomarkers.png" width="240"></a><br><b>Fig 4 — Cox Biomarkers (Forest)</b> | Forest-style plot of top survival-associated genes ranked by Cox effect size (hazard direction / magnitude). |
| <a href="reports/figures/fig5_biomarker_heatmap_risk.png"><img src="reports/figures/fig5_biomarker_heatmap_risk.png" width="240"></a><br><b>Fig 5 — Biomarker Heatmap</b> | Heatmap of top survival biomarkers across patients ordered/grouped by predicted risk to reveal consistent expression patterns. |
| <a href="reports/figures/fig6_volcano_oligo_vs_astro.png"><img src="reports/figures/fig6_volcano_oligo_vs_astro.png" width="240"></a><br><b>Fig 6 — Volcano (DE)</b> | Differential expression volcano plot (Oligodendroglioma vs Astrocytoma), highlighting significant genes by effect size and FDR. |

