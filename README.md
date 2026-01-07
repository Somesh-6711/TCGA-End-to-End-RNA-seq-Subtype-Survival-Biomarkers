# TCGA-End-to-End-RNA-seq-Subtype-Survival-Biomarkers
pipeline that pulls TCGA data and predicts cancer subtype and survival risk, then explains which genes/pathways drive risk.


## Results: Key Figures (Paper-ready)

> Tip: Click any figure to open it full-size.

| Figure | Summary |
|---|---|
| **Fig 1 — Embedding (PCA/UMAP-style)**<br>![Fig 1](reports/figures/fig1_embedding_subtypes.png) | - Shows global structure of RNA-seq samples in low-dim space.<br>- Helps verify whether classes cluster/separate. |
| **Fig 2 — Confusion Matrix**<br>![Fig 2](reports/figures/fig2_confusion_matrix.png) | - Subtype classification performance across classes.<br>- Reveals which classes are commonly confused. |
| **Fig 3 — Kaplan–Meier (Risk Groups)**<br>![Fig 3](reports/figures/fig3_km_risk_groups.png) | - High vs Low risk split using Cox risk score.<br>- Demonstrates survival separation (log-rank). |
| **Fig 4 — Forest Plot (Top Cox Biomarkers)**<br>![Fig 4](reports/figures/fig4_forest_cox_biomarkers.png) | - Effect sizes for top survival-associated genes.<br>- Direction indicates risk-increasing vs protective genes. |
| **Fig 5 — Heatmap (Risk-Stratified Biomarkers)**<br>![Fig 5](reports/figures/fig5_biomarker_heatmap_risk.png) | - Expression patterns of top biomarkers across patients.<br>- Visual confirmation of risk-linked expression shifts. |
| **Fig 6 — Volcano Plot (Differential Expression)**<br>![Fig 6](reports/figures/fig6_volcano_oligo_vs_astro.png) | - DE results between classes (e.g., Oligo vs Astro).<br>- Highlights large fold-change + significant genes. |

---


| Figure | Summary |
|---|---|
| **Fig 1 — Embedding**<br>[![Fig 1](reports/figures/fig1_embedding_subtypes.png)](reports/figures/fig1_embedding_subtypes.png) | - Global sample structure in reduced space.<br>- Quick visual check of cluster separation. |
| **Fig 2 — Confusion Matrix**<br>[![Fig 2](reports/figures/fig2_confusion_matrix.png)](reports/figures/fig2_confusion_matrix.png) | - Class-wise performance overview.<br>- Identifies dominant confusion patterns. |
| **Fig 3 — KM Risk Groups**<br>[![Fig 3](reports/figures/fig3_km_risk_groups.png)](reports/figures/fig3_km_risk_groups.png) | - Survival separation across predicted risk groups.<br>- Evidence of prognostic signal. |
| **Fig 4 — Forest Biomarkers**<br>[![Fig 4](reports/figures/fig4_forest_cox_biomarkers.png)](reports/figures/fig4_forest_cox_biomarkers.png) | - Top genes ranked by Cox effect size.<br>- Interpretable survival biomarkers. |
| **Fig 5 — Biomarker Heatmap**<br>[![Fig 5](reports/figures/fig5_biomarker_heatmap_risk.png)](reports/figures/fig5_biomarker_heatmap_risk.png) | - Biomarker expression vs risk stratification.<br>- Pattern-level biological intuition. |
| **Fig 6 — Volcano DE**<br>[![Fig 6](reports/figures/fig6_volcano_oligo_vs_astro.png)](reports/figures/fig6_volcano_oligo_vs_astro.png) | - DE genes: effect size vs significance.<br>- Highlights candidate subtype drivers. |
