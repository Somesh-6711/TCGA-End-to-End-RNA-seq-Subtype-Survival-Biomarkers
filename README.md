# TCGA LGG: RNA-seq Subtype Classification + Survival Risk + Biomarkers (End-to-End)
pipeline that pulls TCGA data and predicts cancer subtype and survival risk, then explains which genes/pathways drive risk.

End-to-end research pipeline using **TCGA Lower Grade Glioma (LGG)** RNA-seq expression data to:

- Build a **multi-class subtype classifier**
- Train a **Cox proportional hazards survival model**
- Produce **risk stratification (Kaplan–Meier)**
- Generate a **biomarkers table** + **publication-style plots**

Outputs include:
- Clean processed datasets
- Metrics + tables
- 6 research figures (embedding, confusion matrix, KM, forest plot, heatmap, volcano)

---

## Dataset

### Download TCGA LGG from UCSC Xena

This downloads expression + clinical + survival (and PanCan subtype table if enabled).
```
python -m src.data.download_xena --cohort LGG --force
```

### Expected outputs:

data/interim/expression.parquet

data/interim/clinical.parquet

data/interim/labels.parquet

reports/tables/label_summary.csv

---
### Build train/test datasets

Select top-variance genes and split into train/test.
```
python -m src.features.build_dataset --top-genes 5000 --test-size 0.2 --seed 42
```
### Outputs:

data/processed/X_train.parquet

data/processed/X_test.parquet

data/processed/y_class_train.csv

data/processed/y_class_test.csv

data/processed/survival_train.csv

data/processed/survival_test.csv

reports/tables/split_summary.csv

---
## Train subtype classification model

Multi-class logistic regression with CV.
```
python -m src.models.train_classification --cv 5 --seed 42 --top-genes-out 25
```
### Outputs:

reports/metrics/classification_metrics.json

reports/tables/classification_test_report.csv

reports/tables/classification_confusion_matrix.csv

reports/tables/top_genes_coefficients.csv

---
## Results: Key Figures (Paper-ready)

Classification (test): accuracy ~0.58, balanced accuracy ~0.56

Cox survival (patient-agg): test C-index ~0.86

KM risk stratification: significant separation (log-rank p-value reported in km_metrics.json)


---
## Results Summary + Key Figures
> Tip: Click any figure to open it full-size.

| Item | Summary |
|---|---|
| **Overall results (LGG)** | **Subtype classification (LogReg, 5-fold CV):** Acc **0.585**, Bal Acc **0.558**, Macro F1 **0.543**. <br><br> **Survival (CoxPH, patient-aggregated, 500 genes, penalizer=0.5):** C-index **0.902 (train)**, **0.863 (test)**. <br><br> **Risk stratification (KM, median split):** Log-rank p-value **3.52e-4** (significant separation). |
| <a href="reports/figures/fig1_embedding_subtypes.png"><img src="reports/figures/fig1_embedding_subtypes.png" width="240"></a><br><b>Fig 1 — Embedding (PCA/UMAP)</b> | Low-dimensional embedding of RNA-seq profiles using top-variance genes to visually assess class separation and potential outliers/batch effects. |
| <a href="reports/figures/fig2_confusion_matrix.png"><img src="reports/figures/fig2_confusion_matrix.png" width="240"></a><br><b>Fig 2 — Confusion Matrix</b> | Test-set confusion matrix showing class-wise errors and which histological groups are most frequently confused. |
| <a href="reports/figures/fig3_km_risk_groups.png"><img src="reports/figures/fig3_km_risk_groups.png" width="240"></a><br><b>Fig 3 — Kaplan–Meier (Risk Groups)</b> | Kaplan–Meier curves for predicted High vs Low risk patients from Cox risk scores; includes log-rank test for survival separation. |
| <a href="reports/figures/fig4_forest_cox_biomarkers.png"><img src="reports/figures/fig4_forest_cox_biomarkers.png" width="240"></a><br><b>Fig 4 — Cox Biomarkers (Forest)</b> | Forest-style plot of top survival-associated genes ranked by Cox effect size (hazard direction / magnitude). |
| <a href="reports/figures/fig5_biomarker_heatmap_risk.png"><img src="reports/figures/fig5_biomarker_heatmap_risk.png" width="240"></a><br><b>Fig 5 — Biomarker Heatmap</b> | Heatmap of top survival biomarkers across patients ordered/grouped by predicted risk to reveal consistent expression patterns. |
| <a href="reports/figures/fig6_volcano_oligo_vs_astro.png"><img src="reports/figures/fig6_volcano_oligo_vs_astro.png" width="240"></a><br><b>Fig 6 — Volcano (DE)</b> | Differential expression volcano plot (Oligodendroglioma vs Astrocytoma), highlighting significant genes by effect size and FDR. |

---
### Notes / limitations

TCGA clinical labels vary by cohort and hub; some subtype fields may be missing.

Survival models with many features require penalization.

For publication-grade biomarker validation, add:

     external validation cohort

     pathway enrichment (GSEA)

     multivariate Cox including clinical covariates

---
### Reproducibility

#### Run the full pipeline in order:

1. download_xena

2. build_dataset

3. train_classification

4. train_survival_cox

5. km_risk_stratification

6. make_biomarkers_table

7. analysis plots

---
## Future work  

### Data & labels

Use WHO-style molecular labels for LGG (IDH mutation + 1p/19q codeletion; often more meaningful than histology). If Xena doesn’t expose them cleanly, pull them from GDC/clinical biospec or curated glioma datasets and rebuild the labels. (Your pipeline already supports label swapping.)

External validation on independent cohorts (e.g., CGGA / REMBRANDT / GEO). Show that your risk score + biomarkers generalize beyond TCGA.

Batch-effect/confound control: adjust for age, sex, tumor purity, platform/batch; consider ComBat or including covariates in the Cox model.

### Modeling upgrades

Nested CV + patient-level splitting everywhere (you already did patient aggregation for survival—extend that discipline to all analyses).

### Survival modeling beyond CoxPH:

Elastic-net Cox (strong baseline for high-dimensional genomics),

Random Survival Forests,

Deep survival (DeepSurv/DeepHit) if you want a DL angle.

Time-dependent metrics: add time-dependent AUC, calibration curves, and integrated Brier score (C-index alone can look great even when calibration is off).

### Biology/interpretability

Pathway analysis (GSEA/GO/KEGG): instead of “top genes,” report top pathways that explain subtype separation and risk groups.

Stability of biomarkers: bootstrap the training set and report how often each gene appears among top features (biomarkers that “survive” resampling are more credible).

Mechanistic interpretation: link your top risk genes to known glioma biology (IDH pathway, cell cycle, immune microenvironment) and cite supporting literature.

### Clinical utility & communication

Nomogram / risk score report: convert Cox risk score into a clinically interpretable tool (risk groups, predicted survival probabilities at 1/3/5 years).

Ablation studies: compare “expression only” vs “expression + clinical covariates” vs “multi-omics,” and show incremental gains.

### Reproducibility/product polish

Full reproducibility: pin versions (conda env + requirements.txt), add a single “run all” pipeline command, and optionally Docker/Snakemake.

Model cards + data provenance: document cohort, exclusions (missing OS), and exact label definitions.

Data in project is coming from UCSC Xena TCGA hubs (e.g., tcga.xenahubs.net) and the Pan-Cancer Atlas subtype table file you downloaded.

---
## References

1. Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*, 57(1), 289–300.

2. Cox, D. R. (1972). Regression models and life-tables. *Journal of the Royal Statistical Society: Series B*, 34(2), 187–220.

3. Kaplan, E. L., & Meier, P. (1958). Nonparametric estimation from incomplete observations. *Journal of the American Statistical Association*, 53(282), 457–481.

4. The Cancer Genome Atlas (TCGA) Pan-Cancer analysis project. (2013). *Nature Genetics*, 45(10), 1113–1120.

5. Hoadley, K. A., et al. (2018). Cell-of-origin patterns dominate the molecular classification of 10,000 tumors from 33 types of cancer (TCGA Pan-Cancer Atlas). *Cell*.

6. Goldman, M. J., et al. (2020). Visualizing and interpreting cancer genomics data via the Xena platform. *Nature Biotechnology*.

7. UCSC Xena Browser / Xena Hubs (TCGA data distribution). Retrieved from UCSC Xena.

8. scikit-learn developers. (2011). Scikit-learn: Machine learning in Python.

9. Davidson-Pilon, C. (lifelines). Lifelines: survival analysis in Python. Documentation.
---