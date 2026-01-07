# TCGA-End-to-End-RNA-seq-Subtype-Survival-Biomarkers
pipeline that pulls TCGA data and predicts cancer subtype and survival risk, then explains which genes/pathways drive risk.


## Results: Key Figures (Paper-ready)

> Tip: Click any figure to open it full-size.
---
<table>
  <tr>
    <td width="50%" valign="top">
      <b>Fig 1 — Embedding (PCA/UMAP-style)</b><br>
      <a href="reports/figures/fig1_embedding_subtypes.png">
        <img src="reports/figures/fig1_embedding_subtypes.png" width="360">
      </a>
      <ul>
        <li>Low-dimensional embedding of RNA-seq profiles.</li>
        <li>Checks whether samples form coherent groups by label.</li>
        <li>Useful sanity-check for batch effects / outliers.</li>
        <li>Separation suggests transcriptomic signal is learnable.</li>
      </ul>
    </td>
---
   </tr> 
    <td width="50%" valign="top">
      <b>Fig 2 — Confusion Matrix</b><br>
      <a href="reports/figures/fig2_confusion_matrix.png">
        <img src="reports/figures/fig2_confusion_matrix.png" width="360">
      </a>
      <ul>
        <li>Class-wise performance for subtype prediction.</li>
        <li>Highlights dominant misclassification patterns.</li>
        <li>Reveals whether errors are symmetric or biased.</li>
        <li>Helpful for deciding which classes need more features / tuning.</li>
      </ul>
    </td>
  </tr>

  <tr>
    <td width="50%" valign="top">
      <b>Fig 3 — Kaplan–Meier (Risk Groups)</b><br>
      <a href="reports/figures/fig3_km_risk_groups.png">
        <img src="reports/figures/fig3_km_risk_groups.png" width="360">
      </a>
      <ul>
        <li>Survival curves for predicted <i>High</i> vs <i>Low</i> risk patients.</li>
        <li>Risk score computed from Cox model predictions.</li>
        <li>Log-rank test quantifies separation significance.</li>
        <li>Demonstrates clinically meaningful stratification from expression data.</li>
      </ul>
    </td>
   </tr> 

  <tr>
    <td width="50%" valign="top">
      <b>Fig 4 — Forest Plot (Top Cox Biomarkers)</b><br>
      <a href="reports/figures/fig4_forest_cox_biomarkers.png">
        <img src="reports/figures/fig4_forest_cox_biomarkers.png" width="360">
      </a>
      <ul>
        <li>Top survival-associated genes ranked by Cox effect size.</li>
        <li>Positive effects indicate increased hazard (worse prognosis).</li>
        <li>Negative effects indicate protective association (better prognosis).</li>
        <li>Provides an interpretable biomarker shortlist for discussion.</li>
      </ul>
    </td>
  </tr>

  <tr>
    <td width="50%" valign="top">
      <b>Fig 5 — Heatmap (Risk-Stratified Biomarkers)</b><br>
      <a href="reports/figures/fig5_biomarker_heatmap_risk.png">
        <img src="reports/figures/fig5_biomarker_heatmap_risk.png" width="360">
      </a>
      <ul>
        <li>Expression patterns of top biomarkers across patients.</li>
        <li>Columns ordered/grouped by predicted risk group.</li>
        <li>Visualizes concordant shifts between high vs low risk cohorts.</li>
        <li>Helps connect model outputs to biological patterns.</li>
      </ul>
    </td>
   </tr>

  <tr>  
    <td width="50%" valign="top">
      <b>Fig 6 — Volcano Plot (Differential Expression)</b><br>
      <a href="reports/figures/fig6_volcano_oligo_vs_astro.png">
        <img src="reports/figures/fig6_volcano_oligo_vs_astro.png" width="360">
      </a>
      <ul>
        <li>Differential expression between two major classes (e.g., Oligo vs Astro).</li>
        <li>Effect size (log2 fold-change) vs significance (-log10 p/FDR).</li>
        <li>Highlights candidate subtype drivers and marker genes.</li>
        <li>Supports biological interpretation beyond pure prediction metrics.</li>
      </ul>
    </td>
  </tr>
</table>
