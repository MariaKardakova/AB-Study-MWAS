Gut Microbiome and Depression: Analysis Pipelines

This repository provides code and documentation for a multi-layered investigation of the relationship between the human gut microbiome and major depressive disorder (MDD). The study was conducted using 16S rRNA data from the Atlas Biomed cohort and includes discovery-replication design, meta-analysis, risk scoring, and causal inference using Mendelian Randomization (MR).

📁 Repository Contents
1. Primary Microbiome–Depression Analysis
Script: GitHub_depression_V1.R Purpose:
* Identifies bacterial taxa associated with MDD.
* Implements diversity analyses, differential abundance, and meta-analysis.
Main Steps:
* Data cleaning and case/control definition.
* Discovery vs. replication cohort split.
* Alpha and beta diversity (Shannon, Simpson, Bray-Curtis).
* Differential abundance using:
    * MaAsLin2 (relative abundances).
    * ANCOM-BC2 (absolute counts).
* Meta-analysis of differential abundance.
* Benjamini–Hochberg FDR correction.
Outputs:
* maaslin_results_discovery.csv, replication.csv
* ancombc2_results_discovery.csv, replication.csv
* meta_analysis_results.csv, significant_taxa.csv
Visualisations:
* Violin and PCoA plots for diversity.
* Volcano and bar plots for taxa significance.
* Heatmaps of log-fold changes.

2. Microbiome Risk Score (MRS), Mediation, and Causal Inference
Script: GitHub_depression_MRS_MR_V1.R Purpose:
* Explores microbiome-depression links via advanced statistical models.
* Constructs and validates a microbiome risk score (MRS) for depression.
* Performs MR and mediation analysis.
Main Steps:
* Two-sample MR using MiBioGen and GWAS summary statistics:
    * IVW, MR-Egger, Weighted Median, MR-PRESSO.
* Sensitivity tests for pleiotropy and heterogeneity.
* Mediation models (ACME, ADE) using lifestyle → taxa → MDD framework.
* MRS construction and testing.
* Multivariable regression for taxa–trait associations.
Outputs:
* MR_results.csv
* mediation_results.csv
* MRS_associations.csv
* Taxa_trait_associations.csv
Visualisations:
* MR scatter and forest plots.
* Bar plots of mediation effects.
* Boxplots and ROC curves of MRS.
* Dot plots of taxa–trait regressions.

📊 Data Sources
* Microbiome data: 16S rRNA, processed via QIIME2.
* Phenotype data: Self-reported questionnaire (Atlas Biomed).
* GWAS summary statistics: MDD (PGC), microbiome (MiBioGen).

🔧 Software and Packages
* R version ≥ 4.2
* Key packages:
    * phyloseq, vegan, MaAsLin2, ANCOMBC, meta, mediation, TwoSampleMR, MRPRESSO
* For MRS and mediation: glmnet, caret, pROC

📎 Citation
If you use this repository, please cite:
* Kurilshikov et al. Nature Genetics, 2021. DOI: 10.1038/s41588-020-00763-1
* Mallick H, Franzosa EA, Mclver LJ, Banerjee S, Sirota-Madi A, Kos<c AD, et al.   Predictive metabolomic profiling of microbial communites using amplicon or metagenomic sequences. Nat Commun. 2019;10(1):3136.
* Huang Lin, ANCOM-BC2 Tutorial (May, 2025) https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html
* Associated publication (TBD)
