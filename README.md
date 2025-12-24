# Meniere-Blood-Profiles
This repository contains the analysis pipelines, custom scripts, and processed data associated with the study of **Ménière’s Disease (MD)** signatures in peripheral blood. [cite_start]By integrating **serum metabolomics**, **lipidomics**, and **PBMC transcriptomics**, we characterize the systemic immune and metabolic landscape of MD patients.

---

## Research Overview

[cite_start]Ménière’s disease (MD) is a chronic inner ear syndrome characterized by vertigo, fluctuating sensorineural hearing loss, tinnitus, and aural fullness. While its etiology remains unelucidated, our study supports an autoimmune or autoinflammatory origin through systemic molecular profiling.

### **Key Findings**
* **Systemic Metabolic Dysregulation**: We identified significant upregulation of key immune-signaling molecules, specifically **Sphingosine-1-phosphate (S1P)** and **Glutamate**.
    
* **Causal Protective Role**: Two-sample Mendelian Randomization (MR) analysis suggests a protective effect of **total free cholesterol** against MD.
    
* **Th1 Cell Skewing**: Computational deconvolution and independent **flow cytometry validation** in a larger cohort confirmed a significant increase in peripheral **CD4+ Th1 cells** in MD patients.
* **Signaling Mechanism**: Integrative analysis suggests that the **Psychosine-GPR65 pathway** may drive the observed immune features.
    

---

## Repository Structure

```text
.
├── data/
│   ├── metabolomics/          # Annotated metabolite features (n=690)
│   ├── lipidomics/            # Processed lipid species (n=596)
│   ├── transcriptomics/       # PBMC RNA-seq count matrices and metadata
│   └── gwas_summary/          # Summary statistics from FinnGen R11 and BHS
├── scripts/
│   ├── 01_metabolomics/       # PCA, OPLS-DA, and KEGG enrichment
│   ├── 02_lipidomics/         # Lipid class grouping and Mann-Whitney U tests
│   ├── 03_mendelian_rand/     # Bi-directional MR and sensitivity analyses
│   ├── 04_transcriptomics/    # DESeq2 DEGs, GSEA, and GO enrichment
│   ├── 05_deconvolution/      # CIBERSORTx and ABIS cell type estimation
│   └── 06_integration/        # Correlation analysis (Psychosine vs. Th1 genes)
├── figures/                   # Scripts to reproduce manuscript figures 1-7
├── environment.yml            # Conda environment configuration (R v4.2.0)
└── README.md
