# Meniere-Blood-Profiles
[cite_start]This repository contains the analysis pipelines, custom scripts, and processed data associated with the study of **Ménière’s Disease (MD)** signatures in peripheral blood. [cite_start]By integrating **serum metabolomics**, **lipidomics**, and **PBMC transcriptomics**, we characterize the systemic immune and metabolic landscape of MD patients[cite: 12, 13, 209].

---

## 📖 Research Overview

[cite_start]Ménière’s disease (MD) is a chronic inner ear syndrome characterized by vertigo, fluctuating sensorineural hearing loss, tinnitus, and aural fullness[cite: 10, 22, 206, 216]. [cite_start]While its etiology remains unelucidated, our study supports an autoimmune or autoinflammatory origin through systemic molecular profiling[cite: 11, 26, 207, 220].

### **Key Findings**
* [cite_start]**Systemic Metabolic Dysregulation**: We identified significant upregulation of key immune-signaling molecules, specifically **Sphingosine-1-phosphate (S1P)** and **Glutamate**[cite: 14, 45, 79, 210, 250].
    
* [cite_start]**Causal Protective Role**: Two-sample Mendelian Randomization (MR) analysis suggests a protective effect of **total free cholesterol** against MD[cite: 16, 55, 80, 211, 277].
    
* **Th1 Cell Skewing**: Computational deconvolution and independent **flow cytometry validation** in a larger cohort confirmed a significant increase in peripheral **CD4+ Th1 cells** in MD patients[cite: 17, 70, 81, 212, 307].
* **Signaling Mechanism**: Integrative analysis suggests that the **Psychosine-GPR65 pathway** may drive the observed immune features[cite: 19, 76, 213, 316, 338].
    

---

## 📂 Repository Structure

```text
.
├── data/
│   ├── metabolomics/          # Annotated metabolite features (n=690) [cite: 42, 246]
│   ├── lipidomics/            # Processed lipid species (n=596) [cite: 53, 269]
│   ├── transcriptomics/       # PBMC RNA-seq count matrices and metadata [cite: 40, 237]
│   └── gwas_summary/          # Summary statistics from FinnGen R11 and BHS [cite: 53, 268, 270]
├── scripts/
│   ├── 01_metabolomics/       # PCA, OPLS-DA, and KEGG enrichment [cite: 43, 48, 50, 247, 255, 260]
│   ├── 02_lipidomics/         # Lipid class grouping and Mann-Whitney U tests [cite: 283]
│   ├── 03_mendelian_rand/     # Bi-directional MR and sensitivity analyses [cite: 56, 122, 276, 405]
│   ├── 04_transcriptomics/    # DESeq2 DEGs, GSEA, and GO enrichment [cite: 61, 63, 64, 289, 297, 298]
│   ├── 05_deconvolution/      # CIBERSORTx and ABIS cell type estimation [cite: 66, 67, 303]
│   └── 06_integration/        # Correlation analysis (Psychosine vs. Th1 genes) [cite: 75, 173, 316]
├── figures/                   # Scripts to reproduce manuscript figures 1-7 [cite: 136-175]
├── environment.yml            # Conda environment configuration (R v4.2.0) [cite: 134, 403]
└── README.md
