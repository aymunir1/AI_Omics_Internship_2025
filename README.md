# 🧬 AI-Driven Omics Data Analysis Internship (2025)

This repository documents a **comprehensive, hands-on bioinformatics and machine learning training program** focused on **omics data analysis**, with practical implementations in **R**, **microarray transcriptomics**, **machine learning**, and **functional enrichment analysis**.

The project demonstrates an end-to-end analytical pipeline—from raw gene expression data to biological interpretation and predictive modeling—aligned with modern **AI-driven biomedical research** standards.

---

## Project Objectives

- Build strong foundations in **R programming for bioinformatics**
- Perform **microarray differential gene expression analysis**
- Apply **machine learning models** to high-dimensional omics data
- Identify biologically meaningful patterns using **functional enrichment analysis**
- Ensure **reproducibility, interpretability, and robust evaluation**

---

## Repository Structure

AI_Omics_Internship_2025/
├── Codes/
│ ├── Module I Basic of R/
│ │ └── R fundamentals, data structures, functions
│ │
│ ├── Module II Microarray Data Analysis/
│ │ ├── Differential expression analysis
│ │ └── Functional enrichment (GO, KEGG, GSEA)
│ │
│ └── Module III Machine Learning/
│ ├── Data preprocessing
│ ├── Feature selection (Boruta, RFE)
│ └── Model training & evaluation
│
├── Datasets/
│ ├── Raw Datasets/
│ ├── Clean Datasets/
│ └── DEG Results/
│
├── Results/
│ ├── Feature selection outputs
│ ├── Enrichment analysis tables
│ └── Model performance metrics
│
├── Workspace/
│ └── Saved RData workspaces
│
└── README.md

---

## Module Overview

### Module I: Basics of R
- R syntax and data structures
- Data manipulation and visualization
- Script-based analytical thinking
- Preparation for bioinformatics workflows

---

### Module II: Microarray Data Analysis
- Quality-controlled gene expression data handling
- Differential gene expression analysis
- Functional interpretation using:
  - **Gene Ontology (GO):** BP, MF, CC
  - **KEGG pathway analysis**
  - **Gene Set Enrichment Analysis (GSEA)** using MSigDB Hallmark sets
- Biological insight generation from DEG profiles

---

### Module III: Machine Learning for Omics Data
- High-dimensional data preprocessing
- Feature reduction and selection
- Supervised classification models
- Model performance comparison and validation

---

## Methodological Workflow

### 1 Data Preprocessing
- Log₁₀ transformation to stabilize variance
- Data transposition (samples as rows, genes as columns)
- Near-zero variance filtering
- Feature centering and scaling
- KNN-based missing value imputation

---

### 2 Feature Selection
- **Boruta:** Random Forest–based feature importance testing
- **Recursive Feature Elimination (RFE):** Iterative feature pruning
- Extraction of **common informative genes** across methods

---

### 3 Model Training
- 70/30 training–testing split
- 10-fold cross-validation
- Algorithms implemented:
  - Random Forest (RF)
  - Support Vector Machine (SVM – radial kernel)
  - Artificial Neural Network (ANN)

---

### 4 Model Evaluation
- Confusion matrices and accuracy comparison
- ROC curve visualization
- Area Under the Curve (AUC) for model discrimination

---

## Functional Enrichment Analysis

- **Overrepresentation Analysis (ORA):**
  - Gene Ontology (BP, MF, CC)
  - KEGG pathways
- **Gene Set Enrichment Analysis (GSEA):**
  - MSigDB Hallmark gene sets
  - Identification of upregulated and downregulated pathways
- Publication-ready visualizations:
  - Dot plots
  - Bar plots
  - Running enrichment score plots

---

## Outputs & Results

- Selected gene lists (Boruta, RFE, intersected genes)
- GO, KEGG, and GSEA enrichment tables
- Machine learning performance metrics
- ROC and AUC comparisons
- Reproducible `.RData` workspaces
- Exported CSV files for downstream analysis

---

## 🛠 Tools & Technologies

### Programming & Analysis
- **R**

### Key R Packages
- `caret`
- `randomForest`
- `Boruta`
- `kernlab`
- `clusterProfiler`
- `org.Hs.eg.db`
- `msigdbr`
- `pROC`
- `ggplot2`

### Databases
- Gene Ontology (GO)
- KEGG
- MSigDB (Hallmark gene sets)

---

## Reproducibility & Best Practices

- Fixed random seeds
- Modular, well-documented scripts
- Structured directories
- Saved analytical workspaces
- Exported results for transparency

---

## Learning Outcomes

- Practical R programming for omics research
- End-to-end transcriptomics analysis
- Machine learning application in biology
- Feature selection in high-dimensional data
- Biological interpretation of computational results
- Research-ready and reproducible workflows

---

## Intended Use

- Bioinformatics training and education
- Research skill demonstration
- Internship and graduate program portfolio
- Foundation for RNA-seq, proteomics, or multi-omics extensions

---
