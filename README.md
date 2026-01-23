# 🧬 AI_Omics_Internship (2025)
This repository documents a **comprehensive, hands-on bioinformatics and machine learning training program** focused on **omics data analysis**, with practical implementations in **R**, **microarray transcriptomics**, **machine learning**, and **functional enrichment analysis**.
The project demonstrates an end-to-end analytical pipeline from raw gene expression data to biological interpretation and predictive modeling aligned with modern **AI-driven biomedical research** standards.

---
## Project Objectives
- Develop strong foundations in **R for bioinformatics**
- Process and normalize **Affymetrix microarray data**
- Perform **differential gene expression analysis**
- Apply **machine learning models** to high-dimensional omics data
- Identify biologically meaningful patterns via **GO, KEGG, and GSEA**
- Ensure **reproducibility, interpretability, and rigorous evaluation**
---
## 📁 Repository Structure

```
AI_Omics_Internship_2025/
├── Codes/
│   ├── Module I Basic of R/
│   │   └── R fundamentals, data structures, functions
│   │
│   ├── Module II Microarray Data Analysis/
│   │   ├── Differential expression analysis
│   │   └── Functional enrichment (GO, KEGG, GSEA)
│   │
│   └── Module III Machine Learning/
│       ├── Data preprocessing
│       ├── Feature selection (Boruta, RFE)
│       └── Model training & evaluation
│
├── Datasets/
│   ├── Raw Datasets/
│   ├── Clean Datasets/
│   └── DEG Results/
│
├── Results/
│   ├── Feature selection outputs
│   ├── Enrichment analysis tables
│   └── Model performance metrics
│
├── Workspace/
│   └── Saved RData workspaces
│
└── README.md
```


---
## Module Overview
### Module I: Basics of R
1. **Class 1A:** R Program Installation  
2. **Class 1B:** R Basic Operations — [View Script](https://github.com/aymunir1/AI_Omics_Internship_2025/blob/main/Codes/Module%20I%20Basic%20of%20R/Yusuf_Munir_Aliyu_class_1b_Assignment.R)  
3. **Class 1C:** R Basic Syntax — [View Script](https://github.com/aymunir1/AI_Omics_Internship_2025/blob/main/Codes/Module%20I%20Basic%20of%20R/Yusuf_Munir_Aliyu_class_2_Assignment.R)  
4. **Class 2:** Operators, Data Structures, Functions & Loops — [View Script](https://github.com/aymunir1/AI_Omics_Internship_2025/blob/main/Codes/Module%20I%20Basic%20of%20R/Module%20I-Basic_R_Functions-Class_2.R)
  
### Module II: Microarray Data Analysis
- Quality-controlled gene expression data handling [Script](https://github.com/aymunir1/AI_Omics_Internship_2025/blob/main/Codes/Module%20II%20Microarray%20Data%20Analysis/Yusuf_Munir_Aliyu_4_Assignment.R)
- Differential gene expression analysis [Script](https://github.com/aymunir1/AI_Omics_Internship_2025/blob/main/Codes/Module%20II%20Microarray%20Data%20Analysis/Yusuf_Munir_Aliyu_5_Assignment.R)
- Functional interpretation using: [R Script](https://github.com/aymunir1/AI_Omics_Internship_2025/blob/main/Codes/Module%20II%20Microarray%20Data%20Analysis/Functional_Enrichment.R)
  - **Gene Ontology (GO):** BP, MF, CC
  - **KEGG pathway analysis**
  - **Gene Set Enrichment Analysis (GSEA)** using MSigDB Hallmark sets
- Biological insight generation from DEG profiles
### Module III: Machine Learning for Omics Data [R Script](https://github.com/aymunir1/AI_Omics_Internship_2025/tree/main/Codes/Module%20III%20Machine%20Learning)
- High-dimensional data preprocessing
- Feature reduction and selection
- Supervised classification models
- Model performance comparison and validation
---
## Methodological Workflow
##  Probe Annotation & Differential Gene Expression
### Probe-to-Gene Mapping
- Mapped Affymetrix probe IDs to **gene symbols**
- Identified multiple probes per gene
- Collapsed duplicates using **average expression (avereps)**
- Ensured one expression value per gene
### Differential Expression Analysis
- Linear modeling with `limma`
- Design matrix without intercept
- Contrast: **cancer_vs_normal**
- Multiple testing correction (Benjamini–Hochberg)
- Genes classified as:
  - Upregulated (logFC > 1)
  - Downregulated (logFC < −1)
### Visualizations
- Volcano plot (logFC vs −log₁₀ p-value)
- Heatmap of top 25 DEGs
- All plots exported as high-resolution PNGs
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
### Machine learning
- **Data Preprocessing**
- Log₁₀ transformation to stabilize variance
- Data transposition (samples as rows, genes as columns)
- Near-zero variance filtering
- Feature centering and scaling
- KNN-based missing value imputation
- **Feature Selection**
- **Boruta:** Random Forest–based feature importance testing
- **Recursive Feature Elimination (RFE):** Iterative feature pruning
- Extraction of **common informative genes** across methods
- **Model Training**
- 70/30 training–testing split
- 10-fold cross-validation
- Algorithms implemented:
  - Random Forest (RF)
  - Support Vector Machine (SVM – radial kernel)
  - Artificial Neural Network (ANN)
- **Model Evaluation**
- Confusion matrices and accuracy comparison
- ROC curve visualization
- Area Under the Curve (AUC) for model discrimination
- ---  
## Outputs & Results
- Selected gene lists (Boruta, RFE, intersected genes)
- GO, KEGG, and GSEA enrichment tables
- Machine learning performance metrics
- ROC and AUC comparisons
- Reproducible `.RData` workspaces
- Exported CSV files for downstream analysis
---
## 🛠 Software, Packages, and Libraries Used
This project was implemented entirely in **R**, leveraging **CRAN** and **Bioconductor** packages widely used in transcriptomics, machine learning, and functional enrichment analysis.
## 🔹 Core Programming Environment
- **R** – statistical computing and data analysis language
- **RStudio** – integrated development environment (IDE)
## 🔹 Data Acquisition & Microarray Processing
- **GEOquery** – retrieval of GEO series matrices and raw CEL files
- **affy** – preprocessing and RMA normalization of Affymetrix microarray data
- **arrayQualityMetrics** – quality control diagnostics for microarray experiments
## 🔹 Annotation & Probe Mapping
- **AnnotationDbi** – annotation framework for biological databases
- **hgu133plus2.db** – Affymetrix HG-U133 Plus 2.0 probe annotation database
- **org.Hs.eg.db** – human gene annotation (Entrez IDs, gene symbols)
## 🔹 Differential Expression Analysis
- **limma** – linear modeling and empirical Bayes statistics for DEG analysis
- **Biobase** – handling of ExpressionSet objects
- **genefilter** – filtering of low-variance and low-expression features
## 🔹 Data Preprocessing & Feature Engineering
- **caret** – unified framework for data preprocessing and model training
- **DMwR** – KNN-based missing value imputation
- **dplyr** – data manipulation and transformation
- **tidyr** – tidy data reshaping
- **tibble** – modern data frame handling
## 🔹 Feature Selection
- **Boruta** – Random Forest–based all-relevant feature selection
- **randomForest** – ensemble learning and variable importance
- **caret** – Recursive Feature Elimination (RFE) implementation
## 🔹 Machine Learning Models
- **randomForest** – Random Forest classifier
- **kernlab** – Support Vector Machine (SVM) models
- **nnet** – Artificial Neural Network (ANN) implementation
- **e1071** – auxiliary utilities for SVM and classification
## 🔹 Model Evaluation & Performance Metrics
- **pROC** – ROC curve generation and AUC calculation
- **caret** – confusion matrices and cross-validation
- **ROCR** – classifier performance visualization
## Functional Enrichment Analysis
- **clusterProfiler** – GO, KEGG, and GSEA enrichment analysis
- **enrichplot** – visualization of enrichment results
- **msigdbr** – MSigDB Hallmark gene set retrieval
- **DOSE** – enrichment result handling and visualization support
## Data Visualization
- **ggplot2** – statistical graphics and plots
- **pheatmap** – heatmap visualization of gene expression
- **RColorBrewer** – color palettes for plots
- **cowplot** – multi-panel figure assembly
---
### Databases
- NCBI GEO
- Gene Ontology (GO)
- KEGG
- MSigDB (Hallmark gene sets)
## Reproducibility & Best Practices
- Fixed random seeds
- Modular, well-documented scripts
- Structured directories
- Saved analytical workspaces
- Exported results for transparency
## Learning Outcomes
- Practical R programming for omics research
- End-to-end transcriptomics analysis
- Machine learning application in biology
- Feature selection in high-dimensional data
- Biological interpretation of computational results
- Research-ready and reproducible workflows
## Intended Use
- Bioinformatics training and education
- Research skill demonstration
- Internship and graduate program portfolio
- Foundation for RNA-seq, proteomics, or multi-omics extensions
