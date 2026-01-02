# 📘 AI Omics Internship 2025 — Practical Sessions  

This repository contains all practical sessions completed during **Module I (R Programming Basics)** and **Module II (Genomics Data Analysis)**. Each task includes scripts, workflows, and detailed summaries of analytical steps.

---

## 🔷 Module I — Basics of R Programming

### **Class Sessions**
1. **Class 1A:** R Program Installation  
2. **Class 1B:** R Basic Operations — [View Script](https://github.com/aymunir1/AI_Omics_Internship_2025/blob/main/Codes/Module%20I%20Basic%20of%20R/Yusuf_Munir_Aliyu_class_1b_Assignment.R)  
3. **Class 1C:** R Basic Syntax — [View Script](https://github.com/aymunir1/AI_Omics_Internship_2025/blob/main/Codes/Module%20I%20Basic%20of%20R/Yusuf_Munir_Aliyu_class_2_Assignment.R)  
4. **Class 2:** Operators, Data Structures, Functions & Loops — [View Script](https://github.com/aymunir1/AI_Omics_Internship_2025/blob/main/Codes/Module%20I%20Basic%20of%20R/Module%20I-Basic_R_Functions-Class_2.R)

---

## 🧩 Task 1 — Data Cleaning & Project Setup  
📄 [Script:](https://github.com/aymunir1/AI_Omics_Internship_2025/blob/main/Codes/Module%20II%20Microarray%20Data%20Analysis/Yusuf_Munir_Aliyu_class_1b_Assignment.R)

### **Key Activities**
- Created standard project directories (`raw_data`, `clean_data`, `scripts`, `results`, `plot`).  
- Loaded dataset (`patient_info.csv`) and performed structural inspection.  
- Converted categorical variables to factors (`gender`, `diagnosis`, `smoker`, `patient_id`).  
- Created `smoking_status` binary variable (Yes = 1, No = 0).  
- Exported cleaned dataset to `clean_data/patient_info_clean.csv`.
---

## 🧬 Task 2 — Gene Expression Classification Function  
📄 [Script](https://github.com/aymunir1/AI_Omics_Internship_2025/blob/main/Codes/Module%20II%20Microarray%20Data%20Analysis/Yusuf_Munir_Aliyu_class_2_Assignment.R)

### **Summary**
- Implemented custom function `classify_gene()` based on `logFC` and `padj`.  
- **Classification rules:**  
  - Upregulated: `logFC > 1` & `padj < 0.05`  
  - Downregulated: `logFC < -1` & `padj < 0.05`  
  - Not Significant: otherwise  
- Processed two DEG datasets using a loop.  
- Replaced missing adjusted p-values with `1`.  
- Added `status` column and saved results to the `Result` folder.  
- Generated summary frequency tables.

---

# 🔷 Module II — Introduction to Genomics Data Analysis

## 🌐 Task 3 — Dataset Retrieval
- Retrieved microarray datasets from **ArrayExpress** and **NCBI GEO** for analysis.

---

## 🧹 Task 4 — Microarray Preprocessing Workflow  
📄 [Script](https://github.com/aymunir1/AI_Omics_Internship_2025/blob/main/Codes/Module%20II%20Microarray%20Data%20Analysis/Yusuf_Munir_Aliyu_4_Assignment.R)

### **1. Quality Control**
- Performed pre- and post-normalization QC.  
- Identified outliers using boxplots, PCA, and MA plots.  
- Documented outliers before and after normalization.

### **2. Normalization & Probe Filtering**
- Applied normalization methods (RMA, quantile normalization).  
- Filtered low-intensity and non-informative probes.  
- Recorded number of retained transcripts.

### **3. Phenotype Definition**
- Defined biological groups (Normal vs Cancer).  
- Ensured consistent sample labeling.

---

##  Task 5 — Differential Expression Analysis  
📄 [Script](https://github.com/aymunir1/AI_Omics_Internship_2025/blob/main/Codes/Module%20II%20Microarray%20Data%20Analysis/Yusuf_Munir_Aliyu_5_Assignment.R)

### **1. Probe-to-Gene Mapping**
- Used `AnnotationDbi` to map probe IDs to gene symbols.  
- Managed duplicate probes by averaging expression.

### **2. limma Differential Expression**
- Designed contrast: `cancer_vs_normal`.  
- Applied `lmFit()` and `eBayes()` for model fitting.  
- Extracted DEGs using `topTable()` with FDR & log2FC thresholds.

### **3. Visualization**
- Generated volcano plot showing up- & downregulated genes.  
- Created heatmap of top 25 DEGs using `pheatmap`.

---

##  Task 6 — Functonal Enrichment Analysis  
[R Script](https://github.com/aymunir1/AI_Omics_Internship_2025/blob/main/Codes/Module%20II%20Microarray%20Data%20Analysis/Functional_Enrichment.R)

### Enrichment Analysis Setup
- Differentially expressed genes (DEGs) are mapped from **gene symbols to Entrez IDs**
- Invalid, duplicated, or unmapped genes are removed
- Significant genes are separated into **upregulated** and **downregulated** sets
- A background gene universe is defined for overrepresentation testing  


### Overrepresentation Analysis (ORA)

#### Gene Ontology (GO)
- Performs ORA across all three GO domains:
  - **Biological Process (BP)**
  - **Molecular Function (MF)**
  - **Cellular Component (CC)**
- Identifies GO terms significantly enriched in DEGs compared to background genes
- Results visualized using:
  - Dot plots
  - Bar plots
  - Ontology-stratified enrichment plots  

#### KEGG Pathway Analysis
- Identifies significantly enriched **metabolic and signaling pathways**
- Uses Homo sapiens–specific KEGG annotations
- Visualizations include:
  - Dot plots
  - UpSet plots showing gene overlap among pathways  



### Gene Set Enrichment Analysis (GSEA)

- Uses a **ranked gene list** based on log fold change (logFC)
- No hard significance threshold required
- Employs **MSigDB Hallmark Gene Sets**
- Identifies pathways enriched at the top or bottom of the ranked list
- Outputs include:
  - Enrichment score (NES)
  - Direction of regulation (up/down)
  - Core enrichment genes
- Visualization using:
  - Running enrichment score plots
  - Comparative GSEA plots across hallmark pathways  

### Results Summary & Reporting

- Significant results filtered using adjusted p-value < 0.05
- Outputs generated for:
  - GO (BP, MF, CC)
  - KEGG pathways
  - GSEA (upregulated and downregulated pathways)
- Summary tables report:
  - Total number of significant pathways
  - Directionality of enrichment (GSEA)
- All results exported as CSV files for reproducibility and downstream reporting


----

##  Task 7 — Machine Learning
[R Script](https://github.com/aymunir1/AI_Omics_Internship_2025/tree/main/Codes/Module%20III%20Machine%20Learning)
###  Data Preprocessing
- Log₁₀ transformation to stabilize variance  
- Transposition (samples as rows, genes as columns)  
- Near-zero variance filtering  
- Feature scaling and centering  
- KNN imputation for missing values  

### Feature Selection
- **Boruta** identifies statistically significant genes  
- **RFE (Recursive Feature Elimination)** iteratively removes weak predictors  
- Common genes between Boruta and RFE are extracted  

### Model Training
- 70/30 train–test split  
- 10-fold cross-validation  
- **Random Forest (RF)**, **Support Vector Machine (SVM)**, and **Artificial Neural Network (ANN)** trained on selected gene sets  

### Model Evaluation
- Accuracy assessment using confusion matrices  
- ROC curve visualization  
- AUC-based model discrimination comparison  

-----------


## 📁 Repository Structure
AI_Omics_Internship_2025/
├── Codes/
│   ├── Module I Basic of R/
│   ├── Module II Microarray Data Analysis/
│   └── Module III Machine Learning/
│
├── Datasets/
│   ├── Raw datasets
│   ├── Clean Datasets/
│   └── DEG results
│
├── Results/
│   ├── Feature selection
│   ├── Enrichment analysis
│   └── Model evaluation outputs
│
└── README.md
