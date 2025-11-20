# 📘 AI Omics Internship 2025 — Practical Sessions  

This repository contains all practical sessions completed during **Module I (R Programming Basics)** and **Module II (Genomics Data Analysis)**. Each task includes scripts, workflows, and detailed summaries of analytical steps.

---

## 🔷 Module I — Basics of R Programming

### **Class Sessions**
1. **Class 1A:** R Program Installation  
2. **Class 1B:** R Basic Operations — [View Script](Codes/Module%20I-Getting%20Started%20WIth%20R-Class%20Ib.R)  
3. **Class 1C:** R Basic Syntax — [View Script](Codes/Module%20I-Basic_Syntax-CLass%20Ic.R)  
4. **Class 2:** Operators, Data Structures, Functions & Loops — [View Script](Codes/Module%20I-Basic_R_Functions-Class_2.R)

---

## 🧩 Task 1 — Data Cleaning & Project Setup  
📄 **Script:** [Yusuf_Munir_Aliyu_Assignment_1b.R](Codes/Yusuf_Munir_Aliyu_Assignment%201b.R)

### **Key Activities**
- Created standard project directories (`raw_data`, `clean_data`, `scripts`, `results`, `plot`).  
- Loaded dataset (`patient_info.csv`) and performed structural inspection.  
- Converted categorical variables to factors (`gender`, `diagnosis`, `smoker`, `patient_id`).  
- Created `smoking_status` binary variable (Yes = 1, No = 0).  
- Exported cleaned dataset to `clean_data/patient_info_clean.csv`.

---

## 🧬 Task 2 — Gene Expression Classification Function  
📄 **Script:** [Yusuf_Munir_Aliyu_class_2_Assignment.R](Codes/Yusuf_Munir_Aliyu_class_2_Assignment.R)

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
📄 **Script:** [Yusuf_Munir_Aliyu_4_Assignment.R](Codes/Yusuf_Munir_Aliyu_4_Assignment.R)

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

## 📊 Task 5 — Differential Expression Analysis  
📄 **Script:** [Yusuf_Munir_Aliyu_5_Assignment.R](Codes/Yusuf_Munir_Aliyu_5_Assignment.R)

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

## 📁 Repository Structure
