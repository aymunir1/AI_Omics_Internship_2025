## 🧪 Practical Sessions  

**Module I:** Basics of R Programming  
1. **Class 1A:** R Program Installation  
2. **Class 1B:** R Basic Operations — [R Script](https://github.com/aymunir1/AI_Omics_Internship_2025/blob/main/Class_1b.R)  
3. **Class 1C:** R Basic Syntax — [R Script with Practice Exercises](https://github.com/aymunir1/AI_Omics_Internship_2025/blob/main/Class_1c.R)  
4. **Class 2:** Operators in R | Data Structures in R | User-Defined Functions | for-Loop — [R Script](https://github.com/aymunir1/AI_Omics_Internship_2025/edit/main/Class_2.R)
---
**Task 1** — [R Script](https://github.com/aymunir1/AI_Omics_Internship_2025/blob/main/Yusuf_Munir_Aliyu_Assignment%201b.R)
- **Create Subfolders:** Sets up project directories (`raw_data`, `clean_data`, `scripts`, `results`, `plot`) for data management.  
- **Load Dataset:** Imports the dataset (`patient_info.csv`) using `read.csv()` and inspects its contents.  
- **Inspect Data Structure:** Uses `str()` and `class()` to examine variables and data organization.  
- **Convert Variables:** Converts categorical variables (`gender`, `diagnosis`, `smoker`, `patient_id`) to factors.  
- **Create Derived Variable:** Generates a binary variable `smoking_status` (Yes = 1, No = 0).  
- **Save Cleaned Data:** Exports the cleaned dataset to `clean_data/patient_info_clean.csv`.
- 
**Task 2** — [R Script](https://github.com/aymunir1/AI_Omics_Internship_2025/edit/main/Yusuf_Munir_Aliyu_class_2_Assignment.R#L15C0)
- Defines a custom function `classify_gene()` to evaluate each gene based on `logFC` and `padj` values.  
- **Classification Criteria:**  
  - `logFC > 1` and `padj < 0.05` → Upregulated  
  - `logFC < -1` and `padj < 0.05` → Downregulated  
  - Otherwise → Not Significant  
- Processes two datasets (`DEGs_Data_1.csv` and `DEGs_Data_2.csv`) in a loop.  
- Replaces missing `padj` values with 1 to prevent computation errors.  
- Adds a `status` column to each dataset for classification results.  
- Saves processed files into the `Result` folder.  
- Displays summary counts of gene expression status using frequency tables.
---

### **Module II**: Introduction to Genomics Data Analysis  

**Task 3**  
- Retrieve microarray datasets from **ArrayExpress** and **NCBI GEO** for downstream analysis.

**Task 4 — Microarray Data Preprocessing Workflow**  
**1️⃣ Quality Control (Pre- and Post-Normalization)**  
- Perform **quality control (QC)** on raw expression data to assess overall array performance.  
- Identify and flag **outlier arrays** using diagnostic plots (boxplots, MA plots, PCA).  
- Document the **number of outliers detected** before and after normalization.
**2️⃣ Normalization and Probe Filtering**  
- Apply appropriate **normalization** (e.g., RMA, quantile normalization) to correct technical variation.  
- Filter out **low-intensity or non-informative probes** to enhance reliability.  
- Record the **number of transcripts retained** after filtering.
**3️⃣ Phenotype Group Definition**  
- Use **phenotype metadata** to define biological groups (e.g., *Normal* vs *Cancer*).  
- Relabel or encode samples to ensure consistency for **differential expression analysis**.
