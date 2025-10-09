## Practical Sessions  

### Module I: Basics of R Programming  
1. **Class IA:** R Program Installation
3. **Class IB:** R Basic Operations | [R Script ](https://github.com/aymunir1/AI_Omics_Internship_2025/blob/main/Class_1b.R)
5. **Class IC:** R Basic Syntax | [R Script with Practice Excercises](https://github.com/aymunir1/AI_Omics_Internship_2025/blob/main/Class_1c.R)
7. **Class 2:** Operators in R | Data STructures in R  | User Define Functions | for-Loop | [R Script](https://github.com/aymunir1/AI_Omics_Internship_2025/edit/main/Class_2.R)

**Task 1**: [R Script](https://github.com/aymunir1/AI_Omics_Internship_2025/blob/main/Yusuf_Munir_Aliyu_Assignment%201b.R)
Create Subfolders — Sets up project directories (raw_data, clean_data, scripts, results, plot) for data management.
Load Dataset — Imports the dataset (patient_info.csv) using read.csv() and views its contents.
Inspect Data Structure — Uses str() and class() to examine variable types and dataset organization.
Convert Variables — Converts categorical variables (gender, diagnosis, smoker, patient_id) to factors for proper analysis.
Create Derived Variable — Generates a binary variable smoking_status (Yes = 1, No = 0).
Save Cleaned Data — Exports the cleaned dataset to clean_data/patient_info_clean.csv.

**Task 2**: [R Script](https://github.com/aymunir1/AI_Omics_Internship_2025/edit/main/Yusuf_Munir_Aliyu_class_2_Assignment.R#L15C0)
A custom function classify_gene() is defined to evaluate each gene based on logFC and padj values.
Genes with logFC > 1 and padj < 0.05 are labeled Upregulated.
Genes with logFC < -1 and padj < 0.05 are labeled Downregulated.
All others are classified as Not Significant.
Two datasets (DEGs_Data_1.csv and DEGs_Data_2.csv) are processed in a loop.
Missing adjusted p-values are replaced with 1 to prevent computation errors.
Each dataset receives a new column, status, indicating the classification result.
The processed files are saved into the Result folder for further analysis.
Summary counts of upregulated, downregulated, and non-significant genes are displayed using frequency tables.
 
### Module II: Introduction to Genomics Data Analysis 
**Task 3**

**Task 4**
