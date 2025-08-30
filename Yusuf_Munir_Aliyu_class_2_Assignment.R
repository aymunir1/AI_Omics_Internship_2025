



#-------------------------------------------
#             Assignment 2
#-------------------------------------------
# Write a function classify_gene() 

# that takes:
#   - logFC (log2FoldChange)
#   - padj  (adjusted p-value)

# and returns:
#   - "Upregulated" if log2FC > 1 and padj < 0.05
#   - "Downregulated" if log2FC < -1 and padj < 0.05
#   - "Not_Significant" otherwise


# Then:
#   - Apply it in a for-loop to process both datasets (DEGs_data_1.csv, DEGs_data_2.csv)
#   - Replace missing padj values with 1
#   - Add a new column 'status'
#   - Save processed files into Results folder
#   - Print summary counts of significant, upregulated, and downregulated genes
#   - Use table() for summaries





# ---------------------------
# Assignment 2: Differential Gene Expression (DGE) Analysis
# ---------------------------

# 1. Define classify_gene() function
classify_gene <- function(logFC, padj) {
  if (is.na(padj)) padj <- 1   # Replace missing padj with 1
  if (logFC > 1 & padj < 0.05) {
    return("Upregulated")
  } else if (logFC < -1 & padj < 0.05) {
    return("Downregulated")
  } else {
    return("Not_Significant")
  }
}

# 2. List of datasets to process
files <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")

# 3. Loop through datasets
for (file in files) {
  # Read dataset
  DEG <- read.csv(file)
  
  # Replace missing padj with 1
  DEG$padj[is.na(DEG$padj)] <- 1
  
  # Apply classify_gene() and add new column 'status'
  DEG$status <- mapply(classify_gene, DEG$logFC, DEG$padj)
  
  # Save processed file into Results folder
  output_file <- paste0("Result/Processed_", file)
  write.csv(df, output_file, row.names = FALSE)
  
  # Print summary counts
  cat("\nSummary for", file, ":\n")
  print(table(DEG$status))
}

