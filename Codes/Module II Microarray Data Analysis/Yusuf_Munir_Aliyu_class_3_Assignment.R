


# Install and Load Required Packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Bioconductor Packages
BiocManager::install(c("GEOquery", "affy", "arrayQualityMetrics"))




# Load Required Libraries

library(GEOquery)                 # Download GEO data (series, matrix, raw CEL files)    
library(affy)                     # Pre-Processing of Affymetrix microarray data (RMA nor)
library(arrayQualityMetrics)      # QC report for microarray data



normalized_data <- rma(raw_data)



#----------------------------------------------------------------
#            Download Series Matrix Files
#----------------------------------------------------------------


gse_data_3 <- getGEO("GSE79940", GSEMatrix = TRUE)



gse_data2 <- getGEO("GSE79972", GSEMatrix = TRUE)


# Extract expression data

expression_data_a <- exprs(gse_data_3$`GSE79940-GPL96_series_matrix.txt.gz`)
expression_data_b <- exprs(gse_data_3$`GSE79940-GPL97_series_matrix.txt.gz`)

expression_data <- exprs(gse_data2$GSE79972_series_matrix.txt.gz)

# Extract feature data

feature_data_a <- fData(gse_data_3$`GSE79940-GPL96_series_matrix.txt.gz`)
feature_data_b <- fData(gse_data_3$`GSE79940-GPL97_series_matrix.txt.gz`)

# Extract phenotype data
phenotype_data_a <- pData(gse_data_3$`GSE79940-GPL96_series_matrix.txt.gz`)
phenotype_data_b <- pData(gse_data_3$`GSE79940-GPL97_series_matrix.txt.gz`)



# finding missing values 
sum(is.na(phenotype_data_a$source_name_ch1))
sum(is.na(phenotype_data_b$source_name_ch1))

 save.image(file= "GSE79940.RData")



#--------------------------------------------------------------------
#        Download Raw Data (CEL files)
#---------------------------------------------------------------------

getGEOSuppFiles("GSE79940", baseDir = "Raw_Data", makeDirectory = TRUE)


# Untar CEL file
untar("Raw_Data/GSE79940_RAW.tar", exdir = "Raw_Data/CEL_Files_3")




# Read CEL

raw_data_3 <- ReadAffy(celfile.path = "Raw_Data/CEL_files_3")
raw_data_3


#--------------------------------------------------------------
#          QC Before Pre-Processing
#--------------------------------------------------------------

# boxplots 
# MA plots
# Density plots



arrayQualityMetrics(expressionset = raw_data2,
                    outdir = "QC_Report2",
                    force = TRUE,
                    do.logtransform = TRUE) 



arrayQualityMetrics(raw_data2, outdir = "QC_Report2", force = TRUE, do.logtransform = TRUE)



#--------------------------------------------------------------
#              RMA Normalization 
#--------------------------------------------------------------

normalized_data2 <- rma(raw_data)




# QC After normalization 

arrayQualityMetrics(expressionset = normalized_data[, c(1, 8, 9, 11, 20)],
                    outdir = "normalized_data2",
                    force = TRUE)


# Extract normalized expression 

processed_data <- as.data.frame(exprs(normalized_data))
dim(processed_data)


#-------------------------------------------------------------------------------
#              Filter Low Variance Transcription 
#-------------------------------------------------------------------------------

row_median <- rowMedians(as.matrix(processed_data))

row_median

# Plot the distribution of median intensities of probes

hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")


threshold <- 3.5

abline(v = threshold, col = "black", lwd = 2)

index <- row_median > threshold

filtered_data <- processed_data[index, ]

colnames(filtered_data) <- rownames(phenotype_data)

processed_data <- filtered_data

#--------------------------------------------------------------
#              Phenotype Data 
#--------------------------------------------------------------

class(phenotype_data$source_name_ch1)
phenotype_data$source_name_ch

groups <- factor(phenotype_data$source_name_ch1,
                 levels = c("gastric mucosa", "gastric adenocarcinoma"),
                 label = c("normal", "cancer"))


class(groups)
level(groups)


save.image(file = "GSE79973.RData")


#--------------------------------------------------------------
#              Probe ID to Gene Mapping (Collape to ENTREZID) 
#--------------------------------------------------------------


# extract probe ids from processed data


# remove unmapped probes (no ENTREZID)

# collapse multiple probes to single ENTREZID (average expression)