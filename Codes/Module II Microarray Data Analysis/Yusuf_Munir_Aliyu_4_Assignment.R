
# Install and Load Required Packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Bioconductor Packages
BiocManager::install(c("GEOquery", "affy", "arrayQualityMetrics"))

# Load Required Libraries

library(GEOquery)                 # Download GEO data (series, matrix, raw CEL files)    
library(affy)                     # Pre-Processing of Affymetrix microarray data (RMA nor)
library(arrayQualityMetrics)      # QC report for microarray data

#----------------------------------------------------------------
#            Download Series Matrix Files
#----------------------------------------------------------------

gse_data_4 <- getGEO("GSE99807", GSEMatrix = TRUE)

# Extract expression data
expression_data_4 <- exprs(gse_data_4$GSE99807_series_matrix.txt.gz)

# Extract feature data
feature_data_4 <- fData(gse_data_4$GSE99807_series_matrix.txt.gz)

# Extract phenotype data
phenotype_data_4 <- pData(gse_data_4$GSE99807_series_matrix.txt.gz)


# finding missing values 
sum(is.na(phenotype_data_4$source_name_ch1))


#--------------------------------------------------------------------
#        Download Raw Data (CEL files)
#---------------------------------------------------------------------
# Download Raw Data
getGEOSuppFiles("GSE99807", baseDir = "Raw_Data", makeDirectory = TRUE)


# Untar CEL file
untar("Raw_Data/GSE99807_RAW.tar", exdir = "Raw_Data/CEL_Files_4")

# Read CEL
raw_data_4 <- ReadAffy(celfile.path = "Raw_Data/CEL_files_4")
raw_data_4


#--------------------------------------------------------------
#          QC Before Pre-Processing
#--------------------------------------------------------------

# boxplots 
# MA plots
# Density plots

arrayQualityMetrics(expressionset = raw_data_4,
                    outdir = "QC_Report_4a",
                    force = TRUE,
                    do.logtransform = TRUE) 

#--------------------------------------------------------------
#              RMA Normalization 
#--------------------------------------------------------------

normalized_data_4 <- rma(raw_data_4)

# QC After normalization 
arrayQualityMetrics(expressionset = normalized_data_4,
                    outdir = "normalized_data_4a",
                    force = TRUE)


# Extract normalized expression 
processed_data_4 <- as.data.frame(exprs(normalized_data_4))
dim(processed_data_4)


#-------------------------------------------------------------------------------
#              Filter Low Variance Transcription 
#-------------------------------------------------------------------------------

row_median_4 <- rowMedians(as.matrix(processed_data_4))
row_median_4

# Plot the distribution of median intensities of probes
hist(row_median_4,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")


threshold <- 3.5

abline(v = threshold, col = "black", lwd = 2)

index <- row_median_4 > threshold


filtered_data <- processed_data_4[index, ]

colnames(filtered_data) <- rownames(phenotype_data_4)

processed_data_4 <- filtered_data


#--------------------------------------------------------------
#              Phenotype Data 
#--------------------------------------------------------------

class(phenotype_data_4$source_name_ch1)
phenotype_data_4$source_name_ch

groups <- factor(phenotype_data_4$source_name_ch1,
                 levels = c("NASH-HCC liver tissues", "precarcinoma-NASH-HCC liver tissues"),
                 label = c("normal", "cancer"))


class(groups)
levels(groups)


# Save workspace
save.image(file = "GSE79940.RData")

