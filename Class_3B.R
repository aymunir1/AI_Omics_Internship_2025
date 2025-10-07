


# Install and Load Required Packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Bioconductor Packages
BiocManager::install(c("GEOquery", "limma", "affy", "arrayQualityMetrics",
                      "AnnotationDbi", "hgu133plus2.db"))


# CRAN Pacakges
install.packages(dplyr)

# Load Required Libraries

library(GEOquery)                 # Download GEO data (series, matrix, raw CEL files)    
library(affy)                     # Pre-Processing of Affymetrix microarray data (RMA nor)
library(arrayQualityMetrics)      # QC report for microarray data
library(limma)                    
library(AnnotationDbi)           # Interface for annotattion databases
library(hgu133plus2.db)           # Annotation for Affymetrix Human Genome U133 Plus 2.0
library(dplyr)



normalized_data <- rma(raw_data)



#----------------------------------------------------------------
#            Download Series Matrix Files
#----------------------------------------------------------------

gse_data <- getGEO("GSE79973", GSEMatrix = TRUE)
 

# Extract expression data

expression_data <- exprs(gse_data$GSE79973_series_matrix.txt.gz)

# Extract feature data

feature_data <- fData(gse_data$GSE79973_series_matrix.txt.gz)

# Extract phenotype data
phenotype_data <- pData(gse_data$GSE79973_series_matrix.txt.gz)
  


# finding missing values 
sum(is.na(phenotype_data$source_name_ch1))





#--------------------------------------------------------------------
#        Download Raw Data (CEL files)
#---------------------------------------------------------------------

getGEOSuppFiles("GSE79973", baseDir = "Raw_Data", makeDirectory = TRUE)


# Untar CEL file
untar("Raw_Data/GSE79973_RAW.tar", exdir = "Raw_Data/CEL_Files")




# Read CEL

raw_data <- ReadAffy(celfile.path = "Raw_Data/CEL_files")
raw_data


#--------------------------------------------------------------
#          QC Before Pre-Processing
#--------------------------------------------------------------

# boxplots 
# MA plots
# Density plots



arrayQualityMetrics(expressionset = raw_data,
                    outdir = "Results/QC_Raw_Data",
                    force = TRUE,
                    do.logtransform = TRUE) 
 


arrayQualityMetrics(raw_data, outdir = "QC_Report", force = TRUE)



#--------------------------------------------------------------
#              RMA Normalization 
#--------------------------------------------------------------

normalized_data <- rma(raw_data)




# QC After normalization 

arrayQualityMetrics(expressionset = normalized_data[, c(1, 8, 9, 11, 20)],
                    outdir = "normalized_data",
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