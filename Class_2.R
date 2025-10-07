
#            Assignment Operators

# Used to store values inside varibales

# <-  (Most common, rightward assignment operators)
height <- c(1.75, 1.76, 1.82, 1.67)

# -> (Same as above, but leftward assignment operator)
c(68, 78, 85, 75) -> Weight

# = (Also assigns, used in function arguments)
smoking_status = c("Yes", "No", "No", "Yes")


#          Arithmetic operators


# Perform basic math i.e
#   + addition
#   - subtraction
#   * multiplication
#   / division
#   ^ exponential


# Let`s BMI using weight and height

BMI <- Weight/(height^2)
BMI


#         Comparison operators

# Comparison Operators ask logical questions about values.
# They return output as TRUE or FALSE
# They don`t calculate answer, they compare

#    > greater than
BMI > 25
#    < lessthan 
BMI < 18.5



#   >    greater than
#   <    less than
#   >=   greater than equal to
#   <=   lessthan equal to
#   ==   equal to
#   !=   not equal to



# In R
# yes = TRUE
# no = FALSE


#         Logical operators

# Logical operators let us combine conditions:

# & AND (both must be TRUE)
# Is the patient overweight AND a smoker
# BMI cutoff = 25
(BMI > 25) & (smoking_status == "Yes")
(BMI < 25) & (smoking_status == "Yes")



# | OR (at least one must be TRUE)
# Is the patient overweight OR a smoker

(BMI >5) | (smoking_status == "Yes")
BMI
smoking_status



# ! NOT (reverse the condition)

# Is the patient NOT a smoker
smoking_status== "No"
# condition = yes
# output = FALSE


# -------------------------------
# 2.      Data structure in R
# -------------------------------
# Data structure are how R organizes and store information.

# Main structures we will use again and again
#   1. Vectors
#   2. Lists
#   3. Matrices
#   4. Data Frame


# ------------------
# Vectors
# ------------------
# simplest data structure in R
# It stores a sequence of value, but all of the same type.


#   - Numeric vector
num_vec <- c(1, 2, 3, 4)
class(num_vec)

# numeric vectors used to perform mathematical calculation

#   - Character vector
chr_vec <- c("gene1", "gene2", "gene3")

#   - Logical vector
logical_vec <- c(TRUE, FALSE, TRUE)

mix_vec <- c("gene1", 1, "gene2", 2)
mean(mix_vec)

# we can extract values from vectors using indexing with {}
num_vec[2]
num_vec[2:4]

# we can only combine vectors of eqa]ual sequence

# we cam treat vectors as column or as rows

# by column 
vec_col <- cbind(num_vec, chr_vec)
View(vec_col)

#--------------
# Lists
#--------------
# Unlike vectors, a list can hold multiple types together:
# numbers, text, logical even other data frames


all_vec <- list(num_vec, chr_vec, logical_vec)
View(all_vec)

# save your raw_data
#save processed data
# results

# we access element with [[]]
all_vec[[2]]


# ----------------
# Matrices
# ----------------
# A matrice is a 20 structure (rows x columns).
# All values must be the same types (usually numeric).
# Example: gene expression matrix where rows are genes 
# and columbns are samples

my_matrix <- matrix(1:9, nrow = 3, ncol = 3)
my_matrix

# we access elements with [row, column]
# we change using byrow = TRUE


my_matrix <- matrix(1:9, nrow = 3, ncol = 3, byrow = TRUE)

# we can access element with [row, column]
my_matrix[2,3]
my_matrix[2,]
my_matrix[#rows, columns]
  
    
#--------------------- 
# Data Frames
#---------------------  
# A data frame is the the most important structure 
# for real dataset
# Each column can be of a different type: numeric, character, factor

data <- data.frame(
 patient_id = c("P1", "P2", "P3"),
 age = c(65, 78, NA),
 diagnosis = c("cancer", "diabetes", "cancer")
)

print(data)

#---------------------  
# Data Assessment
#---------------------  

# Function like str(), head(), names() help us to explore
# the dataset before analysis.

str(data) # structure of dataset
head(data) # first 6 rows
head(data, n=2) # first 2 rows
tail(data) # last 6 rows
tail(data, n=2) #last 2 rows
dim(data) # rows & columns
names(data) # columns


# Data frame are indexed like matrix with more flexible
# access a column directly
data$patient_id
data[c(1,3), c(2,3)]


# Access data using:
data$patient_id     # extract single column
data[1:2, c(1,3)]   # extract specific rows and columns

# Create new columns:
data$new_column <- c(1, 2, 3)

# ----------------
# Missing Values
# ----------------
# Real data often contains missing values (NA).
# You must check and handle them before analysis.

is.na(data)                # identify missing values
sum(is.na(data))           # total missing values
colSums(is.na(data))       # missing values per column
rowSums(is.na(data))       # missing values per row

# Ways to handle NA:

# remove rows with NA
clean_data1 <- na.omit(data)   

# remove columns with NA
clean_data_2 <- data[, colSums(is.na(data))==0]

# replace NA value with 0
clean_data_3 <- data
clean_data_3[is.na(clean_data_3)] <- 0

# replace NA value with mean
clean_data_4 <- data
clean_data_4[is.na(clean_data_4)] <- mean(data$age, na.rm = TRUE)

# ------------------------------
# Summary of Data Structures:
# ------------------------------

#   - Vectors: simple sequences of same data type
#   - Lists: mix of different data types
#   - Matrices: numeric tables
#   - Data Frames: mixed-type tables 


#--------------------
# 3. Functions in R
#--------------------
# Functions let us wrap code into reusable blocks.

# function is  reusable block of code 
# Why use functions?
#   - Avoid repetition
#   - Organize and simplify code
#   - Reuse across projects (save it for later use)
#   - Share with others

# A function in R has 4 key parts:
#   1. Name         -> the name you give to the function
#   2. Arguments    -> the inputs you provide to the function
#   3. Body         -> the set of operations the function performs
#   4. Return Value -> the output the function gives back

# Example: A function to calculate Body Mass Index (BMI)

# 1. Function Name: Calculate_BMI
# 2. Arguments: weight (in kg), height (in meters)
# 3. Body: performs BMI calculation e.g   # Formula: BMI = weight / (height^2)
# 4. Return Value: the BMI value

Calculate_BMI <- function(weight, height) {
# Perform the BMI calculation
Bmi <- Weight/(height^2)
Bmi
  
Bmi# Return the BMI value
return(Bmi)
}

# Call the function by naming arguments explicitly
calculate(weight = 60, height = 1.5)

# Call the function using variables as arguments
calculate_BMI(weight = weight, height = height)

# If a function expects two arguments, you must provide both
# This will give an error because 'height' is missing
calculate_BMI(60) 

# You can assign default values to function arguments
Calculate_BMI <- function(weight, height = 1.75) {
  # Perform the BMI calculation
  Bmi <- weight / (height ^ 2)
  
  # Return the BMI value
  return(Bmi)
}

# In this case, if you don’t provide height, R automatically uses 1.75 as the default.
Calculate_BMI(weight = 60)

# ----------------------------
# Lazy evaluation in R
# ----------------------------
# If your function has three arguments, but the body only uses two,
# R does not force you to supply the third argument for the calculation.
# Example: 'age' is defined as an argument, but not used in the formula.

# Define the function with three arguments
calculate_BMI <- function(weight, height, age){
  # Perform the BMI calculation
  bmi <- weight / (height^2)
  
  # Return the BMI value
  return(bmi)
}

# Here we pass only 'weight' and 'height'
# Even though 'age' exists as an argument, it is ignored because it is not used
calculate_BMI(60, 1.65)

#---------
# Summary:
#---------
# Functions help us package logic once and apply it to different inputs.

# ----------------------------------
# 4. Automating Workflows with for-Loop
# ----------------------------------
# Suppose you have multiple datasets and you want to:
#   - import them,
#   - check missing values,
#   - clean columns,
#   - compute BMI,
#   - and save results.
#
# Instead of repeating steps for each file, we use loops.

# Typical loop workflow:
#   1. Define input and output folders
input_dir <- "Raw_Data"
output_dir <-  "Result"


# create output folder if not already exist
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

# the Result folder if we be created if already not exists


#   2. List which files to process
file_to_process <- c("raw_data/BMI_ data_1.csv", "raw_data/BMI_ data_1.csv")


# 3. List which files to process
result_list <- list()

#   4. For each file 
#    - import data
#    - handle data
#    - calculate BMI using calculate_BMI function
#    - save results (both as csv inside R list)

for (file_names in file_to_proces) {
  cat("\nProcessing;", file_names, "\n")
  
  input_file_path <- file.path(input_dir, file_names)
  
  
# Import dataset
 data <- read.csv(input_file_path, header = TRUE)
 cat("File imported. checking for missing values...\n")
  

# handling missing values

if("height" %in% names(data)){
  missing_count <- sum(is.na(data$height))
 }
  
  cat("Missing values in 'height':", missing_count, "\n")
  data$height[is.na(data$height)] <- mean(data$height, na.rm = TRUE)
}

if("weight" %in% names(data)){
  missing_count <- sum(is.na(data$weight))
  
  cat("Missing values in 'weight':", missing_count, "\n")
  data$weight[is.na(data$weight)] <- mean(data$weight, na.rm = TRUE)
}
# calculate BMI
data$bmi <- calculate_BMI(data$weight, data$height)
cat("BMI has been calculated successfully.\n")

# save results in R
result_list[[file_names]] <- data 

# save results in Results folder
output_file_path <- file.path(output_dir, paste0("BMI_results", file_names))
write.csv(data, output_file_path, row.names = FALSE)
cat("Results saved to:", output_file_path, "\n")
}

# The loop repeats until all files are processed.

results_1 <- result_list[[1]] 
results_2 <- result_list[[2]]

# --------
# Summary:
# --------
# Loops automate repetitive work 
# making your workflow faster 
# consistent, and reproducible








# 3. For each file:
# 4.   
  
# -----------------------
# Typical loop workflow:
# -----------------------

# Define the input folder (where raw data files are stored) and the output folder (where results will be saved).
# Specify the list of files that need to be processed.
# Prepare an empty list in R to store results for later use.
# For each file in the list:
#          Import the data into R.
#          Check and handle missing values (NA).
#          Calculate BMI using the calculate_BMI function.
#          Save the processed results both as a CSV file and in the R results list.

# ----------------------------------------------------
# Calculate BMI of two dataset within loop

# Define input and output folders
input_dir <- "Raw_Data" 
output_dir <- "Results"


# create output folder if not already exist

if(!dir.exists(output_dir)){
  dir.create(output_dir)
}


# List which files to process
files_to_process <- c("data_1.csv", "data_2.csv") 
# These must match exactly with the names in your working folder,
# otherwise R will not find them.

# For practice, you can instead use the provided datasets:
# BMI_data_1.csv
# BMI_data_2.csv
# (download from the GitHub repository).
https://github.com/AI-Biotechnology-Bioinformatics/AI_and_Omics_Research_Internship_2025/tree/main

# Prepare empty list to store results in R 
result_list <- list()

# For each file with in a loop:
#       - import data
#       - handle NA values
#       - calculate BMI using calculate_BMI function
#       - save results (both as CSV and inside R list)


for (file_names in files_to_process) {
  cat("\nProcessing:", file_names, "\n")
  
  input_file_path <- file.path(input_dir, file_names)
  
  # Import dataset
  data <- read.csv(input_file_path, header = TRUE)
  cat("File imported. Checking for missing values...\n")
  
 



  

  
