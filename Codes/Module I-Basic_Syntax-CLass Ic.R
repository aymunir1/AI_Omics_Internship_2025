#                               CLASS 1C


# variables
# store value in R
# <- assignment operator

# gene name "TPS3"

gene <- "TPS3"

# to retrieve value in console
gene

print(gene)

# 2.3, 4.6, 3.6, 7.2, 4.7

expression_level <- c(2.3, 4.6, 3.6, 7.2, 4.7)

# to import data as variable
raw_data <- read.csv(file.choose())

# Rules

# variable name must be start with letter
1gene <- 25 #variable name can`t stand with number

gene1 <- 25 # add number at the end of the variable


#no space is allowed in variable name
sample id <- "s01"

# instead use underscore (_) or dot (.)
sample_id <- "s01"
sample.id <- "s01"

# R is case sensitive
Glucose_level <- 110
glucose_level <- 110

# R overwrite variable without any warning
glucose_level <- c(110, 90, 120)

data <- raw_data
raw_data$patient_id <- NULL #this code remove patient_id column

# For data and cleaning and transforming create a variable for that data

clean_data <- data[,-1] 
# this code deleted the patient_id column & it assign a new variable

# comments are for our own understanding R doesn`t consider it as code
# data_2 <- 23
data_2 <- 23

# pro tips: turn comments into heading
#### Heading 1 ####
#### Heading 2 ####

# keywords
# These are reserved words in for specific function
#if, else, TRUE, FALSE, for so on ...........

help ("reserved")
help ("mean")
?median
median

# sort values from largest to smallest
sorted_age <- sort(raw_data$age, decreasing = TRUE)
sorted_age
raw_data

# sort values from smallet to largest
sorted_age2 <- sort(raw_data$age, decreasing = FALSE)
sorted_age2

# if & else, which are used for creating logical conditions

gene_expression <- 30

if(gene_expression > 50) {
  print("Gene is highly expressed")
}

# here if is the keyword that check the condition if gene_expressiom
# condition
# in case if the condition is false

if(gene_expression < 50) {
  print ("Gene expression is low")
}else {
  print("Gene is highly expression")
}

# you cannot use keywords as variable name
if <- 28

# for loop: used to repeat same tasks multiple times
# let say we want to convert data type of multiple in the our
str(raw_data)

# gender is categorical data type
# It should be in a factor format
# Gender column from chr to factor
# diagnosis: cancer/normal it is also a categorical variable
# Smoker: chr from factor

#Coverting multiple categorical variables to factor
# Instead of manually conversion we will use this for loop function for all 3 
#column (Gender, Diagnosis, and smoker) with one command

# To convert raw_data column into factor
# I want to ave output in clea_data

# to create a copy of raw_data with name clean_data
clean_data <- raw_data
str(clean_data)

# to convert column automatically into factor
# create for loop function
for (i in 1:ncol(clean_data)) {
  if(is.character(clean_data[[i]])){
    clean_data[[i]] <- as.factor(clean_data[[i]])
  }
  
}
str(clean_data)
str(raw_data)


write.csv(clean_data,file = "results/clean_patient_info_data.csv")
save.image(file = "class_1c_Workspace.RData")
load(class_1c_Workspace.RData"")



