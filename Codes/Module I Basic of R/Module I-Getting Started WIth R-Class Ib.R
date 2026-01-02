#------------------------------------------------------------------
#                  Class 1c   
#------------------------------------------------------------------

# Organize Project folder
# data
# script
# result

# Create folder
dir.create("data")
dir.create("script")
dir.create("result")

# How R code works
# function(your argument/input)

print("Biotechnology")
gene_expression <- c(2.3, 3.8, 3.9, 5.6, 9.4)

# Varible store values
# gene_expression is variable name
# <- assignment operators to store values into variable

meam_result <- mean(gene_expression)
plot(gene_expression)
hist(gene_expression)
barplot(gene_expression)

#check basic stats
summary(gene_expression)

# R Data types

# 1 Numeric witth decimal point positive or negative
x <- 42
class(x)

y <- 40
class(y)

#2 integer/whole number
#R classify decimal& whole num as mumeric

z <- 24L
class(z)

var_1 <- c(28.4, 23.7, 2,9)
class(var_1)
as.integer(var_1)
class(var_1)

# convert numeric into integer
var_2 <- as.integer(var_1)
class(var_2)

# convert numeric into integer
var_3 <- as.numeric(var_2)
class(var_3)

# 3 character/string data type
var_4 <- c("gene1", "gene2", "gene3")
class(var_4)

Var_5 <- c(gene1, "gene2", "gene3")

var_6 <- c("3.8", "4.6", "5.9")
class(var_6)
mean(var_6)

var_7 <- c(3.8, 4.6, 5.9)
class(var_7)
mean(var_7)

# 4 Factor or categorical data
# either written as numerical (0,1) or character ("cancer" or "normal")

disease_status <- c("cancer", "normal", "cancer", "cancer", "normal")
class(disease_status)
disease_status

#convert class into  factor
disease_status <- as.factor(disease_status)
class(disease_status)
disease_status
View(disease_status)

disease_status <- factor(disease_status,
                         levels = c("normal","cancer"))
disease_status

# 5 logical data

age <- c(23,25,28)
var_8 <- age < 25
var_8

# import csv file
data <- read.csv(file.choose())
View(data)
str(data)

# convert height into factor
data$height_fac <- as.factor(data$height)
str(data)

# relevel factor
data$height_fac <- factor(data$height_fac,
                          levels = c("Tall","Medium","Short"))
levels(data$height_fac)

data$gender_fac <- as.factor(data$gender)
str(data)

# Convert factor into numeric factor
data$gender_num <- ifelse(data$gender_fac == "Female", 1, 0)
class(data$gender_num)

data$gender_num <- as.factor(data$gender_num)
class(data$gender_num)

# Save file as csv format
write.csv(disease_status,file = "result/patient_data.csv")

# Save workspace

# save entire workspace
save.image(file = "full_workspace.RData")

# Save specific object
save(gene_expression, disease_status, file = "workspace.RData")

load("workspace.RData")
load("full_workspace.RData")


