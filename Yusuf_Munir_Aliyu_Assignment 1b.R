

############################### ASSIGNED TASK 1   ###########################

# ---------------------------------------------------------------------------
# 1. Set Working Directory
# Create a new folder on your computer "AI_Omics_Internship_2025".

# 2. Create Project Folder
# In RStudio, create a new project named "Module_I" in your "AI_Omics_Internship_2025" folder.
# Inside the project directory, create the following subfolders using R code:
# raw_data, clean_data, scripts, results or Tasks, plots etc

# 3. Download "patient_info.csv" dataset from GitHub repository
# load the dataset into your R environment
# Inspect the structure of the dataset using appropriate R functions
# Identify variables with incorrect or inconsistent data types.
# Convert variables to appropriate data types where needed
# Create a new variable for smoking status as a binary factor:
  # 1 for "Yes", 0 for "No"
# Save the cleaned dataset in your clean_data folder with the name patient_info_clean.csv
# Save your R script in your script folder with name "class_Ib"
# Upload "class_Ib" R script into your GitHub repository
# ---------------------------------------------------------------------------


#  Creating subfolders
dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("plot")

# Loading dataset
data <- read.csv(file.choose())
View(data)

# Inspecting the structure of the dataset 
str(data)

class(data$patient_id)
class(data$age)
class(data$gender)
class(data$diagnosis)
class(data$bmi)
class(data$smoker)


# Coverting to appropriate data types

# Converting gender, diagnosis and smoker varibale to factor 
data$gender <- as.factor(data$gender)
class(data$gender)
str(data)
levels(data$gender)

data$diagnosis <- as.factor(data$diagnosis)
class(data$diagnosis)
str(data)
levels(data$diagnosis)

data$smoker <- as.factor(data$smoker)
class(data$smoker)
levels(data$smoker)
str(data)

data$patient_id <- as.factor(data$patient_id)
levels(data$patient_id)
str(data$patient_id)
View(data$patient_id)

# creating new variable for smoking statusas Yes = 1, No = 0)
data$smoking_status <- ifelse(data$smoker == "Yes", 1, 0)
str(data$smoking_status)
View(data)


write.csv(data, file = "clean_data/patient_info_clean.csv")

save.image(file = "Yusuf_Munir_Aliyu_class_1b_Assignment.RData")



