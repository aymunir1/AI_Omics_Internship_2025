


getwd()
setwd("C:/Users/User/Desktop/AI_Omics_Internship_2025/Module_I")

# Creating subfolders
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


# creating new variable for smoking statusas Yes = 1, No = )
data$smoking_status <- ifelse(data$smoker == "Yes", 1, 0)
str(data$smoking_status)
View(data)



write.csv(data, file = "clean_data/patient_info_clean.csv")

save.image(file = "Yusuf_Munir_Aliyu_class_1b_Assignment.RData")
save(data, file = "Yusuf_Munir_Aliyu_class_1b_Assignment.RData")

