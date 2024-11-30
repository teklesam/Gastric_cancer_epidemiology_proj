###############################
# Epidemiology of gastic cancer
##############################

# 
# Written by : Samuel Tekle,MD
# Copyright (c) - 2024 Samuel Tekle
# 
# Started on : 15.11.2024
# Version control:
# Objectives:
#Tasks (Samuel)

#1. Initial study subject (n= )
#2. Eligible study subjects (n= )
#3. Final study cohort (n= )
#4. Prepare a flow chart for study subject selection
#5. List of all variables in the GC registry
#6. Data cleaning and variable standardization (outliers adjustments, 
#                                               missing values handling, 
#                                              data imputation, 
#                                              data censuring, 
#                                             variable coding, and 
#                                              re-coding)
#7. Frequency table of categorical variables
#8. Descriptive table of continuous variables
#9. Tabulation of study subjects according to their gender, age category, and
#  geographic location
#10. Cross Tabulation of GC with age categories
#11. Cross Tabulation of GC with sex
#12. Cross Tabulation of GC with geo-region
#13. Cross Tabulation of GC with age and sex
#14. Cross Tabulation of GC with geo-region and sex
#15. Preliminary writing and preparation of template for result section of 
#  GC study manuscript


####################################





# Step 1: Load Required Packages

library(pacman)  # Package manager to load multiple libraries at once
pacman::p_load(
  here, 
  rio, 
  tidyverse, 
  janitor, 
  stringr, 
  labelled, 
  dplyr, 
  purrr,
  skimr
)

# Step 2: Import and Clean Data
# Adjust the file path if necessary
gastric_data_original <- "gastric_biopsy_raw.csv"

# Import the data
preliminary_gastric_data <- import(here("1_data", gastric_data_original))

# Clean variable names and de-identify cases
preliminary_gastric_data <- preliminary_gastric_data %>%
  clean_names() %>%
  mutate(case_id = row_number())  # Assign sequential IDs for de-identification


# Step 3: Standardize 'symptomes' Variable
# Replace missing or empty symptom values with "None"
secondary_gastric_data <- preliminary_gastric_data %>%
  mutate(symptomes = if_else(is.na(symptomes) | symptomes == "", 
                             "None documented", symptomes))

# Define symptom mapping
symptom_mapping <- list(
  "epigastric_pain" = c("abd__Dist__EGP", "EGP", "EGP__nausea", "Depression___EGP"),
  "abdominal_pain" = c("Abdo_pain", "abdominal_pain", "pain_abd_", "pain_n_vomiting"),
  "vomiting" = c("vomit", "vomitin", "vomiting", "Vomiting", "vomitinng", "Vomiitng", "bleeding_vomiting"),
  "nausea" = c("nausea"),
  "diarrhea" = c("diaarhea", "diarhea", "diarhhea", "diarrhea", "watery_diarrhea"),
  "weight_loss" = c("weight_loss", "Weight_loss", "DYSPHAGIA_AND_WEIGHT_LOSS"),
  "dysphagia" = c("dysphagea", "Dysphagea", "Dysphagia", "DYSPHAGIA"),
  "melena" = c("melana", "melena"),
  "blood_vomit" = c("blood_vomit", "Blood_vomit", "hemetemesis", "Hemetemisis", "hemetesis"),
  "peptic_ulcer" = c("peptic_ulcer", "ulcer", "stomach_ulcer"),
  "GERD" = c("GERD"),
  "malaise" = c("malaise", "malise"),
  "polyp" = c("polyp"),
  "headache" = c("Headache"),
  "swelling" = c("swelling", "edema"),
  "gastrocutaneous_fistula" = c("gastrocutaneous_fistula"),
  "apetite_loss" = c("apetite_loss"),
  "bleeding" = c("bleeding"),
  "None documented" = "Non-documented"
  
)

# Create a mapping vector
mapping_vector <- unlist(lapply(names(symptom_mapping), function(key) {
  setNames(rep(key, length(symptom_mapping[[key]])), symptom_mapping[[key]])
}))

# Map and standardize symptoms

secondary_gastric_data <- secondary_gastric_data %>%
  mutate(
    symptomes_standardized = str_split(symptomes, ";") %>%
      map(~ str_trim(.x) %>%
            map_chr(~ mapping_vector[.x] %||% .x) %>%
            unique() %>%
            str_c(collapse = ";"))
  )



# Step 5: Create Binary Columns for Symptoms
# Extract unique standardized symptoms
all_symptoms_standardized <- unique(unlist(str_split(secondary_gastric_data$symptomes_standardized, ";")))
all_symptoms_standardized <- all_symptoms_standardized[!is.na(all_symptoms_standardized) & all_symptoms_standardized != ""]



# Ensure column names are clean
clean_col_names <- str_replace_all(all_symptoms_standardized, "[^a-zA-Z0-9]", "_")

# Add binary columns for each symptom
for (symptom in all_symptoms_standardized) {
  col_name <- str_replace_all(symptom, "[^a-zA-Z0-9]", "_")
  secondary_gastric_data[[col_name]] <- str_detect(
    secondary_gastric_data$symptomes_standardized,
    fixed(symptom, ignore_case = TRUE)
  )
}

# Convert binary variables to Yes/No
secondary_gastric_data <- secondary_gastric_data %>%
  mutate(across(
    all_of(clean_col_names),
    ~ if_else(.x, "Yes", "No")
  ))

# Gastric age restructring to population and standard age cuts
secondary_gastric_data <- secondary_gastric_data %>%
  mutate(age_group = cut(age, 
                         breaks = seq(0, 105, by = 5),  # Adjusted to create 21 intervals
                         right = FALSE, 
                         labels = paste(seq(0, 100, by = 5), 
                                        seq(4, 104, by = 5), sep = "-")))

colnames (secondary_gastric_data)

## Defining variable types and labeling

secondary_gastric_data <- secondary_gastric_data %>% 
  mutate(
    case_id                        = as.character (case_id),
    age                            = as.numeric (age),
    age_group                      = as.factor(age_group),
    calender_year                  = as.numeric (calender_year),
    gender                         = as.factor (gender),
    address                        = as.character (address),
    address_zoba                   = as.factor (address_zoba),
    specimen                       = as.factor (specimen),
    symptomes                      = as.character (symptomes),
    duration                       = as.numeric (duration),
    malignant_not                  = as.factor (malignant_not),
    report_diagnosis               = as.character (report_diagnosis),
    grade                          = as.factor (grade),
    diagnosis_category             = as.factor (diagnosis_category),
    non_neoplastic_subtype         = as.factor (non_neoplastic_subtype),
    adenocarcinoma_subtype         = as.factor (adenocarcinoma_subtype),
    other_malignancy_subtype       = as.factor (other_malignancy_subtype),
    symptomes_standardized         = as.character (symptomes_standardized),
    dysphagia                      = as.factor (dysphagia),
    malaise                        = as.factor (malaise),
    headache                       = as.factor (headache),
    diarrhea                       = as.factor (diarrhea),
    epigastric_pain                = as.factor (epigastric_pain),
    blood_vomit                    = as.factor (blood_vomit),
    peptic_ulcer                   = as.factor (peptic_ulcer),
    GERD                           = as.factor (GERD),
    bleeding                       =as.factor (bleeding),
    melena                         = as.factor (melena),
    swelling                       = as.factor (swelling),
    polyp                          = as.factor (polyp),
  )

# Labeling the variables

secondary_gastric_data <- secondary_gastric_data %>% 
  labelled::set_variable_labels(
    case_id                        = "Primary Identifier",
    age                            = "Age in years",
    age_group                      = "Population structure",
    calender_year                  = "Calender year",
    gender                         = "Sex",
    address                        = "Address (town)",
    address_zoba                   = "Address (Zoba)",
    specimen                       = "Site of specimen collection",
    symptomes                      = "Clinical history",
    duration                       = "Duration of illness",
    malignant_not                  = "Malignant or NOn-malignant",
    report_diagnosis               = "Diagnosis on the request form",
    grade                          = "Grade of Cancerous lessions",
    diagnosis_category             = "Reclassified diagnosis",
    non_neoplastic_subtype         = "Classified Non-neoplastic diseases",
    adenocarcinoma_subtype         = "Classified Adenocarcinoma Sub-type",
    other_malignancy_subtype       = "Classified Other Malignancy ",
    symptomes_standardized         = "List of symptomes",
    dysphagia                      = "Difficulty of swallowing",
    malaise                        = "General feeling of sick",
    headache                       = "headache",
    diarrhea                       = "Diarrhea",
    epigastric_pain                = "Epigastric pain",
    blood_vomit                    = "Bloody vomitus",
    peptic_ulcer                   = "Diagnosis of PUD",
    GERD                           = "Diagnosis of GERD",
    bleeding                       = "Unspecified bleeding",
    melena                         = "melena",
    swelling                       = "Extremity swelling",
    polyp                          = "Polyp in the EGD",
  )


# Reorder and make variables flow
secondary_gastric_data <- secondary_gastric_data %>% 
  select(case_id, -c(age_cat), 
         age,
         age_group,
         calender_year,
         gender,
         address,
         address_zoba,
         symptomes,
         duration,
         specimen,
         symptomes_standardized,
         dysphagia,
         malaise,
         headache,
         diarrhea,
         epigastric_pain,
         blood_vomit,
         peptic_ulcer,
         GERD,
         bleeding,
         melena,
         swelling,
         malignant_not,
         report_diagnosis,
         diagnosis_category,
         adenocarcinoma_subtype,
         other_malignancy_subtype,
         grade,
         non_neoplastic_subtype
  )

# Step 6: Review and Export Data
# View the updated dataset
glimpse(secondary_gastric_data)

skimr::skim(secondary_gastric_data)

unique(secondary_gastric_data$diagnosis_category)
# Export the cleaned data
export(secondary_gastric_data, here("3_output", "cleaned_gastric_data.rds"))

export(secondary_gastric_data, here("3_output", "cleaned_gastric_data.csv"))

