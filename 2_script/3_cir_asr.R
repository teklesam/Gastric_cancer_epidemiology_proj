##  NHL gastric  Analysis--------------------------------------------------------

## Title:-----------------------------------------------------------------------

# Analysis of epidemiology of Lymphoma in NHL, Asmara, Eritrea 


# Written by : Samuel Tekle,MD--------------------------------------------------
# Copyright (c) - 2024 Samuel Tekle

# Version control:--------------------------------------------------------------
# Started on : 23.11.2024

# Install required packages
library(pacman)

pacman::p_load(
  tidyr,
  dplyr,
  ggplot2,
  readr,
  epitools,
  apyramid,
  rio,
  here
)


#### Step 1:

# Load datasets

gastric_cancer_data <- import(here("3_output","cleaned_gastric_data.rds")) # Import the cleaned data set

eritrea_population <- import(here("1_data","eritrea_population_2015.csv")) #Import the Eritrean population

who_standard <- import(here("1_data","who_standard.csv"))                 #Import the standard population



# Filter malignant cases and group by age group and gender
malignant_cases_count <- gastric_cancer_data %>%
  filter(malignant_not == "Malignant") %>%
  count(age_group, gender) %>% # Count cases by age group and gender
  rename(cases = n)


malignant_cases <- gastric_cancer_data %>%
  filter(malignant_not == "Malignant") 


apyramid::age_pyramid(data = malignant_cases,
                      age_group = "age_group",
                      split_by = "gender")


apyramid::age_pyramid(data = malignant_cases,
                      age_group = "age_group",
                      split_by = "address_zoba")  




library(dplyr)
library(apyramid)

# Filter the data to include only relevant diagnosis categories and remove NAs
malignant_cases_filtered <- malignant_cases %>%
  filter(diagnosis_category %in% c("Adenocarcinoma", "Lymphoma", "Other Malignancies", 
                                   "Squamous Cell Carcinoma")) %>%
  filter(!is.na(diagnosis_category)) %>%
  # Modify diagnosis_category to retain only the specified levels
  mutate(diagnosis_category = factor(diagnosis_category, 
                                     levels = c("Adenocarcinoma", "Lymphoma", 
                                                "Other Malignancies", "Squamous Cell Carcinoma")))

# Now plot the age pyramid
apyramid::age_pyramid(data = malignant_cases_filtered,
                      age_group = "age_group",
                      split_by = "diagnosis_category")



plot <- apyramid::age_pyramid(data = malignant_cases_filtered,
                              age_group = "age_group",
                              split_by = "diagnosis_category") +
  labs(
    title = "Age Pyramid by Type of Malignancy",  # Modify the title to match the new label
    x = "Age Groups",  # Label for the x-axis
    y = "Counts",  # Label for the y-axis
    fill = "Type of Malignancy",  # Update the legend title for clarity
    color = "Type of Malignancy"  # Update the color legend title
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better visibility
    plot.margin = margin(1, 1, 3, 1, "cm")  # Increase bottom margin for footnote
  ) 

# Ensure 'male' and 'female' data exist in Eritrean population data
# Eritrean population data structure: age_group, male_population, female_population
eritrea_population <- eritrea_population %>%
  mutate(total_population = Male + Female)

# Merge malignant cases with Eritrea population data
age_data <- eritrea_population %>%
  left_join(malignant_cases_count, by = "age_group") %>%
  mutate(
    cases = ifelse(is.na(cases), 0, cases), # Replace NA cases with 0
    crude_rate = (cases / total_population) * 100000 # Crude incidence rate per 100,000
  )
####### Crude incidence rate table

# Ensure 'cases' is numeric and prepare the data
age_gender_table <- age_data %>%
  group_by(age_group, gender) %>%
  summarise(
    cases = sum(cases, na.rm = TRUE),
    total_population = sum(total_population, na.rm = TRUE),
    crude_rate = (cases / total_population) * 100000, # Crude rate calculation
    .groups = "drop"
  )




# Bar Plot of Crude Incidence Rates
ggplot(age_data, aes(x = age_group, y = crude_rate)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "Crude Incidence Rate of Malignant Cancer Cases by Age Group",
    x = "Age Group",
    y = "Crude Incidence Rate (per 100,000)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Line Chart of Crude Incidence Rates
ggplot(age_data, aes(x = age_group, y = crude_rate)) +
  geom_line(color = "darkred", size = 1) +
  geom_point(color = "darkred", size = 3) +
  theme_minimal() +
  labs(
    title = "Crude Incidence Rate of Malignant Cancer Cases by Age Group",
    x = "Age Group",
    y = "Crude Incidence Rate (per 100,000)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Assuming age_data includes a `gender` column for grouped rates
ggplot(age_data, aes(x = age_group, y = crude_rate, fill = gender)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = c("blue", "pink")) +
  theme_minimal() +
  labs(
    title = "Crude Incidence Rate of Malignant Cancer Cases by Age Group and Gender",
    x = "Age Group",
    y = "Crude Incidence Rate (per 100,000)",
    fill = "Gender"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# Merge with WHO standard population
age_data <- age_data %>%
  left_join(who_standard, by = "age_group") %>%
  mutate(
    weighted_rate = crude_rate * weight
  )

# Calculate age-adjusted rate
age_adjusted_rate <- sum(age_data$weighted_rate) / sum(who_standard$standard_population) * 100000

# Results
cat("Crude Incidence Rate (per 100,000):", sum(age_data$cases) / sum(eritrea_population$population) * 100000, "\n")
cat("Age-Adjusted Incidence Rate (per 100,000):", age_adjusted_rate, "\n")



# Incidence rate calculations by years
#R Code for Yearly CIR and ASR

malignant_gastric <- gastric_cancer_data %>% 
  filter(malignant_not == "Malignant")


# Calculate crude incidence rate (CIR) by year
yearly_summary <- malignant_gastric %>%
  group_by(calender_year, gender) %>%
  summarise(
    Count = n(),
    .groups = "drop"
  )

# Add total cases for both sexes
yearly_total <- yearly_summary %>%
  group_by(calender_year) %>%
  summarise(
    Total_Count = sum(Count),
    .groups = "drop"
  )



# Calculate Crude Incidence Rate (CIR) yearly
yearly_cir <- yearly_summary %>%
  left_join(eritrea_population %>% summarise(Population = sum(Male + Female)), by = character()) %>%
  mutate(
    CIR = (Count / Population) * 100000
  )

# Calculate Age-Standardized Rate (ASR)
yearly_asr <- malignant_gastric %>%
  group_by(calender_year, age_group, gender) %>%
  summarise(
    Count = n(),
    .groups = "drop"
  ) %>%
  left_join(eritrea_population, by = "age_group") %>%
  left_join(who_standard, by = "age_group") %>%
  mutate(
    Weighted_Rate = (Count / (Male + Female)) * weight
  ) %>%
  group_by(calender_year, gender) %>%
  summarise(
    ASR = sum(Weighted_Rate, na.rm = TRUE) * 100000,
    .groups = "drop"
  )

# Combine CIR and ASR for final output
final_yearly_output <- yearly_cir %>%
  left_join(yearly_asr, by = c("calender_year", "gender")) %>%
  group_by(calender_year) %>%
  summarise(
    Male_CIR = CIR[gender == "Male"],
    Female_CIR = CIR[gender == "Female"],
    Total_CIR = sum(CIR),
    Male_ASR = ASR[gender == "Male"],
    Female_ASR = ASR[gender == "Female"],
    Total_ASR = sum(ASR),
    Male_to_Female_Ratio = ifelse(
      CIR[gender == "Female"] > 0,
      CIR[gender == "Male"] / CIR[gender == "Female"],
      NA
    )
  )

# Display final output
print(final_yearly_output)


# Summary of Crude incidence rate and ASR by gastric carcinoma type


gastric_carcinoma_type <- malignant_gastric %>% 
  mutate(adeno_other = case_when(
    diagnosis_category == "Adenocarcinoma" ~  "Adenocarcinoma",
    TRUE           ~ "Other Gastric Malignancy"
  )
  )


# Summarize data by classification, year, and sex
gastriccancer_summary_diagnosis <- gastric_carcinoma_type %>%
  group_by(calender_year, adeno_other, gender, age_group) %>%
  summarise(
    Count = n(), .groups = "drop"
  )

# Add population and standard weights

merged_data <- gastriccancer_summary_diagnosis %>%
  left_join(eritrea_population, by = "age_group") %>%
  left_join(who_standard, by = "age_group") %>%
  mutate(
    Population = ifelse(gender == "Male", Male, Female),
    Weighted_Rate = (Count / Population) * weight
  )


# Calculate ASR for each lymphoma type, gender, and year

asr_summary <- merged_data %>%
  group_by(calender_year, adeno_other, gender) %>%
  summarise(
    ASR = sum(Weighted_Rate, na.rm = TRUE) * 100000,
    .groups = "drop"
  )


# Calculate CIR
cir_summary <- gastriccancer_summary_diagnosis %>%
  left_join(
    eritrea_population %>% summarise(Total_Population = sum(Male + Female)),
    by = character()
  ) %>%
  mutate(
    CIR = (Count / Total_Population) * 100000
  ) %>%
  group_by(calender_year, adeno_other, gender) %>%
  summarise(
    CIR = sum(CIR, na.rm = TRUE),
    .groups = "drop"
  )

combined_summary <- asr_summary %>%
  full_join(cir_summary, by = c("calender_year", "adeno_other", "gender")) %>%
  pivot_wider(
    names_from = adeno_other,
    values_from = c(CIR, ASR),
    names_glue = "{.value}_{adeno_other}"
  )

colnames(combined_summary)

# Restructure for the final summary table
final_summary_diagnosis <- combined_summary %>%
  group_by(calender_year) %>%
  summarise(
    Female_ASR_Adenocarcinoma = sum(`ASR_Adenocarcinoma`[gender == "Female"], na.rm = TRUE),
    Male_ASR_Adenocarcinoma = sum(`ASR_Adenocarcinoma`[gender == "Male"], na.rm = TRUE),
    Total_ASR_Adenocarcinoma = sum(`ASR_Adenocarcinoma`, na.rm = TRUE),
    Female_ASR_OGM = sum(`ASR_Other Gastric Malignancy`[gender == "Female"], na.rm = TRUE),
    Male_ASR_OGM = sum(`ASR_Other Gastric Malignancy`[gender == "Male"], na.rm = TRUE),
    Total_ASR_OGM = sum(`ASR_Other Gastric Malignancy`, na.rm = TRUE),
    .groups = "drop"
  )

# Print final summary
print(final_summary_diagnosis)










