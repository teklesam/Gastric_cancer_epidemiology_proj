###############################
# Epidemiology of gastic cancer
##############################

# 
# Written by : Samuel Tekle,MD
# Copyright (c) - 2024 Samuel Tekle
# 
# Started on : 5.11.2024
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

# Load necessary libraries
library(pacman)             # "pacman" package manager - Helps to load multiple packages at a time

pacman::p_load(
  here,
  rio,
  dplyr,
  ggplot2,
  knitr,
  gtsummary,
)

# load data

gastric_cancer_data <- import(here("3_output","cleaned_gastric_data.rds"))


# Create the summary table by diagnosis category
# Create the summary table
# Create the summary table


summary_table_gender <- gastric_cancer_data %>%
  select(
    age, 
    gender,
    calender_year, 
    diagnosis_category,
    address_zoba,
    malignant_not,
  ) %>%
  tbl_summary(
    by = gender, # Stratify by Gender
    missing = "no", # Exclude missing values
    statistic = list(
      age ~ "{median} ({p25} - {p75})", # Median (25th - 75th percentile) for Age
      calender_year ~ "{median} ({p25} - {p75})", # Median (25th - 75th percentile) for Calendar Year
      all_categorical() ~ "{n} ({p}%)" # Counts (%) for categorical variables
    ),
    label = list(
      age ~ "Age (years)",
      calender_year ~ "Calendar Year",
      diagnosis_category ~ "Diagnosis Category",
      address_zoba  ~ "Zonal Address",
      gender ~ "Gender",
      malignant_not ~ "Malignant or Non-Malignant"
    )
  ) %>%
  add_n() %>%   # Add the total number of observations in each group
  add_p()

summary_table_malignant <- gastric_cancer_data %>%
  filter(malignant_not != "Not diagnostic") %>% # Exclude "Not diagnostic"
  mutate(malignant_not = droplevels(malignant_not)) %>% # Drop unused factor levels
  select(
    age, 
    calender_year, 
    address_zoba,
    gender,
    malignant_not
  ) %>%
  tbl_summary(
    by = malignant_not, # Stratify by Malignant or Not
    missing = "no", # Exclude missing values
    statistic = list(
      age ~ "{median} ({p25} - {p75})", # Median (25th - 75th percentile) for Age
      calender_year ~ "{median} ({p25} - {p75})", # Median (25th - 75th percentile) for Calendar Year
      all_categorical() ~ "{n} ({p}%)" # Counts (%) for categorical variables
    ),
    label = list(
      age ~ "Age (years)",
      calender_year ~ "Calendar Year",
      address_zoba ~ "Zonal Address",
      gender ~ "Gender",
      malignant_not ~ "Malignant (Yes/No)" # Correct label for the actual column name
    )
  ) %>%
  add_p() %>% # Add p-values
  add_n() # Add total number of observations in each group


summary_table_pathology <- gastric_cancer_data %>%
  mutate(
    diagnosis_category = factor(
      diagnosis_category,
      levels = c(
        "Adenocarcinoma",
        "Lymphoma",
        "Other Malignancies",
        "Squamous Cell Carcinoma",
        "Non-Neoplastic Diseases",
        "Benign Neoplasms",
        "Other Conditions",
        "Indeterminate or Poor Sample Quality"
      )
    )
  ) %>%
  select(
    age, 
    gender,
    calender_year, 
    diagnosis_category,
    address_zoba
  ) %>%
  tbl_summary(
    by = diagnosis_category, # Stratify by Diagnosis Category
    missing = "no", # Exclude missing values
    statistic = list(
      age ~ "{median} ({p25} - {p75})", # Show median (25th - 75th percentile) for Age
      calender_year ~ "{median} ({p25} - {p75})", # Show median (25th - 75th percentile) for Calendar Year
      all_categorical() ~ "{n} ({p}%)" # Show counts (%) for categorical variables
    ),
    label = list(
      age ~ "Age (years)",
      calender_year ~ "Calendar Year",
      diagnosis_category ~ "Diagnosis Category",
      address_zoba ~ "Zonal Address"
    )
  ) %>%
  add_n()   # Add the total number of observations in each group
# Print the summary table
summary_table_pathology

# Print the summary table
summary_table

# Print the summary table
summary_table



# Visualizations
# 1. Numeric variables: Age
library(ggplot2)

# Create the histogram with a line and borders around the bins
ggplot(gastric_cancer_data, aes(x = age)) +
  geom_histogram(
    binwidth = 5, 
    fill = "skyblue", 
    alpha = 0.8, 
    color = "black",   # Add border around the bins
    size = 0.6         # Border thickness
  ) +
  geom_smooth(
    stat = "bin",
    binwidth = 5,
    aes(y = ..count..), # Align with histogram counts
    color = "red",
    size = 1.2,         # Smooth line thickness
    se = FALSE          # Remove confidence interval shading
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Age Distribution with Peak Trend",
    x = "Age (years)",
    y = "Count"
  ) +
  theme(
    panel.grid = element_blank(),       # Remove grid lines
    panel.border = element_rect(
      color = "black", fill = NA, size = 1.2
    ),                                  # Add border around the plot
    plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold title
    axis.title = element_text(face = "bold")              # Bold axis labels
  )

ggplot(gastric_cancer_data, aes(x = age)) +
  geom_histogram(aes(y = ..density..), binwidth = 5, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_density(color = "red", size = 1.2, adjust = 1) +  # Density line to highlight peaks
  labs(title = "Age Distribution of Lymphoma Cases with Density Curve",
       x = "Age",
       y = "Density") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, max(gastric_cancer_data$age, na.rm = TRUE), by = 5)) +
  scale_y_continuous(labels = scales::percent)  # Optional: shows density



# Create a clearer visualization with separate histograms and transparent density curves
malignant_not_visul <- gastric_cancer_data %>%
  filter(malignant_not != "Not diagnostic") %>% # Exclude "Not diagnostic"
  mutate(malignant_not = droplevels(malignant_not)) 

ggplot(malignant_not_visul, aes(x = age, fill = malignant_not)) +
  # Histograms with transparency and borders, side-by-side
  geom_histogram(
    aes(y = ..density..),  # Use density instead of count for clearer comparison
    binwidth = 5, 
    position = "dodge",    # Place histograms side by side
    alpha = 0.6,           # Slight transparency for histograms
    color = "black",       # Borders around the bins
    size = 0.5
  ) +
  # Transparent density curves
  geom_density(
    aes(y = ..density.., color = malignant_not), # Separate density lines
    size = 1.2,
    adjust = 1.5,            # Smoothing parameter for density line
    alpha = 0.2              # Transparency for density curves
  ) +
  # Custom color palette for fill and density line
  scale_fill_manual(values = c("blue", "red","yellow")) +  
  scale_color_manual(values = c("darkblue", "darkred", "yellow")) +
  theme_minimal() +
  labs(
    title = "Age Distribution by Malignant Status",
    x = "Age (years)",
    y = "Density"
  ) +
  theme(
    panel.grid.major = element_blank(),  # Remove grid lines for clarity
    panel.grid.minor = element_blank(),
    legend.title = element_blank()  # Remove legend title for cleaner look
  )



# Create the histogram with transparent density curves for diagnosis categories
gastric_cancer_cat <- gastric_cancer_data %>%
  mutate(
    diagnosis_category = factor(
      diagnosis_category,
      levels = c(
        "Adenocarcinoma",
        "Lymphoma",
        "Other Malignancies",
        "Squamous Cell Carcinoma",
        "Non-Neoplastic Diseases",
        "Benign Neoplasms",
        "Other Conditions",
        "Indeterminate or Poor Sample Quality"
      )
    )
  ) 


# Filter out levels with fewer than two data points
gastric_cancer_cat_filtered <- gastric_cancer_cat %>%
  group_by(diagnosis_category) %>%
  filter(n() >= 2) %>%
  ungroup()

# Filter out "Indeterminate or Poor Sample Quality" and NA values
gastric_cancer_cat_filtered <- gastric_cancer_cat_filtered %>%
  filter(
    !is.na(diagnosis_category),                # Exclude NA values
    diagnosis_category != "Indeterminate or Poor Sample Quality" # Exclude specific category
  )

# Plot
ggplot(gastric_cancer_cat_filtered, aes(x = age, fill = diagnosis_category)) +
  geom_histogram(
    aes(y = ..density..),  # Use density instead of count for better comparison
    binwidth = 5, 
    position = "dodge",    # Place histograms side by side
    alpha = 0.6,           # Slight transparency for histograms
    color = "black",       # Borders around the bins
    size = 0.5
  ) +
  geom_density(
    aes(y = ..density.., color = diagnosis_category), # Separate density lines by diagnosis category
    size = 1.2,
    adjust = 1.5,            # Smoothing parameter for density line
    alpha = 0.2              # Transparency for density curves
  ) +
  scale_fill_manual(values = c("blue", "red", "green", "purple", "orange", "yellow", "pink", "#00CED1")) +  
  scale_color_manual(values = c("darkblue", "darkred", "darkgreen", "purple4", "darkorange", "gold", "deeppink", "#008B8B")) +
  theme_minimal() +
  labs(
    title = "Age Distribution by Diagnosis Category",
    x = "Age (years)",
    y = "Density"
  ) +
  theme(
    panel.grid.major = element_blank(),  # Remove grid lines for clarity
    panel.grid.minor = element_blank(),
    legend.title = element_blank()  # Remove legend title for cleaner look
  )

# 2. Categorical variables: Gender
ggplot(gastric_cancer_data, aes(x = gender, fill = gender)) +
  geom_bar() +
  theme_minimal() +
  labs(title = "Gender Distribution", x = "Gender", y = "Count")


# Create the histogram with transparent density curves for gender
ggplot(gastric_cancer_data, aes(x = age, fill = gender)) +
  # Histograms with transparency and borders, side-by-side
  geom_histogram(
    aes(y = ..density..),  # Use density instead of count for better comparison
    binwidth = 5, 
    position = "dodge",    # Place histograms side by side
    alpha = 0.6,           # Slight transparency for histograms
    color = "black",       # Borders around the bins
    size = 0.5
  ) +
  # Transparent density curves
  geom_density(
    aes(y = ..density.., color = gender), # Separate density lines by gender
    size = 1.2,
    adjust = 1.5,            # Smoothing parameter for density line
    alpha = 0.2              # Transparency for density curves
  ) +
  # Custom color palette for fill and density line
  scale_fill_manual(values = c("blue", "pink")) +  
  scale_color_manual(values = c("darkblue", "pink")) +
  theme_minimal() +
  labs(
    title = "Age Distribution by Gender",
    x = "Age (years)",
    y = "Density"
  ) +
  theme(
    panel.grid.major = element_blank(),  # Remove grid lines for clarity
    panel.grid.minor = element_blank(),
    legend.title = element_blank()  # Remove legend title for cleaner look
  )

# 3. Population structure

ggplot(gastric_cancer_cat_filtered, aes(x = age_group)) +
  geom_bar(fill = "green", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Population Structure", x = "Age Group", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


library(ggplot2)
library(RColorBrewer)

ggplot(gastric_cancer_cat_filtered, aes(x = age_group, fill = diagnosis_category)) +
  # Bar plot with transparency and black borders
  geom_bar(alpha = 0.7, position = "identity", color = "black") +
  
  # Density curves aligned with counts, with distinct colors
  geom_density(
    aes(y = ..count.., group = diagnosis_category, color = diagnosis_category),
    alpha = 0.4,
    position = "identity"
  ) +
  
  # Minimal theme for clarity
  theme_minimal() +
  
  # Title and axis labels
  labs(
    title = "Population Structure",
    x = "Age Group",
    y = "Count",
    fill = "Diagnosis Category",
    color = "Diagnosis Category"
  ) +
  
  # Rotate x-axis labels for better readability
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  
  # Color palette for fill and density lines with more distinctive colors
  scale_fill_brewer(palette = "Set3") + # Set3 is often more diverse and vibrant
  scale_color_brewer(palette = "Set3") # Matching Set3 palette for density lines






