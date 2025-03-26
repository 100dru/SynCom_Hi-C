# Load necessary libraries

library(tidyverse)
library(dplyr)
library(ggplot2)
library(readxl)
library(readr)
library(ggpubr)
library(RColorBrewer)
library(viridisLite)
library(gplots)
library(patchwork)

# Read in the CSV file containing phage-host contact data
all_phage_host <- read_csv("combined_phage_host_normalized_contact.csv")
View(all_phage_host)

# Summation of chromosomes
# Chromosome summation for 13-15
## Step 1: Group by contig2_name and Mock, then filter rows
filtered_df <- all_phage_host %>%
  group_by(contig2_name, Mock) %>%
  filter(contig1_name %in% c("HOST_PSA-13-15_CHROM1", "HOST_PSA-13-15_CHROM2"))

## Step 2: Calculate the sum of "Score" for each group
sums <- filtered_df %>%
  summarise(total_score = sum(Score), .groups = "drop")

## Step 3: Create new rows with the calculated sums
new_rows <- data.frame(
  Mock = sums$Mock,           # Mock values from the summarized data
  contig1_name = "HOST_PSA-13-15",
  contig2_name = sums$contig2_name,     # contig2_name values from the summarized data
  Score = sums$total_score)

## Append the new rows to the original DataFrame
all_phage_host <- bind_rows(all_phage_host, new_rows)

# Chromosome summation for H71
## Step 1: Group by contig2_name and Mock, then filter rows
filtered_df <- all_phage_host %>%
  group_by(contig2_name, Mock) %>%
  filter(contig1_name %in% c("HOST_PSA-H71_CONTIG1", "HOST_PSA-H71_CONTIG3"))

## Step 2: Calculate the sum of "Score" for each group
sums <- filtered_df %>%
  summarise(total_score = sum(Score), .groups = "drop")

## Step 3: Create new rows with the calculated sums
new_rows <- data.frame(
  Mock = sums$Mock,           # Mock values from the summarized data
  contig1_name = "HOST_PSA-H71",
  contig2_name = sums$contig2_name,     # contig2_name values from the summarized data
  Score = sums$total_score)

## Append the new rows to the original DataFrame
all_phage_host <- bind_rows(all_phage_host, new_rows)

# Chromosome summation for H105
## Step 1: Group by contig2_name and Mock, then filter rows
filtered_df <- all_phage_host %>%
  group_by(contig2_name, Mock) %>%
  filter(contig1_name %in% c("HOST_PSA-H105_CONTIG1", "HOST_PSA-H105_CONTIG2"))

## Step 2: Calculate the sum of "Score" for each group
sums <- filtered_df %>%
  summarise(total_score = sum(Score), .groups = "drop")

## Step 3: Create new rows with the calculated sums
new_rows <- data.frame(
  Mock = sums$Mock,           # Mock values from the summarized data
  contig1_name = "HOST_PSA-H105",
  contig2_name = sums$contig2_name,     # contig2_name values from the summarized data
  Score = sums$total_score)

## Append the new rows to the original DataFrame
all_phage_host <- bind_rows(all_phage_host, new_rows)

# Print the updated DataFrame
print(all_phage_host)

# Remove the redundant split chromosome values and create a new table
values_to_remove <- c(
  "HOST_PSA-13-15_CHROM1",
  "HOST_PSA-13-15_CHROM2",
  "HOST_PSA-H71_CONTIG1",
  "HOST_PSA-H71_CONTIG3",
  "HOST_PSA-H105_CONTIG1",
  "HOST_PSA-H105_CONTIG2")

all_phage_host <- all_phage_host %>%
  filter(!(contig1_name %in% values_to_remove))

# Write the updated data to a new CSV file
write.table(all_phage_host, file = "all_phage_host_chromosome_sum.csv", sep = ',', row.names = FALSE)

# Subsetting all the mocks
mock_1B1 <- all_phage_host %>% filter(Mock == '1-B1')
mock_1B3 <- all_phage_host %>% filter(Mock == '1-B3')
mock_1D1 <- all_phage_host %>% filter(Mock == '1-D1')
mock_1D2 <- all_phage_host %>% filter(Mock == '1-D2')
mock_1D3 <- all_phage_host %>% filter(Mock == '1-D3')
mock_1N1 <- all_phage_host %>% filter(Mock == '1-N1')
mock_1N2 <- all_phage_host %>% filter(Mock == '1-N2')
mock_1N3 <- all_phage_host %>% filter(Mock == '1-N3')

# Calculate Z-scores for each mock subset
mock_1N1$z_score <- (mock_1N1$Score - mean(mock_1N1$Score, na.rm = TRUE)) / sd(mock_1N1$Score, na.rm = TRUE)
mock_1N2$z_score <- (mock_1N2$Score - mean(mock_1N2$Score, na.rm = TRUE)) / sd(mock_1N2$Score, na.rm = TRUE)
mock_1N3$z_score <- (mock_1N3$Score - mean(mock_1N3$Score, na.rm = TRUE)) / sd(mock_1N3$Score, na.rm = TRUE)
mock_1B1$z_score <- (mock_1B1$Score - mean(mock_1B1$Score, na.rm = TRUE)) / sd(mock_1B1$Score, na.rm = TRUE)
mock_1B3$z_score <- (mock_1B3$Score - mean(mock_1B3$Score, na.rm = TRUE)) / sd(mock_1B3$Score, na.rm = TRUE)
mock_1D1$z_score <- (mock_1D1$Score - mean(mock_1D1$Score, na.rm = TRUE)) / sd(mock_1D1$Score, na.rm = TRUE)
mock_1D2$z_score <- (mock_1D2$Score - mean(mock_1D2$Score, na.rm = TRUE)) / sd(mock_1D2$Score, na.rm = TRUE)
mock_1D3$z_score <- (mock_1D3$Score - mean(mock_1D3$Score, na.rm = TRUE)) / sd(mock_1D3$Score, na.rm = TRUE)

# Read in and clean data
scores = all_phage_host %>% 
  rename(sample = Mock, Host = contig1_name, Phage = contig2_name) %>% 
  mutate(experiment = sub("_S", "", sub("[0-9]$", "", sample)),
         Host = sub("HOST_", "", Host),
         Phage = sub("PHAGE_", "", Phage)) %>% 
  relocate(experiment, sample)

scores <- scores %>% 
  group_by(sample) %>% 
  mutate(z_score = scale(Score)) 

scores$z_score <- as.numeric(scores$z_score)

key = read_csv("key.csv") %>% 
  mutate(Host = ifelse(Host == "X", "none", Host)) %>% 
  mutate(Host = sub("HOST_", "", Host),
         Host = fct_relevel(Host, "none", after = Inf),
         Phage = sub("PHAGE_", "", Phage))

# Prepare data frames to join
## Make a matrix of unique samples and hosts, assuming each host was detected in each sample
samps_hosts0 = scores %>% 
  select(experiment, sample, Host) %>% 
  group_by_all() %>% 
  slice(1)

## No host is not an option for the original data, so create it
samps_hosts1 = scores %>% 
  select(experiment, sample) %>% 
  group_by_all() %>% 
  slice(1) %>% 
  mutate(Host = "none")

## One host was entirely missing from one sample
samps_hosts2 = scores %>% 
  select(experiment, sample) %>% 
  group_by_all() %>% 
  slice(1) %>% 
  filter(sample == "2_N3_S3") %>% 
  mutate(Host = "PSA-H71")

## Combine the data frames to make "no host" a possible host
samps_hosts = bind_rows(samps_hosts0, samps_hosts1, samps_hosts2) 

# Create a host-phage contact matrix from the key
key_w = key %>% 
  spread(Phage, Prediction, fill = "Incorrect")

# Join the sample-host matrix with the host-phage matrix
samps_key = full_join(samps_hosts, key_w) %>% 
  gather(Phage, Prediction, -experiment:-Host) %>% 
  mutate(drop = ifelse(Host == "none" & Prediction == "Incorrect", "drop", NA)) %>% 
  filter(is.na(drop)) %>% 
  select(-drop)

# Add the rest of the data
scores_key = full_join(samps_key, scores) %>% 
  mutate(zsc_pass = ifelse(z_score >= 0.5, "pass", "fail"),
         zsc_pass = ifelse(is.na(zsc_pass), "fail", zsc_pass)) %>% 
  mutate(type = ifelse(Host != "none" & Prediction == "Correct", "T-P",
                       ifelse(Host != "none" & Prediction == "Correct", "F-N",
                              ifelse(Host != "none" & Prediction == "Incorrect", "F-P",
                                     ifelse(Host != "none" & Prediction == "Incorrect", "T-N", NA)))) %>% 
  mutate(type = ifelse(Host == "none" & Prediction == "Correct", "T-P",
                       ifelse(Host == "none" & Prediction == "Correct", "T-N",
                              ifelse(Host == "none" & Prediction == "Incorrect", "F-P",
                                     ifelse(Host == "none" & Prediction == "Incorrect", "F-N", type)))) %>% 
  mutate(filt_type = ifelse(Host != "none" & Prediction == "Correct" & zsc_pass == "pass", "T-P",
                            ifelse(Host != "none" & Prediction == "Correct" & zsc_pass == "fail", "F-N",
                                   ifelse(Host != "none" & Prediction == "Incorrect" & zsc_pass == "pass", "F-P",
                                          ifelse(Host != "none" & Prediction == "Incorrect" & zsc_pass == "fail", "T-N", NA)))) %>% 
  mutate(filt_type = ifelse(Host == "none" & Prediction == "Correct" & zsc_pass == "pass", "T-P",
                            ifelse(Host == "none" & Prediction == "Correct" & zsc_pass == "fail", "T-N",
                                   ifelse(Host == "none" & Prediction == "Incorrect" & zsc_pass == "pass", "F-P",
                                          ifelse(Host == "none" & Prediction == "Incorrect" & zsc_pass == "fail", "F-N", filt_type))))

# Summarize to get counts at different levels of grouping
scores_sensspec_sample = scores_key %>% 
  group_by(experiment, sample, filt_type) %>% 
  summarise(count = n())

scores_sensspec_sample_experiment = scores_sensspec_sample %>% 
  group_by(experiment, filt_type) %>% 
  summarise(min = min(count),
            men = mean(count),
            med = median(count),
            max = max(count),
            sum = sum(count))

scores_sensspec_experiment = scores_key %>% 
  group_by(experiment, filt_type) %>% 
  summarise(count = n())

scores_sensspec_host = scores_key %>% 
  group_by(Host, filt_type) %>% 
  summarise(count = n())

scores_sensspec_hostexperiment = scores_key %>% 
  group_by(experiment, Host, filt_type) %>% 
  summarise(count = n())

scores_sensspec_phage = scores_key %>% 
  group_by(Phage, filt_type) %>% 
  summarise(count = n())

scores_sensspec_phageexperiment = scores_key %>% 
  group_by(experiment, Phage, filt_type) %>% 
  summarise(count = n())

# Sensitivity and Specificity calculations

# Function to calculate sensitivity (TP / (TP + FN))
syncom1_sensitivity <- scores_key %>% 
  group_by(sample) %>% 
  filter(Host != 'none') %>% 
  summarise(
    TP = sum(type == "T-P", na.rm = TRUE),  # Count occurrences of "T-P"
    FN = sum(type == "F-N", na.rm = TRUE),  # Count occurrences of "F-N"
    sensitivity = (TP / (TP + (7-TP))) * 100
  )

syncom1_sensitivity_filt <- scores_key %>% 
  group_by(sample) %>% 
  summarise(
    TP = sum(filt_type == "T-P", na.rm = TRUE),  # Count occurrences of "T-P"
    FN = sum(filt_type == "F-N", na.rm = TRUE),  # Count occurrences of "F-N"
    filt_sensitivity = (TP / (TP + FN)) * 100
  )

syncom1_specificity <- scores_key %>% 
  group_by(sample) %>% 
  filter(is.na(Score)) %>% 
  summarise(
    TP = sum(type == "T-P", na.rm = TRUE),
    TN = sum(type == "F-P" & is.na(Score), na.rm =TRUE),  # Count occurrences of "T-P"
    FP = sum(type == "F-P" & !is.na(Score), na.rm =TRUE),  # Count occurrences of "F-N"
    specificity = (TN / 29) * 100
  )

syncom1_specificity_filt <- scores_key %>% 
  group_by(sample) %>% 
  summarise(
    TP = sum(filt_type == "T-P", na.rm = TRUE),
    TN = sum(filt_type == "F-N", na.rm = TRUE),  # Count occurrences of "T-P"
    FP = sum(filt_type == "F-P", na.rm = TRUE),  # Count occurrences of "F-N"
    filt_specificity = (TN / (TN + FP)) * 100
  )

# Join sensitivity and specificity tables
syncom1_combined <- syncom1_sensitivity %>% 
  left_join(syncom1_sensitivity_filt, by = "sample") %>% 
  left_join(syncom1_specificity, by = "sample") %>% 
  left_join(syncom1_specificity_filt, by = "sample")

# Group by sample type and calculate the averages
sens_spec_result <- syncom1_combined %>% 
  mutate(sample_type = str_extract(sample, "^[0-9]-[A-Z]")) %>% 
  group_by(sample_type) %>% 
  summarise(
    Avg_Sensitivity = mean(sensitivity),
    Avg_Specificity = mean(specificity),
    Avg_Filt_Sensitivity = mean(filt_sensitivity),
    Avg_Filt_Specificity = mean(filt_specificity)
  )

sens_spec_long <- sens_spec_result %>% 
  pivot_longer(
    cols = starts_with("Avg"),
    names_to = "Metric",
    values_to = "Value"
  )

