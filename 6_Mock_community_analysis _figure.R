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

         # Filter data and create the main plot
# Filter out rows where Host is "none" or Score is NA
# Reverse the order of Prediction factor levels
# Create a ggplot object with Host on y-axis and Score on x-axis
# Add jittered points with specific fill colors based on Prediction
# Customize fill color manual scale
# Facet the plot by experiment and sample, with free scales
# Apply classic theme and customize axis text and legend position
main_plot <- scores_key %>%
  filter(Host != "none", !is.na(Score)) %>%
  mutate(Prediction = fct_rev(Prediction)) %>%
  ggplot(aes(y = Host, x = Score)) +
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("gray70", "gray30")) +
  facet_grid(experiment ~ sample, scales = "free") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

# Plot specific samples with different faceting and scales
# Filter out rows where Score is NA and specific samples are selected
# Create a ggplot object with Host on y-axis and Score on x-axis
# Add jittered points with specific fill colors based on Prediction
# Customize fill color manual scale
# Facet the plot by sample with 3 rows and free x scales
# Apply classic theme and customize axis text and legend position
p0A1 <- scores_key %>%
  filter(!is.na(Score), sample %in% c("1-N1", "1-N2", "1-N3")) %>%
  ggplot(aes(y = Host, x = Score)) +
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("black", "gray70")) +
  facet_wrap(~sample, nrow = 3, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 0, hjust = 1))

# Repeat the same process for other specific samples
p0A2 <- scores_key %>%
  filter(!is.na(Score), sample %in% c("1-D1", "1-D2", "1-D3")) %>%
  ggplot(aes(y = Host, x = Score)) +
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("black", "gray70")) +
  facet_wrap(~sample, nrow = 3, scales = "free_x", dir = "v") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 0, hjust = 1))

p0A3 <- scores_key %>%
  filter(!is.na(Score), sample %in% c("1-B1", "1-B3")) %>%
  ggplot(aes(y = Host, x = Score)) +
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("black", "gray70")) +
  facet_wrap(~sample, nrow = 3, scales = "free_x", dir = "v") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0A4 <- scores_key %>%
  filter(!is.na(Score), sample %in% c("2-N1", "2-N2", "2-N3")) %>%
  ggplot(aes(y = Host, x = Score)) +
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("black", "gray70")) +
  facet_wrap(~sample, nrow = 3, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0A5 <- scores_key %>%
  filter(!is.na(Score), sample %in% c("3-N1", "3-N3")) %>%
  ggplot(aes(y = Host, x = Score)) +
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("black", "gray70")) +
  facet_wrap(~sample, nrow = 3, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

# Plot with z_score instead of Score and with additional customization
p0B <- scores_key %>%
  filter(!is.na(Score)) %>%
  arrange(filt_type) %>%
  ggplot(aes(y = Host, x = z_score)) +
  geom_vline(xintercept = 0.5, color = "gray", linetype = 2, linewidth = 1) +
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("black", "gray70")) +
  facet_wrap(~sample, nrow = 3, scales = "fixed") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

# Plot specific samples with z_score and additional customization
p0B1 <- scores_key %>%
  filter(!is.na(Score), sample %in% c("1-N1", "1-N2", "1-N3")) %>%
  ggplot(aes(y = Host, x = z_score)) +
  geom_vline(xintercept = 0.5, color = "red", linetype = "dotted", linewidth = 1) +
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("black", "gray70")) +
  xlim(-1, 4) +
  facet_wrap(~sample, nrow = 3, scales = "free") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 0, hjust = 1))

# Repeat the same process for other specific samples with z_score
p0B2 <- scores_key %>%
  filter(Host != "none", !is.na(Score), sample %in% c("1-B1", "1-B3")) %>%
  ggplot(aes(y = Host, x = z_score)) +
  geom_vline(xintercept = 0.5, color = "red", linetype = 'dotted', linewidth = 1) +
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("black", "gray70")) +
  xlim(-1, 4) +
  facet_wrap(~sample, nrow = 3, scales = "free") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0B3 <- scores_key %>%
  filter(Host != "none", !is.na(Score), sample %in% c("1-D1", "1-D2", "1-D3")) %>%
  ggplot(aes(y = Host, x = z_score)) +
  geom_vline(xintercept = 0.5, color = "red", linetype = 'dotted', linewidth = 1) +
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("black", "gray70")) +
  xlim(-1, 4) +
  facet_wrap(~sample, nrow = 3, scales = "free") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 0, hjust = 1))

p0B4 <- scores_key %>%
  filter(Host != 'none', !is.na(Score), sample %in% c("2-N1", "2-N2", "2-N3")) %>%
  ggplot(aes(y = Host, x = z_score)) +
  geom_vline(xintercept = 0.5, color = "red", linetype = 'dotted', linewidth = 1) +
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("black", "gray70")) +
  xlim(-1, 4) +
  facet_wrap(~sample, nrow = 3, scales = "free") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 0, hjust = 1))

p0B5 <- scores_key %>%
  filter(Host != "none", !is.na(Score), sample %in% c("3-N1", "3-N3")) %>%
  ggplot(aes(y = Host, x = z_score)) +
  geom_vline(xintercept = 0.5, color = "red", linetype = 'dotted', linewidth = 1) +
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("black", "gray70")) +
  xlim(-1, 4) +
  facet_wrap(~sample, nrow = 3, scales = "free") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 0, hjust = 1))

# Create heatmap plots for specific samples
p0C1 <- scores_key %>%
  filter(Host != "none", !is.na(Score), sample %in% c("1-N1", "1-N2", "1-N3")) %>%
  ggplot(aes(y = Host, x = Phage)) +
  geom_tile(aes(fill = type), color = "white") +
  scale_fill_manual(values = c("T-P" = "black", "T-N" = "white", "F-P" = "gray60", "F-N" = "gray60"),
                    labels = c("T-P" = "True Positive", "T-N" = "True Negative", "F-P" = "False Positive", "F-N" = "False Negative"),
                    na.value = NA) +
  facet_wrap(~sample, nrow = 3, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

# Repeat the same process for other specific samples
p0C2 <- scores_key %>%
  filter(Host != "none", !is.na(Score), sample %in% c("1-D1", "1-D2", "1-D3")) %>%
  ggplot(aes(y = Host, x = Phage)) +
  geom_tile(aes(fill = type), color = "white") +
  scale_fill_manual(values = c("T-P" = "black", "T-N" = "gray90", "F-P" = "gray60", "F-N" = "gray60"),
                    labels = c("T-P" = "True Positive", "T-N" = "True Negative", "F-P" = "False Positive", "F-N" = "False Negative"),
                    na.value = NA) +
  facet_wrap(~sample, nrow = 3, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0C3 <- scores_key %>%
  filter(Host != "none", !is.na(Score), sample %in% c("1-B1", "1-B3")) %>%
  ggplot(aes(y = Host, x = Phage)) +
  geom_tile(aes(fill = type), color = "white") +
  scale_fill_manual(values = c("T-P" = "black", "T-N" = "gray90", "F-P" = "gray60", "F-N" = "gray60"),
                    labels = c("T-P" = "True Positive", "T-N" = "True Negative", "F-P" = "False Positive", "F-N" = "False Negative"),
                    na.value = NA) +
  facet_wrap(~sample, nrow = 3, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0C4 <- scores_key %>%
  filter(Host != "none", !is.na(Score), sample %in% c("2-N1", "2-N2", "2-N3")) %>%
  ggplot(aes(y = Host, x = Phage)) +
  geom_tile(aes(fill = type), color = "white") +
  scale_fill_manual(values = c("T-P" = "black", "T-N" = "gray90", "F-P" = "gray60", "F-N" = "gray60"),
                    labels = c("T-P" = "True Positive", "T-N" = "True Negative", "F-P" = "False Positive", "F-N" = "False Negative"),
                    na.value = NA) +
  facet_wrap(~sample, nrow = 3, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0C5 <- scores_key %>%
  filter(Host != "none", !is.na(Score), sample %in% c("3-N1", "3-N3")) %>%
  ggplot(aes(y = Host, x = Phage)) +
  geom_tile(aes(fill = type), color = "white") +
  scale_fill_manual(values = c("T-P" = "black", "T-N" = "gray90", "F-P" = "gray60", "F-N" = "gray60"),
                    labels = c("T-P" = "True Positive", "T-N" = "True Negative", "F-P" = "False Positive", "F-N" = "False Negative"),
                    na.value = NA) +
  facet_wrap(~sample, nrow = 3, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

# Create heatmap plots with conditional fill
p0D <- scores_key %>%
  filter(Host != "none", !is.na(Score)) %>%
  ggplot(aes(y = Host, x = Phage)) +
  geom_tile(aes(fill = ifelse(z_score > 0.5, type, NA)), color = "white") +
  scale_fill_brewer(palette = "Greys", direction = 1, na.value = NA) +
  facet_wrap(~sample, nrow = 3, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

         p0D1 = scores_key %>% 
  filter(Host != "none",sample %in% c("1-N1", "1-N2", "1-N3")) %>% 
  ggplot(aes(y = Host, x = Phage)) + 
  geom_tile(aes(fill = ifelse(z_score>0.5, type, NA)), color = "white") +
  scale_fill_manual(values = c("T-P" = "black", "T-N" = "gray90", "F-P" = "gray60", "F-N" = "gray60"), 
                    labels = c("T-P" = "True Positive", "T-N" = "True Negative", "F-P" = "False Positive", "F-N" = "False Negative"), 
                    na.value = NA) + 
  facet_wrap(sample~., nrow=3, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0D1

p0D2 = scores_key %>% 
  filter(Host != "none",sample %in% c( "1-D1", "1-D2", "1-D3")) %>% 
  ggplot(aes(y = Host, x = Phage)) + 
  geom_tile(aes(fill = ifelse(z_score>0.5, type, NA)), color = "white") +
  scale_fill_manual(values = c("T-P" = "black", "T-N" = "gray90", "F-P" = "gray60", "F-N" = "gray60"), 
                    labels = c("T-P" = "True Positive", "T-N" = "True Negative", "F-P" = "False Positive", "F-N" = "False Negative"), 
                    na.value = NA) + 
  facet_wrap(sample~., nrow=3, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0D2

p0D3 = scores_key %>% 
  filter(Host != "none",sample %in% c("1-B1", "1-B3")) %>% 
  ggplot(aes(y = Host, x = Phage)) + 
  geom_tile(aes(fill = ifelse(z_score>0.5, type, NA)), color = "white") +
  scale_fill_manual(values = c("T-P" = "black", "T-N" = "gray90", "F-P" = "gray60", "F-N" = "gray60"), 
                    labels = c("T-P" = "True Positive", "T-N" = "True Negative", "F-P" = "False Positive", "F-N" = "False Negative"), 
                    na.value = NA) + 
  facet_wrap(sample~., nrow=3, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0D3

p0D4 = scores_key %>% 
  filter(Host != "none",sample %in% c("2-N1", "2-N2", "2-N3")) %>% 
  ggplot(aes(y = Host, x = Phage)) + 
  geom_tile(aes(fill = ifelse(z_score>0.5, 
                              type, NA)), color = "white") +
  scale_fill_manual(values = c("T-P" = "black", "T-N" = "gray90", "F-P" = "gray60", "F-N" = "gray60"), 
                    labels = c("T-P" = "True Positive", "T-N" = "True Negative", "F-P" = "False Positive", "F-N" = "False Negative"), 
                    na.value = NA) + 
  facet_wrap(sample~., nrow=3, scales = "free_x")+
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0D4

p0D5 = scores_key %>% 
  filter(Host != "none",sample %in% c("3-N1", "3-N3")) %>% 
  ggplot(aes(y = Host, x = Phage)) + 
  geom_tile(aes(fill = ifelse(z_score>0.5, type, NA)), color = "white") +
  scale_fill_manual(values = c("T-P" = "black", "T-N" = "gray90", "F-P" = "gray60", "F-N" = "gray60"), 
                    labels = c("T-P" = "True Positive", "T-N" = "True Negative", "F-P" = "False Positive", "F-N" = "False Negative"), 
                    na.value = NA) + 
  facet_wrap(sample~., nrow=3, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0D5

         # Create bar plots for sensitivity and specificity metrics
# Filter data for specific sample types and metrics
# Create a ggplot object with sample_type on x-axis and Value on y-axis
# Add bar geom with identity stat and dodge position
# Customize fill color manual scale and labels
# Apply classic theme and customize axis text
p0E <- sens_spec_long %>%
  filter(sample_type %in% c("1-D", "1-B", "1-N"), Metric %in% c("Avg_Filt_Sensitivity", "Avg_Filt_Specificity")) %>%
  ggplot(aes(x = sample_type, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Avg_Filt_Specificity" = "#32373AFF", "Avg_Filt_Sensitivity" = "#9F9994FF"),
                    labels = c("Avg_Filt_Specificity" = "Specificity", "Avg_Filt_Sensitivity" = "Sensitivity")) +
  theme_classic() +
  labs(x = "Synthetic Communities", y = "Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Repeat the same process for other sample types
p0F <- sens_spec_long %>%
  filter(sample_type %in% c("1-N", "2-N", "3-N"), Metric %in% c("Avg_Filt_Sensitivity", "Avg_Filt_Specificity")) %>%
  ggplot(aes(x = sample_type, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Avg_Filt_Specificity" = "#32373AFF", "Avg_Filt_Sensitivity" = "#9F9994FF"),
                    labels = c("Avg_Filt_Specificity" = "Specificity", "Avg_Filt_Sensitivity" = "Sensitivity")) +
  theme_classic() +
  labs(x = "Synthetic Communities", y = "Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine plots into figures
# Arrange plots p0A1, p0B1, p0C1, and p0D1 in a single row for Figure 1
Fig_1 <- ggarrange(p0A1, p0B1, p0C1, p0D1, ncol = 4, align = "h")

# Arrange plots p0B3 and p0D2 in a single row for Figure 2A
Fig_2_A <- ggarrange(p0B3, p0D2, ncol = 2, align = "h")

# Arrange plots p0B2 and p0D3 in a single row for Figure 2B
Fig_2_B <- ggarrange(p0B2, p0D3, ncol = 2, align = 'h')

# Arrange plots p0B4 and p0D4 in a single row for Figure 3A
Fig_3_A <- ggarrange(p0B4, p0D4, ncol = 2, align = 'h')

# Arrange plots p0B5 and p0D5 in a single row for Figure 3B
Fig_3_B <- ggarrange(p0B5, p0D5, ncol = 2, align = 'h')

# Create heatmap plot for specific samples with customized facet labels
# Define custom labels for samples
syncom1_names <- c(
  `1-N1` = "SynCom 1 replicate 1",
  `1-N2` = "SynCom 1 replicate 2",
  `1-N3` = "SynCom 1 replicate 3"
)

# Filter data for specific samples and create a heatmap plot
# Customize fill color manual scale and labels
# Facet the plot by sample with custom labels
# Apply classic theme and customize axis text and legend position
p0E1 <- scores_key %>%
  filter(Host != "none", !is.na(Score), sample %in% c("1-N1", "1-N2", "1-N3")) %>%
  ggplot(aes(y = Host, x = Phage)) +
  geom_tile(aes(fill = ifelse(z_score > 0.5, filt_type, NA)), color = "white") +
  scale_fill_manual(values = c("T-P" = "black", "F-P" = "gray30", "F-N" = "gray60", "T-N" = "gray90")) +
  facet_grid(~ sample, labeller = as_labeller(syncom1_names)) +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

# Create scatter plot for specific samples with customized facet labels
# Filter data for specific samples and create a scatter plot
# Customize color manual scale and labels
# Apply classic theme and customize axis text
SOM1 <- scores_key %>%
  filter(sample %in% c("1-N1", "1-N2", "1-N3")) %>%
  ggplot(aes(x = Phage, y = Score, color = Prediction)) +
  geom_point() +
  scale_colour_manual(values = setNames(c('black', 'gray'), c("Correct", "Incorrect"))) +
  theme_classic() +
  facet_wrap(~sample, labeller = as_labeller(syncom1_names)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# Summarize data by z-score pass and fail counts
# Drop NA values and create binary columns for z-score pass and fail
# Summarize counts for z-score pass and fail
count <- scores_key %>%
  drop_na() %>%
  mutate(
    z_pass = zsc_pass == "pass",
    z_fail = zsc_pass == "fail"
  ) %>%
  summarize(
    z_pass_count = sum(z_pass, na.rm = TRUE),
    z_fail_count = sum(z_fail, na.rm = TRUE)
  )

# Pivot summarized counts to long format and create bar plot
# Filter data for specific samples and create a bar plot
# Customize fill color manual scale and labels
# Apply classic theme and customize axis text
count_long <- count %>%
  pivot_longer(
    cols = c(z_pass_count, z_fail_count), names_to = 'zpass', values_to = 'count'
  )

som2 <- count_long %>%
  filter(sample %in% c("1-N1", "1-N2", "1-N3")) %>%
  ggplot(aes(x = count, y = sample, fill = zpass)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c("z_pass_count" = "black", "z_fail_count" = "grey70"),
                    labels = c("z_pass_count" = "Z-Score >= 0.54", "z_fail_count" = "Z-Score <= 0.5")) +
  labs(x = "Percentage of Linkages", y = "Synthetic Community") +
  theme_classic()

# Create ROC curve plots for z_score and Score thresholds
# Convert Prediction to binary values: 1 for "Correct" and 0 for "Incorrect"
scores_key$Prediction_binary <- ifelse(scores_key$Prediction == "Correct", 1, 0)

# Filter data for specific samples and calculate ROC curve for z_score
syncom1 <- scores_key %>%
  filter(sample %in% c("1-N1", "1-N2", "1-N3"))
roc_obj <- roc(syncom1$Prediction_binary, syncom1$z_score, direction = "<")

# Plot the ROC curve for z_score
plot(roc_obj, col = "blue", main = "ROC Curve for z_score Thresholds")

# Calculate ROC curve for Score and plot it
roc_score <- roc(syncom1$Prediction_binary, syncom1$Score, direction = "<")
plot(roc_score, col = "red", main = "ROC Curve for Score Thresholds")

# Add the z_score ROC curve to the same plot for comparison
plot(roc_obj, col = "blue", add = TRUE)

# Add a legend to the plot
legend("bottomright", legend = c("z_score", "Score"), col = c("blue", "red"), lwd = 2)

# Calculate and print AUC for z_score
auc_z_score <- auc(roc_obj)
print(paste("AUC for z_score:", auc_z_score))

# Calculate and print AUC for Score
auc_score <- auc(roc_score)
print(paste("AUC for Score:", auc_score))

# Save the ROC curves to a PNG file
png("roc_curve_mock1.png", width = 400, height = 400)

# Plot both ROC curves on the same plot and add a legend with AUC values
plot(roc_score, col = "red", main = "ROC Curve for SynCom 1 Replicates")
plot(roc_obj, col = "blue", add = TRUE)
legend("bottomright",
       legend = c(paste("z_score (AUC =", round(auc_z_score, 2), ")"),
                  paste("Score (AUC =", round(auc_score, 2), ")")),
       col = c("blue", "red"),
       lwd = 2)
dev.off()

# Calculate and print the threshold for z_score with the highest specificity
highest_specificity_threshold <- coords(roc_obj, "local maximas", ret = "threshold", input = "specificity", transpose = TRUE)
print(paste("Threshold for z_score with highest specificity:", highest_specificity_threshold))
