setwd("C:/Users/shatadru.2/OneDrive - The Ohio State University/LAB/Hi-C/paper_unit/final_analysis/mock_comm")

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

all_phage_host <- read_csv("C:/Users/shatadru.2/OneDrive - The Ohio State University/LAB/Hi-C/paper_unit/final_analysis/mock_comm/combined_phage_host_normalized_contact.csv")
View(all_phage_host)

#Summation of chromosomes
## chromosome summation 13-15

### Step 1: Group by contig2_name and Mock, then filter rows
filtered_df <- all_phage_host %>%
  group_by(contig2_name, Mock) %>%
  filter(contig1_name %in% c("HOST_PSA-13-15_CHROM1", "HOST_PSA-13-15_CHROM2"))

### Step 2: Calculate the average of "Score" for each group
sums <- filtered_df %>%
  summarise(total_score = sum(Score), .groups = "drop")

### Step 3: Create new rows with the calculated sums
new_rows <- data.frame(
  Mock = sums$Mock,           # Mock values from the summarized data
  contig1_name = "HOST_PSA-13-15",
  contig2_name = sums$contig2_name,     # contig2_name values from the summarized data
  Score = sums$total_score)

### Append the new rows to the original DataFrame
all_phage_host <- bind_rows(all_phage_host, new_rows)

## chromosome summation H71

### Step 1: Group by contig2_name and Mock, then filter rows
filtered_df <- all_phage_host %>%
  group_by(contig2_name, Mock) %>%
  filter(contig1_name %in% c("HOST_PSA-H71_CONTIG1", "HOST_PSA-H71_CONTIG3"))

### Step 2: Calculate the average of "Score" for each group
sums <- filtered_df %>%
  summarise(total_score = sum(Score), .groups = "drop")

### Step 3: Create new rows with the calculated sums
new_rows <- data.frame(
  Mock = sums$Mock,           # Mock values from the summarized data
  contig1_name = "HOST_PSA-H71",
  contig2_name = sums$contig2_name,     # contig2_name values from the summarized data
  Score = sums$total_score)

### Append the new rows to the original DataFrame
all_phage_host <- bind_rows(all_phage_host, new_rows)

## chromosome summation H105

### Step 1: Group by contig2_name and Mock, then filter rows
filtered_df <- all_phage_host %>%
  group_by(contig2_name, Mock) %>%
  filter(contig1_name %in% c("HOST_PSA-H105_CONTIG1", "HOST_PSA-H105_CONTIG2"))

### Step 2: Calculate the average of "Score" for each group
sums <- filtered_df %>%
  summarise(total_score = sum(Score), .groups = "drop")

### Step 3: Create new rows with the calculated sums
new_rows <- data.frame(
  Mock = sums$Mock,           # Mock values from the summarized data
  contig1_name = "HOST_PSA-H105",
  contig2_name = sums$contig2_name,     # contig2_name values from the summarized data
  Score = sums$total_score)

### Append the new rows to the original DataFrame
all_phage_host <- bind_rows(all_phage_host, new_rows)

print(all_phage_host)

## remove the redundant split chromosome values and make new table
values_to_remove <- c(
  "HOST_PSA-13-15_CHROM1",
  "HOST_PSA-13-15_CHROM2",
  "HOST_PSA-H71_CONTIG1",
  "HOST_PSA-H71_CONTIG3",
  "HOST_PSA-H105_CONTIG1",
  "HOST_PSA-H105_CONTIG2")

all_phage_host <- all_phage_host %>%
  filter(!(contig1_name %in% values_to_remove))

write.table(all_phage_host, file = "all_phage_host_chromosome_sum.csv", sep = ',', row.names = FALSE)


### subsetting all the mocks

mock_1B1 <- all_phage_host %>% filter(Mock == '1-B1')
mock_1B3 <- all_phage_host %>% filter(Mock == '1-B3')
mock_1D1 <- all_phage_host %>% filter(Mock == '1-D1')
mock_1D2 <- all_phage_host %>% filter(Mock == '1-D2')
mock_1D3 <- all_phage_host %>% filter(Mock == '1-D3')
mock_1N1 <- all_phage_host %>% filter(Mock == '1-N1')
mock_1N2 <- all_phage_host %>% filter(Mock == '1-N2')
mock_1N3 <- all_phage_host %>% filter(Mock == '1-N3')


# Statistics Test

mock_1N1$z_score <- (mock_1N1$Score - mean(mock_1N1$Score, na.rm = TRUE)) / sd(mock_1N1$Score, na.rm = TRUE)
mock_1N2$z_score <- (mock_1N2$Score - mean(mock_1N2$Score, na.rm = TRUE)) / sd(mock_1N2$Score, na.rm = TRUE)
mock_1N3$z_score <- (mock_1N3$Score - mean(mock_1N3$Score, na.rm = TRUE)) / sd(mock_1N3$Score, na.rm = TRUE)
mock_1B1$z_score <- (mock_1B1$Score - mean(mock_1B1$Score, na.rm = TRUE)) / sd(mock_1B1$Score, na.rm = TRUE)
mock_1B3$z_score <- (mock_1B3$Score - mean(mock_1B3$Score, na.rm = TRUE)) / sd(mock_1B3$Score, na.rm = TRUE)
mock_1D1$z_score <- (mock_1D1$Score - mean(mock_1D1$Score, na.rm = TRUE)) / sd(mock_1D1$Score, na.rm = TRUE)
mock_1D2$z_score <- (mock_1D2$Score - mean(mock_1D2$Score, na.rm = TRUE)) / sd(mock_1D2$Score, na.rm = TRUE)
mock_1D3$z_score <- (mock_1D3$Score - mean(mock_1D3$Score, na.rm = TRUE)) / sd(mock_1D3$Score, na.rm = TRUE)

#read in and clean data
scores = all_phage_host %>% 
  rename(sample = Mock, Host = contig1_name, Phage = contig2_name) %>% 
  mutate(experiment = sub("_S", "", sub("[0-9]$", "", sample)),
         #replicate = sub(".*[A-Z]", "", sub("_S[0-9]", "", sample)),
         Host = sub("HOST_", "", Host),
         #Host = fct_relevel(Host, "none", after = Inf),
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

#

#prep dfs to join
#first make a matrix of unique samples and hosts, assuming that each host was detected in each sample
samps_hosts0 = scores %>% 
  select(experiment, sample, Host) %>% 
  group_by_all() %>% 
  slice(1)

#no host is not an option for the original data, so create it
samps_hosts1 = scores %>% 
  select(experiment,  sample) %>% 
  group_by_all() %>% 
  slice(1) %>% 
  mutate(Host = "none")

#One host was entirely missing from one sample
samps_hosts2 = scores %>% 
  select(experiment,  sample) %>% 
  group_by_all() %>% 
  slice(1) %>% 
  filter(sample == "2_N3_S3") %>% 
  mutate(Host = "PSA-H71")

#put the dfs together to make no host a possible host
samps_hosts = bind_rows(samps_hosts0, samps_hosts1, samps_hosts2) 

#make a host-phage contact matrix from the key
key_w = key %>% 
  #filter(Host != "X") %>% 
  spread(Phage, Prediction, fill= "Incorrect")

#join the sample-host matrix with the host-phage matrix
samps_key = full_join(samps_hosts, key_w) %>% 
  gather(Phage, Prediction, -experiment:-Host) %>% 
  mutate(drop = ifelse(Host == "none" & Prediction == "Incorrect", "drop", NA)) %>% 
  filter(is.na(drop)) %>% 
  select(-drop)

#add the rest of the data
scores_key = full_join(samps_key, scores) %>% 
  mutate(zsc_pass = ifelse(z_score >= 0.5, "pass", "fail"),
         zsc_pass = ifelse(is.na(zsc_pass), "fail", zsc_pass)) %>%
  #create the T/F-neg/pos infomration in new column
  mutate(type = ifelse(Host != "none" & Prediction == "Correct", "T-P",
                       ifelse(Host != "none" & Prediction == "Correct", "F-N",
                              ifelse(Host != "none" & Prediction == "Incorrect", "F-P",
                                     ifelse(Host != "none" & Prediction == "Incorrect", "T-N", NA)))),
         #adjust for hosts that aren't possible
         type = ifelse(Host == "none" & Prediction == "Correct", "T-P",
                       ifelse(Host == "none" & Prediction == "Correct", "T-N",
                              ifelse(Host == "none" & Prediction == "Incorrect", "F-P",
                                     ifelse(Host == "none" & Prediction == "Incorrect", "F-N", type)))),
         filt_type = ifelse(Host != "none" & Prediction == "Correct" & zsc_pass == "pass", "T-P",
                            ifelse(Host != "none" & Prediction == "Correct" & zsc_pass == "fail", "F-N",
                                   ifelse(Host != "none" & Prediction == "Incorrect" & zsc_pass == "pass", "F-P",
                                          ifelse(Host != "none" & Prediction == "Incorrect" & zsc_pass == "fail", "T-N", NA)))),
         #adjust for hosts that aren't possible
         filt_type = ifelse(Host == "none" & Prediction == "Correct" & zsc_pass == "pass", "T-P",
                            ifelse(Host == "none" & Prediction == "Correct" & zsc_pass == "fail", "T-N",
                                   ifelse(Host == "none" & Prediction == "Incorrect" & zsc_pass == "pass", "F-P",
                                          ifelse(Host == "none" & Prediction == "Incorrect" & zsc_pass == "fail", "F-N", filt_type)))))#,
#         Host = fct_relevel(Host, "none", after = Inf))

#

#summarise to get counts at different levels of grouping
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

#Sensitivity Specificity calculations

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

#make plots
  filter(Host != "none",
         !is.na(Score)) %>% 
  mutate(Prediction = fct_rev(Prediction)) %>% 
  ggplot(aes(y = Host, x = Score)) + 
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("gray70", "gray30")) + 
  #scale_x_log10() +
  facet_grid(experiment ~ sample, scales = "free") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0A

p0A1 = scores_key %>% 
  filter(
         !is.na(Score),
         sample %in% c("1-N1", "1-N2", "1-N3")) %>% 
  ggplot(aes(y = Host, x = Score)) + 
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("black","gray70")) + 
  #scale_x_log10() +
  facet_wrap(sample~., nrow=3, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle =0, hjust = 1))

p0A1

p0A2 = scores_key %>% 
  filter(
         !is.na(Score),
         sample %in% c("1-D1", "1-D2", "1-D3")) %>% 
  ggplot(aes(y = Host, x = Score)) + 
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("black","gray70")) + 
  #scale_x_log10() +
  facet_wrap(sample~., nrow=3, scales = "free_x", dir = "v") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 0, hjust = 1))

p0A2

p0A3 = scores_key %>% 
  filter(
    !is.na(Score),
    sample %in% c("1-B1", "1-B3")) %>% 
  ggplot(aes(y = Host, x = Score)) + 
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("black","gray70")) + 
  #scale_x_log10() +
  facet_wrap(sample~., nrow=3, scales = "free_x", dir = "v") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0A3

p0A4 = scores_key %>% 
  filter(
         !is.na(Score),
         sample %in% c("2-N1", "2-N2", "2-N3")) %>% 
  ggplot(aes(y = Host, x = Score)) + 
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("black","gray70")) + 
  #scale_x_log10() +
  facet_wrap(sample~., nrow=3, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0A4

p0A5 = scores_key %>% 
  filter(
    !is.na(Score),
    sample %in% c("3-N1", "3-N3")) %>% 
  ggplot(aes(y = Host, x = Score)) + 
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("black","gray70")) + 
  #scale_x_log10() +
  facet_wrap(sample~., nrow=3, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0A5

p0B = scores_key %>% 
  filter(
         !is.na(Score)) %>% 
  arrange(filt_type) %>% 
  ggplot(aes(y = Host, x = z_score)) + 
  geom_vline(xintercept = 0.5, color = "gray", linefilt_type = 2, linewidth = 1) +
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("black","gray70")) + 
  facet_wrap(sample~., nrow=3, scales = "fixed")+
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0B

p0B1 = scores_key %>% 
  filter(
         !is.na(Score),
  sample %in% c("1-N1", "1-N2", "1-N3")) %>% 
  ggplot(aes(y = Host, x = z_score)) + 
  geom_vline(xintercept = 0.5, color = "red", linetype = "dotted", linewidth = 1) +
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("black","gray70")) + 
  xlim(-1, 4)+
  facet_wrap(sample~., nrow=3, scales = "free") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 0, hjust = 1))

p0B1

p0B2 = scores_key %>% 
  filter(Host != "none",
         !is.na(Score),
         sample %in% c("1-B1", "1-B3")) %>% 
  ggplot(aes(y = Host, x = z_score)) + 
  geom_vline(xintercept = 0.5, color = "red", linetype = 'dotted', linewidth = 1) +
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("black","gray70")) + 
  xlim(-1, 4)+
  facet_wrap(sample~., nrow=3, scales = "free")+
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0B2

p0B3 = scores_key %>% 
  filter(Host != "none",
         !is.na(Score),
         sample %in% c("1-D1", "1-D2", "1-D3")) %>% 
  ggplot(aes(y = Host, x = z_score)) + 
  geom_vline(xintercept = 0.5, color = "red", linetype = 'dotted', linewidth = 1) +
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("black","gray70")) + 
  xlim(-1, 4)+
  facet_wrap(sample~., nrow=3, scales = "free")+
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 0, hjust = 1))

p0B3

p0B4 = scores_key %>% 
  filter(Host != 'none',
         !is.na(Score),
         sample %in% c("2-N1", "2-N2", "2-N3")) %>% 
  ggplot(aes(y = Host, x = z_score)) + 
  geom_vline(xintercept = 0.5, color = "red", linetype = 'dotted', linewidth = 1) +
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("black","gray70")) + 
  xlim(-1, 4)+
  facet_wrap(sample~., nrow=3, scales = "free") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 0, hjust = 1))

p0B4

p0B5 = scores_key %>% 
  filter(Host != "none",
         !is.na(Score),
         sample %in% c("3-N1", "3-N3")) %>% 
  ggplot(aes(y = Host, x = z_score)) + 
  geom_vline(xintercept = 0.5, color = "red", linetype = 'dotted', linewidth = 1) +
  geom_point(aes(fill = Prediction), shape = 21, color = "black", size = 1.5, alpha = 0.9, position = position_jitter(height = 0.3, width = 0)) +
  scale_fill_manual(values = c("black","gray70")) +
  xlim(-1, 4)+
  facet_wrap(sample~., nrow=3, scales = "free") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 0, hjust = 1))

p0B5

p0C1 = scores_key %>% 
  filter(Host != "none",!is.na(Score), sample %in% c("1-N1", "1-N2", "1-N3")) %>% 
  ggplot(aes(y = Host, x = Phage)) + 
  geom_tile(aes(fill = type), color = "white") +
  scale_fill_manual(values = c("T-P" = "black", "T-N" = "white", "F-P" = "gray60", "F-N" = "gray60"), 
                    labels = c("T-P" = "True Positive", "T-N" = "True Negative", "F-P" = "False Positive", "F-N" = "False Negative"), 
                    na.value = NA) + 
  facet_wrap(sample~., nrow=3, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0C1

p0C2 = scores_key %>% 
  filter(Host != "none",!is.na(Score),sample %in% c("1-D1", "1-D2", "1-D3")) %>% 
  ggplot(aes(y = Host, x = Phage)) + 
  geom_tile(aes(fill = type), color = "white") +
  scale_fill_manual(values = c("T-P" = "black", "T-N" = "gray90", "F-P" = "gray60", "F-N" = "gray60"), 
                    labels = c("T-P" = "True Positive", "T-N" = "True Negative", "F-P" = "False Positive", "F-N" = "False Negative"), 
                    na.value = NA) + 
  facet_wrap(sample~., nrow=3, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0C2

p0C3 = scores_key %>% 
  filter(Host != "none",!is.na(Score),sample %in% c("1-B1", "1-B3")) %>% 
  ggplot(aes(y = Host, x = Phage)) + 
  geom_tile(aes(fill = type), color = "white") +
  scale_fill_manual(values = c("T-P" = "black", "T-N" = "gray90", "F-P" = "gray60", "F-N" = "gray60"), 
                    labels = c("T-P" = "True Positive", "T-N" = "True Negative", "F-P" = "False Positive", "F-N" = "False Negative"), 
                    na.value = NA) + 
  facet_wrap(sample~., nrow=3, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0C3

p0C4 = scores_key %>% 
  filter(Host != "none",!is.na(Score),sample %in% c("2-N1", "2-N2", "2-N3")) %>% 
  ggplot(aes(y = Host, x = Phage)) + 
  geom_tile(aes(fill = type), color = "white") +
  scale_fill_manual(values = c("T-P" = "black", "T-N" = "gray90", "F-P" = "gray60", "F-N" = "gray60"), 
                    labels = c("T-P" = "True Positive", "T-N" = "True Negative", "F-P" = "False Positive", "F-N" = "False Negative"), 
                    na.value = NA) + 
  facet_wrap(sample~., nrow=3, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0C4

p0C5 = scores_key %>% 
  filter(Host != "none",!is.na(Score),sample %in% c("3-N1", "3-N3")) %>% 
  ggplot(aes(y = Host, x = Phage)) + 
  geom_tile(aes(fill = type), color = "white") +
  scale_fill_manual(values = c("T-P" = "black", "T-N" = "gray90", "F-P" = "gray60", "F-N" = "gray60"), 
                    labels = c("T-P" = "True Positive", "T-N" = "True Negative", "F-P" = "False Positive", "F-N" = "False Negative"), 
                    na.value = NA) + 
  facet_wrap(sample~., nrow=3, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0C5

p0D = scores_key %>% 
  filter(Host != "none",!is.na(Score)) %>% 
  ggplot(aes(y = Host, x = Phage)) + 
  geom_tile(aes(fill = ifelse(z_score>0.5, type, NA)), color = "white") +
  scale_fill_brewer(palette = "Greys", direction = 1, na.value = NA)+ 
  facet_wrap(sample~., nrow=3, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0D

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

#sensitivity_Specificity_plots
p0E <- sens_spec_long %>%
  filter(sample_type %in% c("1-D", "1-B", "1-N"),
         Metric %in% c("Avg_Filt_Sensitivity", "Avg_Filt_Specificity")) %>%
  ggplot(aes(x = sample_type, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Avg_Filt_Specificity" = "#32373AFF", "Avg_Filt_Sensitivity" = "#9F9994FF"),
                    labels = c("Avg_Filt_Specificity" = "Specificity", "Avg_Filt_Sensitivity" = "Sensitivity")) +
  theme_classic() +
  labs(, x = "Synthetic Communities", y = "Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p0E

p0F <- sens_spec_long %>%
  filter(sample_type %in% c("1-N", "2-N", "3-N"),
         Metric %in% c("Avg_Filt_Sensitivity", "Avg_Filt_Specificity")) %>%
  ggplot(aes(x = sample_type, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Avg_Filt_Specificity" = "#32373AFF", "Avg_Filt_Sensitivity" = "#9F9994FF"),
                    labels = c("Avg_Filt_Specificity" = "Specificity", "Avg_Filt_Sensitivity" = "Sensitivity")) +
  theme_classic() +
  labs(, x = "Synthetic Communities", y = "Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p0F

#Figure_1
Fig_1 <- ggarrange(p0A1, p0B1, p0C1, p0D1,  ncol = 4, align = "h")
Fig_1

#Figure_2 
Fig_2_A<- ggarrange(p0B3, p0D2, ncol = 2, align = "h")
Fig_2_A

Fig_2_B <- ggarrange( p0B2, p0D3, ncol = 2, align ='h')
Fig_2_B


#Figure_3
Fig_3_A <- ggarrange(p0B4, p0D4, ncol = 2, align = 'h')
Fig_3_A

Fig_3_B <- ggarrange(p0B5, p0D5, ncol = 2, align = 'h')
Fig_3_B

p0F

syncom1_names <- c(
  `1-N1` = "SynCom 1 replicate 1",
  `1-N2` = "SynCom 1 replicate 2",
  `1-N3` = "SynCom 1 replicate 3"
)

p0E1 = scores_key %>% 
  filter(Host != "none",
         !is.na(Score),
         sample %in% c("1-N1", "1-N2", "1-N3"))%>% 
  ggplot(aes(y = Host, x = Phage)) + 
  geom_tile(aes(fill = ifelse(z_score>0.5, filt_type, NA)), color = "white") +
  scale_fill_manual(values = c("T-P" = "black", "F-P" = "gray30", "F-N" = "gray60", "T-N" = "gray90")) + 
  facet_grid(~ sample, labeller = as_labeller(syncom1_names)) +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1))

p0E1

#### SOM

SOM1 <- scores_key %>%  
  filter(sample %in% c("1-N1", "1-N2", "1-N3")) %>%
  ggplot(aes(x=Phage, y=Score, color=Prediction))+
  geom_point()+
  scale_colour_manual(values = setNames(c('black', 'gray'), c("Correct", "Incorrect"))) +
  theme_classic()+
  facet_wrap(~sample, labeller = as_labeller(syncom1_names))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
SOM1

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

count_long <- count %>%
  pivot_longer(
    cols = c(z_pass_count, z_fail_count), names_to = 'zpass', values_to = 'count')

som2 <- count_long %>%
  filter(sample %in% c("1-N1", "1-N2", "1-N3")) %>%
  ggplot( aes( x = count, y = sample, fill = zpass)) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c("z_pass_count" = "black", "z_fail_count" = "grey70"),
                    labels = c("z_pass_count" = "Z-Score >= 0.54", "z_fail_count" = "Z-Score =< 0.5")) +
  labs(x = "Percentage of Linkages", y = "Synthetic Community") +
  theme_classic()
som2



#roc curve
library(pROC)

# Convert Prediction to binary values: 1 for "Correct" and 0 for "Incorrect"
scores_key$Prediction_binary <- ifelse(scores_key$Prediction == "Correct", 1, 0)
syncom1 <- scores_key %>%
  filter(sample %in% c("1-N1", "1-N2", "1-N3"))
roc_obj <- roc(syncom1$Prediction_binary, syncom1$z_score, direction = "<")

# Plot the ROC curve
plot(roc_obj, col = "blue", main = "ROC Curve for z_score Thresholds")

# Generate the ROC curve using Score as the predictor
roc_score <- roc(syncom1$Prediction_binary, syncom1$Score, direction = "<")

# Plot the ROC curve
plot(roc_score, col = "red", main = "ROC Curve for Score Thresholds")

# Add the z_score ROC curve to the same plot for comparison
plot(roc_obj, col = "blue", add = TRUE)

# Add a legend
legend("bottomright", legend = c("z_score", "Score"), col = c("blue", "red"), lwd = 2)

# Calculate AUC for z_score
auc_z_score <- auc(roc_obj)
print(paste("AUC for z_score:", auc_z_score))

# Calculate AUC for Score
auc_score <- auc(roc_score)
print(paste("AUC for Score:", auc_score))

png("roc_curve_mock1.png", width = 400, height = 400)

# Plot both ROC curves on the same plot
plot(roc_score, col = "red", main = "ROC Curve for SynCom 1 Replicates")
plot(roc_obj, col = "blue", add = TRUE)

# Add a legend with AUC values
legend("bottomright", 
       legend = c(paste("z_score (AUC =", round(auc_z_score, 2), ")"), 
                  paste("Score (AUC =", round(auc_score, 2), ")")),
       col = c("blue", "red"), 
       lwd = 2)
dev.off()
highest_specificity_threshold <- coords(roc_obj, "local maximas", ret = "threshold", input = "specificity", transpose = TRUE)
print(paste("Threshold for z_score with highest specificity:", highest_specificity_threshold))


