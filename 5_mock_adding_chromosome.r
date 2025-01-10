setwd("C:/Users/shatadru.2/OneDrive - The Ohio State University/LAB/Hi-C/informatic_pipeline/6th_test/metacc")

library(tidyverse)
library(dplyr)
library(ggplot2)
library(readxl)
library(readr)
library(ggpubr)
library(RColorBrewer)
library(viridisLite)
library(gplots)

all_phage_host <- read_csv("all_phage_host.csv")
View(all_phage_host)

#Summation of chromosomes
## chromosome summation 13-15

### Step 1: Group by contig2_name and Mock, then filter rows
filtered_df <- all_phage_host %>%
  group_by(contig2_name, Mock, Prediction) %>%
  filter(contig1_name %in% c("HOST_PSA-13-15_CHROM1", "HOST_PSA-13-15_CHROM2"))

### Step 2: Calculate the average of "Score" for each group
sums <- filtered_df %>%
  summarise(total_score = sum(Score), .groups = "drop")

### Step 3: Create new rows with the calculated sums
new_rows <- data.frame(
  Mock = sums$Mock,           # Mock values from the summarized data
  contig1_name = "HOST_PSA-13-15",
  contig2_name = sums$contig2_name,     # contig2_name values from the summarized data
  Score = sums$total_score,
  Prediction = sums$Prediction)

### Append the new rows to the original DataFrame
all_phage_host <- bind_rows(all_phage_host, new_rows)

## chromosome summation H71

### Step 1: Group by contig2_name and Mock, then filter rows
filtered_df <- all_phage_host %>%
  group_by(contig2_name, Mock, Prediction) %>%
  filter(contig1_name %in% c("HOST_PSA-H71_CONTIG1", "HOST_PSA-H71_CONTIG3"))

### Step 2: Calculate the average of "Score" for each group
sums <- filtered_df %>%
  summarise(total_score = sum(Score), .groups = "drop")

### Step 3: Create new rows with the calculated sums
new_rows <- data.frame(
  Mock = sums$Mock,           # Mock values from the summarized data
  contig1_name = "HOST_PSA-H71",
  contig2_name = sums$contig2_name,     # contig2_name values from the summarized data
  Score = sums$total_score,
  Prediction = sums$Prediction)

### Append the new rows to the original DataFrame
all_phage_host <- bind_rows(all_phage_host, new_rows)

## chromosome summation H105

### Step 1: Group by contig2_name and Mock, then filter rows
filtered_df <- all_phage_host %>%
  group_by(contig2_name, Mock, Prediction) %>%
  filter(contig1_name %in% c("HOST_PSA-H105_CONTIG1", "HOST_PSA-H105_CONTIG2"))

### Step 2: Calculate the average of "Score" for each group
sums <- filtered_df %>%
  summarise(total_score = sum(Score), .groups = "drop")

### Step 3: Create new rows with the calculated sums
new_rows <- data.frame(
  Mock = sums$Mock,           # Mock values from the summarized data
  contig1_name = "HOST_PSA-H105",
  contig2_name = sums$contig2_name,     # contig2_name values from the summarized data
  Score = sums$total_score,
  Prediction = sums$Prediction)

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

write.table(all_phage_host, file = "all_phage_host_chromosome_sum.csv", sep = ',')

## calculating z score
mock_1N1$z_score <- (mock_1N1$Score - mean(mock_1N1$Score, na.rm = TRUE)) / sd(mock_1N1$Score, na.rm = TRUE)
mock_1N2$z_score <- (mock_1N2$Score - mean(mock_1N2$Score, na.rm = TRUE)) / sd(mock_1N2$Score, na.rm = TRUE)
mock_1N3$z_score <- (mock_1N3$Score - mean(mock_1N3$Score, na.rm = TRUE)) / sd(mock_1N3$Score, na.rm = TRUE)
mock_1B1$z_score <- (mock_1B1$Score - mean(mock_1B1$Score, na.rm = TRUE)) / sd(mock_1B1$Score, na.rm = TRUE)
mock_1B3$z_score <- (mock_1B3$Score - mean(mock_1B3$Score, na.rm = TRUE)) / sd(mock_1B3$Score, na.rm = TRUE)
mock_1D1$z_score <- (mock_1D1$Score - mean(mock_1D1$Score, na.rm = TRUE)) / sd(mock_1D1$Score, na.rm = TRUE)
mock_1D2$z_score <- (mock_1D2$Score - mean(mock_1D2$Score, na.rm = TRUE)) / sd(mock_1D2$Score, na.rm = TRUE)
mock_1D3$z_score <- (mock_1D3$Score - mean(mock_1D3$Score, na.rm = TRUE)) / sd(mock_1D3$Score, na.rm = TRUE)


