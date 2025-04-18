# Load required libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(readr)

# Define file paths for datasets
file_paths <- list(
  palsa = "Palsa_hic_gtdb_iphop_grouped.csv",
  bog = "Bog_hic_gtdb_iphop_grouped.csv",
  fen = "Fen_hic_gtdb_iphop_grouped.csv"
)

# Function to process links
process_links <- function(file_path, sample_name) {
  read_csv(file_path) %>%
    separate(classification,
             into = c("hic_Domain", "hic_Phylum", "hic_Class", "hic_Order", "hic_Family", "hic_Genus", "hic_Species"), sep = ";") %>%
    separate(`Host taxonomy`,
             into = c("ip_Domain", "ip_Phylum", "ip_Class", "ip_Order", "ip_Family", "ip_Genus", "ip_Species"), sep = ";") %>%
    filter(!is.na(ip_Domain)) %>%
    mutate(across(c(hic_Species, ip_Species), ~ na_if(., "s__"))) %>%
    group_by(contig2_name) %>%
    slice_max(Score, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(sample = sample_name)
}

# Process datasets
all_links <- bind_rows(
  lapply(names(file_paths), function(sample) process_links(file_paths[[sample]], sample))
)

# Add comparison columns
all_links <- all_links %>%
  mutate(
    Phylum = ifelse(hic_Phylum == ip_Phylum, 'same', 'different'),
    Domain = ifelse(hic_Domain == ip_Domain, 'same', 'different'),
    Class = ifelse(hic_Class == ip_Class, 'same', 'different'),
    Order = ifelse(hic_Order == ip_Order, 'same', 'different'),
    Family = ifelse(hic_Family == ip_Family, 'same', 'different'),
    Genus = ifelse(hic_Genus == ip_Genus, 'same', 'different'),
    Species = ifelse(hic_Species == ip_Species, 'same', 'different')
  )

# Function to calculate counts for each subset
calculate_counts <- function(data, z_score_threshold, threshold_label) {
  data %>%
    filter(Z_Score >= z_score_threshold) %>%
    select(Phylum, Domain, Class, Order, Family, Genus, Species) %>%
    pivot_longer(cols = everything(), names_to = "match_type", values_to = "match_status") %>%
    group_by(match_type, match_status) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(Z_Score_threshold = threshold_label)
}

# Define Z_Score thresholds
thresholds <- list(
  list(threshold = -0.6269769306, label = "None"),
  list(threshold = 0.5, label = "0.5")
)

# Calculate counts for each sample and threshold
combined_counts <- bind_rows(
  lapply(thresholds, function(t) {
    bind_rows(lapply(unique(all_links$sample), function(sample) {
      calculate_counts(all_links %>% filter(sample == sample), t$threshold, t$label)
    }))
  })
)

# Set factor levels
desired_order <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
combined_counts$match_type <- factor(combined_counts$match_type, levels = desired_order)
combined_counts$Z_Score_threshold <- factor(combined_counts$Z_Score_threshold, levels = c('None', '0.5'))
combined_counts <- combined_counts %>%
  mutate(match_status = ifelse(is.na(match_status), "NA", match_status))
combined_counts$match_status <- factor(combined_counts$match_status, levels = c('same', 'different', 'NA'))

# Plot the counts
bar_plot <- combined_counts %>%
  ggplot(aes(x = match_type, y = count, fill = match_status)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = c('black', 'grey80', 'white')) +
  facet_wrap(~ Z_Score_threshold, ncol = 2) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  xlab("Match Type") +
  ylab("Count") +
  ggtitle("Count of 'Same' vs 'Different' for Each Match Type Across Z_Score Thresholds and Samples")

# Display the plot
bar_plot
