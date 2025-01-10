
setwd("C:/Users/shatadru.2/OneDrive - The Ohio State University/LAB/Hi-C/informatic_pipeline/7th_test/metacc_mire")

library(ggplot2)
library(tidyr)
library (dplyr)
library(tibble)
library(pheatmap)
library(likert)
library(paletteer)
library(ggpubr)
library(ggpattern)
# Read the CSV file
data <- read.csv("C:/Users/shatadru.2/OneDrive - The Ohio State University/LAB/Hi-C/informatic_pipeline/7th_test/metacc_mire/bog_hic_gtdb_iphop_grouped_mag.csv")

#drop the hic predictions that do not match with iphop
data_clean <- drop_na(data)

subset_df <- data_clean %>%
  group_by(contig2_name) %>%
  slice(which.max(Z_Score)) %>%
  ungroup()

same_phylum <- subset(subset_df, hic_phyla == iphop_phyla)
same_class <- subset(subset_df, hic_class == iphop_class)
same_order <- subset(subset_df, hic_order == iphop_order)
same_family <- subset(subset_df, hic_family == iphop_family)
same_genus <- subset(subset_df, hic_genus == iphop_genus)

not_common <- anti_join(subset_df, same_phylum)

com <- data.frame(Taxonomy = c("Phylum", "Class", "Order", "Family", "Genus"),
                  match =  c(48,44,42,32,27),
                  no_match = c(17,21,23,33,38),
                  z_score_above = c(15,15,14,10,10))

com$Taxonomy <- factor(com$Taxonomy, 
                       levels = c("Phylum", "Class", "Order", "Family", "Genus"))


## the no match count above and below z score is calculated manually

com <- com %>%
  mutate(total = match + no_match,
         match_percent = match / total * 100,
         no_match_belowz_percent = (no_match -2) / total * 100,
         no_match_abovez_percent = 2 / total * 100,
         above_z_score_percent = ifelse(match > z_score_above, z_score_above / total * 100, 0),
         below_z_score_percent = match_percent - above_z_score_percent)
# Reshape the data for plotting
# Reshape the data for plotting
com_long <- com %>%
  select(Taxonomy, no_match_belowz_percent, no_match_abovez_percent, above_z_score_percent, below_z_score_percent) %>%
  pivot_longer(cols = c(no_match_belowz_percent, no_match_abovez_percent, above_z_score_percent, below_z_score_percent),
               names_to = "type",
               values_to = "percent")
com_long$type = factor(com_long$type, levels = c("no_match_belowz_percent", "no_match_abovez_percent", "below_z_score_percent", "above_z_score_percent"))

# Plot the stacked bar chart
new_stack <- ggplot(com_long, aes(y = Taxonomy, x = percent, fill = type, pattern = type)) +
  geom_bar_pattern(stat = "identity", pattern_color = "black", color = "black") +
  scale_fill_manual(values = c("above_z_score_percent" = "black", 
                               "below_z_score_percent" = "grey", 
                               "no_match_abovez_percent" = "white", 
                               "no_match_belowz_percent" = "white"),
                    labels = c("above_z_score_percent" = "Same Prediction High Z Score", 
                               "below_z_score_percent" = "Same Prediction Low Z Score", 
                               "no_match_abovez_percent" = "Different Predictions High Z Score",
                               "no_match_belowz_percent" = "Different Predictions Low Z Score")) +
  scale_pattern_manual(values = c("above_z_score_percent" = "none", 
                                  "below_z_score_percent" = "none", 
                                  "no_match_abovez_percent" = "circle", 
                                  "no_match_belowz_percent" = "crosshatch")) +
  labs(y = "Taxonomy", fill = " ", pattern = " ") +
  ggtitle("Iphop vs Hi-C predictions before applying threshold") +
  theme_classic() +
  theme(legend.position = "bottom")
new_stack

ggsave("new_stack_iphop_hic.png", new_stack)


## now we apply z_score > 0.5 [subject to change]

z_com <- data.frame(
  Taxonomy = c("Phylum", "Class", "Order", "Family", "Genus"),
  match = c(15, 15, 14, 10, 10),
  no_match = c(2, 2, 3, 7, 7)
)

# Reverse the order of Taxonomy levels
z_com$Taxonomy <- factor(z_com$Taxonomy, 
                         levels = c("Phylum", "Class", "Order", "Family", "Genus"))

# Pivot to long format
z_com_long <- pivot_longer(
  z_com, 
  cols = c(match, no_match), 
  names_to = "type", 
  values_to = "value")
z_com_long$type <- factor(z_com_long$type, 
                         levels = c("no_match", "match"))

# Plot
new_stack_z <- ggplot(z_com_long, aes(y = Taxonomy, x = value, fill = type, pattern = type)) +
  geom_bar_pattern(stat = "identity", pattern_color = "black", color = "black") +
  scale_fill_manual(values = c("match" = "black", 
                               "no_match" = "white"),
                    labels = c("match" = "Same Prediction High Z Score", 
                               "no_match" = "Different Predictions High Z Score")) +
  scale_pattern_manual(values = c("match" = "none", 
                                  "no_match" = "circle")) +
  labs(y = "Taxonomy", fill = " ", pattern = " ") +
  ggtitle("Iphop vs Hi-C predictions after applying threshold") +
  theme_classic() +
  theme(legend.position = "bottom")

new_stack_z

hic_iphop <- ggarrange(new_stack, new_stack_z, nrow = 2, ncol = 1)
hic_iphop
