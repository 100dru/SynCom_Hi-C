library(tidyverse)
library(ggpubr)
library(readxl)

# Read all data files
palsa         <- read_csv("palsa.csv") %>% mutate(sample = "palsa") #Hi-C links from Palsa
bog           <- read_csv("bog.csv")   %>% mutate(sample = "bog") #Hi-C links from Bog
fen           <- read_csv("fen.csv")   %>% mutate(sample = "fen") #Hi-C links from Fen
peat          <- bind_rows(palsa, bog, fen)
ip_test       <- read_excel("pbio.3002083.s015.xlsx") #training dataset from iPHoP
iphop         <- read_csv("iphop.csv") #iPHoP Host_prediction_to_genus_m90.csv,
iphop_genome  <- read.csv("iphop_genome.csv") # iPHoP Host_prediction_to_genome_mXX.csv
gtdb          <- read_tsv("gtdb_combined_summary.tsv") # GTDB outputs for all the MAGs

# =========================
# Section 2: Data Curation
# =========================

# Clean peat contig names
peat <- peat %>%
  mutate(contig1_name = str_remove(contig1_name, "~.*$"))

# Prepare GTDB metadata
gtdb <- gtdb %>%
  mutate(
    extracted_user_genome = str_extract(user_genome, "EMERGEV2_([^\\s]+)"),
    user_genome = str_replace_all(user_genome, "EMERGEV2", "MAG")
  )

# Join GTDB to peat data and split taxonomy
peat_gtdb <- peat %>%
  left_join(gtdb, by = c("contig1_name" = "user_genome")) %>%
  separate(classification, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";", remove = FALSE)

# Separate taxonomy for iPHOP and iPHOP genome
iphop <- iphop %>%
  separate(`Host genus`, into = c("domain", "phylum", "class", "order", "family", "genus"), sep = ";", remove = FALSE)

iphop_genome <- iphop_genome %>%
  separate(`Host.taxonomy`, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";", remove = FALSE)

# Get top hits
peat_gtdb_high <- peat_gtdb %>%
  group_by(contig2_name) %>%
  slice_max(order_by = Z_Score, n = 1, with_ties = FALSE) %>%
  ungroup()

iphop_high <- iphop %>%
  group_by(Virus) %>%
  slice_max(order_by = `Confidence score`, n = 1, with_ties = FALSE) %>%
  ungroup()

iphop_genome_high <- iphop_genome %>%
  group_by(Virus) %>%
  slice_max(order_by = Confidence.score, n = 1, with_ties = FALSE) %>%
  ungroup()

# Calculate phylum and family counts for merging/regression
peat_counts         <- peat_gtdb %>% count(phylum, name = "peat_count")
iphop_counts        <- iphop %>% count(phylum, name = "iphop_count")
ip_test_counts      <- ip_test %>% count(phylum, name = "ip_test_count")
peat_counts_family  <- peat_gtdb %>% count(family, name = "peat_counts_family")
ip_test_counts_family <- ip_test %>% count(family, name = "ip_test_counts_family")
peat_counts_high    <- peat_gtdb_high %>% count(phylum, name = "peat_count")
iphop_counts_high   <- iphop_high %>% count(phylum, name = "iphop_count")
iphop_genome_counts_high   <- iphop_genome_high %>% count(phylum, name = "iphop_count")

# Merge counts as needed
merged_counts          <- inner_join(peat_counts, iphop_counts, by = "phylum")
merged_counts_high     <- inner_join(peat_counts_high, iphop_counts_high, by = "phylum")
merged_hic_genome_counts_high     <- inner_join(peat_counts_high, iphop_genome_counts_high, by = "phylum")
merge_hic_test_fam     <- left_join(peat_counts_family, ip_test_counts_family, by = "family")
merge_iphop_ip_test    <- left_join(iphop_counts, ip_test_counts, by = "phylum")

# Prepare family merges for regression (replace NA with 0)
merge_hic_test_fam <- merge_hic_test_fam %>%
  mutate(
    ip_test_counts_family = ifelse(is.na(ip_test_counts_family), 0, ip_test_counts_family),
    peat_counts_family = ifelse(is.na(peat_counts_family), 0, peat_counts_family)
  )

# =========================
# Section 3: Plotting
# =========================


# Whole phylum counts
lm_model <- lm(peat_count ~ iphop_count, data = merged_counts)
summary_stats <- summary(lm_model)
r_squared <- summary_stats$r.squared
p_value <- summary_stats$coefficients[2, 4]

ggplot(merged_counts, aes(x = iphop_count, y = peat_count)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(
    title = "Linear Regression between iphop and peat_gtdb phylum counts",
    x = "iphop phylum count",
    y = "peat_gtdb phylum count"
  ) +
  annotate("text", x = Inf, y = Inf, label = paste0("R² = ", round(r_squared, 3), "\nP = ", signif(p_value, 3)),
           hjust = 1.1, vjust = 1.5, size = 5, color = "black") +
  theme_classic2()

# Whole phylum counts - GENOME
lm_model3 <- lm(peat_count ~ iphop_count, data = merged_hic_genome_counts_high)
summary_stats3 <- summary(lm_model)
r_squared3 <- summary_stats$r.squared
p_value3 <- summary_stats$coefficients[2, 4]

ggplot(merged_hic_genome_counts_high, aes(x = iphop_count, y = peat_count)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(
    title = "Linear Regression between iphop and peat_gtdb phylum counts",
    x = "iphop phylum count",
    y = "peat_gtdb phylum count"
  ) +
  annotate("text", x = Inf, y = Inf, label = paste0("R² = ", round(r_squared3, 3), "\nP = ", signif(p_value3, 3)),
           hjust = 1.1, vjust = 1.5, size = 5, color = "black") +
  theme_classic2()


# Top hits phylum counts
lm_model_high <- lm(peat_count ~ iphop_count, data = merged_counts_high)
summary_stats_high <- summary(lm_model_high)
r_squared_h <- summary_stats_high$r.squared
p_value_h <- summary_stats_high$coefficients[2, 4]

plot_line_phyla <- ggplot(merged_counts_high, aes(x = iphop_count, y = peat_count)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(
    title = "Linear Regression between iphop and peat_gtdb phylum counts TOP HITS",
    x = "iphop phylum count",
    y = "peat_gtdb phylum count"
  ) +
  annotate("text", x = Inf, y = Inf, label = paste0("R² = ", round(r_squared_h, 3), "\nP = ", signif(p_value_h, 3)),
           hjust = 1.1, vjust = 1.5, size = 5, color = "black") +
  theme_classic2()
plot_line_phyla
ggsave("phylum_regression_plot.pdf", plot = plot_line_phyla, width = 8, height = 6)

# Family counts regression
lm_model2 <- lm(peat_counts_family ~ ip_test_counts_family, data = merge_hic_test_fam)
summary_stats2 <- summary(lm_model2)
r_squared2 <- summary_stats2$r.squared
p_value2 <- summary_stats2$coefficients[2, 4]

plot_line_phyla2 <- ggplot(merge_hic_test_fam, aes(x = ip_test_counts_family, y = peat_counts_family)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(
    title = "Comparison between the natural community and iPHoP training data",
    x = "Counts of families in iPHoP training data",
    y = "Counts of families in Hi-C inferred linkages in soil"
  ) +
  annotate("text", x = Inf, y = Inf, label = paste0("R² = ", round(r_squared2, 3), "\nP = ", signif(p_value2, 3)),
           hjust = 1.1, vjust = 1.5, size = 5, color = "black") +
  theme_classic2()
plot_line_phyla2

# =========================
# Section: Regression - iPHOP Test vs GTDB
# =========================

# Count phylum occurrences
ip_test_counts <- ip_test %>%
  count(phylum, name = "ip_test_count")

gtdb_counts <- gtdb %>%
  separate(classification, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";", remove = FALSE) %>%
  count(phylum, name = "gtdb_count")

# Merge counts
merged_iptest_gtdb <- full_join(gtdb_counts, ip_test_counts, by = "phylum")
merged_iptest_gtdb <- merged_iptest_gtdb %>%
  mutate(
    ip_test_count = replace_na(ip_test_count, 0),
    gtdb_count = replace_na(gtdb_count, 0)
  )


# Linear regression
lm_model_iptest_gtdb <- lm(gtdb_count ~ ip_test_count, data = merged_iptest_gtdb)
summary_stats_iptest_gtdb <- summary(lm_model_iptest_gtdb)
r_squared_iptest_gtdb <- summary_stats_iptest_gtdb$r.squared
p_value_iptest_gtdb <- summary_stats_iptest_gtdb$coefficients[2, 4]

# Plot

lot_iptest_gtdb <- ggplot(merged_iptest_gtdb, aes(x = ip_test_count, y = gtdb_count)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(
    x = "count of phylum in iPHOP training data",
    y = "count of phylum in all MAGs from peat sample"
  ) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("R² = ", round(r_squared_iptest_gtdb, 3), "\nP = ", signif(p_value_iptest_gtdb, 3)),
           hjust = 1.1, vjust = 1.5, size = 5, color = "black") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_classic2()

lot_iptest_gtdb
