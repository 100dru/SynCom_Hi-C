import pandas as pd
import re

# Load the input CSV file
df = pd.read_csv("phage_host_nomalized_contact.csv")

# Function to remove trailing numbers only from contig names that do not start with "PHAGE"


# Calculate the Z-score for the 'Score' column
mean_score = df['Score'].mean()
std_dev_score = df['Score'].std()
df['Z_Score'] = (df['Score'] - mean_score) / std_dev_score


# Save the result to a new CSV file
df.to_csv("phage_host_zscore.csv", index=False)
