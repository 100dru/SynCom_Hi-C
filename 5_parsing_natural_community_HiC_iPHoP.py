import pandas as pd
import re

# Load the input CSV file from the extracted matrix from MetaCC output into a DataFrame
df = pd.read_csv('phage_host_nomalized_contact_S.csv')

# Extract the MAG identifier from the "contig1_name" column [this is for the naming convention of our dataset]
df["MAG"] = df["contig1_name"].str.extract(r"MAG_(.*?)~")

# Aggregate the linkages by each MAG
#Group the DataFrame by the "MAG" column and sum the "Score" column
grouped = df.groupby("MAG", as_index=False)["Score"].sum()

# Save the grouped DataFrame to a new CSV file
output_file = "grouped_scores_by_mag_S.csv"
grouped.to_csv(output_file, index=False)

print(f"Grouped scores have been saved to {output_file}.")

# Calculate the mean and standard deviation of the "Score" column
mean_score = df['Score'].mean()
std_dev_score = df['Score'].std()

# Calculate the Z-scores for the "Score" column and add them as a new column
df['Z_Score'] = (df['Score'] - mean_score) / std_dev_score

# Load additional tables from CSV files
#Load GTDB predictions for the MAGS
#Load iPHOP predictions for the viruses

table1 = df  # First table (comma-separated)
table2 = pd.read_csv('gtdb_combined_summary.tsv', sep='\t')  # Second table (tab-separated)
iphop = pd.read_csv('iphop.csv')  # Third table (comma-separated)

# Keep the row with the highest "Confidence score" for each Virus in the iphop table
iphop = iphop.loc[iphop.groupby('Virus')['Confidence score'].idxmax()]

# Extract the Virus identifier from the "contig2_name" column
df["Virus"] = df["contig2_name"]

# Create a list to store rows with additional information from table2 and iphop
matches = []

# Loop over each row in table1 and find matching rows in table2 and iphop
for _, row1 in table1.iterrows():
    mag = row1['MAG']
    virus = row1['Virus']
    
    # Find matching rows in table2 where the user_genome matches the MAG
    matching_rows_table2 = table2[table2['user_genome'].str.extract(r'EMERGEV2_(.*?)(?:\s|$)')[0] == mag]
    
    # Find matching rows in iphop where the Virus matches
    matching_rows_iphop = iphop[iphop['Virus'] == virus]
    
    # Merge row1 with the first matching row in table2 and iphop if matches are found
    combined_row = row1.to_dict()  # Convert row1 to a dictionary
    if not matching_rows_table2.empty:
        row2 = matching_rows_table2.iloc[0]  # Take the first matching row from table2
        for col in table2.columns:
            combined_row[col] = row2[col]  # Add columns from row2
    
    if not matching_rows_iphop.empty:
        row3 = matching_rows_iphop.iloc[0]  # Take the first matching row from iphop
        for col in iphop.columns:
            combined_row[col] = row3[col]  # Add columns from row3
    
    matches.append(combined_row)  # Append the combined row to the matches list

# Convert the list of matches to a DataFrame
result = pd.DataFrame(matches)

# Save the resulting DataFrame to a CSV file
result.to_csv('Bog_hic_gtdb_iphop_grouped.csv', index=False)
print("Merged table has been saved to Bog_hic_gtdb_iphop_grouped.csv.")
