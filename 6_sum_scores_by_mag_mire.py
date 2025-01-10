import pandas as pd
import re

df = pd.read_csv('phage_host_nomalized_contact.csv')

# Extract the MAG from the six letters after "HOST_ENA|"
df["MAG"] = df["contig1_name"].str.extract(r"HOST_ENA\|(\w{6})")

# Group by MAG and contig2_name, summing the scores
grouped = df.groupby(["MAG", "contig2_name"], as_index=False)["Score"].sum()

# Save the result to a new CSV file
output_file = "grouped_scores_by_mag.csv"
grouped.to_csv(output_file, index=False)

print(f"Grouped scores have been saved to {output_file}.")

mean_score = df['Score'].mean()
std_dev_score = df['Score'].std()
df['Z_Score'] = (df['Score'] - mean_score) / std_dev_score

table1 = df  # First table (comma-separated)
table2 = pd.read_csv('gtdbtk.bac120.summary.tsv', sep='\t')  # Second table (tab-separated)

# Create a list to store rows with additional info from table2
matches = []

# Loop over each row in table1 and search for substring matches with each user_genome in table2
for _, row1 in table1.iterrows():
    contig1_name = row1['contig1_name']
    
    # Find matching rows in table2 where user_genome is a substring in contig1_name
    matching_rows = table2[table2['user_genome'].apply(lambda x: x in contig1_name)]
    
    # If matches are found, merge row1 with each matching row in table2
    if not matching_rows.empty:
        for _, row2 in matching_rows.iterrows():
            combined_row = row1.to_dict()  # Convert row1 to dict
            for col in table2.columns:
                combined_row[col] = row2[col]  # Add columns from row2
            matches.append(combined_row)  # Append to matches list
    else:
        # If no match, append row1 with NaNs for table2 columns
        combined_row = row1.to_dict()
        for col in table2.columns:
            combined_row[col] = pd.NA  # Assign NaN for unmatched rows
        matches.append(combined_row)

# Convert the list of matches to a DataFrame
result = pd.DataFrame(matches)

merged_table =result  # Adjust filename if needed
new_table = pd.read_csv('iphop.csv')  # Replace with the actual filename for the new table

# List to store the combined rows
combined_rows = []

# Loop over each row in merged_table and check for substring matches with Virus in new_table
for _, row1 in merged_table.iterrows():
    contig2_name = row1['contig2_name']
    
    # Find matching rows in new_table where Virus is a substring in contig2_name
    matching_rows = new_table[new_table['Virus'].apply(lambda x: x in contig2_name)]
    
    # If matches are found, merge row1 with each matching row in new_table
    if not matching_rows.empty:
        for _, row2 in matching_rows.iterrows():
            combined_row = row1.to_dict()  # Convert row1 to dict
            for col in new_table.columns:
                combined_row[col] = row2[col]  # Add columns from row2
            combined_rows.append(combined_row)  # Append to combined_rows list
    else:
        # If no match, add row1 with NaNs for columns from new_table
        combined_row = row1.to_dict()
        for col in new_table.columns:
            combined_row[col] = pd.NA  # Assign NaN for unmatched rows
        combined_rows.append(combined_row)

# Convert the list of combined rows to a DataFrame
final_result = pd.DataFrame(combined_rows)

# Save or display the result
final_result.to_csv('bog_hic_gtdb_iphop_grouped_mag.csv', index=False)
