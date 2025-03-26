import pandas as pd
import re

# Load the CSV file
df = pd.read_csv('phage_host_nomalized_contact_S.csv')

# Extract the MAG from the string after "MAG_" and before "~"
df["MAG"] = df["contig1_name"].str.extract(r"MAG_(.*?)~")

# Group by MAG and sum the scores
grouped = df.groupby("MAG", as_index=False)["Score"].sum()

# Save the result to a new CSV file
output_file = "grouped_scores_by_mag_S.csv"
grouped.to_csv(output_file, index=False)

print(f"Grouped scores have been saved to {output_file}.")

mean_score = df['Score'].mean()
std_dev_score = df['Score'].std()
df['Z_Score'] = (df['Score'] - mean_score) / std_dev_score

table1 = df  # First table (comma-separated)
table2 = pd.read_csv('gtdb_combined_summary.tsv', sep='\t')  # Second table (tab-separated)
iphop = pd.read_csv('iphop.csv')  # Third table (comma-separated)

iphop = iphop.loc[iphop.groupby('Virus')['Confidence score'].idxmax()]

# Extract Virus identifier from contig2_name
df["Virus"] = df["contig2_name"]
 
# Create a list to store rows with additional info from table2 and iphop
matches = []

# Loop over each row in table1 and search for substring matches with each user_genome in table2
for _, row1 in table1.iterrows():
    mag = row1['MAG']
    virus = row1['Virus']
    
    # Find matching rows in table2 where the string after 'EMERGEV2_' matches the MAG
    matching_rows_table2 = table2[table2['user_genome'].str.extract(r'EMERGEV2_(.*?)(?:\s|$)')[0] == mag]
    
    # Find matching rows in iphop where the Virus matches
    matching_rows_iphop = iphop[iphop['Virus'] == virus]
    
    # If matches are found, merge row1 with the first matching row in table2 and iphop
    combined_row = row1.to_dict()  # Convert row1 to dict
    if not matching_rows_table2.empty:
        row2 = matching_rows_table2.iloc[0]  # Take the first matching row
        for col in table2.columns:
            combined_row[col] = row2[col]  # Add columns from row2
    
    if not matching_rows_iphop.empty:
        row3 = matching_rows_iphop.iloc[0]  # Take the first matching row
        for col in iphop.columns:
            combined_row[col] = row3[col]  # Add columns from row3
    
    matches.append(combined_row)  # Append to matches list

# Convert the list of matches to a DataFrame
result = pd.DataFrame(matches)

# Save or display the result
result.to_csv('Bog_hic_gtdb_iphop_grouped.csv', index=False)
print("Merged table has been saved to Bog_hic_gtdb_iphop_grouped.csv.")



