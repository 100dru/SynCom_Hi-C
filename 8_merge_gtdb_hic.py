import pandas as pd

# Load the tables
# Assuming table1.csv and table2.tsv are the filenames
table1 = pd.read_csv('phage_host_zscore.csv')  # First table (comma-separated)
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

# Save or display the result
result.to_csv('hic_gtdb.csv', index=False)
