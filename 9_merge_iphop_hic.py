import pandas as pd

# Load the merged table and the new table
merged_table = pd.read_csv('hic_gtdb.csv')  # Adjust filename if needed
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
final_result.to_csv('hic_gtdb_iphop.csv', index=False)
