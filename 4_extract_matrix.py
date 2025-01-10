import numpy as np
import pandas as pd
import scipy.sparse as sp

# Load the sparse matrix from the npz file
matrix = sp.load_npz('Normalized_contact_matrix.npz')

#read contig data
contig_data = pd.read_csv('contig_info.csv')

# Convert the sparse matrix to a dense NumPy array
dense_matrix = matrix.toarray()

# Create a Pandas DataFrame from the dense array
df = pd.DataFrame(dense_matrix)

# You can optionally set column and index names if needed
df.columns = contig_data["Contig name"]
df.index = contig_data["Contig name"]

df.index.name = 'contig1_name'
df.columns.name = 'contig2_name'

# Reset the index to make 'contig1_name' a regular column
df.reset_index(inplace=True)
df.rename(columns={'index': 'contig1_name'}, inplace=True)

#print(df)

# Now, 'df' contains the data in table form with the desired row and column names
#print(df)

df_longer = pd.melt(df, id_vars=['contig1_name'], var_name='contig2_name', value_name='Score')

# Filter out rows with '0.00000' values in the 'Score' column
df_longer = df_longer[df_longer['Score'] != 0.00000]
#print(df_longer)

# Assuming you have already created and populated the 'df' DataFrame
df_longer.to_csv('all_nomalized_contact.csv', index=False)

# Assuming df_longer contains your data

# Filter and keep only rows where 'contig1_name' starts with "HOST" and 'contig2_name' starts with "PHAGE"
filtered_df = df_longer[df_longer['contig1_name'].str.startswith("HOST") & df_longer['contig2_name'].str.startswith("PHAGE")]

# Now, filtered_df contains the desired rows
#print(filtered_df.reset_index())
filtered_df.to_csv('phage_host_nomalized_contact.csv', index=False)
