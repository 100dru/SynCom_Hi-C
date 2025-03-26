import numpy as np
import pandas as pd
import scipy.sparse as sp

# Load the sparse matrix from the npz file
# BB: You've loaded a sparse matrix - it's sparse because it doesn't hold all the 0s
matrix = sp.load_npz('Normalized_contact_matrix.npz')
sp_matrix_df = pd.DataFrame.sparse.from_spmatrix(matrix)
# BB: You're taking a sparse matrix and immediately making it dense, don't do that because you've lost
# all of the performance gains and memory-saving features of the sparse structure
# BB: What I've done here is use a pandas sparse dataframe to have all the features of a pandas dataframe with the
# performance of a sparse matrix

# BB: No
# # Convert the sparse matrix to a dense NumPy array
# dense_matrix = matrix.toarray()
# # Create a Pandas DataFrame from the dense array
# df = pd.DataFrame(dense_matrix)

#read contig data
# BB: Ensure you're handling headers and indices if available
contig_data = pd.read_csv('contig_info.csv', header=0, index_col=False)

# BB: None of this needs to be done until the end. Know that your columns and indices are identical.
# You can optionally set column and index names if needed
# df.columns = contig_data["Contig name"]
# df.index = contig_data["Contig name"]
#
# df.index.name = 'contig1_name'
# df.columns.name = 'contig2_name'
#
# # Reset the index to make 'contig1_name' a regular column
# df.reset_index(inplace=True)
# df.rename(columns={'index': 'contig1_name'}, inplace=True)

# Now, 'df' contains the data in table form with the desired row and column names

# BB: Instead, we take our pandas sparse datframe and convert it to a CSR (Compressed Sparse Row) matrix
# This is literally designed to work on these "pairwise" comparisons, all without ever using more memory
csr_mat = sp_matrix_df.sparse.to_coo().tocsr()
print((csr_mat.data.nbytes + csr_mat.indptr.nbytes + csr_mat.indices.nbytes) / (1024 * 1024))

# BB: Technically, there shouldn't be 0s, but we're extracting the non-zero values (i.e. your data)
# BB: Because of the CSR format, the rows and column positions (i.e. indices) are stored in a paired fashion
row_idx, col_idx = csr_mat.nonzero()  # BB: Extract nonzero indices and values
scores = csr_mat.data  # BB: The CSR matrix stores the values in its matrix in its .data attribute

# BB: NOW we have all the pairs that have scores AND their data, now we want to convert from the index values (0, 1, 2)
# into names that we humans can use
name1 = contig_data["Contig name"].iloc[row_idx].values
name2 = contig_data["Contig name"].iloc[col_idx].values

# BB: With our names for both rows and indices, and the scores from the data, we can create a new dataframe that
# only has the values you want
df_longer = pd.DataFrame({'contig1_name': name1, 'contig2_name': name2, 'Score': scores})

# BB: NONE of this is necessary anymore. With sparse data structures, you - by definition - never need to filter out 0s
# BB: When you melt such a huge dataframe, you'll be getting n^2 rows.
# BB: There are 2,653,164 values, but with a dense dataframe and melt you're looking at 7,477,233,841 rows!
# df_longer = pd.melt(df, id_vars=['contig1_name'], var_name='contig2_name', value_name='Score')
#
# # Filter out rows with '0.00000' values in the 'Score' column
# df_longer = df_longer[df_longer['Score'] != 0.00000]
# # #print(df_longer)

# BB: Rest of code is fine

# Assuming you have already created and populated the 'df' DataFrame
df_longer.to_csv('all_nomalized_contact.csv', index=False)

# Filter and keep only rows where 'contig1_name' starts with "HOST" and 'contig2_name' starts with "PHAGE"
filtered_df = df_longer[df_longer['contig1_name'].str.startswith("MAG") & df_longer['contig2_name'].str.startswith("Virus")]

# Now, filtered_df contains the desired rows
#print(filtered_df.reset_index())
filtered_df.to_csv('phage_host_nomalized_contact.csv', index=False)