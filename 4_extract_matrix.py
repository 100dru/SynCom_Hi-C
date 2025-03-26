import numpy as np  # Importing numpy library for numerical operations
import pandas as pd  # Importing pandas library for data manipulation and analysis
import scipy.sparse as sp  # Importing scipy.sparse for sparse matrix operations

matrix = sp.load_npz('Normalized_contact_matrix.npz')  # Loading a sparse matrix from a .npz file
sp_matrix_df = pd.DataFrame.sparse.from_spmatrix(matrix)  # Converting the sparse matrix to a sparse DataFrame

contig_data = pd.read_csv('contig_info.csv', header=0, index_col=False)  # Reading contig information from a CSV file

csr_mat = sp_matrix_df.sparse.to_coo().tocsr()  # Converting the sparse DataFrame to a CSR (Compressed Sparse Row) matrix
print((csr_mat.data.nbytes + csr_mat.indptr.nbytes + csr_mat.indices.nbytes) / (1024 * 1024))  # Printing the memory usage of the CSR matrix in MB

row_idx, col_idx = csr_mat.nonzero()  # Getting the row and column indices of non-zero elements in the CSR matrix
scores = csr_mat.data  # Extracting the data (scores) from the CSR matrix

name1 = contig_data["Contig name"].iloc[row_idx].values  # Mapping row indices to contig names
name2 = contig_data["Contig name"].iloc[col_idx].values  # Mapping column indices to contig names

df_longer = pd.DataFrame({'contig1_name': name1, 'contig2_name': name2, 'Score': scores})  # Creating a DataFrame with contig names and scores

df_longer.to_csv('all_nomalized_contact.csv', index=False)  # Saving the DataFrame to a CSV file

filtered_df = df_longer[df_longer['contig1_name'].str.startswith("MAG") & df_longer['contig2_name'].str.startswith("Virus")]  # Filtering the DataFrame for specific contig name patterns

filtered_df.to_csv('phage_host_nomalized_contact.csv', index=False)  # Saving the filtered DataFrame to a CSV file
