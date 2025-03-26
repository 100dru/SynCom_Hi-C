import numpy as np
import pandas as pd
import scipy.sparse as sp


matrix = sp.load_npz('Normalized_contact_matrix.npz')
sp_matrix_df = pd.DataFrame.sparse.from_spmatrix(matrix)

contig_data = pd.read_csv('contig_info.csv', header=0, index_col=False)

csr_mat = sp_matrix_df.sparse.to_coo().tocsr()
print((csr_mat.data.nbytes + csr_mat.indptr.nbytes + csr_mat.indices.nbytes) / (1024 * 1024))

row_idx, col_idx = csr_mat.nonzero() 
scores = csr_mat.data

name1 = contig_data["Contig name"].iloc[row_idx].values
name2 = contig_data["Contig name"].iloc[col_idx].values

df_longer = pd.DataFrame({'contig1_name': name1, 'contig2_name': name2, 'Score': scores})

df_longer.to_csv('all_nomalized_contact.csv', index=False)


filtered_df = df_longer[df_longer['contig1_name'].str.startswith("MAG") & df_longer['contig2_name'].str.startswith("Virus")]


filtered_df.to_csv('phage_host_nomalized_contact.csv', index=False)
