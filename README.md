# Hi-C_draft

## Step 1
The following scripts are applicable for all the samples, mock community and natural community. These scripts are command lines as prescribed by the tool and literature. 

1. Clean and filter the raw reads (```1_clean_raw_sequence.sh```)
2. Align the Hi-C reads to a reference file (```2_align_reference.sh```)
3. Normalization of the Hi-C contacts using MetaCC (```3_normalize_metacc.sh```)

The directory is in OSC:
Mock Community directory: ``` /fs/ess/PAS1117/HiC/6thRoundTests/HiC ```

Natural Community directory: ``` /fs/ess/PAS1117/HiC/7thRoundTests/HiC ```

## Step 2
The output of metaCC is in the Sparse Matrix format. The extraction was done using the ```4_extract_matrix.py``` script. This is done in the output directory of MetaCC ususally named as ```${filename}_out```. 
All the python scripts here are ran on python3. 

Because there are multiple contigs for each genome, the scores of multiple contigs were summed together for each genome.

For the mock community, ```5_mock_adding_chromosome.r``` was used locally on the raw ```phage_host_nomalized_contact.csv```. 
For the natural community, ```6_sum_scores_by_mag_mire.py``` was used on the output folder. 

## Step 3 

For further analysis, the rest of the scripts are generated for merging tables, Z-score calculation and figure analysis. 

For the mock community, ```11_mock_figure.r``` was used to generate the figures locally.

For the natural community the following were done:
1. generate z-score (```7_z_score_mire.py```)
2. merge GTDB prediction based on MAG (```8_merge_gtdb_hic.py```)
3. merge iPHoP results and Hi-C results based on viruses (```9_merge_iphop_hic.py```)
4. Generate figure locally (```10_hic_iphop_analysis.r```)
