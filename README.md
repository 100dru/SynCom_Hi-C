# Hi-C_draft

For our analysis, we follow an alignment based approach; where we align the Hi-C reads to a reference genome and count aligned reads that have chimeras consisting of viruses and hosts.

## The Hi-C Analysis Pipeline
---

### Step 1: Clean Hi-C reads

We use the BBduk script from the BBmap suite to remove the adapters and the phix from the raw reads.
Then we use the same tool to trimming and filtering the reads.

``` bash 1_clean_raw_sequence.sh ```

### Step 2: Align the  Hi-C reads

We use BWA mem to align the Hi-C reads to the reference genomes NewHiCRef.fasta. Then we use samtools to convert the aligned sam file to a bam file and sort it. 

```bash 2_align_reference.sh ```

### Step 5: Run MetaCC

Used the mapping file from step 2, the coverage file from step 4 and the reference genome as the input for HiCZin.

```bash 3_normalize_metacc.sh ```

### Analysis of the HiCZin output

- **contig_info.csv**: information of assembled contigs with three columns (contig name, the number of restriction sites on contigs, and contig length).
- **Normalized_contact_matrix.npz**: a sparse matrix of normalized Hi-C contact maps in csr format and can be reloaded using Python command *'scipy.sparse.load_npz('Normalized_contact_matrix.npz')'*.
- **NormCC_normalized_contact.gz**: Compressed format of the normalized contacts and contig information by pickle. This file can further serve as the input of MetaCC binning module.
- **MetaCC.log**: the specific implementation information of NormCC normalization module.

To read the Normalized_contact_matrix.npz file and have it in  a readable format, I wrote a python script that loads the sparse matrix and adds the contig names from the contig_info.csv. It would produce two different tables, one with ALL the contacts and one with only Virus-Host contacts for further analysis. 

``` python3 4_extract_matrix.py ```

## Further Analysis with Mock Community

For Further analysis, the mock community Hi-C linkages were joined by each bacterial strain. Their sensitivity and Specificity was calculated and figures were generated.

``` r 6_Mock_community_analysis _figure.R ```


## Further Analysis with Natural Community

For further analysis, the rest of the scripts are generated for merging tables, Z-score calculation and figure analysis. 

The Hi-C Linkages were aggregated by MAGs and then they were joined to their GTDB prediction by MAG and the iPHoP Prediction by virus.

```python3 5_parsing_natural_community_HiC_iPHoP.py ```

The similar was done for VirMatcher Prediction.

```python3 5_parsing_natural_community_HiC_VirMatcher.py ```


