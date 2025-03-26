#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --time=08:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@example.com # Replace with your email

## Define variables
REFERENCE="NewHiCRef.fasta"
READS_DIR="/path/to/clean_reads" # Replace with your reads directory
OUTPUT_DIR="/path/to/output"     # Replace with your output directory

# Load necessary modules
module load bwa/0.7.17-r1198
module load samtools/1.10

# Change to the working directory
cd $OUTPUT_DIR

# Index the reference genome
bwa index $REFERENCE

# Align the reads to the reference genome
bwa mem -5SP -T0 -t16 $REFERENCE ${READS_DIR}/${filename}_R1_clean.fastq.gz ${READS_DIR}/${filename}_R2_clean.fastq.gz > ${filename}_aligned.sam

# Convert SAM to BAM and sort
samtools view -F 0x904 -bS ${filename}_aligned.sam > ${filename}_MAP_UNSORTED.bam
samtools sort -n ${filename}_MAP_UNSORTED.bam -o ${filename}_MAP_SORTED.bam

# Optional: Remove intermediate files to save space
# rm ${filename}_aligned.sam ${filename}_MAP_UNSORTED.bam
