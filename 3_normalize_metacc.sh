#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --account=YOUR_ACCOUNT_NAME
#SBATCH --job-name=metacc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=YOUR_EMAIL@example.com

# Load necessary modules (adjust the module path as needed)
module use /path/to/your/modulefiles
module load MetaCC

# Change to the working directory (adjust the path as needed)
cd /path/to/your/working/directory

# Normalize MetaCC data
# Replace 'Sau3AI', 'MluCI', 'NewHiCRef.fasta', and paths with your actual enzyme names, reference file, and output paths
MetaCC norm -e Sau3AI -e MluCI NewHiCRef.fasta ${filename}_MAP_SORTED.bam /path/to/your/output/directory/${filename}_out
