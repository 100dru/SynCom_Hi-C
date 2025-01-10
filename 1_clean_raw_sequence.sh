#!/bin/bash
#SBATCH -t 8:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -A PAS1117
#SBATCH -J cleaning_hic
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shatadru.2@osu.edu

## here filename=${i}; for i in `ls -1 *_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//'`

/fs/ess/PAS1117/bioinformatic_tools/bbmap_38.51/bbduk.sh -Xmx48G threads=12 overwrite=t in=${filename}_R1_001.fastq.gz in2=${filename}_R2_001.fastq.gz ref=/fs/ess/PAS1117/bioinformatic_tools/bbmap_38.51/resources/adapters.fa,/fs/ess/PAS1117/bioinformatic_tools/bbmap_38.51/resources/phix174_ill.ref.fa.gz ktrim=r minlength=30 k=23 mink=11 hdist=1 hdist2=1 out=/fs/ess/PAS1117/HiC/6thRoundTests/HiC/clean_reads/${filename}_R1_no_adapter_phix.fastq.gz out2=/fs/ess/PAS1117/HiC/6thRoundTests/HiC/clean_reads/${filename}_R2_no_adapter_phix.fastq.gz refstats=/fs/ess/PAS1117/HiC/6thRoundTests/HiC/clean_reads/${filename}_no_adapter_phix_Trimming.stats statscolumns=5

/fs/ess/PAS1117/bioinformatic_tools/bbmap_38.51/bbduk.sh -Xmx48G threads=12 overwrite=t in=/fs/ess/PAS1117/HiC/6thRoundTests/HiC/clean_reads/${filename}_R1_no_adapter_phix.fastq.gz in2=/fs/ess/PAS1117/HiC/6thRoundTests/HiC/clean_reads/${filename}_R2_no_adapter_phix.fastq.gz  qtrim=rl maq=20 maxns=0 minlength=30 trimq=20 out=/fs/ess/PAS1117/HiC/6thRoundTests/HiC/clean_reads/${filename}_R1_clean.fastq.gz  out2=/fs/ess/PAS1117/HiC/6thRoundTests/HiC/clean_reads/${filename}_R2_clean.fastq.gz outs=/fs/ess/PAS1117/HiC/6thRoundTests/HiC/clean_reads/${filename}_S_clean.fastq.gz
