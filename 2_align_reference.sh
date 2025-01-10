#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH --account=PAS1117
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --time=08:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shatadru.2@osu.edu

### here filename=${i}; for i in `ls -1 *_R1_clean.fastq.gz | sed 's/_R1_clean.fastq.gz//'`

module use /fs/project/PAS1117/modulefiles
module load bwa/0.7.17-r1198
module load samtools/1.10

cd /fs/ess/PAS1117/HiC/6thRoundTests/HiC 

bwa index NewHiCRef.fasta

bwa mem -5SP -T0 -t16 NewHiCRef.fasta /fs/ess/PAS1117/HiC/6thRoundTests/HiC/clean_reads/${filename}_R1_clean.fastq.gz /fs/ess/PAS1117/HiC/6thRoundTests/HiC/clean_reads/${filename}_R2_clean.fastq.gz -o ${filename}_aligned.sam


samtools view -F 0x904 -bS ${filename}_aligned.sam > ${filename}_MAP_UNSORTED.bam
samtools sort -n ${filename}_MAP_UNSORTED.bam -o ${filename}_MAP_SORTED.bam

