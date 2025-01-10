#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --account=PAS1117
#SBATCH --job-name=metacc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=shatadru.2@osu.edu

module use /fs/ess/PAS1117/modulefiles
module load MetaCC

cd /fs/ess/PAS1117/HiC/6thRoundTests/HiC


## here filename=${i}. the terminal code is:

#for i in `ls -1 *_MAP_SORTED.bam | sed 's/_MAP_SORTED.bam//'`
#do sbatch --export=filename=${i} bwa_mock.sh
#done
MetaCC norm -e Sau3AI -e MluCI NewHiCRef.fasta ${filename}_MAP_SORTED.bam /fs/ess/PAS1117/HiC/6thRoundTests/HiC/metacc/${filename}_out
