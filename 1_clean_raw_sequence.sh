#!/bin/bash
#SBATCH -t 8:00:00
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -A <account_name>
#SBATCH -J cleaning_hic
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<your_email>

## Here filename=${i}. The terminal code is:

#for i in `ls -1 *_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//'`
#do sbatch --export=filename=${i} clean.sh
#done

bbduk_path="<path_to_bbduk>"
bbduk_ref_path="<path_to_reference_files>"

$bbduk_path/bbduk.sh -Xmx48G threads=12 overwrite=t in=${filename}_R1_001.fastq.gz in2=${filename}_R2_001.fastq.gz ref=$bbduk_ref_path/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=10 minlen=20 out=${filename}_R1_clean.fastq.gz out2=${filename}_R2_clean.fastq.gz stats=${filename}_clean_stats.txt

$bbduk_path/bbduk.sh -Xmx48G threads=12 overwrite=t in=${filename}_R1_clean.fastq.gz in2=${filename}_R2_clean.fastq.gz ref=$bbduk_ref_path/phix174_ill.ref.fa k=31 hdist=1 stats=${filename}_phix_stats.txt out=${filename}_R1_no_adapter_phix.fastq.gz out2=${filename}_R2_no_adapter_phix.fastq.gz
