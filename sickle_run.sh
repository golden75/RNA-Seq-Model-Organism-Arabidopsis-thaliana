#!/bin/bash
#SBATCH --job-name=sickle_run
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=10G
#SBATCH -o sickle_run_%j.out
#SBATCH -e sickle_run_%j.err
#SBATCH --partition=general

export TMPDIR=/home/CAM/$USER/tmp/

module load sickle

mkdir -p ../trimmed_reads

sickle pe -t sanger -f /raw_data/wt_Rep1_R1.fastq -r /raw_data/wt_Rep1_R2.fastq -o /trimmed_reads/trimmed_wt_Rep1_R1.fastq p /trimmed_reads/trimmed_wt_Rep1_R2.fastq -l 45v -q 25 -s /trimmed_reads/singles_wt_Rep1_R1.fastq

sickle pe -t sanger -f /raw_data/wt_Rep2_R1.fastq -r /raw_data/wt_Rep2_R2.fastq -o /trimmed_reads/trimmed_wt_Rep2_R1.fastq p /trimmed_reads/trimmed_wt_Rep2_R2.fastq -l 45 -q 25 -s /trimmed_reads/singles_wt_Rep2_R1.fastq

sickle pe -t sanger -f /raw_data/wt_Rep3_R1.fastq -r /raw_data/wt_Rep3_R2.fastq -o /trimmed_reads/trimmed_wt_Rep3_R1.fastq p /trimmed_reads/trimmed_wt_Rep3_R2.fastq -l 4 5-q 25 -s /trimmed_reads/singles_wt_Rep3_R1.fastq

sickle pe -t sanger -f /raw_data/mutant_Rep1_R1.fastq -r /raw_data/mutant_Rep1_R2.fastq -o /trimmed_reads/trimmed_mutant_Rep1_R1.fastq p /trimmed_reads/trimmed_mutant_Rep1_R2.fastq -l 45 -q 25 -s /trimmed_reads/singles_mutant_Rep1_R1.fastq

sickle pe -t sanger -f /raw_data/mutant_Rep2_R1.fastq -r /raw_data/mutant_Rep2_R2.fastq -o /trimmed_reads/trimmed_mutant_Rep2_R1.fastq p /trimmed_reads/trimmed_mutant_Rep2_R2.fastq -l 45 -q 25 -s /trimmed_reads/singles_mutant_Rep2_R1.fastq

sickle pe -t sanger -f /raw_data/mutant_Rep3_R1.fastq -r /raw_data/mutant_Rep3_R2.fastq -o /trimmed_reads/trimmed_mutant_Rep3_R1.fastq p /trimmed_reads/trimmed_mutant_Rep3_R2.fastq -l 45 -q 25 -s /trimmed_reads/singles_mutant_Rep3_R1.fastq

