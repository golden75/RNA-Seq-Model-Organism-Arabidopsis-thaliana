#!/bin/bash
#SBATCH --job-name=sam_sort_bam
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=40G
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --partition=general

export TMPDIR=/home/CAM/$USER/tmp/

module load samtools

samtools view -@ 8 -bhS ../mapping/wt_Rep1.sam -o ../mapping/wt_Rep1.bam
samtools sort -@ 8 ../mapping/wt_Rep1.bam -o ../mapping/wt_Rep1_sort.bam

samtools view -@ 8 -bhS ../mapping/wt_Rep2.sam -o ../mapping/wt_Rep2.bam
samtools sort -@ 8 ../mapping/wt_Rep2.bam -o ../mapping/wt_Rep2_sort.bam

samtools view -@ 8 -bhS ../mapping/wt_Rep3.sam -o ../mapping/wt_Rep3.bam
samtools sort -@ 8 ../mapping/wt_Rep3.bam -o ../mapping/wt_Rep3_sort.bam

samtools view -@ 8 -bhS ../mapping/mutant_Rep1.sam -o ../mapping/mutant_Rep1.bam
samtools sort -@ 8 ../mapping/mutant_Rep1.bam -o ../mapping/mutant_Rep1_sort.bam

samtools view -@ 8 -bhS ../mapping/mutant_Rep2.sam -o ../mapping/mutant_Rep2.bam
samtools sort -@ 8 ../mapping/mutant_Rep2.bam -o ../mapping/mutant_Rep2_sort.bam

samtools view -@ 8 -bhS ../mapping/mutant_Rep3.sam -o ../mapping/mutant_Rep3.bam
samtools sort -@ 8 ../mapping/mutant_Rep3.bam -o ../mapping/mutant_Rep3_sort.bam

