#!/usr/bin/env bash

#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=pibu_el8
#SBATCH --time=24:00:00
#SBATCH --output=/data/users/dbassi/internship/outputs/02_count_WT_pAEC_WT_37_%j.o
#SBATCH --mail-user=dario.bassi@students.unibe.ch
#SBATCH --mail-type=begin,end,fail
#SBATCH --job-name=cellranger_count


SAMPLES=/data/users/dbassi/internship/data
REFERENCE=/data/users/dbassi/internship/results
TRANSCRIPTOME=/data/users/dbassi/internship/results/Sus_Sscrofa_and_h1n1
THREADS=$SLURM_CPUS_PER_TASK
# load the module cellranger
module load CellRanger/7.1.0

cellranger count --id=WT_pAEC_WT_37_count \
        --transcriptome=$TRANSCRIPTOME \
        --nosecondary \
        --mempercore=$THREADS \
        --disable-ui \
        --fastqs=$SAMPLES/WT_pAEC_WT_37 
        


