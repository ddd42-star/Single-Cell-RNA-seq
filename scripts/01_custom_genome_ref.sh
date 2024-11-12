#!/usr/bin/env bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=pibu_el8
#SBATCH --time=06:00:00
#SBATCH --output=/data/users/dbassi/internship/outputs/01_custom_genome_ref_%j.o


REFERENCE=/data/users/dbassi/internship/data/genome_ref
THREADS=$SLURM_CPUS_PER_TASK
# load the module cellranger
module load CellRanger/7.1.0

# filter gtf files using only protein_coding
cellranger mkgtf $REFERENCE/pdm09H1N1.gtf $REFERENCE/pdm09H1N1_filtered.gtf --attribute=gene_biotype:protein_coding

cellranger mkgtf $REFERENCE/Sus_scrofa.Sscrofa11.1.112.gtf $REFERENCE/Sus_scrofa.Sscrofa11.1.112_filtered.gtf --attribute=gene_biotype:protein_coding


# make custom genome reference
cellranger mkref --genome=Sus_Sscrofa --fasta=$REFERENCE/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa --genes=$REFERENCE/Sus_scrofa.Sscrofa11.1.112_filtered.gtf \
    --genome=h1n1 --fasta=$REFERENCE/pdm09H1N1.fasta --genes=$REFERENCE/pdm09H1N1_filtered.gtf --nthreads=$THREADS




