#!/bin/bash


# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=4:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # 36 processor core(s) per node 
#SBATCH --mem=8G   # maximum memory per node
#SBATCH --job-name="multiqc"
#SBATCH --array=1-9
#SBATCH --output="logs/out.addRG.%A_%a.txt" # job standard output file (%j replaced by job id)
#SBATCH --error="logs/err.addRG.%A_%a.txt" # job standard error file (%j replaced by job id)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE


i=$(($SLURM_ARRAY_TASK_ID - 1))

gtf=Sus_scrofa.Sscrofa11.1.92.gtf
out=Sus_scrofa.genome
mkdir -p $out

cellranger mkgtf  Sus_scrofa.Sscrofa11.1.92.gtf   Sus_scrofa.Sscrofa11.1.92.gtf.filtered.gtf \
                 --attribute=gene_biotype:IG_C_gene \
                 --attribute=gene_biotype:IG_V_gene \
                 --attribute=gene_biotype:lincRNA \
                 --attribute=gene_biotype:miRNA \
                 --attribute=gene_biotype:misc_RNA \
                 --attribute=gene_biotype:Mt_rRNA \
                 --attribute=gene_biotype:Mt_tRNA \
                 --attribute=gene_biotype:protein_coding \
                 --attribute=gene_biotype:pseudogene \
                 --attribute=gene_biotype:ribozyme \
                 --attribute=gene_biotype:rRNA \
                 --attribute=gene_biotype:scaRNA \
                 --attribute=gene_biotype:snoRNA \
                 --attribute=gene_biotype:snRNA \
                 --attribute=gene_biotype:sRNA \
                 --attribute=gene_biotype:TR_J_gene \
                 --attribute=gene_biotype:TR_V_gene
