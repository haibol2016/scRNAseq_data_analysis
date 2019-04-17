#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=6:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # 36 processor core(s) per node 
#SBATCH --mem=96G   # maximum memory per node
#SBATCH --job-name="trim_AT"
#SBATCH --array=1-4
#SBATCH --output="logs/out.rcorrect.%A_%a.txt" # job standard output file (%j replaced by job id)
#SBATCH --error="logs/err.rcorrect.%A_%a.txt" # job standard error file (%j replaced by job id)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

i=$(($SLURM_ARRAY_TASK_ID - 1))

R2=(`ls  Rcorrected.fastq/*_R2_001.cor.fq.gz` )
name=(`ls  Rcorrected.fastq/*_R2_001.cor.fq.gz | perl -p -e 's/cor.fq.gz/fastq.gz/'`)

zcat ${R2[$i]} | ./trimming.AT.pl RNA-seq | gzip > ${name[$i]}
