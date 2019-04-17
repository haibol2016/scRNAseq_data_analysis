#!/bin/bash


# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=4:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=8   # 36 processor core(s) per node 
#SBATCH --mem=8G   # maximum memory per node
#SBATCH --job-name="multiqc"
#SBATCH --array=1-9
#SBATCH --output="logs/out.addRG.%A_%a.txt" # job standard output file (%j replaced by job id)
#SBATCH --error="logs/err.addRG.%A_%a.txt" # job standard error file (%j replaced by job id)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE


i=$(($SLURM_ARRAY_TASK_ID - 1))


cellranger mkref --genome=ssc3 \
                 --fasta=Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
                 --genes=Sus_scrofa.Sscrofa11.1.92.gtf.filtered.gtf \
                 --ref-version=11.1.92
