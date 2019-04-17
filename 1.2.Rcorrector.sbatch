#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=72:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 36 processor core(s) per node 
#SBATCH --mem=96G   # maximum memory per node
#SBATCH --job-name="rcorrector"
#SBATCH --array=4
#SBATCH --output="logs/out.rcorrect.%A_%a.txt" # job standard output file (%j replaced by job id)
#SBATCH --error="logs/err.rcorrect.%A_%a.txt" # job standard error file (%j replaced by job id)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

i=$(($SLURM_ARRAY_TASK_ID - 1))
module load java/1.8.0_171
module load zlib/1.2.11-vhzh5cf
mkdir -p logs


libs=(A B C CT2-1NOV  CT2-30OCT)

mkdir -p Rcorrected.fastq

R2=(`ls  -1 fastq/${libs[$i]}*_R2_001.fastq.gz |perl -p -e 's/\n/,/' |perl -p -e  's/,$//'` )
run_rcorrector.pl -s  ${R2[0]}  -k 25  -od ./Rcorrected.fastq  -t 36
