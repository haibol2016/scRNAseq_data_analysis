#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=6:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=1   # 36 processor core(s) per node 
#SBATCH --mem=16G   # maximum memory per node
#SBATCH --job-name="repair"
#SBATCH --array=13-16
#SBATCH --output="logs/out.rcorrect.%A_%a.txt" # job standard output file (%j replaced by job id)
#SBATCH --error="logs/err.rcorrect.%A_%a.txt" # job standard error file (%j replaced by job id)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load jdk/8u141-b15-xf726oe

i=$(($SLURM_ARRAY_TASK_ID - 1))

R1=(`ls Rcorrected.fastq/2.Ecor.trimAT.fastq/*_R1_001.fastq.gz`)
R2=(`ls Rcorrected.fastq/2.Ecor.trimAT.fastq/*_R2_001.fastq.gz`)
name=(`ls Rcorrected.fastq/2.Ecor.trimAT.fastq/*_R2_001.fastq.gz | perl -p -e 's{.+/(.+?)_R2_001.fastq.gz}{$1}'`)
out=repaired.fastq
mkdir -p ${out}


## using BBmap to re-pair the trimmed reads
repair.sh in1=${R1[$i]} in2=${R2[$i]} \
          out1=${out}/${R1[$i]}  out2=${out}/${R2[$i]} \
          outs=${out}/${name[$i]}.SE.fastq.gz  \
          repair   ziplevel=2 -Xmx16G 
