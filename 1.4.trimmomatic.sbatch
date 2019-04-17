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

module load java/1.8.0_171

R1=(`ls fastq/*_R1_001.fastq.gz`)
R2=(`ls fastq/*_R2_001.fastq.gz`)

name=(`ls fastq/*_R1_001.fastq.gz | perl -p -e   's/_R1_001.fastq.gz//g' `)

java -jar ~/project/umw_mccb/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 8 \
          -phred33 ${R1[$i]} ${R2[$i]} \
          ${name[$i]}.R1Ptrim.fastq.gz  ${name[$i]}.R1Utrim.fastq.gz \
          ${name[$i]}.R2Ptrim.fastq.gz  ${name[$i]}.R2Utrim.fastq.gz  \
          ILLUMINACLIP:adapter.fa:2:30:10:1:true LEADING:0 TRAILING:3 \
          SLIDINGWINDOW:4:20 LEADING:0 TRAILING:3 MINLEN:25                 

