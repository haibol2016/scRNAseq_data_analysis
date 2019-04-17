#!/bin/bash
# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # 36 processor core(s) per node 
#SBATCH --mem=128G   # maximum memory per node
#SBATCH --job-name="cellranger"
#SBATCH --array=4
#SBATCH --output="logs/out.cellranger.%A_%a.txt" # job standard output file (%j replaced by job id)
#SBATCH --error="logs/err.cellranger.%A_%a.txt" # job standard error file (%j replaced by job id)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load jdk/8u141-b15-xf726oe

i=$(($SLURM_ARRAY_TASK_ID - 1))

ids=(A B C CT2-1NOV  CT2-30OCT)
dir=(A B C CT2-1NOV  CT2-30OCT)
fastq_in=repaired.fastq

cellranger count --id=${ids[$i]} \
                   --transcriptome=ssc3 \
                   --fastqs=${fastq_in}/${dir[$i]}\
                   --sample=${dir[$i]}-1,${dir[$i]}-2,${dir[$i]}-3,${dir[$i]}-4 \
                   --chemistry=SC3Pv2 \
                   --localcores=36  \
                   --localmem=252 
