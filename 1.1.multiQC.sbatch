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

in_dir=(`ls results/fastqc.out/* -d  | uniq`)
project=(`ls results/fastqc.out/* -d | uniq | perl -p -e 's{.+/(.+)}{$1}'`)
out_dir=results/MultiQC.out/${project[$i]}

mkdir -p ${out_dir}

## using virtual environment
source ~/bin/system.py3.6.5_env/bin/activate

multiqc --filename  ${project[$i]}.mapping.multiQC  --outdir  ${out_dir}   ${in_dir[$i]}
