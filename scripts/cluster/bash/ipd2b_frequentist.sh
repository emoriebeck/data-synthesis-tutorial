#!/bin/bash
#SBATCH --account=p31385  ## YOUR ACCOUNT p31385 or bXXXX
#SBATCH --partition=normal  ### PARTITION (buyin, short, normal, w10001, etc)
#SBATCH --array=0-26 ## number of jobs to run "in parallel" 
#SBATCH --nodes=1 ## how many computers do you need
#SBATCH --ntasks-per-node=4 ## how many cpus or processors do you need on each computer
#SBATCH --time=24:00:00 ## how long does this need to run (remember different partitions have restrictions on this param)
#SBATCH --mem-per-cpu=4G ## how much RAM do you need per CPU (this effects your FairShare score so be careful to not ask for more than you need))
#SBATCH --job-name="ipd2b_frequentist_\${SLURM_ARRAY_TASK_ID}" ## use the task id in the name of the job
#SBATCH --output=ipd2b_frequentist.%A_%a.out ## use the jobid (A) and the specific job index (a) to name your log file
#SBATCH --mail-type=all ## you can receive e-mail alerts from SLURM when your job begins and when your job finishes (completed, failed, etc)
#SBATCH --mail-user=emorie_beck@northwestern.edu


# Load the environmental variables necessary for running R
module purge all
module load R/4.0.3

## Run R script based on the array number. 
Rscript ipd2b.R
