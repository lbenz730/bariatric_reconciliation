#!/usr/bin/bash

#SBATCH -c 2 ## number of cores
#SBATCH -t 0-1:00 ## amount of time in D-HH:MM
#SBATCH -p serial_requeue ## Partition to submit to
#SBATCH --mem=180000 ## memory pool for all cores
#SBATCH -o logs/bootstrap/log.stdout_%a ## STDOUT
#SBATCH -e logs/bootstrap/log.stderr_%a ## STDERR
#SBATCH --account=haneuse_lab
#SBATCH --array=1-1000

module load R/4.2.2-fasrc01
export R_LIBS_USER=$HOME/apps/R_4.2.2:$R_LIBS_USER

cd $HOME/bariatric_tte

Rscript scripts/analysis/bootstrap_CI.R  $SLURM_ARRAY_TASK_ID