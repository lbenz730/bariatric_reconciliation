#!/usr/bin/bash

#SBATCH -c 4 ## number of cores
#SBATCH -t 0-00:20 ## amount of time in D-HH:MM
#SBATCH -p serial_requeue ## Partition to submit to
#SBATCH --mem=170000 ## memory pool for all cores
#SBATCH -o logs/build/log.stdout_%a ## STDOUT
#SBATCH -e logs/build/log.stderr_%a ## STDERR
#SBATCH --account=haneuse_lab
#SBATCH --array=1-96

module load R/4.2.2-fasrc01
export R_LIBS_USER=$HOME/apps/R_4.2.2:$R_LIBS_USER

cd $HOME/bariatric_tte

Rscript scripts/data/build_bariatric_CVD_trial_datasets.R  $SLURM_ARRAY_TASK_ID
