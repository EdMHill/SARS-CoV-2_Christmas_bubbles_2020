#!/bin/bash
#SBATCH --job-name=household_bubble_IBM
#SBATCH --time=010:00:00
#SBATCH --mem=2g
#SBATCH --export=ALL
#SBATCH --partition=ntd

# Submit a job array
#SBATCH --array=1-20
#SBATCH --ntasks=1

# ${SLURM_SUBMIT_DIR} points to the path where this script was submitted from
cd ${SLURM_SUBMIT_DIR}

module load julia/1.5

# ARGS list
# ARGS[1] job_ID
# ARGS[2] RNG_seed
# ARGS[3] n_simns
# ARGS[4] n_households
# ARGS[5] time_duration_from_christmas_bubble_start
julia run_stochastic_householdIBM.jl ${SLURM_ARRAY_TASK_ID} 1234 100 100000 15

echo "Finished"
