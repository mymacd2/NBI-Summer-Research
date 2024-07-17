#!/bin/bash
#SBATCH --job-name=2nu_decay_calcs
#SBATCH --output=2nuNO/outputs/array_output_%A_%a.txt
#SBATCH --error=2nuNO/errors/array_error_%A_%a.txt
#SBATCH --ntasks=1
#SBATCH --partition=astro2_short
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=11:59:59
#SBATCH --array=1-189

module load julia

# Read the tuple for this task
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" input_tuples_2nu.txt)
channel=$(echo $line | cut -d' ' -f1)
prog=$(echo $line | cut -d' ' -f2)
a=$(echo $line | cut -d' ' -f3)

# Run the Julia script with the (x, y) pair
julia-1.10.4/bin/julia dsnb_2nu_cluster_calcs.jl --channel $channel --prog $prog --alpha $a