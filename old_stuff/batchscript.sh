#!/bin/bash
#SBATCH --job-name=julia_array_job_first_pass
#SBATCH --output=array_output_%A_%a.txt
#SBATCH --error=array_error_%A_%a.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=70G
#SBATCH --time=00:15:00
#SBATCH --array=1-726

module load julia

# Read the tuple for this task
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" input_pairs.txt)
bar=$(echo $line | cut -d' ' -f1)
prog=$(echo $line | cut -d' ' -f2)
a1=$(echo $line | cut -d' ' -f3)
a2=$(echo $line | cut -d' ' -f4)

# Run the Julia script with the (x, y) pair
julia-1.10.4/bin/julia dsnb_3nu_IO_cluster_calcs_funcevals.jl --prog $prog --alpha1 $a1 --alpha2 $a2