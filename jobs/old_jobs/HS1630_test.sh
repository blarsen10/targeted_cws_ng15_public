#!/bin/bash
#SBATCH --job-name=HS1630_20240211_test
#SBATCH --partition=day
#SBATCH --output=/home/bbl29/targeted_cws_ng15/logs/HS1630_20240211.txt
#SBATCH --error=/home/bbl29/targeted_cws_ng15/logs/error/HS1630_20240211.txt
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1:00:00
#SBATCH --array=0-1
#SBATCH --mail-type=END
#SBATCH --mail-user=bjorn.larsen@yale.edu

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load miniconda
module load OpenMPI
conda activate PTA_env

source_name=HS1630
noisedict=15yr_wn_dict
psrdists=pulsar_distances_15yr
path=/home/bbl29/targeted_cws_ng15
human=blarsen
ed_name=HS1630_init_emp_dist
Niter=250000
T=5

srun -n $SLURM_NTASKS python3 ./scripts/run_targeted_search.py --source_name $source_name --noisedict $noisedict --psrdists $psrdists --project_path $path --human $human --emp_dist_name $ed_name -N $Niter -T $T --write_hot_chains
