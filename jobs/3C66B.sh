#!/bin/bash
#SBATCH --job-name=3C66B_20240220
#SBATCH --partition=day
#SBATCH --output=/home/bbl29/targeted_cws_ng15/logs/3C66B_20240220.txt
#SBATCH --error=/home/bbl29/targeted_cws_ng15/logs/error/3C66B_20240220.txt
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-
#SBATCH --array=0-9
#SBATCH --mail-type=END
#SBATCH --mail-user=bjorn.larsen@yale.edu

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load miniconda
module load OpenMPI
conda activate PTA_env

source_name=3C66B
noisedict=15yr_wn_dict
psrdists=pulsar_distances_15yr
path=/home/bbl29/targeted_cws_ng15
outdir=/home/bbl29/project/targeted_cws_ng15
human=blarsen
ed_name=3C66B_emp_dist
Niter=500000
T=10

srun -n $SLURM_NTASKS python3 ./scripts/run_targeted_search.py --source_name $source_name --noisedict $noisedict --psrdists $psrdists --project_path $path --outdir_path $outdir --human $human --emp_dist_name $ed_name -N $Niter -T $T
