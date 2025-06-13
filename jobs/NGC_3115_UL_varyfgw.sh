#!/bin/bash
#SBATCH --job-name=NGC_3115_UL_varyfgw_20240412
#SBATCH --partition=pi_mingarelli
#SBATCH --exclude=r806u23n04
#SBATCH --output=/home/bbl29/targeted_cws_ng15/logs/NGC_3115_UL_varyfgw_20240412.txt
#SBATCH --error=/home/bbl29/targeted_cws_ng15/logs/error/NGC_3115_UL_varyfgw_20240412.txt
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=2-
#SBATCH --array=0-9
#SBATCH --mail-type=END
#SBATCH --mail-user=bjorn.larsen@yale.edu

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load miniconda
module load OpenMPI
conda activate PTA_env

source_name=NGC_3115
dataset=ng15_v1p1
noisedict=15yr_wn_dict
ra_dec=
orf=
psrdists=pulsar_distances_15yr
path=/home/bbl29/targeted_cws_ng15
outdir=/home/bbl29/project/targeted_cws_ng15
human=blarsen
ed_name=15yr_emp_distr
Niter=500000
T=10

srun -n $SLURM_NTASKS python3 ./scripts/run_targeted_search.py --source_name $source_name --dataset $dataset --upper_limit --vary_fgw --ra_dec $ra_dec --orf $orf --noisedict $noisedict --psrdists $psrdists --project_path $path --outdir_path $outdir --human $human --emp_dist_name $ed_name -N $Niter -T $T --write_hot_chains
