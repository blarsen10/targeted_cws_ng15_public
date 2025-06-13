#!/bin/bash
#SBATCH --job-name=simnoise_NGC_3115_det_varyfgw_20240527
#SBATCH --partition=pi_mingarelli
#SBATCH --exclude=r806u23n04
#SBATCH --output=/home/bbl29/targeted_cws_ng15/logs/simnoise_NGC_3115_det_varyfgw_20240527.txt
#SBATCH --error=/home/bbl29/targeted_cws_ng15/logs/error/simnoise_NGC_3115_det_varyfgw_20240527.txt
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=2-
#SBATCH --array=0-15
#SBATCH --mail-type=END
#SBATCH --mail-user=bjorn.larsen@yale.edu

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load miniconda
module load OpenMPI
conda activate PTA_env

source_name=NGC_3115
dataset=ng15_sim_noiseonly
noisedict=sim_PTA_noisedict
orf=
psrdists=pulsar_distances_15yr
path=/home/bbl29/targeted_cws_ng15
outdir=/home/bbl29/project/targeted_cws_ng15
human=blarsen
ed_name=15yr_emp_distr
Niter=800000
T=10

srun -n $SLURM_NTASKS python3 ./scripts/run_targeted_search.py --source_name $source_name --dataset $dataset --detection --vary_fgw --no_gwb --orf $orf --noisedict $noisedict --psrdists $psrdists --project_path $path --outdir_path $outdir --human $human --emp_dist_name $ed_name -N $Niter -T $T --write_hot_chains
