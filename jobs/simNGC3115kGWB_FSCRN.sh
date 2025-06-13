#!/bin/bash
#SBATCH --job-name=simNGC3115kGWB_FSCRN_20240731
#SBATCH --partition=scavenge
#SBATCH --exclude=r806u23n04
#SBATCH --output=/home/bbl29/targeted_cws_ng15/logs/simNGC3115kGWB_FSCRN_20240731.txt
#SBATCH --error=/home/bbl29/targeted_cws_ng15/logs/error/simNGC3115kGWB_FSCRN_20240731.txt
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-
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
dataset=ng15_sim_kneeGWB_NGC3115
noisedict=sim_NGC3115_kneeGWB
orf=
log10_fgw_prior=uniform
psrdists=pulsar_distances_15yr
path=/home/bbl29/targeted_cws_ng15
outdir=/vast/palmer/scratch/mingarelli/bbl29/targeted_cws_ng15
human=blarsen
ed_name=15yr_emp_distr
Niter=500000
T=10

srun -n $SLURM_NTASKS python3 ./scripts/run_targeted_search.py --source_name $source_name --dataset $dataset --no_cw --vary_crn --freespec --orf $orf --noisedict $noisedict --psrdists $psrdists --project_path $path --outdir_path $outdir --human $human --emp_dist_name $ed_name -N $Niter -T $T --write_hot_chains
