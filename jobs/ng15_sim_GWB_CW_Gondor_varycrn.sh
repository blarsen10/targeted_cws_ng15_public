#!/bin/bash
#SBATCH --job-name=ng15_sim_GWB_CW_Gondor_varycrn_20250606
#SBATCH --partition=day
#SBATCH --exclude=r806u23n04
#SBATCH --output=/home/bbl29/targeted_cws_ng15/logs/ng15_sim_GWB_CW_Gondor_varycrn_20250606.txt
#SBATCH --error=/home/bbl29/targeted_cws_ng15/logs/error/ng15_sim_GWB_CW_Gondor_varycrn_20250606.txt
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-
#SBATCH --requeue
#SBATCH --array=0-15
#SBATCH --mail-type=END
#SBATCH --mail-user=bjorn.larsen@yale.edu

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load miniconda
module load OpenMPI
conda activate PTA_env

source_name=SDSS_J072908.71+400836.6
gw_components=14
dataset=ng15_sim_GWB_Rondor
noisedict=15yr_wn_dict
orf=
log10_fgw_prior=uniform
psrdists=pulsar_distances_15yr
path=/home/bbl29/targeted_cws_ng15
pkldir=/gpfs/gibbs/project/mingarelli/bbl29/targeted_cws_ng15
outdir=/vast/palmer/scratch/mingarelli/bbl29/targeted_cws_ng15
human=blarsen
ed_name=
Niter=750000
T=10

srun -n $SLURM_NTASKS python3 ./scripts/run_targeted_search.py --source_name $source_name --dataset $dataset --vary_crn --orf $orf --gw_components $gw_components --detection --log10_fgw_prior $log10_fgw_prior --noisedict $noisedict --psrdists $psrdists --project_path $path --outdir_path $outdir --pkldir_path $pkldir --human $human --emp_dist_name $ed_name -N $Niter -T $T
