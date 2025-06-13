#!/bin/bash
#SBATCH --job-name=sim3_HS_1630+2355_UL_varycrn_20240302
#SBATCH --partition=pi_mingarelli
#SBATCH --output=/home/bbl29/targeted_cws_ng15/logs/sim3_HS_1630+2355_UL_varycrn_20240302.txt
#SBATCH --error=/home/bbl29/targeted_cws_ng15/logs/error/sim3_HS_1630+2355_UL_varycrn_20240302.txt
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

source_name=HS_1630+2355_UL
dataset=ng15_sim_GWB_HS1630
noisedict=15yr_wn_dict
psrdists=pulsar_distances_15yr
path=/home/bbl29/targeted_cws_ng15
outdir=/home/bbl29/project/targeted_cws_ng15
human=blarsen
ed_name=HS1630_init_emp_dist
Niter=500000
T=10

srun -n $SLURM_NTASKS python3 ./scripts/run_targeted_search.py --source_name $source_name --dataset $dataset --upper_limit --vary_crn --noisedict $noisedict --psrdists $psrdists --project_path $path --outdir_path $outdir --human $human --emp_dist_name $ed_name -N $Niter -T $T
