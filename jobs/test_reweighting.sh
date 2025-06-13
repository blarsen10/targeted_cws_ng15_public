#!/bin/bash
#SBATCH --job-name=test_reweighting_1cpu
#SBATCH --partition=pi_mingarelli
#SBATCH --exclude=r806u23n04
#SBATCH --output=/home/bbl29/targeted_cws_ng15/logs/test_reweighting_1cpu.txt
#SBATCH --error=/home/bbl29/targeted_cws_ng15/logs/error/test_reweighting_1cpu.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=30:00
#SBATCH --array=0-7
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
psrdists=pulsar_distances_15yr
path=/home/bbl29/targeted_cws_ng15
coredir=/home/bbl29/project/targeted_cws_ng15/data/chains/ng15_sim_kneeGWB_NGC3115/no_cw_varycrn_BPL
Niter=1000

python3 ./scripts/run_reweighting.py --coredir $coredir --source_name $source_name --dataset $dataset --no_cw --vary_crn --bpl --noisedict $noisedict --psrdists $psrdists --project_path $path -N $Niter --restart
