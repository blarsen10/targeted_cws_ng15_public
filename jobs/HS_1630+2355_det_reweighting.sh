#!/bin/bash
#SBATCH --job-name=HS_1630+2355_det_HD_20240923
#SBATCH --partition=pi_mingarelli
#SBATCH --exclude=r806u23n04
#SBATCH --output=/home/bbl29/targeted_cws_ng15/logs/HS_1630+2355_det_HD_20240923.txt
#SBATCH --error=/home/bbl29/targeted_cws_ng15/logs/error/HS_1630+2355_det_HD_20240923.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=2-
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

source_name=HS_1630+2355
dataset=ng15_v1p1
noisedict=15yr_wn_dict
psrdists=pulsar_distances_15yr
path=/home/bbl29/targeted_cws_ng15
coredir=/gpfs/gibbs/project/mingarelli/bbl29/targeted_cws_ng15/data/chains/ng15_v1p1/HS_1630+2355_det
Niter=100000

python3 ./scripts/run_reweighting.py --coredir $coredir --source_name $source_name --dataset $dataset --detection --noisedict $noisedict --psrdists $psrdists --project_path $path -N $Niter
