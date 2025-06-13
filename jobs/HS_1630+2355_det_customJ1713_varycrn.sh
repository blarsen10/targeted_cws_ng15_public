#!/bin/bash
#SBATCH --job-name=HS_1630+2355_det_customJ1713_varycrn_20240923
#SBATCH --partition=pi_mingarelli
#SBATCH --exclude=r806u23n04
#SBATCH --output=/home/bbl29/targeted_cws_ng15/logs/HS_1630+2355_det_customJ1713_varycrn_20240923.txt
#SBATCH --error=/home/bbl29/targeted_cws_ng15/logs/error/HS_1630+2355_det_customJ1713_varycrn_20240923.txt
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=2-
#SBATCH --requeue
#SBATCH --array=0-9
#SBATCH --mail-type=END
#SBATCH --mail-user=bjorn.larsen@yale.edu

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

module load miniconda
module load OpenMPI
conda activate PTA_env

source_name=HS_1630+2355
dataset=ng15_v1p1_customJ1713
noisedict=15yr_customJ1713_noisedict
orf=
psrdists=pulsar_distances_15yr
path=/home/bbl29/targeted_cws_ng15
outdir=/vast/palmer/scratch/mingarelli/bbl29/targeted_cws_ng15
human=blarsen
ed_name=15yr_custom6_emp_distr
Niter=800000
T=10

srun -n $SLURM_NTASKS python3 ./scripts/run_targeted_search.py --source_name $source_name --dataset $dataset --detection --vary_crn --orf $orf --noisedict $noisedict --psrdists $psrdists --project_path $path --outdir_path $outdir --human $human --emp_dist_name $ed_name -N $Niter -T $T
