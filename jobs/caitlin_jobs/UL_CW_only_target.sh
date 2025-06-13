#!/bin/sh

## Array
##PBS -t 0-4

### Specify job name
#PBS -N 12p5_UL_CW_only_target


# Specify the resources need for the job
# Walltime is specified as hh:mm:ss (hours:minutes:seconds)
#PBS -l nodes=1:ppn=1,pvmem=10gb

#email preferences
#PBS -m ae
#PBS -M caw0057@mix.wvu.edu

#which queue
#PBS -q sbs0016

cd /scratch/caw0057/nano12p5_gwb/cw_12p5/target_scripts_e_rn_emp_nd/
source /users/caw0057/.bashrc
conda activate enterprise_dev
python 12p5_UL_CW_only_target_check.py -N_iter 2e6
