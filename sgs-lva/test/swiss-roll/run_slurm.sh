#!/bin/bash
###sbatch --constraint=node-2017 run_slurm.sh
###scontrol show jobid -dd 4429
###scancel 4432

#SBATCH --chdir=/scratch/nas/4/operedo/gslib-lva/sgs-lva/test/swiss-roll/
##SBATCH --exclusive
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=5GB
#SBATCH --time 3:00:00

source /Soft/intelf/composer_xe_2013.3.163/bin/compilervars.sh intel64

#cd ../../src
#make clean;make levels
#cp SGS_LVA_levels ../test/swiss-roll/
#cd -

#export KMP_AFFINITY=granularity=fine,scatter
export OMP_NUM_THREADS=20
/usr/bin/time ./SGS_LVA_levels sgs_lva.par_v6_360_10_1 > err27.txt 2>&1
#/usr/bin/time ./SGS_LVA sgs_lva.par_v6_360_10_1 > err28.txt 2>&1


