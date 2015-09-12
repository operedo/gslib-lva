#!/bin/sh
#SBATCH --job-name="MPI-BOOST TEST"
#SBATCH --nodes=16
#SBATCH --exclusive
#SBATCH --mail-user=operedo@alges.cl
#SBATCH --mail-type=end

/bin/hostname

echo "SLURM_JOBID: $SLURM_JOBID"
echo "SLURM_JOB_NODELIST: $SLURM_JOB_NODELIST"
echo "SLURM_NODELIST: $SLURM_NODELIST"
echo "SLURM_NODENAMES: $SLURM_NODENAMES"
echo "SLURM_JOB_NUM_NODES: $SLURM_JOB_NUM_NODES"
echo "SLURM_NTASKS_PER_NODE: $SLURM_NTASKS_PER_NODE"
echo "SLURM_NTASKS: $SLURM_NTASKS"

module load gcc/4.8.2
module load openmpi/1.8.3
rm nodes.txt
scontrol show hostname $SLURM_NODELIST > nodes.txt

#NPROCS=$SLURM_JOB_NUM_NODES
NPROCS=96
NPERNODE=6
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/home/apps/lib/boost/1.57"
mpirun -np $NPROCS --npernode $NPERNODE /home/operedo/Projects/dgraph/dgraph 
echo "done :)"
