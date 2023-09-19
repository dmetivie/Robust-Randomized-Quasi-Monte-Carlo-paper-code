#!/bin/bash

## BEGIN SBATCH directives
#SBATCH --partition=mistea #! change that

#SBATCH --job-name="rqmc"
#SBATCH --nodes=1                # node count
#SBATCH --ntasks-per-node=28     # total number of tasks across each nodes (>1 if distributed)
## SBATCH --time=80:00:00        # I though this was mendatory but apparently it is not on all cluster...

#SBATCH --mail-user=your.emmail@something.edu #! change that!!!
#SBATCH --mail-type=END
#SBATCH --output=output.txt
## END SBATCH directives

echo "SLURM_JOB_ID " $SLURM_JOB_ID # id of the job useful to track/cancel it
echo "SLURM_NODELIST "$SLURM_NODELIST # all the node id you requested
echo "SLURM_NTASKS_PER_NODE " $SLURM_NTASKS_PER_NODE # number of task per node 
echo "SLURM_NPROCS " $SLURM_NPROCS # all the cpu requested on the different nodes (total number of CPUs allocated)

## To clean and load modules defined at the compile and link phases
module purge
module load julia/1.9.1

## Execution
julia run_g1.jl