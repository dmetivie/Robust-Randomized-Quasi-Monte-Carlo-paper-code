#!/bin/bash
# Shell to use for running the job
#$ -S /bin/bash

#Job name
#$ -N rqmc

# Queue name
#$ -q short.q #! change that

# Identifier of a node of the chosen queue
# -l hostname=n12 #! you might want to pick

# Export of all environment variables
#$-V

# Output standard
#$ -o sge_julia.out

# Error output
#$ -e sge_julia.err

# Run the command from the directory where the script is launched
#$ -cwd

# Use n CPU
#$ -pe wire 64 #! you can change that

# Example for using Julia:

conda activate julia-1.9.3
julia run_g1.jl
conda deactivate