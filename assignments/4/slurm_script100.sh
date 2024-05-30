#!/bin/bash
#SBATCH --job-name=mpi_job              # Job name
#SBATCH --nodes=8                       # Number of nodes
#SBATCH --ntasks=33                     # Total number of tasks
#SBATCH --ntasks-per-node=5             # Number of tasks per node
#SBATCH --time=01:00:00                 # Time limit
#SBATCH --output=output_%j.txt
#SBATCH --error=error_%j.txt

# Load necessary modules (adjust to your environment)
# module load mpi

# Run the MPI program
mpirun -np $SLURM_NTASKS ./4 100 1000000 1