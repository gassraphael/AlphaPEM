#!/bin/bash
#SBATCH --job-name=AlphaPEM
#SBATCH --partition=cpu
#SBATCH --time=167:00:00
#SBATCH --output=stdout.slurm_EH-31_2.0&2.25
#SBATCH --error=stderr.slurm_EH-31_2.0&2.25
#SBATCH --cpus-per-task=80
srun -N1 -n1 python 'parameter_calibration.py'
