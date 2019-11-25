#!/bin/bash

#SBATCH -N 1
#SBATCH -p RM-small
#SBATCH --ntasks-per-node 2
#SBATCH -t 00:30:00
# echo commands to stdout 
set -x

# move to your appropriate pylon5 directory
# this job assumes:
#  - all input data is stored in this directory
#  - all output should be stored in this directory
cd /pylon5/eg560mp/angranl/NeuronTransportOptimization/src_heat/

make all
# run OpenMP program
# export OMP_NUM_THREADS=28
mpiexec -n 2 ./neuron_opt -f ../io/heat2/ -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -ksp_error_if_not_converged
