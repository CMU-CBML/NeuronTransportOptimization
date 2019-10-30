#!/bin/bash

# make all

# mpiexec -n 1 ./neuron_opt -f ../io/single_pipe/ -ksp_type fgmres -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_schur_fact_type lower -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type jacobi -fieldsplit_1_pc_type jacobi -fieldsplit_1_inner_ksp_type preonly -fieldsplit_1_inner_pc_type jacobi -fieldsplit_1_upper_ksp_type preonly -fieldsplit_1_upper_pc_type jacobi
#  -ksp_error_if_not_converged

# mpiexec -n 1 ./neuron_opt -f ../io/burger/ -ksp_type fgmres -pc_type fieldsplit -pc_fieldsplit_type schur -fieldsplit_0_ksp_type cg -fieldsplit_0_pc_type jacobi -ksp_error_if_not_converged

# mpiexec -n 1 ./neuron_opt -f ../io/burger/ -ksp_type tfqmr -pc_type fieldsplit -pc_fieldsplit_type schur -fieldsplit_0_ksp_type gmres -fieldsplit_0_pc_type jacobi -fieldsplit_1_pc_type jacobi -fieldsplit_1_inner_ksp_type preonly -fieldsplit_1_inner_pc_type jacobi -fieldsplit_1_upper_ksp_type preonly -fieldsplit_1_upper_pc_type jacobi -ksp_error_if_not_converged

# mpiexec -n 1 ./neuron_opt -f ../io/burger/ -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -ksp_error_if_not_converged

# mpiexec -n 1 ./neuron_opt -f ../io/burger2/ -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -ksp_error_if_not_converged

# mpiexec -n 1 ./neuron_opt -f ../io/burger3/ -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -ksp_error_if_not_converged

# mpiexec -n 1 ./neuron_opt -f ../io/burger4/ -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -ksp_error_if_not_converged

# mpiexec -n 1 ./neuron_opt -f ../io/heat3/ -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -ksp_error_if_not_converged

mpiexec -n 1 ./neuron_opt -f ../io/heat4/ -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type mumps -ksp_error_if_not_converged
