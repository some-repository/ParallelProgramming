#!/bin/bash

#PBS -l walltime=00:10:00,nodes=1:ppn=4
#PBS -N parallel_solver
#PBS -q batch

cd $PBS_O_WORKDIR
mpirun --hostfile $PBS_NODEFILE -np 4 ./parallel_solver