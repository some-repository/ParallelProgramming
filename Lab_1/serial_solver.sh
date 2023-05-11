#!/bin/bash

#PBS -l walltime=00:10:00,nodes=1:ppn=1
#PBS -N serial_solver
#PBS -q batch

cd $PBS_O_WORKDIR
mpirun --hostfile $PBS_NODEFILE -np 1 ./serial_solver