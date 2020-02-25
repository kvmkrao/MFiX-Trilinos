#!/bin/bash
#SBATCH -J MPIj           # job name
#SBATCH -o output%j.txt     # output and error file name (%j expands to jobID)
#SBATCH -n 16              # total number of mpi tasks requested
#SBATCH -N 1              # total number of nodes 
#SBATCH -p flat-quadrant  # normal  
#SBATCH -t 00:60:00       # run time (hh:mm:ss) - 1.5 hours

export OMP_PROC_BIND=spread  #"spread" is also good for Intel and CCE compilers
export OMP_PLACES=threads    #to run with 8 threads per task (unpacked)
export OMP_NUM_THREADS=1
export OMP_DISPLAY_ENV=verbose
#export I_MPI_PIN_DOMAIN=auto

# Unset any MPI Affinities
export MV2_USE_AFFINITY=0
export MV2_ENABLE_AFFINITY=0
export VIADEV_USE_AFFINITY=0
export VIADEV_ENABLE_AFFINITY=0

#Run the job
#time mpirun -np 60  /work/04051/vkotteda/stampede2/trilinos_3rd_openmp/run/precond_ifpack2/modified_src_matrixmarket/example
time ibrun tacc_affinity /work/04051/vkotteda/stampede2/trilinos_3rd_openmp/run/precond_ifpack2/modified_src_matrixmarket/example
#time ibrun  numa.sh /work/04051/vkotteda/stampede2/trilinos_3rd_openmp/run/precond_ifpack2/modified_src_matrixmarket/example

#mpirun  -n 2 -device=p100 -count=2  ./clean4.exe 
#./get_local_rank ./osu_latency D D

