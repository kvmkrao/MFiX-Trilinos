#!/bin/bash
#SBATCH -J MPIj           # job name
#SBATCH -o output%j.txt     # output and error file name (%j expands to jobID)
#SBATCH -n 960              # total number of mpi tasks requested
#SBATCH -N 16              # total number of nodes 
#SBATCH -p flat-quadrant  # normal  
#SBATCH -t 00:10:00       # run time (hh:mm:ss) - 1.5 hours

export OMP_PROC_BIND=spread  #"spread" is also good for Intel and CCE compilers
export OMP_PLACES=threads    #to run with 8 threads per task (unpacked)
export OMP_NUM_THREADS=2
export OMP_DISPLAY_ENV=verbose
#export I_MPI_PIN_DOMAIN=auto

#Run the job
time mpirun -np 960  /work/04051/vkotteda/stampede2/trilinos_3rd_openmp/run/precond_ifpack2/solve_ifpack.exe   
#mpirun  -n 2 -device=p100 -count=2  ./clean4.exe 
#./get_local_rank ./osu_latency D D

