#!/bin/bash
#SBATCH -A mp107d
#SBATCH --qos=debug
#SBATCH -C cpu
#SBATCH -t 0:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -o /pscratch/sd/x/xzackli/joboutput/%x.o%j
#SBATCH -J no768
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zackli@berkeley.edu

cd /pscratch/sd/x/xzackli/websky_convert/LagrangianPerturbationTheory.jl/run/
module load cray-mpich
module load cray-hdf5-parallel
export JULIA_NUM_THREADS=16

julia --project=. -e \
    'using Pkg; using InteractiveUtils;
     Pkg.instantiate(); Pkg.precompile(); Pkg.status(); versioninfo();
     using MPI; println("MPI: ", MPI.identify_implementation());'

/global/homes/x/xzackli/.julia/bin/mpiexecjl --project=. -n 16 -c 16 julia test_fft_no768.jl
