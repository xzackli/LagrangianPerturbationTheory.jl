#!/bin/bash
#SBATCH -A mp107d
#SBATCH --qos=debug
#SBATCH -C cpu
#SBATCH -t 0:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1
#SBATCH -o /pscratch/sd/x/xzackli/joboutput/%x.o%j
#SBATCH -J no768
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zackli@berkeley.edu

# go to your working directory containing the batch script, code and data

cd /pscratch/sd/x/xzackli/websky_convert/LagrangianPerturbationTheory.jl/run/
module load cray-mpich
module load cray-hdf5-parallel
export JULIA_NUM_THREADS=1

# julia --project -e 'using Pkg; Pkg.instantiate()'
# julia --project -e 'using Pkg; Pkg.precompile()'

/global/homes/x/xzackli/.julia/bin/mpiexecjl --project=. -n 3 julia test_fft_no768.jl
