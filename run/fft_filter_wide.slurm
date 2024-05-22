#!/bin/bash
#SBATCH -A mp107d
#SBATCH --qos=preempt
#SBATCH -C cpu
#SBATCH -t 2:00:00
#SBATCH --nodes=32
#SBATCH --ntasks-per-node=128
#SBATCH -o /pscratch/sd/x/xzackli/joboutput/%x.o%j
#SBATCH -J fullresfilter_wide
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zackli@berkeley.edu

# go to your working directory containing the batch script, code and data

cd /pscratch/sd/x/xzackli/websky_convert/LagrangianPerturbationTheory.jl/run/
module load cray-mpich
module load cray-hdf5-parallel
export JULIA_NUM_THREADS=1

# julia --project -e 'using Pkg; Pkg.instantiate()'
# julia --project -e 'using Pkg; Pkg.precompile()'

/global/homes/x/xzackli/.julia/bin/mpiexecjl --project=. -n 4096 --cpu-bind=cores -c 2 julia fft_filter_6144.jl
