#!/bin/bash -l
#PBS -l nodes=3:ppn=128
#PBS -l mem=64gb
#PBS -l walltime=1:00:00
#PBS -r n
#PBS -j oe
#PBS -q starq

# go to your working directory containing the batch script, code and data
cd /fs/lustre/cita/zack/jl/dev/LagrangianPerturbationTheory/projectamd
module load openmpi

julia --project -e 'using Pkg; Pkg.instantiate()'
julia --project -e 'using Pkg; Pkg.precompile()'

/fs/lustre/cita/zack/jl/bin/mpiexecjl --project -n 192 julia fft_filter_6144.jl
