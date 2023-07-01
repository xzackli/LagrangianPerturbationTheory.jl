using Base.Iterators: partition
using Printf

# set up parameters which we'll be looping over
ΔlogM = 0.01f0
masses = 11:ΔlogM:12
masses_per_job = 10

# set up directories (change for your situation)
outdir = "/fs/lustre/scratch/zack/ICs/low_mass_halos/"
workdir = "/fs/lustre/cita/zack/jl/dev/LagrangianPerturbationTheory/examples/PBS/"
mkpath(outdir)
mkpath(joinpath(outdir, "scripts"))

# loop over different parameters and create a bunch of PBS files
for job_masses in partition(masses, masses_per_job)
    logMmin, logMmax = minimum(job_masses), maximum(job_masses) + ΔlogM  # inclusive
    PBS_script = """
#!/bin/bash -l
#PBS -l nodes=1:ppn=32
#PBS -l mem=128gb
#PBS -l walltime=1:00:00
#PBS -r n
#PBS -j oe
#PBS -q starq

# go to your working directory containing the batch script, code and data
cd $workdir

julia --project=. -t 32 generate_tracers.jl $logMmin $logMmax $ΔlogM $outdir
    """
    mass_string = @sprintf("%.2f",logMmin) * "_" * @sprintf("%.2f",logMmax)
    scriptfile = joinpath(outdir, "scripts", "job_$(mass_string).pbs")
    open(scriptfile, "w") do file write(file, PBS_script) end
    run(`qsub $scriptfile`)  # submit PBS file to cluster
end
