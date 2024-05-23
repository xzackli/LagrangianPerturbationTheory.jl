using MPI
using HDF5
using PencilArrays
using PencilFFTs
using AbstractFFTs
using LinearAlgebra  # for mul!, ldiv!

MPI.Init()
comm = MPI.COMM_WORLD

# Input data dimensions (Nx × Ny × Nz)
dims = (6144, 6144, 6144)
fname = "/pscratch/sd/x/xzackli/websky_convert/displacementprocessing/Fvec_7700Mpc_n6144_nb30_nt16_p.h5"

pen = Pencil(dims, comm)
transform = Transforms.RFFT()
plan = PencilFFTPlan(pen, transform, Float32)

function setup_grid(fname, plan, comm)
    u = allocate_input(plan)
    v = allocate_output(plan)
    open(PencilArrays.PencilIO.PHDF5Driver(), fname, comm, read=true) do ff
        read!(ff, u, "data")
    end
    println("FINISHED READING, I am $(MPI.Comm_rank(comm)) of $(MPI.Comm_size(comm))")
    return u, v
end

function filter_grid!(u::AbstractArray{T}, v, plan, dims) where T
    mul!(v, plan, u)

    kx = fftfreq(dims[1])
    ky = fftfreq(dims[2])
    kz = fftfreq(dims[3])

    kvec = (kx, ky, kz)
    grid_fourier = localgrid(v, kvec)

    σ = 2.5
    filter_factor = T(-2 * π^2 * σ^2)
    @. v *= exp( filter_factor * (grid_fourier[1]^2 + grid_fourier[2]^2 +  grid_fourier[3]^2)  )

    ldiv!(u, plan, v)  # now w ≈ u
    return u
end

u, v = setup_grid(fname, plan, comm)
Base.GC.gc()
filter_grid!(u, v, plan, dims)
println("FINISHED FILTERING, I am $(MPI.Comm_rank(comm)) of $(MPI.Comm_size(comm))")

open(PencilArrays.PencilIO.PHDF5Driver(), "filtered_$(fname)", MPI.COMM_WORLD; write=true) do ff
    ff["data"] = u
end
