using MPI
using HDF5
using PencilArrays
using PencilFFTs
using AbstractFFTs
using LinearAlgebra  # for mul!, ldiv!

MPI.Init()
comm = MPI.COMM_WORLD

# Input data dimensions (Nx × Ny × Nz)
dims = (768, 768, 768)
comm = MPI.COMM_WORLD
fname = "/pscratch/sd/x/xzackli/websky_convert/data/Fvec_7700Mpc_n6144_nb30_nt16_no768_p.h5"



function filter_grid(fname, dims, comm)
    pen = Pencil(dims, comm)
    transform = Transforms.RFFT()
    # Create plan
    plan = PencilFFTPlan(pen, transform)

    # In our example, this returns a 3D PencilArray of real data (Float64).
    u = allocate_input(plan)

    # read in the HDF5 handle
    # read local slice
    # set the local pencil to it

    global_u = global_view(u)
    global_axes = map( ax -> first(ax):last(ax), axes(global_u) )
    u .= h5read(fname, "data", global_axes)

    println("$(MPI.Comm_rank(comm)): $(global_axes)")

    # In our example, this returns a 3D PencilArray of complex data (Complex{Float64}).
    v = allocate_output(plan)

    # Apply plan on `u` with `v` as an output
    mul!(v, plan, u)

    kx = fftfreq(dims[1])
    ky = fftfreq(dims[2])
    kz = fftfreq(dims[3])

    kvec = (kx, ky, kz)
    grid_fourier = localgrid(v, kvec)

    σ = 2.5
    filter_factor = -2 * π^2 * σ^2
    @. v *= exp( filter_factor * (grid_fourier[1]^2 + grid_fourier[2]^2 +  grid_fourier[3]^2)  )

    ldiv!(u, plan, v)  # now w ≈ u
    return u
end

u = filter_grid(fname, dims, comm)

println("WRITING, I am $(MPI.Comm_rank(comm)) of $(MPI.Comm_size(comm))")

open(PencilArrays.PencilIO.PHDF5Driver(), "filtered_output.h5", MPI.COMM_WORLD; write=true) do ff
    ff["data"] = u
end
