using MPI
using HDF5
using PencilArrays
using PencilFFTs
using AbstractFFTs
using Logging
using LinearAlgebra  # for mul!, ldiv!

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)

if rank == 1
    io = open("log.txt", "w+")
    logger = SimpleLogger(io)
end

function rank1log(s, rank)
    if rank == 1
        with_logger(logger) do
            @info(s)
        end
        flush(io)
    end
end

# Input data dimensions (Nx × Ny × Nz)
dims = (6144, 6144, 6144)
comm = MPI.COMM_WORLD
pen = Pencil(dims, comm)


# Apply a 3D real-to-complex (r2c) FFT.
transform = Transforms.RFFT()

# Note that, for more control, one can instead separately specify the transforms along each dimension:
# transform = (Transforms.RFFT(), Transforms.FFT(), Transforms.FFT())

# Create plan
plan = PencilFFTPlan(pen, transform)


rank1log("plan made, about to allocate input", rank)

# In our example, this returns a 3D PencilArray of real data (Float64).
u = allocate_input(plan)

# read in the HDF5 handle
# read local slice
# set the local pencil to it

fname = "/fs/lustre/project/act/zack/ICs/Fvec_7700Mpc_n6144_nb30_nt16_p.h5"
global_u = global_view(u)
global_axes = map( ax -> first(ax):last(ax), axes(global_u) )


rank1log("$(MPI.Comm_rank(comm)): $(global_axes)", rank)

u .= h5read(fname, "data", global_axes)

rank1log("data read", rank)

# In our example, this returns a 3D PencilArray of complex data (Complex{Float64}).
v = allocate_output(plan)


rank1log("output allocated", rank)


# Apply plan on `u` with `v` as an output
mul!(v, plan, u)

rank1log("fft applied", rank)

# Filter with a Gaussian 

# In our case (Lx = 2π and Nx even), this gives kx = [0, 1, 2, ..., Nx/2].
kx = rfftfreq(dims[1])

# In our case (Ly = 2π and Ny even), this gives
# ky = [0, 1, 2, ..., Ny/2-1, -Ny/2, -Ny/2+1, ..., -1] (and similarly for kz).
ky = fftfreq(dims[2])
kz = fftfreq(dims[3])

kvec = (kx, ky, kz)
grid_fourier = localgrid(v, kvec)

σ = 2.5
filter_factor = -2 * π^2 * σ^2
@. v *= exp( filter_factor * (grid_fourier[1]^2 +  grid_fourier[2]^2 +  grid_fourier[3]^2)  )


rank1log("filter applied", rank)

# Apply backward plan on `v` with `w` as an output
# w = similar(u)
ldiv!(u, plan, v)  # now w ≈ u

rank1log("everything done, saving to disk", rank)

open(PencilArrays.PencilIO.PHDF5Driver(), "filtered_$(fname).h5", MPI.COMM_WORLD; write=true) do ff
    ff["data"] = u
end


rank == 1 && close(io)
