using MPI
using HDF5
using PencilArrays
using PencilFFTs
using AbstractFFTs
using LinearAlgebra  # for mul!, ldiv!
PencilFFTs.FFTW.set_num_threads(Threads.nthreads())  # match thread count
MPI.Init()
comm = MPI.COMM_WORLD

function setup_grid(fname, plan, comm)
    u = allocate_input(plan)
    open(PencilArrays.PencilIO.PHDF5Driver(), fname, comm, read=true) do ff
        read!(ff, u, "data")
    end
    return u
end

function apply_gaussian!(v::AbstractArray{Complex{T}}, dims) where T
    kvec = map(d -> T.(fftfreq(d)), dims)
    grid_fourier = localgrid(v, kvec)
    kx, ky, kz = grid_fourier[1], grid_fourier[2], grid_fourier[3]
    σ = 2.5
    filter_factor = T(-2 * π^2 * σ^2)
    Threads.@threads for k in axes(v, 3)
        for j in axes(v, 2), i in axes(v, 1)
            v[i,j,k] *= exp( filter_factor * (kx[i]^2 + ky[j]^2 + kz[k]^2) )
        end
    end
end

function filter_grid!(u, v, plan, dims)
    mul!(v, plan, u)
    apply_gaussian!(v, dims)
    ldiv!(u, plan, v)
end

function run(T, comm)
    # Input data dimensions (Nx × Ny × Nz)
    dims = (6144, 6144, 6144)
    fname = "/pscratch/sd/x/xzackli/websky_convert/displacementprocessing/Fvec_7700Mpc_n6144_nb30_nt16_p.h5"

    pen = Pencil(dims, comm)
    transform = Transforms.RFFT()
    plan = PencilFFTPlan(pen, transform, T)
    if MPI.Comm_rank(comm) == 1
        print(plan)
    end

    u = setup_grid(fname, plan, comm)
    println("FINISHED READING, I am $(MPI.Comm_rank(comm)) of $(MPI.Comm_size(comm))")
    GC.gc(true); GC.gc(false); GC.gc(true); GC.gc(false); GC.gc(true); GC.gc(false);
    GC.gc(true); GC.gc(false); GC.gc(true); GC.gc(false); GC.gc(true); GC.gc(false);
    
    v = allocate_output(plan)
    println("FINISHED ALLOCATING OUTPUT, I am $(MPI.Comm_rank(comm)) of $(MPI.Comm_size(comm))")
    filter_grid!(u, v, plan, dims)
    println("FINISHED FILTERING, I am $(MPI.Comm_rank(comm)) of $(MPI.Comm_size(comm))")

    open(PencilArrays.PencilIO.PHDF5Driver(), "filtered_$(fname)", MPI.COMM_WORLD; write=true) do ff
        ff["data"] = u
    end
end

run(Float32, comm)
