using MPI
using HDF5
using PencilArrays
using PencilFFTs
using AbstractFFTs
using LinearAlgebra  # for mul!, ldiv!
PencilFFTs.FFTW.set_num_threads(Threads.nthreads())  # match thread count


function read_into_pencil_array!(u, fname, comm, info)
    ff = h5open(fname, "r", comm, info)
    slices = range_local(u)
    dset = ff["data"]
    u .= dset[slices...]
    close(ff)
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

function write_grid(u::AbstractArray{Complex{T}}, fname, comm, info) where T
    out_fname = replace(fname, (".h5" =>  "_filtered.h5"))
    slices = range_local(u)
    ff = h5open(out_fname, "w", comm, info)
    dset = create_dataset(ff, "data", datatype(T), dataspace(size_global(u)))
    local_array = parent(u)
    dset[slices...] = real.(local_array)
    close(ff)
end

function run(T, dims, fname)
    MPI.Init()
    comm = MPI.COMM_WORLD       # MPI communicator
    rank = MPI.Comm_rank(comm)  # rank of local process
    info = MPI.Info()

    pen = Pencil(dims, comm)
    transform = Transforms.FFT!()
    plan = PencilFFTPlan(pen, transform, T)
    rank == 1 && print(plan)

    us = allocate_input(plan)
    u = first(us)
    û = last(us)

    read_into_pencil_array!(u, fname, comm, info)
    println("FINISHED READING, I am $(MPI.Comm_rank(comm)) of $(MPI.Comm_size(comm)) with " 
        * string(range_local(first(us))))
    MPI.Barrier(comm)
    GC.gc()
    
    plan * us  # apply FFT
    apply_gaussian!(û, dims)  # filter û
    plan \ us  # apply IFFT
    println("FINISHED FILTERING, I am $(MPI.Comm_rank(comm)) of $(MPI.Comm_size(comm))")
    GC.gc()

    write_grid(u, fname, comm, info)  # write the real-space array
    MPI.Comm_rank(comm) == 1 && println("FINISHED WRITING")
end


dims = (6144, 6144, 6144)
fname = "/pscratch/sd/x/xzackli/websky_convert/displacementprocessing/Fvec_7700Mpc_n6144_nb30_nt16_p.h5"
run(Float32, dims, fname)
