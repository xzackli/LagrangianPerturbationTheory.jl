using MPI
using PencilFFTs
using Random
# using LinearAlgebra  # for mul!, ldiv!

function main()
    MPI.Init()
    comm = MPI.COMM_WORLD

    # Input data dimensions (Nx × Ny × Nz)
    dims = (4800, 4800, 4800)
    comm = MPI.COMM_WORLD
    pen = Pencil(dims, comm)

    # u_re_io = PencilArray{Float32}(undef, pen)  # allocate a real array for reading from disk

    # Apply a 3D real-to-complex (r2c) FFT.
    transform = Transforms.FFT!()

    # Create plan
    plan = PencilFFTPlan(pen, transform, Float32)

    # In our example, this returns a 3D PencilArray of real data (Float64).
    us = allocate_input(plan)
    u = first(us); û = last(us)

    # Fill the array with some (random) data
    randn!(u)
    plan * us
    plan \ us

    println("us: $(Base.summarysize(us) / 2^20), plan: $(Base.summarysize(plan) / 2^20), pen: $(Base.summarysize(pen) / 2^20) $(MPI.Comm_rank(comm)) of $(MPI.Comm_size(comm))")
end
main()
