

function example_ics_path()
    rootpath = artifact"density_example"
    return joinpath(rootpath, "Fvec_7700Mpc_n6144_nb30_nt16_ks192.h5")
end

function load_example_ics()
    path = example_ics_path()
    delta = h5open(path, "r") do file
        read(file, "delta")
    end
    # uploaded data was reversed from what it should be
    return permutedims(delta, (3,2,1))  
end

# todo: read the sliced array
function read_websky_ics(file_den, file_sx, file_sy, file_sz, cosmo, boxsize=7700.f0u"Mpc")
    
    delta_array = h5open(file_den, "r") do f read(f, "data") end
    sx_array = h5open(file_sx, "r") do f read(f, "data") end
    sy_array = h5open(file_sy, "r") do f read(f, "data") end
    sz_array = h5open(file_sz, "r") do f read(f, "data") end

    # assuming cube
    @assert (
        size(delta_array,1) == size(delta_array,2) == size(delta_array,3) ==
        size(sx_array,1) == size(sx_array,2) == size(sx_array,3) ==
        size(sy_array,1) == size(sy_array,2) == size(sy_array,3) ==
        size(sz_array,1) == size(sz_array,2) == size(sz_array,3))

    grid_spacing = boxsize / size(delta_array,1)
    box_sizes = (boxsize, boxsize, boxsize)
    offset = grid_spacing / 2
    nx, ny, nz = size(delta_array)  # make a grid for each Lagrangian coordinate
    q_axes = (LinRange(offset, offset + (nx - 1) * grid_spacing, nx),
              LinRange(offset, offset + (ny - 1) * grid_spacing, ny),
              LinRange(offset, offset + (nz - 1) * grid_spacing, nz))
    grid = LagrangianGridWebsky(cosmo, grid_spacing, box_sizes, q_axes)
    δ₀ = ICFieldWebsky(FirstOrderLPT, grid, delta_array)
    sx = ICFieldWebsky(FirstOrderLPT, grid, sx_array * u"Mpc")
    sy = ICFieldWebsky(FirstOrderLPT, grid, sy_array * u"Mpc")
    sz = ICFieldWebsky(FirstOrderLPT, grid, sz_array * u"Mpc")
    𝚿₀ =  DisplacementICFieldWebsky(sx, sy, sz)
    return (; δ₀, 𝚿₀)
end
