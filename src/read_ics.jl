

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

# redshift = 1 / scale_factor - 1

# a::T      # scale factor
# chi::LT   # comoving radial distance


# """
#     lattice_0(ic, i, j, k, oi, oj, ok)

# Returns the lattice value ``f_0(\\vec{q})`` at ``z=0`` for some field ``f``,
# where `i`, `j`, and `k` are the Lagrangian lattice indices so that 
# ``\\vec{q} = (q_i, q_j, q_k)``, and `oi`, `oj`, and `ok` refer to the octant 
# index (Websky is repeated twice in each dimension with the observer at the 
# center, to form eight subvolumes). Octant indices are either -1 or 0.

# **Do not access the array directly, use this function instead.** Python's array
# format flips the dimensions of arrays in storage.
# """
# function lattice_0(ic::InitialConditionsWebsky{T}, 
#                   i::Int, j::Int, k::Int, oi::Int, oj::Int, ok::Int) where T
#     loc = lattice_location(ic, i, j, k, oi, oj, ok)
#     f₀ = ic.field[loc.k+1,loc.j+1,loc.i+1]  # (py → jl) reverses dimensions
#     return f₀
# end

# function lattice_0(ic::InitialConditionsWebsky{T}, 
#                    loc::LatticeLocation) where T
#     f₀ = ic.field[loc.k+1,loc.j+1,loc.i+1]  # (py → jl) reverses dimensions
#     return f₀
# end



# function lattice_location(ic::InitialConditionsWebsky{T}, 
#         i::Int, j::Int, k::Int, oi::Int, oj::Int, ok::Int) where T
#     half = T(1//2)
#     x = (i + half) * ic.grid_spacing + oi * ic.boxsize_x
#     y = (j + half) * ic.grid_spacing + oj * ic.boxsize_y
#     z = (k + half) * ic.grid_spacing + ok * ic.boxsize_z
#     chi = sqrt(x^2+y^2+z^2)
#     scale_factor = scale_factor_of_chi(ic.cosmo, chi)
#     return LatticeLocation(scale_factor, chi, x, y, z, i, j, k)
# end


# function index_to_pos(ic::InitialConditionsWebsky{T}, 𝐢ₕ, oi::Int, oj::Int, ok::Int) where T
#     half = T(1//2)
#     i,j,k = 𝐢ₕ
#     x = (i + half) * ic.grid_spacing + oi * ic.boxsize_x
#     y = (j + half) * ic.grid_spacing + oj * ic.boxsize_y
#     z = (k + half) * ic.grid_spacing + ok * ic.boxsize_z
#     return SA[x, y, z]
# end

# function pos_to_scale_factor(ic, 𝐪::SVector{3,T}) where T
#     chi = sqrt(𝐪[1]^2+𝐪[2]^2+𝐪[3]^2)
#     scale_factor = scale_factor_of_chi(ic.cosmo, chi)
#     return scale_factor
# end

# function getindex(ic::InitialConditionsWebsky, q_lagrangian::LatticeLocation)
#     return lattice_0(ic, q_lagrangian)
# end

# function getindex(ic::InitialConditionsWebsky, ijk::SVector)
#     i, j, k = ijk
#     f₀ = ic.field(k+1, j+1, i+1)  # (py → jl) reverses dimensions
#     return f₀
# end


# function read_websky_ics(file_den, file_sx, file_sy, file_sz, boxsize, cosmo)
#     delta_array = h5open(file_den, "r") do f read(f, "data") end
#     sx_array = h5open(file_sx, "r") do f read(f, "data") end
#     sy_array = h5open(file_sy, "r") do f read(f, "data") end
#     sz_array = h5open(file_sz, "r") do f read(f, "data") end

#     grid_spacing = boxsize / size(delta_array,1)
#     δ₀ = InitialConditionsWebsky(FirstOrderLPT, grid_spacing, cosmo, 
#         interpolate(delta_array, BSpline(Linear())))
#     sx = InitialConditionsWebsky(FirstOrderLPT, grid_spacing, cosmo, 
#         interpolate(sx_array, BSpline(Linear())))
#     sy = InitialConditionsWebsky(FirstOrderLPT, grid_spacing, cosmo, 
#         interpolate(sy_array, BSpline(Linear())))
#     sz = InitialConditionsWebsky(FirstOrderLPT, grid_spacing, cosmo, 
#         interpolate(sz_array, BSpline(Linear())))

#     return (; δ₀, sx, sy, sz)
# end

