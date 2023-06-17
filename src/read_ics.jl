
function load_example_ics()
    rootpath = artifact"density_example"
    path = joinpath(rootpath, "Fvec_7700Mpc_n6144_nb30_nt16_ks192.h5")
    delta = h5open(path, "r") do file
        read(file, "delta")
    end
    return delta
end

# redshift = 1 / scale_factor - 1

"""
octants are indexed -1 and 0
"""
function lattice_value(lpt::InitialConditionsWebsky{T}, i, j, k, octant_i, octant_j, octant_k) where T

    x = (i+T(0.5)) * lpt.grid_spacing + octant_i * lpt.boxsize_x
    y = (j+T(0.5)) * lpt.grid_spacing + octant_j * lpt.boxsize_y
    z = (k+T(0.5)) * lpt.grid_spacing + octant_k * lpt.boxsize_z
    
    chi = sqrt(x^2+y^2+z^2)
    scale_factor = scale_factor_of_chi(lpt.cosmology, chi)
    δ₀ = lpt.field[k+1,j+1,i+1]  # (py → jl) reverses dimensions, index+1
    field_linear = δ₀ * growth_factor(lpt.cosmology, scale_factor)

    return (x=x, y, z, scale_factor, field_linear)
end
