
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
    lattice0(ic, i, j, k, oi, oj, ok)

Returns the lattice value ``f_0(\\vec{q})`` at ``z=0`` for some field ``f``,
where `i`, `j`, and `k` are the Lagrangian lattice indices so that 
``\\vec{q} = (q_i, q_j, q_k)``, and `oi`, `oj`, and `ok` refer to the octant 
index (Websky is repeated twice in each dimension with the observer at the 
center, to form eight subvolumes). Octant indices are either -1 or 0.  
"""
function lattice0(ic::InitialConditionsWebsky{T}, i, j, k, oi, oj, ok) where T
    half = T(1//2)
    x = (i + half) * ic.grid_spacing + oi * ic.boxsize_x
    y = (j + half) * ic.grid_spacing + oj * ic.boxsize_y
    z = (k + half) * ic.grid_spacing + ok * ic.boxsize_z
    
    chi = sqrt(x^2+y^2+z^2)
    scale_factor = scale_factor_of_chi(ic.cosmo, chi)
    δ₀ = ic.field[k+1,j+1,i+1]  # (py → jl) reverses dimensions, index+1
    field_linear = δ₀ * growth_factor(ic.cosmo, scale_factor)

    return (x=x, y, z, scale_factor, field_linear)
end
