
abstract type AbstractLPT end

struct FirstOrderLPT <: AbstractLPT end

struct InitialConditionsWebsky{T, LPT, C, AA, TL} <: AbstractLPT 
    cosmo::C
    field::AA
    grid_spacing::TL
    boxsize_x::TL
    boxsize_y::TL
    boxsize_z::TL
end

"""
    InitialConditionsWebsky(::Type{LPT}, grid_spacing, cosmo, field)

Returns a wrapper around an array `field` in Lagrangian coordinates, defined by 
a kind of Lagrangian perturbation theory `LPT`, grid spacing, and cosmo. For 
Websky, this kind of field is typically something like ``\\delta_0(\\vec{q})``, 
the Lagrangian density at ``z=0``.
"""
function InitialConditionsWebsky(::Type{LPT}, grid_spacing::TL, cosmo::C, field::AA
        ) where {T, LPT, C, AA<:AbstractArray{T}, TL}
    boxsize_x = size(field, 1) * grid_spacing
    boxsize_y = size(field, 2) * grid_spacing
    boxsize_z = size(field, 3) * grid_spacing
    return InitialConditionsWebsky{T, LPT, C, AA, TL}(
        cosmo, field, grid_spacing, boxsize_x, boxsize_y, boxsize_z)
end
