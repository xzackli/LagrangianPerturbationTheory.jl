
abstract type AbstractLPT end

struct FirstOrderLPT <: AbstractLPT end

struct InitialConditionsWebsky{T, LPT, C, AA} <: AbstractLPT 
    grid_spacing::T
    cosmology::C
    field::AA

    boxsize_x::T
    boxsize_y::T
    boxsize_z::T
end

function InitialConditionsWebsky(::Type{LPT}, grid_spacing, cosmology::C, 
                                 field::AA) where {T, LPT, C, AA<:AbstractArray{T}}
    boxsize_x = size(field, 1) * grid_spacing
    boxsize_y = size(field, 2) * grid_spacing
    boxsize_z = size(field, 3) * grid_spacing
    return InitialConditionsWebsky{T, LPT, C, AA}(
        grid_spacing, cosmology, field, boxsize_x, boxsize_y, boxsize_z)
end
