
abstract type AbstractLPT end
struct FirstOrderLPT <: AbstractLPT end

abstract type AbstractLagrangianGrid end
struct LagrangianGridWebsky{T,TL,C,R} <: AbstractLagrangianGrid 
    cosmo::C
    grid_spacing::TL
    box_sizes::NTuple{3,TL}
    q_axes::NTuple{3,R}  # axis grids for each Lagrangian coordinate
end

function LagrangianGridWebsky(cosmo::C, grid_spacing::TL, 
        box_sizes, q_axes::NTuple{3,R}) where {T, C, R, TL<:Quantity{T}}
    LagrangianGridWebsky{T,TL,C,R}(cosmo, grid_spacing, box_sizes, q_axes)
end

struct LagrangianCoordinate{T}
    x::T
    y::T
    z::T
end

function lagrangian_coordinate(grid, i, j, k, oi, oj, ok)
    Δq = grid.grid_spacing
    offset = Δq / 2
    box_size_x, box_size_y, box_size_z = grid.box_sizes
    x = i * Δq + offset + oi * box_size_x
    y = j * Δq + offset + oj * box_size_y
    z = k * Δq + offset + ok * box_size_z
    return LagrangianCoordinate(x, y, z)
end

# websky grid is normalized so that |𝐪| is comoving distance
function scale_factor(grid::LagrangianGridWebsky, 𝐪::LagrangianCoordinate)
    chi = √(𝐪.x^2 + 𝐪.y^2 + 𝐪.z^2)
    return scale_factor_of_chi(grid.cosmo, chi)
end

random_number_in_cell(x, grid_spacing::Q) where {T, Q<:Quantity{T}} =
    x + (rand() - T(1//2)) * grid_spacing
random_number_in_cell(x, grid_spacing::T) where {T} =
    x + (rand() - T(1//2)) * grid_spacing
function random_position_in_cell(grid::AbstractLagrangianGrid, 𝐪)
    Δq = grid.grid_spacing
    return LagrangianCoordinate(
        random_number_in_cell(𝐪.x, Δq), 
        random_number_in_cell(𝐪.y, Δq), 
        random_number_in_cell(𝐪.z, Δq))
end


# chi = sqrt(x^2+y^2+z^2)
# scale_factor = scale_factor_of_chi(ic.cosmo, chi)
# return LatticeLocation(scale_factor, chi, x, y, z, i, j, k)


# struct LatticeLocation{T,LT}
#     x::LT     # first lagrangian coordinate
#     y::LT     # second lagrangian coordinate
#     z::LT     # third lagrangian coordinate
#     i::Int    # array index for parent
#     j::Int    # array index for parent
#     k::Int    # array index for parent
# end


struct ICFieldWebsky{T, LPT, AA, G, QITP}
    grid::G
    field::AA
    lagrangian_interp::QITP  # access field with Lagrangian coordinate 𝐪
end



"""
    ICFieldWebsky(::Type{LPT}, grid_spacing, cosmo, field)

Returns a wrapper around an array `field` in Lagrangian coordinates, defined by 
a kind of Lagrangian perturbation theory `LPT`, grid spacing, and cosmo. For 
Websky, this kind of field is typically something like ``\\delta_0(\\vec{q})``, 
the Lagrangian density at ``z=0``.
"""
function ICFieldWebsky(::Type{LPT}, grid::G, 
        field::AA, interp_type=BSpline(Quadratic(Periodic(OnCell())))
        ) where {T, TL, LPT, AA, G<:LagrangianGridWebsky{T,TL}}

    itp = interpolate(field, interp_type)
    sitp = scale(itp, grid.q_axes...)
    
    return ICFieldWebsky{T, LPT, AA, G, typeof(sitp)}(
        grid, field, sitp)
end

function getindex(ic::ICFieldWebsky{T}, 𝐪::LagrangianCoordinate{LT}) where {T, LT}
    x = 𝐪.x ≥ zero(LT) ? 𝐪.x : (𝐪.x + ic.grid.box_sizes[1])
    y = 𝐪.y ≥ zero(LT) ? 𝐪.y : (𝐪.y + ic.grid.box_sizes[2])
    z = 𝐪.z ≥ zero(LT) ? 𝐪.z : (𝐪.z + ic.grid.box_sizes[3])
    return T(ic.lagrangian_interp(x, y, z))
end
