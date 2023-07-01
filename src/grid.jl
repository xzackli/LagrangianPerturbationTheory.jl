abstract type AbstractLagrangianGrid end
struct LagrangianGridWebsky{T,TL,C,R,TV} <: AbstractLagrangianGrid 
    cosmo::C
    grid_spacing::TL
    box_sizes::NTuple{3,TL}
    q_axes::NTuple{3,R}  # axis grids for each Lagrangian coordinate
    dV::TV
end

function LagrangianGridWebsky(cosmo::C, grid_spacing::TL, 
        box_sizes, q_axes::NTuple{3,R}) where {T, C, R, TL<:Quantity{T}}
    dV = grid_spacing^3
    LagrangianGridWebsky{T,TL,C,R,typeof(dV)}(cosmo, grid_spacing, box_sizes, q_axes, dV)
end


function lagrangian_coordinate(grid, i, j, k, oi, oj, ok)
    Δq = grid.grid_spacing
    offset = Δq / 2
    box_size_x, box_size_y, box_size_z = grid.box_sizes
    x = (i-1) * Δq + offset + oi * box_size_x
    y = (j-1) * Δq + offset + oj * box_size_y
    z = (k-1) * Δq + offset + ok * box_size_z
    return SVector(x, y, z)
end

# websky grid is normalized so that |𝐪| is comoving distance
function scale_factor(grid::LagrangianGridWebsky, 𝐪::SVector{3})
    chi = √(𝐪.x^2 + 𝐪.y^2 + 𝐪.z^2)
    return scale_factor_of_chi(grid.cosmo, chi)
end

random_number_in_cell(x, grid_spacing::Q) where {T, Q<:Quantity{T}} =
    x + (rand(T) - T(1//2)) * grid_spacing
random_number_in_cell(x, grid_spacing::T) where {T} =
    x + (rand(T) - T(1//2)) * grid_spacing

function random_position_in_cell(grid::AbstractLagrangianGrid, 𝐪::SVector{3,T}) where T
    Δq = grid.grid_spacing
    return SVector(
        random_number_in_cell(𝐪.x, Δq), 
        random_number_in_cell(𝐪.y, Δq), 
        random_number_in_cell(𝐪.z, Δq))
end


struct ICFieldWebsky{T, LPT, TL, AA, G, QITP}
    grid::G
    field::AA
    lagrangian_interp::QITP  # access field with Lagrangian coordinate 𝐪
end

"""
    ICFieldWebsky

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
    
    return ICFieldWebsky{T, LPT, TL, AA, G, typeof(sitp)}(
        grid, field, sitp)
end


ensure_periodic_positive(x, L) = x < zero(x) ? (x + L) : x
box_size(g::AbstractLagrangianGrid, i) = g.box_sizes[i]

function getindex(ic::ICFieldWebsky{T}, 𝐪::SVector{3}) where T
    x = ensure_periodic_positive(𝐪.x, box_size(ic.grid,1))
    y = ensure_periodic_positive(𝐪.y, box_size(ic.grid,2))
    z = ensure_periodic_positive(𝐪.z, box_size(ic.grid,3))
    return T(ic.lagrangian_interp(x, y, z))
end

function Base.show(io::IO, ic::ICFieldWebsky{T, LPT, TL, AA, G, QITP}) where {T, LPT, TL, AA, G, QITP}
    expr = "ICFieldWebsky{$T, $LPT, $TL, $AA, ...}"
    print(io, expr)
end


# for three at once
struct DisplacementICFieldWebsky{IC}
    sx::IC
    sy::IC
    sz::IC
end

function getindex(ics::DisplacementICFieldWebsky{IC}, 
                  𝐪::SVector{3,TL}) where {T,LPT,TL,IC<:ICFieldWebsky{T, LPT, TL}}
    x = ensure_periodic_positive(𝐪.x, box_size(ics.sx.grid,1))
    y = ensure_periodic_positive(𝐪.y, box_size(ics.sy.grid,2))
    z = ensure_periodic_positive(𝐪.z, box_size(ics.sz.grid,3))
    return SVector(
        T(ics.sx.lagrangian_interp(x, y, z)), 
        T(ics.sy.lagrangian_interp(x, y, z)), 
        T(ics.sz.lagrangian_interp(x, y, z)))
end


function Base.show(io::IO, ic::DisplacementICFieldWebsky{IC}) where {T, LPT, TL, AA, G, QITP, IC<:ICFieldWebsky{T, LPT, TL, AA, G, QITP}}
    expr = "DisplacementICFieldWebsky containing \n ⋅ 3x ICFieldWebsky{$T, $LPT, $TL, $AA, ...}"
    print(io, expr)
end
