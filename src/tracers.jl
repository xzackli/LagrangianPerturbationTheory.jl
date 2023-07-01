
"""
Octants are indexed -1 and 0.
"""
function draw_tracer!(halo_positions, δ₀::ICFieldWebsky{T, LPT, TL}, 𝚿⁽¹⁾₀, tracer, 
                      i_range, j_range, k_range, octant) where {T, LPT, TL}
    grid, cosmo = δ₀.grid, δ₀.grid.cosmo
    oi, oj, ok = octant
    for k in k_range, j in j_range, i in i_range
        𝐪 = lagrangian_coordinate(grid, i, j, k, oi, oj, ok)  # Lagrangian 𝐪 of cell
        a = scale_factor(grid, 𝐪)           # grid normalization gives you a(𝐪)
        if tracer.a_min ≤ a ≤ tracer.a_max
            D = growth_factor(cosmo, a)         # growth factor at grid cell
            n̄ = mean_density(tracer, a)         # mean tracer density at cell center
            b⁽¹⁾ᴸ = bias_lagrangian(tracer, a)  # Lagrangian bias at cell center
            δ⁽¹⁾ᴸ = D * δ₀[𝐪]                   # Lagrangian density at grid cell
            N = pois_rand(n̄ * (1 + b⁽¹⁾ᴸ * δ⁽¹⁾ᴸ) * grid.dV)  # find N to Poisson draw
            for _ in 1:N                        # now generate N halos
                𝐪ₕ = random_position_in_cell(grid, 𝐪) # randomly distribute halo in cell
                𝚿⁽¹⁾ₕ = 𝚿⁽¹⁾₀[𝐪ₕ]               # interpolate the displacement
                aₕ = scale_factor(grid, 𝐪ₕ)     # a(𝐪) depends on grid normalization
                Dₕ = growth_factor(cosmo, aₕ)   # growth factor at the random 𝐪ₕ in cell
                𝐱ₕ = 𝐪ₕ + Dₕ * 𝚿⁽¹⁾ₕ            # halo Eulerian position
                push!(halo_positions, SVector(𝐱ₕ.x, 𝐱ₕ.y, 𝐱ₕ.z))
            end
        end
    end
end

const FULL_WEBSKY_OCTANTS = (
    (-1,-1,-1), (-1,-1,0), (-1,0,-1), (-1,0,0), (0,-1,-1), (0,-1,0), (0,0,-1), (0,0,0))

# draw_tracer!(halo_positions, δ₀, 𝚿⁽¹⁾₀, tracer) = draw_tracer!(
#     halo_positions, δ₀, 𝚿⁽¹⁾₀, tracer, 
#     axes(δ₀.field,1), axes(δ₀.field,2), axes(δ₀.field,3), FULL_WEBSKY_OCTANTS)

struct TopHatMassBinTracer{T, MT, ITP1 <: AbstractInterpolation, ITP2 <: AbstractInterpolation}
    a_min::T
    a_max::T
    M_min::MT              # top-hat bin mass minimum
    M_max::MT              # top-hat bin mass maximum
    density::ITP1             # object density
    bias_lagrangian::ITP2  # lagrangian bias
end

function build_massfunc_interpolator(hmf::MassFunc, a_grid, M_min, M_max, rtol=1e-6,
                                    interp_type=BSpline(Cubic(Line(OnGrid()))))
    nbar_grid = [
        quadgk(m -> dndm(hmf, m, a), M_min, M_max, rtol=rtol)[1] for a in a_grid]
    zero_val = zero(nbar_grid[begin])
    return extrapolate(scale(interpolate(nbar_grid, interp_type), a_grid), zero_val)
end

function build_lagrangian_bias_interpolator(hb::HaloBias, a_grid, M_min, M_max, 
                                    interp_type=BSpline(Cubic(Line(OnGrid()))))
    M̄ = sqrt(M_min * M_max)
    nbar_grid = [halo_bias(hb, M̄, a) - 1 for a in a_grid]
    return extrapolate(scale(interpolate(nbar_grid, interp_type), a_grid), Flat())
end


function TopHatMassBinTracer(M_min, M_max, hmf::MassFunc, hb::HaloBias, a_grid)
    nbar = build_massfunc_interpolator(hmf, a_grid, M_min, M_max)
    lag_bias = build_lagrangian_bias_interpolator(hb, a_grid, M_min, M_max)
    return TopHatMassBinTracer(
        minimum(a_grid), maximum(a_grid), M_min, M_max, nbar, lag_bias)
end

mean_density(tracer::TopHatMassBinTracer, a::Real) = tracer.density(a)
bias_lagrangian(tracer::TopHatMassBinTracer, a::Real) = tracer.bias_lagrangian(a)


const VECLIKE = Base.AbstractVecOrTuple

"""
    package_tracer_positions(vector_of_vectors_of_svector_positions)

Convert a triply-nested linearly-indexed of SVector{3,Quantity{T}} 
(unitful positions) into an Array{T,2} of (3, total_elements) to be stored to disk, 
without units. Lengths are converted to Mpc.
"""
function package_tracer_positions(
        vvv::V) where {T, Q<:Quantity{T}, V<:VECLIKE{VECLIKE{VECLIKE{SVector{3,Q}}}}}
    total_number_of_elements = 0
    for vv in vvv, v in vv
        total_number_of_elements += length(v)
    end
    result = zeros(T, (3, total_number_of_elements))
    result_index = 1
    for vv in vvv, v in vv, p in v
        result[:, result_index] .= ustrip.(u"Mpc", p)
        result_index += 1
    end
    return result
end
