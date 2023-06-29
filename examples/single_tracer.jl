using LagrangianPerturbationTheory, Unitful, UnitfulAstro, HDF5, PoissonRandom
import PythonPlot; const plt = PythonPlot
using QuadGK
icdir = "/fs/lustre/scratch/zack/ICs/"
files = (den=joinpath(icdir, "Fvec_7700Mpc_n6144_nb30_nt16_no768.h5"),
         sx=joinpath(icdir, "sx1_7700Mpc_n6144_nb30_nt16_no768.h5"),
         sy=joinpath(icdir, "sy1_7700Mpc_n6144_nb30_nt16_no768.h5"),
         sz=joinpath(icdir, "sz1_7700Mpc_n6144_nb30_nt16_no768.h5"))

ccl = CCLCosmology(Float32;
    Omega_c=0.2589, Omega_b=0.0486, h=0.6774, sigma8=0.8159,
    n_s=0.9667, transfer_function="boltzmann_camb")
cosmo = InterpolatedCosmology(ccl)
δ₀, Ψ₀_x, Ψ₀_y, Ψ₀_z = read_websky_ics(
    files.den, files.sx, files.sy, files.sz, 7700.0f0u"Mpc", cosmo)

dV = δ₀.grid_spacing^3

##
𝐪 = lattice_location(δ₀, 191, 191, 191, -1, -1, -1)
𝐪.x, 𝐪.y, 𝐪.z

##

massdef = CCLMassDef(200, "matter")
hmf = CCLMassFuncTinker08(ccl, massdef)

struct TopHatMassBinTracer{T, MT, BT}
    nbar::T              # object density
    M_min::MT            # top-hat bin mass minimum
    M_max::MT            # top-hat bin mass maximum
    bias_lagrangian::BT  # lagrangian bias
end

M_min, M_max = (10.0^11.0)u"Msun", (10.0^11.1)u"Msun"
n̄, n̄_err = quadgk(m -> dndm(hmf, m, 1.0), M_min, M_max, rtol=1e-8)
b⁽¹⁾ᴱ = halo_bias(tinker10, 1e7u"Msun", 1.0)  # tinker
b⁽¹⁾ᴸ = b⁽¹⁾ᴱ - 1
tb = TopHatMassBinTracer(n̄, M_min, M_max, b⁽¹⁾ᴸ)

##
function draw_tracer(δ₀, 𝚿⁽¹⁾₀, n̄, b⁽¹⁾ᴸ, octants=(-1,0))
    grid, cosmo = δ₀.grid, δ₀.grid.cosmo
    nx, ny, nz = size(δ₀.field)
    halos = SVector{3, TL}[]
    for k in 0:(nz-1), j in 0:(ny-1), i in 0:(nx-1)
        for oi in octants, oj in octants, ok in octants
            𝐪 = lagrangian_coordinate(grid, i, j, k, oi, oj, ok)  # grid cell 𝐪
            a = scale_factor(grid, 𝐪)  # a(𝐪) depends on normalization of grid
            D = growth_factor(cosmo, a)
            δ⁽¹⁾ᴸ = D * δ₀[𝐪]  # Lagrangian density at grid cell
            N = pois_rand(n̄ * (1 + b⁽¹⁾ᴸ * δ⁽¹⁾ᴸ) * grid.dV)  # draw halo number

            for halo_index in 1:N
                𝐪ₕ = random_position_in_cell(grid, 𝐪) # distribute halos in cell
                𝚿⁽¹⁾ₕ = 𝚿⁽¹⁾₀[𝐪ₕ]              # interpolate displacement
                aₕ = scale_factor(grid, 𝐪ₕ)    # a(𝐪) depends on grid norm
                Dₕ = growth_factor(cosmo, aₕ)  # D(a) at the random 𝐪ₕ in cell
                𝐱ₕ = 𝐪ₕ + Dₕ * 𝚿⁽¹⁾ₕ           # halo Eulerian position
                push!(halos, 𝐱ₕ)
            end
        end
    end
    return halos
end


##

##
