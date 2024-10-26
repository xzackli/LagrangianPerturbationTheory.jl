using LagrangianPerturbationTheory, Unitful, UnitfulAstro, HDF5, PoissonRandom
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
import PythonPlot; const plt = PythonPlot
using NPZ
plt.clf()
# plt.imshow(δ₀.field[20:60, 20:60,  4])
# xv = Ψ₀_x.field
xv = h5open(joinpath(icdir, "sx1_7700Mpc_n6144_nb30_nt16_no768.h5"), "r") do f read(f, "data") end

N = 768

##

# 𝐪 = lattice_location(δ₀, 3, 4, 5, 0, 0, 0)  # coordinate and inds
# D = growth_factor(cosmo, 𝐪.a)
# 𝚿⁽¹⁾₀ = (Ψ₀_x[𝐪], Ψ₀_y[𝐪], Ψ₀_z[𝐪]) .* u"Mpc"
# δ⁽¹⁾₀ = δ₀[𝐪]

# 1LPT solution, 𝐱 = 𝐪 .+ D .* 𝚿⁽¹⁾₀
# δ⁽¹⁾ᴸ = D * δ⁽¹⁾₀


# M₁, M₂ = 12.0f0, 12.1f0
# integrate hmf to get n̄ of redshift
# n̄ = 0.0001f0 / 1u"Mpc^3"
# b⁽¹⁾ᴱ = 3.0  # tinker
# b⁽¹⁾ᴸ = b⁽¹⁾ᴱ - 1

N = pois_rand(n̄ * (1 + b⁽¹⁾ᴸ * δ⁽¹⁾ᴸ) * dV)

##
massdef = CCLMassDef(200, "matter")
hmf = CCLMassFuncTinker08(ccl, massdef)

##
import PythonPlot
plt.close()
plt.hist([pois_rand(n̄ * (1 + b⁽¹⁾ᴸ * δ⁽¹⁾ᴸ) * dV) for i in 1:100])
plt.gcf()

##
masses = 10 .^ (10:0.1:15) .* u"Msun"

plt.close()
plt.plot(ustrip.((u"Msun",), masses),
    [ustrip(u"Mpc^(-3)* Msun^(-1)", dndm(hmf, m, 1.0)) for m in masses])
plt.xlabel(raw"Mass [$M_{\odot}$]")
plt.ylabel(raw"$dn/dM_h$ $[Mpc^{-3} M^{-1}_{\odot}]$")
plt.xscale("log"); plt.yscale("log")
plt.gcf()


##
using QuadGK
integral, err = quadgk(m -> dndlogm(hmf, m, 1.0), 11.0, 11.1, rtol=1e-8)

##
struct TophatTracerBin{T, MT, BT}
    nbar::T              # object density
    M_min::MT            # top-hat bin mass minimum
    M_max::MT            # top-hat bin mass maximum
    bias_lagrangian::BT  # lagrangian bias
end

##
M_min, M_max = (10.0^11.0)u"Msun", (10.0^11.1)u"Msun"
n̄, n̄_err = quadgk(m -> dndm(hmf, m, 1.0), M_min, M_max, rtol=1e-8)
b⁽¹⁾ᴱ = halo_bias(tinker10, 1e7u"Msun", 1.0)  # tinker
b⁽¹⁾ᴸ = b⁽¹⁾ᴱ - 1
ttb = TophatTracerBin(n̄, M_min, M_max)


##

# LagrangianPerturbationTheory.pyccl.halos.hbias.HaloBiasTinker10(ccl.p, massdef.p)
tinker10 = CCLHaloBiasTinker10(ccl, massdef)

plt.close()
plt.plot(ustrip.((u"Msun",), masses),
    [halo_bias(tinker10, m, 1.0) for m in masses], label="Tinker10")
plt.xlabel(raw"Mass [$M_{\odot}$]")
plt.ylabel(raw"$b_h(M)$")
plt.xscale("log")
plt.legend()
plt.gcf()


##
plt.close()
plt.plot(ustrip.((u"Msun",), masses),
    [ustrip(u"Mpc^(-3)* Msun^(-1)", dndm(hmf, m, 1.0)) for m in masses],)
plt.xlabel(raw"Mass [$M_{\odot}$]")
plt.ylabel(raw"$dn/dM_h$ $[Mpc^{-3} M^{-1}_{\odot}]$")
plt.xscale("log"); plt.yscale("log")
plt.gcf()

##
using StaticArrays
# i=3; j=4; k=5

# integrate hmf to get n̄ of redshift\
b⁽¹⁾ᴱ = halo_bias(tinker10, 1e7u"Msun", 1.0)  # tinker
b⁽¹⁾ᴸ = b⁽¹⁾ᴱ - 1

##
xx = Ref{Any}()


using LagrangianPerturbationTheory: LatticeLocation
randomcellpos(x, grid_spacing::Q) where {T, Q<:Quantity{T}} =
    x + (rand() - T(1//2)) * grid_spacing
randomcellpos(x, grid_spacing::T) where {T} =
    x + (rand() - T(1//2)) * grid_spacing
random_cell_position(qx, qy, qz, Δx) where T = 
    SA[randomcellpos(qx, Δx), randomcellpos(qy, Δx), randomcellpos(qz, Δx)]

##
# plt.clf()
# plt.hist([random_coord_in_cell(1.0, 1.0) for rr in 1:100],
#     range=(-1.5,4.5), bins=6)
# plt.gcf()


##
function draw_tracer(δ₀::IC, Ψ₀_x::IC, Ψ₀_y::IC, Ψ₀_z::IC, dV, n̄, b⁽¹⁾ᴸ) where {
        T, LPT, C, AA, TL, IC<:InitialConditionsWebsky{T, LPT, C, AA, TL}}

    cosmo = δ₀.cosmo
    nx, ny, nz = size(δ₀.field, 1), size(δ₀.field, 2), size(δ₀.field, 3)
    octants = (-1, 0)
    halos = SVector{3, TL}[]
    for k in 0:(nz-1), j in 0:(ny-1), i in 0:(nx-1)
        for oi in octants, oj in octants, ok in octants
            𝐪 = lattice_location(δ₀, i, j, k, oi, oj, ok)
            D = growth_factor(δ₀.cosmo, 𝐪.a)
            δ⁽¹⁾ᴸ = D * δ₀[𝐪]
            # 𝚿⁽¹⁾₀ = SA[Ψ₀_x[𝐪]*u"Mpc", Ψ₀_y[𝐪]*u"Mpc", Ψ₀_z[𝐪]*u"Mpc"] 
            N = pois_rand(n̄ * (1 + b⁽¹⁾ᴸ * δ⁽¹⁾ᴸ) * dV) 

            for halo_i in 1:N
                𝐢ₕ = random_cell_position(i, j, k, one(T)) # rand pos in index units
                𝚿⁽¹⁾ₕ = SA[Ψ₀_x[𝐢ₕ]*u"Mpc", Ψ₀_y[𝐢ₕ]*u"Mpc", Ψ₀_z[𝐢ₕ]*u"Mpc"]
                a = LagrangianPerturbationTheory.pos_to_scale_factor(δ₀, 𝚿⁽¹⁾ₕ)
                D = growth_factor(δ₀.cosmo, a)
                𝐪ₕ = LagrangianPerturbationTheory.index_to_pos(δ₀, 𝐢ₕ, oi, oj, ok)
                push!(halos, 𝐪ₕ + D * 𝚿⁽¹⁾ₕ)
            end
        end
    end
    return halos
end

##
using BenchmarkTools
# @btime draw_tracer($δ₀, $Ψ₀_x, $Ψ₀_y, $Ψ₀_z, $dV, $n̄, $b⁽¹⁾ᴸ)
halos = draw_tracer(δ₀, Ψ₀_x, Ψ₀_y, Ψ₀_z, dV, n̄, b⁽¹⁾ᴸ)

##
xs = [ustrip(u"Mpc", h[2])  for h in halos] 
ys = [ustrip(u"Mpc", h[3]) for h in halos]

##
# import PythonPlot; const plt = PythonPlot
plt.clf()
plt.hexbin(xs, ys)

# for i in 20:40
#     plt.axvline(i * ustrip(u"Mpc", δ₀.grid_spacing), lw=0.5)
#     plt.axhline(i * ustrip(u"Mpc", δ₀.grid_spacing), lw=0.5)
# end

plt.colorbar()
plt.gca().set_aspect("equal")
plt.xlabel("x [Mpc]"); plt.ylabel("y [Mpc]")
plt.gcf()

##
using FFTW
plt.clf()
# plt.imshow(δ₀.field[20:60, 20:60,  4])
# plt.imshow(Ψ₀_z.field[10:60, 10:60,  4])
# plt.imshow(Ψ₀_y.field[20:60, 20:60,  4])
plt.imshow(abs.(ifft(Ψ₀_z.field[:, :, 21]))[20:80,20:80])
plt.colorbar()
plt.gcf()

##
