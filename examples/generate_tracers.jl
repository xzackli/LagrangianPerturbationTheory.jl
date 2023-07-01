using LagrangianPerturbationTheory, Unitful, UnitfulAstro, HDF5, PoissonRandom, StaticArrays
using BenchmarkTools, ThreadPools, JLD2, FileIO
import PythonPlot; const plt = PythonPlot
icdir = "/fs/lustre/scratch/zack/ICs/"
files = (den=joinpath(icdir, "Fvec_7700Mpc_n6144_nb30_nt16_no768.h5"),
         sx=joinpath(icdir, "sx1_7700Mpc_n6144_nb30_nt16_no768.h5"),
         sy=joinpath(icdir, "sy1_7700Mpc_n6144_nb30_nt16_no768.h5"),
         sz=joinpath(icdir, "sz1_7700Mpc_n6144_nb30_nt16_no768.h5"))

ccl = CCLCosmology(Float32;
    Omega_c=0.2589, Omega_b=0.0486, h=0.6774, sigma8=0.8159,
    n_s=0.9667, transfer_function="boltzmann_camb")
cosmo = InterpolatedCosmology(ccl)
massdef = CCLMassDef(200, "matter")
hmf_tinker08 = CCLMassFuncTinker08(ccl, massdef)
bias_tinker10 = CCLHaloBiasTinker10(ccl, massdef)
δ₀, 𝚿₀ = read_websky_ics(files.den, files.sx, files.sy, files.sz, cosmo, 7700.f0u"Mpc")

##
function draw_tracer_threaded(δ₀::ICFieldWebsky{T, LPT, TL}, 𝚿₀, tracer, octant) where {T, LPT, TL}
    tasks = map(axes(δ₀.field,1)) do i
        Threads.@spawn begin
            halo_positions = SVector{3, TL}[]
            draw_tracer!(halo_positions, δ₀, 𝚿₀, tracer, 
                i:i, axes(δ₀.field,2), axes(δ₀.field,3), octant)
            return halo_positions    # return one octant at a time
        end
    end
    return fetch.(tasks)
end

function save_positions(filename, δ₀, 𝚿₀, tracer)
    tasks = map(LagrangianPerturbationTheory.FULL_WEBSKY_OCTANTS) do octant
        Threads.@spawn draw_tracer_threaded(δ₀, 𝚿₀, tracer, octant)
    end
    positions = package_tracer_positions(fetch.(tasks))
    save(filename, "positions", positions)
end

function save_multiple_masses(logmasses, δ₀::ICFieldWebsky{T, LPT, TL}, 𝚿₀, hmf, bias) where {T, LPT, TL}
    mass_indices = 1:(length(logmasses)-1)
    tracers = [TopHatMassBinTracer((10^logmasses[iₘ])u"Msun", (10^logmasses[iₘ+1])u"Msun", 
        hmf, bias, 0.2f0:0.01f0:1f0) for iₘ in mass_indices]
    for i_m in mass_indices
        tracer = tracers[i_m]
        save_positions("/fs/lustre/scratch/zack/ICs/low_mass_halos/lowres_positions_" * 
            "$(logmasses[i_m])_$(logmasses[i_m+1]).jld2", δ₀, 𝚿₀, tracer)
        GC.gc(true); GC.gc(false)  # force full GC collection
    end
end
@time save_multiple_masses(11:0.01f0:12, δ₀, 𝚿₀, hmf_tinker08, bias_tinker10)

##
# mass_indices = 1:(length(11:0.05f0:12)-1)
# tracer = TopHatMassBinTracer((10^12.4f0)u"Msun", (10^12.5f0)u"Msun", 
#     hmf_tinker08, bias_tinker10, 0.2f0:0.01f0:1f0)
# ts = [draw_tracer_threaded(δ₀, 𝚿₀, tracer, (0,0,0))]
# package_tracer_positions(ts)

##

# @time save_multiple_masses(12:0.1f0:14, δ₀, 𝚿₀, hmf_tinker08, bias_tinker10)

        # for octant in LagrangianPerturbationTheory.FULL_WEBSKY_OCTANTS
        #     append!(positions, draw_tracer_threaded(δ₀, 𝚿₀, tracer, octant))
        # end

##

##
# halo_positions = SVector{3, typeof(δ₀.grid.grid_spacing)}[]
# @time draw_tracer!(halo_positions, δ₀, 𝚿₀, tracer, 
#     1:16, axes(δ₀.field,2), axes(δ₀.field,3))


# halo_positions = SVector{3, typeof(δ₀.grid.grid_spacing)}[]
# M_min, M_max = (10.0^12.0)u"Msun", (10.0^12.1)u"Msun"
# tracer = TopHatMassBinTracer(M_min, M_max, hmf_tinker08, bias_tinker10, 0.2f0:0.01f0:1.0f0)
# @time draw_tracer!(halo_positions, δ₀, 𝚿₀, tracer, 
#     1:2, axes(δ₀.field,2), axes(δ₀.field,3), (0,0,0))


# # @profview draw_tracer!(halo_positions, δ₀, 𝚿₀, tracer, (-1, 0))


# ##
# plt.clf()
# a_grid = 0.1:0.01:1.0
# plt.plot(a_grid, [ustrip(bias_lagrangian(tracer, a)) for a in a_grid] )
# plt.gcf()
