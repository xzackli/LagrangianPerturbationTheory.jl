using LagrangianPerturbationTheory, Unitful, UnitfulAstro, StaticArrays, JLD2, FileIO, Printf
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
    save(filename, "positions", positions)  # WebSky reverses the dims: saves as z, y, x
end

function save_multiple_masses(logmasses, δ₀::ICFieldWebsky{T, LPT, TL}, 𝚿₀, hmf, bias) where {T, LPT, TL}
    mass_indices = 1:(length(logmasses)-1)
    scale_factor_grid = 0.2f0:0.01f0:1f0
    tracers = [TopHatMassBinTracer((10^logmasses[iₘ])u"Msun", (10^logmasses[iₘ+1])u"Msun", 
        hmf, bias, scale_factor_grid) for iₘ in mass_indices]
    for i_m in mass_indices
        tracer = tracers[i_m]
        m_str = @sprintf("%.2f",logmasses[i_m]) * "_" * @sprintf("%.2f",logmasses[i_m+1])
        save_positions("/fs/lustre/scratch/zack/ICs/low_mass_halos/" * 
            "lowres_positions_$(m_str).jld2", δ₀, 𝚿₀, tracer)
        GC.gc(true); GC.gc(false)  # force full GC collection
    end
end
@time save_multiple_masses(11:0.01f0:12, δ₀, 𝚿₀, hmf_tinker08, bias_tinker10)
