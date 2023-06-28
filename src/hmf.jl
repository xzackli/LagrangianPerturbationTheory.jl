

struct CCLMassDef
    p::Py
end

function CCLMassDef(Δ, rho_type)
    massdef = pyccl.halos.massdef.MassDef(Δ, rho_type)
    return CCLMassDef(massdef)
end

struct MassFunc{T}
    hmf::Py
    cosmopy::Py
end


function CCLMassFuncTinker08(cosmo::CCLCosmology{T}, 
        mass_def::CCLMassDef=CCLMassDef(200, "matter"), 
        mass_def_strict=true) where T
    return MassFunc{T}(
        pyccl.halos.hmfunc.MassFuncTinker08(
            cosmo.p, mass_def.p, mass_def_strict), cosmo.p)
end


"""
as in pyccl, returns number density per Mpc³ per log₁₀(Mₕ/Msun)
"""
function dndlogm(hmf::MassFunc{T}, log₁₀M, scale_factor) where T
    pyconvert(T, hmf.hmf.get_mass_function(
        hmf.cosmopy, T(10)^log₁₀M, scale_factor)) * u"Mpc^(-3)"
end

function dndm(hmf::MassFunc{T}, mass, scale_factor) where T
    pyconvert(T, hmf.hmf.get_mass_function(
        hmf.cosmopy, ustrip(u"Msun", mass), scale_factor)) * u"Mpc^(-3)" / (
            T(mass) * log(T(10))) 
end


struct HaloBias{T}
    hb::Py
    cosmopy::Py
end

function CCLHaloBiasTinker10(cosmo::CCLCosmology{T}, 
        mass_def::CCLMassDef=CCLMassDef(200, "matter"), 
        mass_def_strict=true) where T
    return HaloBias{T}(pyccl.halos.hbias.HaloBiasTinker10(
        cosmo.p, mass_def.p, mass_def_strict), cosmo.p)
end

function halo_bias(hb::HaloBias{T}, mass, scale_factor) where T
    pyconvert(T, hb.hb.get_halo_bias(
        hb.cosmopy, ustrip(u"Msun", mass), scale_factor))
end


# function mean_number(hmf, z, logM₁, logM₂)
    



