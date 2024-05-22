using PythonCall

const pyccl = PythonCall.pynew()

function __init__()
    PythonCall.pycopy!(pyccl, pyimport("pyccl"))
end


abstract type AbstractCosmologyBackground{T} end


"""
    growth_factor(cosmo, a)

Returns the linear growth factor ``D(a)`` as a function of cosmo `cosmo`
and scale factor ``a``.
"""
function growth_factor end

"""
    scale_factor_of_chi(cosmo, chi)

Returns the scale factor ``a(\\chi)`` as a function of comoving radial
distance ``\\chi`` (Mpc).
"""
function scale_factor_of_chi end


struct CCLCosmology{T} <: AbstractCosmologyBackground{T}
    p::Py
end

"""
    CCLCosmology(T; kwargs...)

Returns a wrapped [pyccl cosmology][1] object with output type `T` and 
cosmological parameters as keyword arguments. For example, a standard Planck
2018 cosmology can be configured with
```julia
cclcosmo = CCLCosmology(Float32; 
    Omega_c=0.2589, Omega_b=0.0486, h=0.6774, sigma8=0.8159,
    n_s=0.9667, transfer_function="boltzmann_camb")
```
[1]: https://ccl.readthedocs.io/en/latest/api/pyccl.core.html#pyccl.core.Cosmology
"""
function CCLCosmology(::Type{T}; kwargs...) where T
    c = pyccl.Cosmology(;kwargs...)
    return CCLCosmology{T}(c)
end


growth_factor(c::CCLCosmology{T}, a::Real) where T = 
    pyconvert(T, pyccl.background.growth_factor(c.p, a))
scale_factor_of_chi(c::CCLCosmology{T}, chi::Quantity) where T = 
    pyconvert(T, pyccl.background.scale_factor_of_chi(c.p, ustrip(u"Mpc", chi)))

angular_diameter_distance(c::CCLCosmology{T}, a::Real) where T = 
    pyconvert(T, pyccl.background.angular_diameter_distance(c.p, a))
luminosity_distance(c::CCLCosmology{T}, a::Real) where T = 
    pyconvert(T, pyccl.background.luminosity_distance(c.p, a))
comoving_radial_distance(c::CCLCosmology{T}, a::Real) where T = 
    pyconvert(T, pyccl.background.comoving_radial_distance(c.p, a)) * 1u"Mpc"
h_over_h0(c::CCLCosmology{T}, a::Real) where T = 
    pyconvert(T, pyccl.background.h_over_h0(c.p, a))

struct InterpolatedCosmology{T, ITP1, ITP2}  <: AbstractCosmologyBackground{T}
    growth_factor::ITP1
    scale_factor_of_chi::ITP2
end


"""
    InterpolatedCosmology(cosmo)

Given a different cosmology object (at the moment, only [`CCLCosmology`](@ref)), 
construct a set of interpolators that provide the quantities needed for LPT.
```julia
cclcosmo = CCLCosmology(Float32; 
    Omega_c=0.2589, Omega_b=0.0486, h=0.6774, sigma8=0.8159,
    n_s=0.9667, transfer_function="boltzmann_camb")
cosmo = InterpolatedCosmology(cclcosmo)
```
"""
function InterpolatedCosmology(cosmo::AbstractCosmologyBackground{T}; 
                               amin=0.004, amax=1.0, N_grid=128) where T
    agrid = LinRange(T(amin), T(amax), N_grid)
    growth_factor_grid = T[growth_factor(cosmo, a) for a in agrid]

    chi_min = comoving_radial_distance(cosmo, amax)
    chi_max = comoving_radial_distance(cosmo, amin)
    chi_grid = LinRange(chi_min, chi_max, N_grid)
    scale_factor_of_chi_grid = [
        scale_factor_of_chi(cosmo, chi) for chi in chi_grid]

    growth_factor_itp = cubic_spline_interpolation(agrid, growth_factor_grid)
    scale_factor_of_chi_itp = cubic_spline_interpolation(
        chi_grid, scale_factor_of_chi_grid)
    
    return InterpolatedCosmology{
        T, typeof(growth_factor_itp), typeof(scale_factor_of_chi_itp)}(
        growth_factor_itp, scale_factor_of_chi_itp)
end

growth_factor(c::InterpolatedCosmology, a::Real) = c.growth_factor(a)
scale_factor_of_chi(c::InterpolatedCosmology, chi::Quantity) = c.scale_factor_of_chi(chi)


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
            mass_def=mass_def.p, mass_def_strict=mass_def_strict), cosmo.p)
end

"""
as in pyccl, returns number density per Mpc³ per log₁₀(Mₕ/Msun)
"""
function dndlogm(hmf::MassFunc{T}, log₁₀M, scale_factor) where T
    pyconvert(T, hmf.hmf(
        hmf.cosmopy, T(10)^log₁₀M, scale_factor)) * u"Mpc^(-3)"
end

function dndm(hmf::MassFunc{T}, mass, scale_factor) where T
    pyconvert(T, hmf.hmf(
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
        mass_def=mass_def.p, mass_def_strict=mass_def_strict), cosmo.p)
end

function halo_bias(hb::HaloBias{T}, mass, scale_factor) where T
    pyconvert(T, hb.hb(
        hb.cosmopy, ustrip(u"Msun", mass), scale_factor))
end
