using PythonCall

const pyccl = PythonCall.pynew()

function __init__()
    PythonCall.pycopy!(pyccl, pyimport("pyccl"))
end

struct CCLCosmology{T}
    p::Py
end

"""
Omega_c=0.2589, Omega_b=0.0486, h=0.6774, sigma8=0.8159,
        n_s=0.9667, transfer_function="boltzmann_camb"
"""
function CCLCosmology(T; kwargs...)
    c = pyccl.Cosmology(;kwargs...)
    return CCLCosmology{T}(c)
end

angular_diameter_distance(c::CCLCosmology{T}, a) where T = 
    pyconvert(T, pyccl.background.angular_diameter_distance(c.p, a))
luminosity_distance(c::CCLCosmology{T}, a) where T = 
    pyconvert(T, pyccl.background.luminosity_distance(c.p, a))
comoving_radial_distance(c::CCLCosmology{T}, a) where T = 
    pyconvert(T, pyccl.background.comoving_radial_distance(c.p, a))
growth_factor(c::CCLCosmology{T}, a) where T = 
    pyconvert(T, pyccl.background.growth_factor(c.p, a))
h_over_h0(c::CCLCosmology{T}, a) where T = 
    pyconvert(T, pyccl.background.h_over_h0(c.p, a))
scale_factor_of_chi(c::CCLCosmology{T}, chi) where T = 
    pyconvert(T, pyccl.background.scale_factor_of_chi(c.p, chi))


# linear_interpolation(xs, A)
struct InterpolatedCosmology{T, ITP}
    growth_factor::ITP
    scale_factor_of_chi::ITP
end

function InterpolatedCosmology(T, cosmology; amin=0.004, amax=1.0, N_grid=2048)
    agrid = LinRange(T(amin), T(amax), N_grid)
    growth_factor_grid = T[growth_factor(cosmology, a) for a in agrid]

    chi_min = comoving_radial_distance(cosmology, amax)
    chi_max = comoving_radial_distance(cosmology, amin)
    chi_grid = LinRange(T(chi_min), T(chi_max), N_grid)
    scale_factor_of_chi_grid = T[
        scale_factor_of_chi(cosmology, chi) for chi in chi_grid]

    growth_factor_itp = cubic_spline_interpolation(agrid, growth_factor_grid)
    scale_factor_of_chi_itp = cubic_spline_interpolation(
        chi_grid, scale_factor_of_chi_grid)
    
    return InterpolatedCosmology{T, typeof(growth_factor_itp)}(
        growth_factor_itp, scale_factor_of_chi_itp)
end

growth_factor(c::InterpolatedCosmology, a) = c.growth_factor(a)
scale_factor_of_chi(c::InterpolatedCosmology, chi) = c.scale_factor_of_chi(chi)
