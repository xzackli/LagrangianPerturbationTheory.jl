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
