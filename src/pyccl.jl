using PythonCall

const pyccl = PythonCall.pynew()

function __init__()
    PythonCall.pycopy!(pyccl, pyimport("pyccl"))
end

struct CCLCosmology
    p::Py
end

"""
Omega_c=0.2589, Omega_b=0.0486, h=0.6774, sigma8=0.8159,
        n_s=0.9667, transfer_function="boltzmann_camb"
"""
function CCLCosmology(;kwargs...)
    c = pyccl.Cosmology(;kwargs...)
    return CCLCosmology(c)
end

growth_factor(c::CCLCosmology, a) = pyccl.background.growth_factor(c.p, a)
h_over_h0(c::CCLCosmology, a) = pyccl.background.h_over_h0(c.p, a)

