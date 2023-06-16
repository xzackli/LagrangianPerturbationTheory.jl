using LagrangianPerturbationTheory
const LPT = LagrangianPerturbationTheory
using HDF5

cosmology = CCLCosmology(Float64;
    Omega_c=0.2589, Omega_b=0.0486, h=0.6774, sigma8=0.8159,
    n_s=0.9667, transfer_function="boltzmann_camb")

ic_dir = "/fs/lustre/scratch/zack/ICs/"
boxsize = 7700. # comoving Mpc
N=192
grid_spacing = boxsize / N
icsfname="Fvec_7700Mpc_n6144_nb30_nt16_ks192"
delta = h5open(ic_dir * icsfname * ".h5", "r") do file
    read(file, "delta")
end

##

delta2 = LPT.load_example_ics()

sum(abs.(delta .- delta2))

##

i=3; j=4; k=5

for xd in (-1,0)
    x = (i+0.5)*grid_spacing + xd * boxsize
    for yd in (-1,0)
        y = (j+0.5)*grid_spacing + yd * boxsize
        for zd in (-1,0)
            z = (k+0.5)*grid_spacing + zd * boxsize

            chi = sqrt(x^2+y^2+z^2)
            scale_factor = LPT.scale_factor_of_chi(cosmology, chi)
            redshift = 1 / scale_factor - 1
            δ₀ = delta[k+1,j+1,i+1]  # (py → jl) reverses dimensions, index+1
            delta_linear = δ₀ * LPT.growth_factor(cosmology, scale_factor)
            println((x, y, z,xd,yd,zd,delta_linear))
        end
    end
end
