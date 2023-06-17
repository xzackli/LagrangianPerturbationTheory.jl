using LagrangianPerturbationTheory, Unitful, UnitfulAstro, HDF5

ccl = CCLCosmology(Float32;
    Omega_c=0.2589, Omega_b=0.0486, h=0.6774, sigma8=0.8159,
    n_s=0.9667, transfer_function="boltzmann_camb")
cosmo = InterpolatedCosmology(ccl)
icfile = "/fs/lustre/scratch/zack/ICs/Fvec_7700Mpc_n6144_nb30_nt16_ks192.h5"
delta_array = h5open(icfile, "r") do file read(file, "delta") end

grid_spacing = 7700.0f0u"Mpc" / size(delta_array,1)
δ₀ =  InitialConditionsWebsky(FirstOrderLPT, grid_spacing, cosmo, delta_array)

##

δ₀.grid_spacing

##

loc = lattice_location(δ₀, 3, 4, 5, 0, 0, 0)
δ₁ = lattice_0(δ₀, 3, 4, 5, 0, 0, 0)
