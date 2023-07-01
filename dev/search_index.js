var documenterSearchIndex = {"docs":
[{"location":"api/#Index","page":"API","title":"Index","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"","category":"page"},{"location":"api/","page":"API","title":"API","text":"Modules = [LagrangianPerturbationTheory]","category":"page"},{"location":"api/#LagrangianPerturbationTheory.CCLCosmology-Union{Tuple{Type{T}}, Tuple{T}} where T","page":"API","title":"LagrangianPerturbationTheory.CCLCosmology","text":"CCLCosmology(T; kwargs...)\n\nReturns a wrapped [pyccl cosmology][1] object with output type T and  cosmological parameters as keyword arguments. For example, a standard Planck 2018 cosmology can be configured with\n\ncclcosmo = CCLCosmology(Float32; \n    Omega_c=0.2589, Omega_b=0.0486, h=0.6774, sigma8=0.8159,\n    n_s=0.9667, transfer_function=\"boltzmann_camb\")\n\n[1]: https://ccl.readthedocs.io/en/latest/api/pyccl.core.html#pyccl.core.Cosmology\n\n\n\n\n\n","category":"method"},{"location":"api/#LagrangianPerturbationTheory.ICFieldWebsky-Union{Tuple{G}, Tuple{AA}, Tuple{LPT}, Tuple{TL}, Tuple{T}, Tuple{Type{LPT}, G, AA}, Tuple{Type{LPT}, G, AA, Any}} where {T, TL, LPT, AA, G<:(LagrangianGridWebsky{T, TL})}","page":"API","title":"LagrangianPerturbationTheory.ICFieldWebsky","text":"ICFieldWebsky\n\nReturns a wrapper around an array field in Lagrangian coordinates, defined by  a kind of Lagrangian perturbation theory LPT, grid spacing, and cosmo. For  Websky, this kind of field is typically something like delta_0(vecq),  the Lagrangian density at z=0.\n\n\n\n\n\n","category":"method"},{"location":"api/#LagrangianPerturbationTheory.InterpolatedCosmology-Union{Tuple{LagrangianPerturbationTheory.AbstractCosmologyBackground{T}}, Tuple{T}} where T","page":"API","title":"LagrangianPerturbationTheory.InterpolatedCosmology","text":"InterpolatedCosmology(cosmo)\n\nGiven a different cosmology object (at the moment, only CCLCosmology),  construct a set of interpolators that provide the quantities needed for LPT.\n\ncclcosmo = CCLCosmology(Float32; \n    Omega_c=0.2589, Omega_b=0.0486, h=0.6774, sigma8=0.8159,\n    n_s=0.9667, transfer_function=\"boltzmann_camb\")\ncosmo = InterpolatedCosmology(cclcosmo)\n\n\n\n\n\n","category":"method"},{"location":"api/#LagrangianPerturbationTheory.dndlogm-Union{Tuple{T}, Tuple{LagrangianPerturbationTheory.MassFunc{T}, Any, Any}} where T","page":"API","title":"LagrangianPerturbationTheory.dndlogm","text":"as in pyccl, returns number density per Mpc³ per log₁₀(Mₕ/Msun)\n\n\n\n\n\n","category":"method"},{"location":"api/#LagrangianPerturbationTheory.draw_tracer!-Union{Tuple{TL}, Tuple{LPT}, Tuple{T}, Tuple{Any, ICFieldWebsky{T, LPT, TL}, Vararg{Any, 6}}} where {T, LPT, TL}","page":"API","title":"LagrangianPerturbationTheory.draw_tracer!","text":"Octants are indexed -1 and 0.\n\n\n\n\n\n","category":"method"},{"location":"api/#LagrangianPerturbationTheory.growth_factor","page":"API","title":"LagrangianPerturbationTheory.growth_factor","text":"growth_factor(cosmo, a)\n\nReturns the linear growth factor D(a) as a function of cosmo cosmo and scale factor a.\n\n\n\n\n\n","category":"function"},{"location":"api/#LagrangianPerturbationTheory.package_tracer_positions-Union{Tuple{V}, Tuple{Q}, Tuple{T}} where {T, Q<:(Unitful.Quantity{T}), V<:Union{Tuple{Vararg{Union{Tuple{Vararg{Union{Tuple{Vararg{StaticArraysCore.SVector{3, Q}}}, AbstractVector{<:StaticArraysCore.SVector{3, Q}}}}}, AbstractVector{<:Union{Tuple{Vararg{StaticArraysCore.SVector{3, Q}}}, AbstractVector{<:StaticArraysCore.SVector{3, Q}}}}}}}, AbstractVector{<:Union{Tuple{Vararg{Union{Tuple{Vararg{StaticArraysCore.SVector{3, Q}}}, AbstractVector{<:StaticArraysCore.SVector{3, Q}}}}}, AbstractVector{<:Union{Tuple{Vararg{StaticArraysCore.SVector{3, Q}}}, AbstractVector{<:StaticArraysCore.SVector{3, Q}}}}}}}}","page":"API","title":"LagrangianPerturbationTheory.package_tracer_positions","text":"package_tracer_positions(vector_of_vectors_of_svector_positions)\n\nConvert a triply-nested linearly-indexed of SVector{3,Quantity{T}}  (unitful positions) into an Array{T,2} of (3, total_elements) to be stored to disk,  without units. Lengths are converted to Mpc.\n\n\n\n\n\n","category":"method"},{"location":"api/#LagrangianPerturbationTheory.scale_factor_of_chi","page":"API","title":"LagrangianPerturbationTheory.scale_factor_of_chi","text":"scale_factor_of_chi(cosmo, chi)\n\nReturns the scale factor a(chi) as a function of comoving radial distance chi (Mpc).\n\n\n\n\n\n","category":"function"},{"location":"background_cosmology/","page":"Background Cosmology","title":"Background Cosmology","text":"CurrentModule = LagrangianPerturbationTheory","category":"page"},{"location":"background_cosmology/#Background-Cosmology","page":"Background Cosmology","title":"Background Cosmology","text":"","category":"section"},{"location":"background_cosmology/","page":"Background Cosmology","title":"Background Cosmology","text":"using LagrangianPerturbationTheory, Unitful, UnitfulAstro\n\n# example , ra and dec in radians, halo mass in M200c (Msun)\nra, dec, redshift, halo_mass = XGPaint.load_example_halos()\nprint(\"Number of halos: \", length(halo_mass))\n\ncosmo_ccl = CCLCosmology(Float32;\n    Omega_c=0.2589, Omega_b=0.0486, h=0.6774, sigma8=0.8159,\n    n_s=0.9667, transfer_function=\"boltzmann_camb\")\n\nmassdef = CCLMassDef(200, \"critical\")\nhmf = CCLMassFuncTinker08(cosmo_ccl, massdef)\n\nimport PythonPlot; const plt = PythonPlot\nmasses = 10 .^ (10:0.1:15) .* u\"Msun\"\nplt.close()\nplt.plot(ustrip.((u\"Msun\",), masses),\n    [ustrip(u\"Mpc^(-3)* Msun^(-1)\", dndm(hmf, m, 1.0)) for m in masses])\nplt.xlabel(raw\"Mass [$M_{\\odot}$]\")\nplt.ylabel(raw\"$dn/dM_h$ $[Mpc^{-3} M^{-1}_{\\odot}]$\")\nplt.xscale(\"log\"); plt.yscale(\"log\")\nplt.gcf()","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = LagrangianPerturbationTheory","category":"page"},{"location":"#LagrangianPerturbationTheory","page":"Home","title":"LagrangianPerturbationTheory","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for LagrangianPerturbationTheory.","category":"page"}]
}