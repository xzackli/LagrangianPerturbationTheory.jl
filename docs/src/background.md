```@meta
CurrentModule = LagrangianPerturbationTheory
```

# Background Cosmology

```@example bg
using LagrangianPerturbationTheory, Unitful, UnitfulAstro

# example , ra and dec in radians, halo mass in M200c (Msun)
ra, dec, redshift, halo_mass = XGPaint.load_example_halos()
print("Number of halos: ", length(halo_mass))

cosmo_ccl = CCLCosmology(Float32;
    Omega_c=0.2589, Omega_b=0.0486, h=0.6774, sigma8=0.8159,
    n_s=0.9667, transfer_function="boltzmann_camb")

massdef = CCLMassDef(200, "critical")
hmf = CCLMassFuncTinker08(cosmo_ccl, massdef)

import PythonPlot; const plt = PythonPlot
masses = 10 .^ (10:0.1:15) .* u"Msun"
plt.close()
plt.plot(ustrip.((u"Msun",), masses),
    [ustrip(u"Mpc^(-3)* Msun^(-1)", dndm(hmf, m, 1.0)) for m in masses])
plt.xlabel(raw"Mass [$M_{\odot}$]")
plt.ylabel(raw"$dn/dM_h$ $[Mpc^{-3} M^{-1}_{\odot}]$")
plt.xscale("log"); plt.yscale("log")
plt.gcf()
```

