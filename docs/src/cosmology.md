```@meta
CurrentModule = LagrangianPerturbationTheory
```

# Cosmology Calculations

This package wraps the [Core Cosmology Library (CCL)](https://ccl.readthedocs.io/en/latest/) 
for calculations involving the background cosmology, halo mass functions, and halo bias.
The wrapped functions provide (zero-cost, compile-time) physical units support. 

To start, you'll need a basic CCL cosmology object.

```@example bg
using LagrangianPerturbationTheory, Unitful, UnitfulAstro

cosmo_ccl = CCLCosmology(Float32;
    Omega_c=0.2589, Omega_b=0.0486, h=0.6774, sigma8=0.8159,
    n_s=0.9667, transfer_function="boltzmann_camb")
```

**Wrapped CCL cannot be called in multithreaded mode.** Do not use CCL objects in multithreaded
code, as your program will likely immediately segfault. For high-performance code, you will
want to use interpolators anyway:

```@example bg
cosmo = InterpolatedCosmology(cosmo_ccl)
```

Similarly, we provide wrapped CCL mass definitions and halo bias.

```@example bg
massdef = CCLMassDef(200, "critical")

dndm(hmf, 1e12u"Msun", 1.0)  # scale factor = 1
```

Let's take a look.

```@example bg
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


```@example bg
tinker10 = CCLHaloBiasTinker10(cosmo_ccl, massdef)

halo_bias(tinker10, 1e12u"Msun", 1.0)
```


```@example bg
plt.close()
plt.plot(ustrip.((u"Msun",), masses),
    [halo_bias(tinker10, m, 1.0) for m in masses], label="Tinker10")
plt.xlabel(raw"Mass [$M_{\odot}$]")
plt.ylabel(raw"$b_h(M)$")
plt.xscale("log")
plt.legend()
plt.gcf()
```
