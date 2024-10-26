using LagrangianPerturbationTheory, Unitful, UnitfulAstro, HDF5, PoissonRandom

icdir = "/fs/lustre/project/act/zack/ICs/"
files = (den=joinpath(icdir, "Fvec_7700Mpc_n6144_nb30_nt16_no768_p.h5"),
         sx=joinpath(icdir, "sx1_7700Mpc_n6144_nb30_nt16_no768_p.h5"),
         sy=joinpath(icdir, "sy1_7700Mpc_n6144_nb30_nt16_no768_p.h5"),
         sz=joinpath(icdir, "sz1_7700Mpc_n6144_nb30_nt16_no768_p.h5"))

##
# δ₀ = LagrangianPerturbationTheory.load_example_ics()
# δ₀ = permutedims(δ₀, (3,2,1))

δ₀ = h5open(files.den, "r") do f read(f, "data") end
sx_ref = h5open(files.sz, "r") do f read(f, "data") end
sy_ref = h5open(files.sy, "r") do f read(f, "data") end
sz_ref = h5open(files.sx, "r") do f read(f, "data") end

##
box_size = 7700
kv = LagrangianPerturbationTheory.wavevectors3D(Float32, size(δ₀), box_size)

##
@time ϕ⁽¹⁾ᵢₙᵢ = lpt(FirstOrderLPT(), δ₀, kv)

##
# sx = permutedims(ϕ⁽¹⁾ᵢₙᵢ[1], (3,2,1))
# sy = permutedims(ϕ⁽¹⁾ᵢₙᵢ[2], (3,2,1))
# sz = permutedims(ϕ⁽¹⁾ᵢₙᵢ[3], (3,2,1));

sx = ϕ⁽¹⁾ᵢₙᵢ[1]
sy = ϕ⁽¹⁾ᵢₙᵢ[2]
sz = ϕ⁽¹⁾ᵢₙᵢ[3]

##

f = Figure(size = (800, 400))
ax1 = Axis(f[1, 1], title = "DIY")
ax2 = Axis(f[1, 2], title = "REF")

hm = heatmap!(ax1, real(sx[2, 1:50, 1:50]))

heatmap!(ax2, real(sx_ref[2, 1:50, 1:50]))

Colorbar(f[1,3], hm)

f

##
f, ax, hm = heatmap( real(sx[1:50, 1:50, 2] .- sx_ref[1:50, 1:50, 2]))
Colorbar(f[1,2], hm)
f

##

# using PythonPlot; const plt = PythonPlot
# plt.close()
# fig, axes = plt.subplots(1,2, figsize=(10,5))
# axes[0].imshow(real.(sx[1:50, 1:50, 1]))
# axes[0].text(2,47,"DIY x", color="w", fontweight="medium")
# axes[1].imshow(real.(sx_ref[1:50, 1:50, 1]))
# axes[1].text(2,47, "REF x", color="w", fontweight="medium")
# plt.gcf()

# ##
# plt.close()
# fig, axes = plt.subplots(1,2, figsize=(10,5))
# axes[0].imshow(real.(sy[1:50, 1:50, 1]))
# axes[0].text(2,47,"DIY y", color="w", fontweight="medium")
# axes[1].imshow(real.(sy_ref[1:50, 1:50, 1]))
# axes[1].text(2,47, "REF y", color="w", fontweight="medium")
# plt.gcf()

# ##
# plt.close()
# fig, axes = plt.subplots(1,2, figsize=(10,5))
# axes[0].imshow(real.(sz[1:50, 1:50, 1]))
# axes[0].text(2,47,"DIY z", color="w", fontweight="medium")
# axes[1].imshow(real.(sz_ref[1:50, 1:50, 1]))
# axes[1].text(2,47, "REF z", color="w", fontweight="medium")
# plt.gcf()
