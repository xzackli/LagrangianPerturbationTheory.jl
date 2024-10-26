using XGPaint
websky_pos, websky_mass = read_halo_catalog_hdf5("/fs/lustre/cita/zack/projects/websky/websky_halos-light.hdf5")

using LagrangianPerturbationTheory

##
# positions = jldopen("/fs/lustre/scratch/zack/ICs/halos/lowres_positions_12.4_12.5.jld2", "r")["positions"]
# function strip(positions)
#     N_halos = length(positions)
#     mirror = zeros(Float32, (3, N_halos))
#     for (i_halo, pos) in enumerate(positions)
#         mirror[1, i_halo] = ustrip(u"Mpc", pos[1])
#         mirror[2, i_halo] = ustrip(u"Mpc", pos[2])
#         mirror[3, i_halo] = ustrip(u"Mpc", pos[3])
#     end
#     return mirror
# end

# fullbox = strip(positions)

# flipbox = deepcopy(positions)
# flipbox[1,:] .= flipbox[3,:]
# flipbox[3,:] .= flipbox[1,:]


##
using Makie, CairoMakie
# using PythonPlot
# plt = PythonPlot
# websky_ref = websky_pos[:, 10^12.4 .< websky_mass .< 10^12.5]  # real websky halo mass resolution is 12.2
websky_ref = websky_pos

XY_LENGTH = 500
websky_ref = websky_ref[:, 
    (500 .< websky_ref[1,:] .< 700) .& 
    (-XY_LENGTH .< websky_ref[2,:] .< XY_LENGTH) .& 
    (-XY_LENGTH .< websky_ref[3,:] .< XY_LENGTH)]

##
f = Figure(size=(580,500)); ax = Axis(f[1,1])
hb = hexbin!( websky_ref[2,:], websky_ref[3,:], bins=200, threshold=0, colorrange=(0,50))
limits!(ax, (-XY_LENGTH, XY_LENGTH), (-XY_LENGTH, XY_LENGTH))
Colorbar(f[1,2], hb)
text!(-480, -480; text="Websky Halo Catalog", color="white")
save("examples/plots/websky_halos.png", f)
f

# plt.clf()
# plt.hist2D( websky_ref[2,:], websky_ref[3,:], bins=100)
# plt.gca().set_aspect(1)
# # plt.colorbar()
# plt.text(-470,-470,"Websky Halos", color="w", fontweight="medium")
# plt.gcf()


##
using JLD2
positions = jldopen("/fs/lustre/scratch/zack/ICs/low_mass_halos/lowres_positions_11.00_11.10.jld2", "r")["positions"]
# positions = jldopen("/fs/lustre/scratch/zack/ICs/low_mass_halos/xyz.jld2", "r")["positions"]


lpt_ref = positions[:, (500 .< positions[1,:] .< 700)]
lpt_ref = lpt_ref[:, (-XY_LENGTH .< lpt_ref[2,:] .< XY_LENGTH)]
lpt_ref = lpt_ref[:, (-XY_LENGTH .< lpt_ref[3,:] .< XY_LENGTH)]

##
f = Figure(size=(580,500)); ax = Axis(f[1,1])
hb = hexbin!( lpt_ref[2,:], lpt_ref[3,:], bins=200, threshold=0, colorrange=(0,120))
limits!(ax, (-XY_LENGTH, XY_LENGTH), (-XY_LENGTH, XY_LENGTH))
Colorbar(f[1,2], hb)
text!(-480, -480; text="LPT with density (3,2,1)", color="white")
save("examples/plots/zyx.png", f)
f
##
# using JLD2
# positions_xyz = jldopen("/fs/lustre/scratch/zack/ICs/low_mass_halos/xyz.jld2", "r")["positions"]

# lpt_ref_xyz = positions_xyz[:, (500 .< positions_xyz[3,:] .< 700)]
# lpt_ref_xyz = lpt_ref_xyz[:, (-XY_LENGTH .< lpt_ref_xyz[2,:] .< XY_LENGTH)]
# lpt_ref_xyz = lpt_ref_xyz[:, (-XY_LENGTH .< lpt_ref_xyz[1,:] .< XY_LENGTH)]

# ##
# f = Figure(size=(580,500)); ax = Axis(f[1,1])
# hb = hexbin!( lpt_ref_xyz[2,:], lpt_ref_xyz[1,:], bins=200, threshold=0, colorrange=(0,90))
# limits!(ax, (-XY_LENGTH, XY_LENGTH), (-XY_LENGTH, XY_LENGTH))
# Colorbar(f[1,2], hb)
# text!(-480, -480; text="LPT with density(3,2,1)", color="white")
# save("examples/plots/xyz.png", f)
# f

# plt.clf()
# plt.hist2D( lpt_ref[2,:], lpt_ref[3,:], bins=100)
# # plt.colorbar()
# plt.gca().set_aspect(1)
# plt.text(-470,-470,"LPT", color="w", fontweight="medium")
# plt.gcf()

##
# lpt_ref = fullbox[:, (1000 .< fullbox[3,:] .< 1200)]
# lpt_ref = lpt_ref[:, (-XY_LENGTH .< lpt_ref[1,:] .< XY_LENGTH)]
# lpt_ref = lpt_ref[:, (-XY_LENGTH .< lpt_ref[2,:] .< XY_LENGTH)]

# plt.clf()
# plt.hist2D( lpt_ref[1,:], lpt_ref[2,:], bins=100)
# # plt.colorbar()
# plt.gca().set_aspect(1)
# plt.text(-470,-470,"LPT", color="w", fontweight="medium")
# plt.gcf()

# ##

