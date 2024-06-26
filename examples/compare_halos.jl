using XGPaint
websky_pos, websky_mass = read_halo_catalog_hdf5("/fs/lustre/cita/zack/projects/websky/websky_halos-light.hdf5")

##
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
using JLD2
positions = jldopen("/fs/lustre/scratch/zack/ICs/test_halos/lowres_positions_12.00_12.10.jld2", "r")["positions"]

fullbox = deepcopy(positions)
fullbox[1,:] .= positions[3,:]
fullbox[3,:] .= positions[1,:]


##
using PythonPlot
plt = PythonPlot
websky_ref = websky_pos[:, 10^12.4 .< websky_mass .< 10^12.5]  # real websky halo mass resolution is 12.2

XY_LENGTH = 500
websky_ref = websky_ref[:, 
    (1000 .< websky_ref[1,:] .< 1200) .& 
    (-XY_LENGTH .< websky_ref[2,:] .< XY_LENGTH) .& 
    (-XY_LENGTH .< websky_ref[3,:] .< XY_LENGTH)]

plt.clf()
plt.hist2D( websky_ref[2,:], websky_ref[3,:], bins=100)
plt.gca().set_aspect(1)
# plt.colorbar()
plt.text(-470,-470,"Websky Halos", color="w", fontweight="medium")
plt.gcf()


##
lpt_ref = fullbox[:, (1000 .< fullbox[1,:] .< 1200)]
lpt_ref = lpt_ref[:, (-XY_LENGTH .< lpt_ref[2,:] .< XY_LENGTH)]
lpt_ref = lpt_ref[:, (-XY_LENGTH .< lpt_ref[3,:] .< XY_LENGTH)]

plt.clf()
plt.hist2D( lpt_ref[2,:], lpt_ref[3,:], bins=100)
# plt.colorbar()
plt.gca().set_aspect(1)
plt.text(-470,-470,"LPT", color="w", fontweight="medium")
plt.gcf()

##
lpt_ref = fullbox[:, (1000 .< fullbox[3,:] .< 1200)]
lpt_ref = lpt_ref[:, (-XY_LENGTH .< lpt_ref[1,:] .< XY_LENGTH)]
lpt_ref = lpt_ref[:, (-XY_LENGTH .< lpt_ref[2,:] .< XY_LENGTH)]

plt.clf()
plt.hist2D( lpt_ref[1,:], lpt_ref[2,:], bins=100)
# plt.colorbar()
plt.gca().set_aspect(1)
plt.text(-470,-470,"LPT", color="w", fontweight="medium")
plt.gcf()

##

