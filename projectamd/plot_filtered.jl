
using Makie, CairoMakie, HDF5

fname_orig = "/fs/lustre/project/act/zack/ICs/Fvec_7700Mpc_n6144_nb30_nt16_no768_p.h5"
fname_filt = "filtered_output.h5"

plot_slice = (1:100, 1:100, 1)
m_original = h5read(fname_orig, "data", plot_slice)
m_filtered = h5read(fname_filt, "data", plot_slice)

##
heatmap(m_original)

##
heatmap(m_filtered)

##
