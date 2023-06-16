
function load_example_ics()
    rootpath = artifact"density_example"
    path = joinpath(rootpath, "Fvec_7700Mpc_n6144_nb30_nt16_ks192.h5")
    delta = h5open(path, "r") do file
        read(file, "delta")
    end
    return delta
end
