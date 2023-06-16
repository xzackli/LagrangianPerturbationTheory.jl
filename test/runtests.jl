using LagrangianPerturbationTheory
using Test

const LPT = LagrangianPerturbationTheory

@testset "LagrangianPerturbationTheory.jl" begin

    delta_array = LPT.load_example_ics()
    boxsize = 7700. # comoving Mpc
    N=192
    grid_spacing = boxsize / N

    cosmology = CCLCosmology(Float64;
        Omega_c=0.2589, Omega_b=0.0486, h=0.6774, sigma8=0.8159,
        n_s=0.9667, transfer_function="boltzmann_camb")

    refs = [[-7.559635416666667e+03,-7.519531250000000e+03,-7.479427083333333e+03,-1,-1,-1,+2.887501840140354e-03],
        [-7.559635416666667e+03,-7.519531250000000e+03,+2.205729166666667e+02,-1,-1,+0,+1.996219688295767e-02],
        [-7.559635416666667e+03,+1.804687500000000e+02,-7.479427083333333e+03,-1,+0,-1,+2.026693311264019e-02],
        [-7.559635416666667e+03,+1.804687500000000e+02,+2.205729166666667e+02,-1,+0,+0,+6.555612890825477e-02],
        [+1.403645833333333e+02,-7.519531250000000e+03,-7.479427083333333e+03,+0,-1,-1,+2.057481531635379e-02],
        [+1.403645833333333e+02,-7.519531250000000e+03,+2.205729166666667e+02,+0,-1,+0,+6.633302590566920e-02],
        [+1.403645833333333e+02,+1.804687500000000e+02,-7.479427083333333e+03,+0,+0,-1,+6.711872580568620e-02],
        [+1.403645833333333e+02,+1.804687500000000e+02,+2.205729166666667e+02,+0,+0,+0,+2.634667773933499e-01]]
    
    i=3; j=4; k=5

    delta = LPT.FirstOrderLPTWebsky(grid_spacing, cosmology, delta_array)

    counter = 1
    for i_oct_x in (-1,0)
        for i_oct_y in (-1,0)
            for i_oct_z in (-1,0)
                x, y, z, δ₁ = LPT.lptcell(delta, i, j, k, i_oct_x, i_oct_y, i_oct_z)
                @test [x, y, z, i_oct_x, i_oct_y, i_oct_z, δ₁] ≈ refs[counter]
                counter += 1
            end
        end
    end
end
