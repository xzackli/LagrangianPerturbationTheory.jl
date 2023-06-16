using LagrangianPerturbationTheory
using Test

const LPT = LagrangianPerturbationTheory

@testset "LagrangianPerturbationTheory.jl" begin

    delta = LPT.load_example_ics()
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

    counter = 1
    for xd in (-1,0)
        x = (i+0.5)*grid_spacing + xd * boxsize
        for yd in (-1,0)
            y = (j+0.5)*grid_spacing + yd * boxsize
            for zd in (-1,0)
                z = (k+0.5)*grid_spacing + zd * boxsize
    
                chi = sqrt(x^2+y^2+z^2)
                scale_factor = LPT.scale_factor_of_chi(cosmology, chi)
                redshift = 1 / scale_factor - 1
                δ₀ = delta[k+1,j+1,i+1]  # (py → jl) reverses dimensions, index+1
                delta_linear = δ₀ * LPT.growth_factor(cosmology, scale_factor)
                @test [x, y, z,xd,yd,zd,delta_linear] ≈ refs[counter]
                counter += 1
            end
        end
    end
end
