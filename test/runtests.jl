using LagrangianPerturbationTheory
using Test

const LPT = LagrangianPerturbationTheory

const refs1 = [[-7.559635416666667e+03,-7.519531250000000e+03,-7.479427083333333e+03,-1,-1,-1,+2.887501840140354e-03,+1.229087635970462e+02],
    [-7.559635416666667e+03,-7.519531250000000e+03,+2.205729166666667e+02,-1,-1,+0,+1.996219688295767e-02,+1.655453509493413e+01],
    [-7.559635416666667e+03,+1.804687500000000e+02,-7.479427083333333e+03,-1,+0,-1,+2.026693311264019e-02,+1.628962079675666e+01],
    [-7.559635416666667e+03,+1.804687500000000e+02,+2.205729166666667e+02,-1,+0,+0,+6.555612890825477e-02,+4.318618236109571e+00],
    [+1.403645833333333e+02,-7.519531250000000e+03,-7.479427083333333e+03,+0,-1,-1,+2.057481531635379e-02,+1.602996576205677e+01],
    [+1.403645833333333e+02,-7.519531250000000e+03,+2.205729166666667e+02,+0,-1,+0,+6.633302590566920e-02,+4.255748743727389e+00],
    [+1.403645833333333e+02,+1.804687500000000e+02,-7.479427083333333e+03,+0,+0,-1,+6.711872580568620e-02,+4.193636431562039e+00],
    [+1.403645833333333e+02,+1.804687500000000e+02,+2.205729166666667e+02,+0,+0,+0,+2.634667773933499e-01,+7.303013583997586e-02]]

const refs2 = [[-7.359114583333333e+03,-7.599739583333333e+03,-7.599739583333333e+03,-1,-1,-1,+9.619948761213303e-04,+1.231481275491521e+02],
    [-7.359114583333333e+03,-7.599739583333333e+03,+1.002604166666667e+02,-1,-1,+0,+6.965441787313151e-03,+1.578987998679561e+01],
    [-7.359114583333333e+03,+1.002604166666667e+02,-7.599739583333333e+03,-1,+0,-1,+6.965441787313144e-03,+1.578987998679563e+01],
    [-7.359114583333333e+03,+1.002604166666667e+02,+1.002604166666667e+02,-1,+0,+0,+2.319256693055088e-02,+4.015121448826588e+00],
    [+3.408854166666666e+02,-7.599739583333333e+03,-7.599739583333333e+03,+0,-1,-1,+6.358556122253649e-03,+1.739857528156851e+01],
    [+3.408854166666666e+02,-7.599739583333333e+03,+1.002604166666667e+02,+0,-1,+0,+2.161071180536011e-02,+4.385928502315838e+00],
    [+3.408854166666666e+02,+1.002604166666667e+02,-7.599739583333333e+03,+0,+0,-1,+2.161071180536012e-02,+4.385928502315837e+00],
    [+3.408854166666666e+02,+1.002604166666667e+02,+1.002604166666667e+02,+0,+0,+0,+8.738518980451039e-02,+8.512036581035365e-02]]


@testset "lptcell" begin
    delta_array = LPT.load_example_ics()
    boxsize = 7700. # comoving Mpc
    grid_spacing = boxsize / size(delta_array,1)
    cosmology = CCLCosmology(Float64;
        Omega_c=0.2589, Omega_b=0.0486, h=0.6774, sigma8=0.8159,
        n_s=0.9667, transfer_function="boltzmann_camb")
    delta = LPT.FirstOrderLPTWebsky(grid_spacing, cosmology, delta_array)

    i=3; j=4; k=5
    counter = 1
    for i_oct_x in (-1,0), i_oct_y in (-1,0), i_oct_z in (-1,0)
        x, y, z, a, δ₁ = LPT.lptcell(delta, i, j, k, i_oct_x, i_oct_y, i_oct_z)
        @test [x, y, z, i_oct_x, i_oct_y, i_oct_z, δ₁, 1/a-1] ≈ refs1[counter]
        counter += 1
    end

    i=8; j=2; k=2
    counter = 1
    for i_oct_x in (-1,0), i_oct_y in (-1,0), i_oct_z in (-1,0)
        x, y, z, a, δ₁ = LPT.lptcell(delta, i, j, k, i_oct_x, i_oct_y, i_oct_z)
        @test [x, y, z, i_oct_x, i_oct_y, i_oct_z, δ₁, 1/a-1] ≈ refs2[counter]
        counter += 1
    end
end

@testset "lptcell interpolated" begin
    delta_array = LPT.load_example_ics()
    boxsize = 7700. # comoving Mpc
    grid_spacing = boxsize / size(delta_array,1)
    cosmo_pyccl = CCLCosmology(Float64;
        Omega_c=0.2589, Omega_b=0.0486, h=0.6774, sigma8=0.8159,
        n_s=0.9667, transfer_function="boltzmann_camb")
    amax = LPT.scale_factor_of_chi(cosmo_pyccl, √(3) * boxsize)
    cosmology = LPT.InterpolatedCosmology(Float64, cosmo_pyccl)
    delta = LPT.FirstOrderLPTWebsky(grid_spacing, cosmology, delta_array)

    i=3; j=4; k=5
    counter = 1
    for i_oct_x in (-1,0), i_oct_y in (-1,0), i_oct_z in (-1,0)
        x, y, z, a, δ₁ = LPT.lptcell(delta, i, j, k, i_oct_x, i_oct_y, i_oct_z)
        @test [x, y, z, i_oct_x, i_oct_y, i_oct_z, δ₁, 1/a-1] ≈ refs1[counter]
        counter += 1
    end

    i=8; j=2; k=2
    counter = 1
    for i_oct_x in (-1,0), i_oct_y in (-1,0), i_oct_z in (-1,0)
        x, y, z, a, δ₁ = LPT.lptcell(delta, i, j, k, i_oct_x, i_oct_y, i_oct_z)
        @test [x, y, z, i_oct_x, i_oct_y, i_oct_z, δ₁, 1/a-1] ≈ refs2[counter] atol=1e-4
        counter += 1
    end
end
