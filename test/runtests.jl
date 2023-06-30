using LagrangianPerturbationTheory
using Test, Unitful, UnitfulAstro

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

##
@testset "lattice_0" begin
    delta_array = (LagrangianPerturbationTheory.load_example_ics())
    grid_spacing = 7700.f0u"Mpc" / size(delta_array,1)
    box_sizes = (7700.f0u"Mpc", 7700.f0u"Mpc", 7700.f0u"Mpc")
    offset = grid_spacing / 2
    nx, ny, nz = size(delta_array)  # make a grid for each Lagrangian coordinate
    q_axes = (LinRange(offset, offset + (nx - 1) * grid_spacing, nx),
              LinRange(offset, offset + (ny - 1) * grid_spacing, ny),
              LinRange(offset, offset + (nz - 1) * grid_spacing, nz))
    cosmo = CCLCosmology(Float32;
        Omega_c=0.2589, Omega_b=0.0486, h=0.6774, sigma8=0.8159,
        n_s=0.9667, transfer_function="boltzmann_camb")
    grid = LagrangianGridWebsky(cosmo, grid_spacing, box_sizes, q_axes)
    δ₀ = ICFieldWebsky(FirstOrderLPT, grid, delta_array)

    octants = (-1, 0)
    i=3+1; j=4+1; k=5+1
    counter = 1
    for oi in octants, oj in octants, ok in octants
        𝐪 = lagrangian_coordinate(grid, i, j, k, oi, oj, ok)
        a = scale_factor(grid, 𝐪)  # a(𝐪) depends on normalization of grid
        D = growth_factor(cosmo, a)
        δ⁽¹⁾ᴸ = D * δ₀[𝐪]
        @test [
            ustrip(u"Mpc", 𝐪.x), ustrip(u"Mpc", 𝐪.y), ustrip(u"Mpc", 𝐪.z), 
            oi, oj, ok, δ⁽¹⁾ᴸ, 1/a-1] ≈ refs1[counter]
        counter += 1
    end

    i=8+1; j=2+1; k=2+1
    counter = 1
    for oi in octants, oj in octants, ok in octants
        𝐪 = lagrangian_coordinate(grid, i, j, k, oi, oj, ok)
        a = scale_factor(grid, 𝐪)  # a(𝐪) depends on normalization of grid
        D = growth_factor(cosmo, a)
        δ⁽¹⁾ᴸ = D * δ₀[𝐪]
        @test [
            ustrip(u"Mpc", 𝐪.x), ustrip(u"Mpc", 𝐪.y), ustrip(u"Mpc", 𝐪.z), 
            oi, oj, ok, δ⁽¹⁾ᴸ, 1/a-1] ≈ refs2[counter]
        counter += 1
    end
end


## todo turn into test: periodicity

# ##
# delta_array = LagrangianPerturbationTheory.load_example_ics()
# grid_spacing = 7700.0f0u"Mpc" / size(delta_array,1)
# box_sizes = (7700.0f0u"Mpc", 7700.0f0u"Mpc", 7700.0f0u"Mpc")
# offset = grid_spacing / 2

# # make a grid for each Lagrangian coordinate
# q_axes = (LinRange(offset, box_sizes[1] + offset, size(delta_array, 1)),
#            LinRange(offset, box_sizes[2] + offset, size(delta_array, 2)),
#            LinRange(offset, box_sizes[3] + offset, size(delta_array, 3)))
# cosmo = CCLCosmology(Float32;
#     Omega_c=0.2589, Omega_b=0.0486, h=0.6774, sigma8=0.8159,
#     n_s=0.9667, transfer_function="boltzmann_camb")
# lgrid = LagrangianGridWebsky(cosmo, grid_spacing, box_sizes, q_axes)

# δ₀ = ICFieldWebsky(FirstOrderLPT, lgrid, delta_array)

# ys = [δ₀[ SVector(offset, offset, offset + box_sizes[1] - i * grid_spacing) ] for i in 0:10]

# ##
# plt.clf()
# plt.plot( delta_array[1,1,192:-1:182] )
# plt.plot( ys )
# plt.gcf()


# ##


# ys = [δ₀[ SVector(offset, offset, offset + (i) * grid_spacing) ] for i in -2:0.1f0:2]


# ##
# plt.clf()
# plt.plot( ys)
# # plt.plot( ys )
# plt.gcf()


