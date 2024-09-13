using LegendHDF5IO

sim = Simulation(SSD_examples[:InvertedCoax])

@testset "Conversion Simulation <--> NamedTuple" begin
    timed_calculate_electric_potential!(sim, verbose = false)
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)

    timed_calculate_electric_field!(sim)
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)

    timed_calculate_weighting_potential!(sim, 1, verbose = false)
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)

    for i in findall(ismissing.(sim.weighting_potentials)) timed_calculate_weighting_potential!(sim, i, verbose = false) end
    nt = NamedTuple(sim)
    @test sim == Simulation(nt)
end

@testset "Test LegendHDF5IO" begin
    @test_nowarn ssd_write("legendhdf5io_test.lh5", sim)
    @test_logs (:warn, "Destination `legendhdf5io_test.lh5` already exists. Overwriting...") ssd_write("legendhdf5io_test.lh5", sim)
    @test sim == ssd_read("legendhdf5io_test.lh5", Simulation)
end