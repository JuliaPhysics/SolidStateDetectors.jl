using Plots; pyplot();
using SolidStateDetectors

T = Float32
detector = SolidStateDetector{T}(SSD_examples[:InvertedCoax])

simulation = Simulation(detector)

simulate!(simulation)

p_ep = plot(simulation.electric_potential)

