# # Example 1: Inverted Coax Detector 

# ToDo... explanation

using Plots; pyplot();
using SolidStateDetectors

T = Float32
detector = SolidStateDetector{T}(SSD_examples[:InvertedCoax])

simulation = Simulation(detector)

simulate!(simulation, max_refinements = 3)

p_ep = plot(simulation.electric_potential)



# # Example 2: ToDo... 

# ToDo... 