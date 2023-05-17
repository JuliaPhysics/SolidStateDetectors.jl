# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

module SolidStateDetectorsGeant4Ext

@static if isdefined(Base, :get_extension)
    using Geant4
else
    using ..Geant4
end

import SolidStateDetectors


function _ssdsim2gdml(io::IO, sim::SolidStateDetectors.Simulation)
    # Dummy implementation:

    println(io, """<?xml version="1.0" encoding="UTF-8" standalone="no" ?>""")
    println(io, """<gdml>""")
    println(io, """</gdml>""")
end


function Geant4.G4JLDetector(sim::SolidStateDetectors.Simulation)
    mktemp() do gdml_filename, io
        _ssdsim2gdml(io, sim)
        close(io)
        G4JLDetectorGDML(gdml_filename)
    end
end





function G4JLApplication(
    sim::SolidStateDetectors.Simulation,
    # ... source type/position, ...
)

    detector = G4JLDetector(sim)

    throw(ErrorException("Not implemented yet"))

    # ToDo:
    #=
    generator = G4JLGunGenerator(...)

    G4JLApplication(
        detector = detector,               # detector with parameters
        simdata = TestEm3SimData(),                 # simulation data structure
        generator = generator,                    # primary particle generator 
        runmanager_type = G4RunManager,             # what RunManager to instantiate
        physics_type = FTFP_BERT,                   # what physics list to instantiate
        stepaction_method = stepaction,             # step action method
        pretrackaction_method = pretrackaction,     # pre-tracking action
        posttrackaction_method = posttrackaction,   # post-tracking action
        beginrunaction_method = beginrun,           # begin-run action (initialize counters and histograms)
        endrunaction_method = endrun,               # end-run action (print summary)               
        begineventaction_method = beginevent,       # begin-event action (initialize per-event data)
        endeventaction_method = endevent            # end-event action (fill histogram per event data)
    )
    =#
end


end # module Geant4
