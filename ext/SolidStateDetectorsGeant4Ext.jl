# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

module SolidStateDetectorsGeant4Ext


@static if isdefined(Base, :get_extension)
    using Geant4
else
    using ..Geant4
end

include("io_gdml.jl")
using Suppressor
using LightXML

# Given an SSD simulation object, create corresponding GDML file to desired location
function Geant4.G4JLDetector(sim::SolidStateDetectors.Simulation, output_filename::String = "tmp.gdml"; verbose::Bool = true)
    # Create basis for GDML file
    x_doc = XMLDocument()
    x_root = create_root(x_doc, "gdml")
    set_attributes(x_root, 
        OrderedDict(
            "xmlns:xsi" => "https://www.w3.org/2001/XMLSchema-instance",
            "xsi:noNamespaceSchemaLocation" => "https://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd"
        )
    )

    x_define = new_child(x_root, "define")
    x_materials = parse_file(joinpath(@__DIR__, "materials.xml"))
    add_child(x_root, root(x_materials))
    x_solids = new_child(x_root, "solids")
    x_structure = new_child(x_root, "structure")
    x_setup = new_child(x_root, "setup")
    set_attributes(x_setup, OrderedDict("name" => "Default", "version" => "1.0"))
    x_setup_world = new_child(x_setup, "world")
    set_attribute(x_setup_world, "ref", "wd")

    # Build <solids> section from detector geometry, add volumes to <structure>
    # Parse semiconductor
    parse_geometry(sim.detector.semiconductor.geometry, x_solids, x_define, 1, "sc_", verbose)
    create_volume(x_structure, "sc", parse_material(sim.detector.semiconductor.material.name))
    # Parse contacts
    for (i, contact) in enumerate(sim.detector.contacts)
        if has_volume(contact.geometry, x_solids, x_define, 1, "ct$(i)_", verbose, parse = false)
            parse_geometry(contact.geometry, x_solids, x_define, 1, "ct$(i)_", verbose)
            create_volume(x_structure, "ct$(i)", parse_material(contact.material.name))
        end
    end
    # Parse passives
    if !ismissing(sim.detector.passives)
        for (i, passive) in enumerate(sim.detector.passives)
            if has_volume(passive.geometry, x_solids, x_define, 1, "pv$(i)_", verbose, parse = false)
                parse_geometry(passive.geometry, x_solids, x_define, 1, "pv$(i)_", verbose)
                create_volume(x_structure, "pv$(i)", parse_material(passive.material.name))
            end
        end
    end
    # Parse world (given by grid)
    parse_geometry(sim.world, x_solids, x_define, 1, "wd_", verbose)
    x_wd_vol = create_volume(x_structure, "wd", parse_material(sim.medium.name))


    # "physvol" contains the list of all children volumes inside the world
    z = new_child(x_wd_vol, "physvol")
    # Append semiconductor to physical volume
    add_to_world(z, x_define, sim.detector.semiconductor.geometry, "sc")
    # Append contacts to physical volume
    for (i, contact) in enumerate(sim.detector.contacts)
        if has_volume(contact.geometry, x_solids, x_define, 1, "ct$(i)_", verbose, parse = false)
            z = new_child(x_wd_vol, "physvol")
            add_to_world(z, x_define, contact.geometry, "ct$(i)")
        end
    end
    # Append passives to physical volume
    if !ismissing(sim.detector.passives)
        for (i, passive) in enumerate(sim.detector.passives)
            if has_volume(passive.geometry, x_solids, x_define, 1, "pv$(i)_", verbose, parse = false)
                z = new_child(x_wd_vol, "physvol")
                add_to_world(z, x_define, passive.geometry, "pv$(i)")
            end
        end
    end

    save_file(x_doc, output_filename)
    @suppress_out Geant4.G4JLDetectorGDML(output_filename)
end

function Geant4.G4JLDetector(input_filename::String, output_filename::String = "tmp.gdml"; verbose::Bool = true)
    sim = SolidStateDetectors.Simulation(input_filename)
    Geant4.G4JLDetector(sim, output_filename, verbose = verbose)
end


function Geant4.G4JLApplication(
    sim::SolidStateDetectors.Simulation,
    # ... source type/position, ...
)

    detector = Geant4.G4JLDetector(sim)

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
