# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

module SolidStateDetectorsGeant4Ext


@static if isdefined(Base, :get_extension)
    using Geant4
    using Geant4.SystemOfUnits
else
    using ..Geant4
    using ..Geant4.SystemOfUnits
end

using LightXML
using Parameters
using RadiationDetectorSignals
using StaticArrays
using Suppressor
using TypedTables
using Unitful

include(joinpath(@__DIR__, "Geant4", "io_gdml.jl"))
include(joinpath(@__DIR__, "Geant4", "g4jl_application.jl"))

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
    x_materials = parse_file(joinpath(@__DIR__, "Geant4", "materials.xml"))
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
        if has_volume(contact.geometry, verbose)
            parse_geometry(contact.geometry, x_solids, x_define, 1, "ct$(i)_", verbose)
            create_volume(x_structure, "ct$(i)", parse_material(contact.material.name))
        end
    end
    # Parse passives
    if !ismissing(sim.detector.passives)
        for (i, passive) in enumerate(sim.detector.passives)
            if has_volume(passive.geometry, verbose)
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
        if has_volume(contact.geometry, verbose)
            z = new_child(x_wd_vol, "physvol")
            add_to_world(z, x_define, contact.geometry, "ct$(i)")
        end
    end
    # Append passives to physical volume
    if !ismissing(sim.detector.passives)
        for (i, passive) in enumerate(sim.detector.passives)
            if has_volume(passive.geometry, verbose)
                z = new_child(x_wd_vol, "physvol")
                add_to_world(z, x_define, passive.geometry, "pv$(i)")
            end
        end
    end

    save_file(x_doc, output_filename)
    detector = @suppress_out Geant4.G4JLDetectorGDML(output_filename)
    
    verbose && @warn "Temporary file $(output_filename) will be deleted."
    rm(output_filename)
    detector 
end

function Geant4.G4JLDetector(input_filename::String, output_filename::String = "tmp.gdml"; verbose::Bool = true)
    if endswith(input_filename, ".gdml")
        @suppress_out Geant4.G4JLDetectorGDML(input_filename)
    else
        sim = SolidStateDetectors.Simulation(input_filename)
        Geant4.G4JLDetector(sim, output_filename, verbose = verbose)
    end
end


function Geant4.G4JLApplication(
    sim::Union{SolidStateDetectors.Simulation, String},
    source::SolidStateDetectors.AbstractParticleSource;
    physics_type = SSDPhysics,
    endeventaction_method = endeventaction,
    verbose = true,
    kwargs...
)

    SensitiveDetector = G4JLSensitiveDetector(
        "SensitiveDetector", SDData{Float32}();           # SD name and associated data are mandatory
        processhits_method=_processHits,    # process hist method (also mandatory)
        initialize_method=_initialize,      # intialize method
        endofevent_method=_endOfEvent       # end of event method
    );

    # @warn "Please never ever re-run this code"
    app = @suppress_out G4JLApplication(; 
        detector = Geant4.G4JLDetector(sim, verbose = verbose),
        sdetectors = ["sc" => SensitiveDetector],
        generator = SSDGenerator(source),
        physics_type = physics_type,
        endeventaction_method = endeventaction_method,
        kwargs...
    )     

    configure(app)
    initialize(app)
    app

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


function SolidStateDetectors.run_geant4_simulation(app::G4JLApplication, number_of_events::Int)
    evts = DetectorHit[]
    
    evtno = 1
    while evtno <= number_of_events
        out = missing
        while ismissing(out) || isempty(out)      
            @suppress_out beamOn(app,1)
            out = app.sdetectors["SensitiveDetector"][1].data.detectorHits
            out = out[findall(!iszero, out.edep)]
            out.evtno .= evtno
	   
        end
        append!(evts, out)
        evtno += 1
    end
    return RadiationDetectorSignals.group_by_evtno(Table(evts))
end 


end # module Geant4
