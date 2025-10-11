struct SSDPhysics <: G4VUserPhysicsList
    function SSDPhysics(verbose)
        pl = G4VModularPhysicsList()
        # pl = FTFP_BERT(verbose)
        # lp = G4StepLimiterPhysics()
        # SetApplyToAll(lp, true)            # Apply to all particles
        # RegisterPhysics(pl, move!(lp))     # Register to the physics list
        RegisterPhysics(pl, move!(G4DecayPhysics(verbose)))              # Default physics
        RegisterPhysics(pl, move!(G4EmStandardPhysics_option4(verbose))) # EM physics
        RegisterPhysics(pl, move!(G4RadioactiveDecayPhysics(verbose)))   # Radioactive decay
        return pl
    end
end

struct SDData{T} <: G4JLSDData
    detectorHits::DetectorHits
    SDData{T}() where {T} = new(DetectorHits((evtno = Int32[], detno = Int32[], thit = typeof(zero(T)*u"s")[], edep = typeof(zero(T)*u"keV")[], pos = SVector{3,typeof(zero(T)*u"mm")}[])))
end

function _initialize(::G4HCofThisEvent, data::SDData)::Nothing
    empty!(data.detectorHits)
    return
end

function _endOfEvent(::G4HCofThisEvent, data::SDData)::Nothing
    return
end

function _processHits(step::G4Step, ::G4TouchableHistory, data::SDData{T})::Bool where {T}
    edep = step |> GetTotalEnergyDeposit
    iszero(edep) && return false
    pos = step |> GetPostStepPoint |> GetPosition
    push!(data.detectorHits, 
        (
            # evtno
            evtno = Int32(0),
            # detno
            detno = Int32(1), #step |> GetPostStepPoint |> GetPhysicalVolume, # get ID?
            # thit
            thit = (step |> GetPostStepPoint |> GetGlobalTime) / Geant4.SystemOfUnits.nanosecond * u"ns",
            # edep
            edep = T(edep * 1000u"keV"),
            # pos
            pos = SVector{3}(T(x(pos))*u"mm",T(y(pos))*u"mm",T(z(pos))*u"mm")
        )
    )
    return true
end

function endeventaction(evt::G4Event, app::G4JLApplication)
    # direction = evt |> (evt -> GetPrimaryVertex(evt, 0)) |> (vertex -> GetPrimary(vertex, 0)) |> GetMomentumDirection
    # pos = evt |> (evt -> GetPrimaryVertex(evt, 0)) |> GetPosition
    # d = [x(direction), y(direction), z(direction)]
    # @info [x(pos), y(pos), z(pos)], normalize(d)
    return
end


# We are using G4SingleParticleSource instead of G4ParticleGun for uniformity 
# across MonoenergeticSource and IsotopeSource, see discussion
# https://github.com/JuliaPhysics/SolidStateDetectors.jl/pull/526
@with_kw mutable struct GeneratorData <: G4JLGeneratorData
    gun::Union{Nothing, CxxPtr{G4SingleParticleSource}} = nothing
    particle::Union{Nothing, CxxPtr{G4ParticleDefinition}} = nothing
    position::G4ThreeVector = G4ThreeVector(0, 0, 0)
    direction::G4ThreeVector = G4ThreeVector(0, 0, 0)
end

function SSDGenerator(source::SolidStateDetectors.MonoenergeticSource;  kwargs...)

    iszero(source.direction) && throw(ArgumentError("The direction of the source cannot be zero."))

    data = GeneratorData(;kwargs...)

    function _init(data::GeneratorData, ::Any)
        data.gun = move!(G4SingleParticleSource())
        particle = data.particle = FindParticle(source.particle_type)
        SetParticleDefinition(data.gun, particle)

        # Place the SingleParticleSource as point source
        data.position = G4ThreeVector(
            source.position.x * Geant4.SystemOfUnits.meter,
            source.position.y * Geant4.SystemOfUnits.meter,
            source.position.z * Geant4.SystemOfUnits.meter
        )
        posdist = GetPosDist(data.gun)
        SetPosDisType(posdist, "Point")
        SetCentreCoords(posdist, data.position)

        # Define the source emission
        d::CartesianVector = normalize(source.direction)
        a::CartesianVector = normalize(d × (abs(d.x) == 1 ? CartesianVector(0,1,0) : CartesianVector(1,0,0)))
        b::CartesianVector = normalize(a × d)
        data.direction = G4ThreeVector(d.x, d.y, d.z)

        angdist = GetAngDist(data.gun)
        SetAngDistType(angdist, "iso")
        DefineAngRefAxes(angdist, "angref1", G4ThreeVector(a...))
        DefineAngRefAxes(angdist, "angref2", G4ThreeVector(b...))
        SetMinTheta(angdist, 0.0)
        SetMaxTheta(angdist, ustrip(NoUnits, source.opening_angle))  # convert to radians
        # SetParticleMomentumDirection(angdist, data.direction)

        # set energy
        SetMonoEnergy(GetEneDist(data.gun), ustrip(u"MeV", source.energy))
    end
    
    function _gen(evt::G4Event, data::GeneratorData)::Nothing
        GeneratePrimaryVertex(data.gun, CxxPtr(evt))
    end

    G4JLPrimaryGenerator("SSDGenerator", data; init_method=_init, generate_method=_gen)
end

function SSDGenerator(source::SolidStateDetectors.IsotopeSource; kwargs...)

    iszero(source.direction) && throw(ArgumentError("The direction of the source cannot be zero."))

    data = GeneratorData(;kwargs...)
    function _init(data::GeneratorData, ::Any)
        data.gun = move!(G4SingleParticleSource())

        # Place the SingleParticleSource as point source
        data.position = G4ThreeVector(
            source.position.x * Geant4.SystemOfUnits.meter,
            source.position.y * Geant4.SystemOfUnits.meter,
            source.position.z * Geant4.SystemOfUnits.meter
        )
        posdist = GetPosDist(data.gun)
        SetPosDisType(posdist, "Point")
        SetCentreCoords(posdist, data.position)

        # Define the source emission
        d::CartesianVector = normalize(source.direction)
        a::CartesianVector = normalize(d × (abs(d.x) == 1 ? CartesianVector(0,1,0) : CartesianVector(1,0,0)))
        b::CartesianVector = normalize(a × d)
        data.direction = G4ThreeVector(d.x, d.y, d.z)

        angdist = GetAngDist(data.gun)
        SetAngDistType(angdist, "iso")
        DefineAngRefAxes(angdist, "angref1", G4ThreeVector(a...))
        DefineAngRefAxes(angdist, "angref2", G4ThreeVector(b...))
        SetMinTheta(angdist, 0.0)
        SetMaxTheta(angdist, ustrip(NoUnits, source.opening_angle))  # convert to radians

        # Set kinetic energy of the decaying isotope to zero ("decay at rest")
        SetMonoEnergy(GetEneDist(data.gun), 0.0 * Geant4.SystemOfUnits.MeV)
    end

    function _gen(evt::G4Event, data::GeneratorData)::Nothing
        if isnothing(data.particle)  # late initialize (after physics processes)
            data.particle = GetIon(source.Z, source.A, Float64(source.excitEnergy))
            SetParticleDefinition(data.gun, data.particle)
            SetParticleCharge(data.gun, source.ionCharge)
        end
        GeneratePrimaryVertex(data.gun, CxxPtr(evt))
    end
    
    G4JLPrimaryGenerator("SSDGenerator", data; init_method=_init, generate_method=_gen)
end
