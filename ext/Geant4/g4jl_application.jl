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
            thit = T(0)*u"s", #step |> GetPostStepPoint |> GetGlobalTime,
            # edep
            edep = T(edep * 1000u"keV"),
            # pos
            pos = SVector{3}(T(x(pos))*u"mm",T(y(pos))*u"mm",T(z(pos))*u"mm")
        )
    )
    return true
end

function endeventaction(evt::G4Event, app::G4JLApplication)
    return
end


@with_kw mutable struct GeneratorData <: G4JLGeneratorData
    gun::Union{Nothing, CxxPtr{G4ParticleGun}, CxxPtr{G4SingleParticleSource}} = nothing
    particle::Union{Nothing, CxxPtr{G4ParticleDefinition}} = nothing
    position::G4ThreeVector = G4ThreeVector(0, 0, 0)
    direction::G4ThreeVector = G4ThreeVector(0, 0, 0)
end

function SSDGenerator(source::SolidStateDetectors.MonoenergeticSource;  kwargs...)

    iszero(source.direction) && throw(ArgumentError("The direction of the source cannot be zero."))

    data = GeneratorData(;kwargs...)

    function _init(data::GeneratorData, ::Any)
        gun = data.gun = move!(G4ParticleGun())
        particle = data.particle = FindParticle(source.particle_type)
        data.position = G4ThreeVector(
            source.position.x * Geant4.SystemOfUnits.meter,
            source.position.y * Geant4.SystemOfUnits.meter,
            source.position.z * Geant4.SystemOfUnits.meter
        )
        SetParticlePosition(gun, data.position)
        data.direction = G4ThreeVector(source.direction.x, source.direction.y, source.direction.z)
        SetParticleMomentumDirection(gun, data.direction)
        SetParticleEnergy(gun, ustrip(u"MeV", source.energy))
        SetParticleDefinition(gun, particle)
    end
    
    function _gen(evt::G4Event, data::GeneratorData)::Nothing
        if !iszero(source.opening_angle)
            d::CartesianVector = normalize(source.direction)
            a::CartesianVector = normalize(d × (abs(d.x) == 1 ? CartesianVector(0,1,0) : CartesianVector(1,0,0)))
            b::CartesianVector = normalize(a × d)
            ϕ = rand()*2π
            θ = acos(1 - (1 - cos(source.opening_angle))*rand())
            v = (cos(θ) * d + sin(θ) * (cos(ϕ) * a + sin(ϕ) * b))
            direction = G4ThreeVector(v.x, v.y, v.z)
            SetParticleMomentumDirection(data.gun, direction)
        end
        GeneratePrimaryVertex(data.gun, CxxPtr(evt))
    end

    G4JLPrimaryGenerator("SSDGenerator", data; init_method=_init, generate_method=_gen)
end

function SSDGenerator(source::SolidStateDetectors.IsotopeSource; kwargs...)

    iszero(source.direction) && throw(ArgumentError("The direction of the source cannot be zero."))

    data = GeneratorData(;kwargs...)
    function _init(data::GeneratorData, ::Any)
        gun = data.gun = move!(G4SingleParticleSource())
        data.position = G4ThreeVector(
            source.position.x * Geant4.SystemOfUnits.meter,
            source.position.y * Geant4.SystemOfUnits.meter,
            source.position.z * Geant4.SystemOfUnits.meter
        )
        SetParticlePosition(gun, data.position)
        data.direction = G4ThreeVector(source.direction.x, source.direction.y, source.direction.z)
        SetParticleMomentumDirection(GetAngDist(gun), data.direction)
        SetMonoEnergy(gun |> GetEneDist, 0.0 * Geant4.SystemOfUnits.MeV)
    end

    function _gen(evt::G4Event, data::GeneratorData)::Nothing
        if isnothing(data.particle)  # late initialize (after physics processes)
            data.particle = GetIon(source.Z, source.A, Float64(source.excitEnergy))
            SetParticleDefinition(data.gun, data.particle)
            SetParticleCharge(data.gun, source.ionCharge)
        end
        if !iszero(source.opening_angle)
            d::CartesianVector = normalize(source.direction)
            a::CartesianVector = normalize(d × (abs(d.x) == 1 ? CartesianVector(0,1,0) : CartesianVector(1,0,0)))
            b::CartesianVector = normalize(a × d)
            ϕ = rand()*2π
            θ = acos(1 - (1 - cos(source.opening_angle))*rand())
            v = (cos(θ) * d + sin(θ) * (cos(ϕ) * a + sin(ϕ) * b))
            direction = G4ThreeVector(v.x, v.y, v.z)
            SetParticleMomentumDirection(GetAngDist(data.gun), direction)
        end
        GeneratePrimaryVertex(data.gun, CxxPtr(evt))
    end
    
    G4JLPrimaryGenerator("SSDGenerator", data; init_method=_init, generate_method=_gen)
end
