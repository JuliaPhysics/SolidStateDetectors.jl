struct SSDPhysics <: G4VUserPhysicsList
    function SSDPhysics(verbose)
        pl = G4VModularPhysicsList()
        RegisterPhysics(pl, move!(G4DecayPhysics(verbose)))             # Default physics
        RegisterPhysics(pl, move!(G4EmStandardPhysics(verbose)))        # EM physics
        RegisterPhysics(pl, move!(G4RadioactiveDecayPhysics(verbose)))  # Radioactive decay
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
    # hits = getSDdata(app, "SensitiveDetector").detectorHits
    # eventID = evt |> GetEventID
    return
end

struct PhysicsList <: G4VUserPhysicsList
  function PhysicsList(verbose)
      pl = FTFP_BERT(verbose)
      lp = G4StepLimiterPhysics()
      SetApplyToAll(lp, true)            # Apply to all particles
      RegisterPhysics(pl, move!(lp))     # Register to the physics list
      return pl
  end 
end

@with_kw mutable struct GeneratorData <: G4JLGeneratorData
    gun::Union{Nothing, CxxPtr{G4ParticleGun}} = nothing
    particle::Union{Nothing, CxxPtr{G4ParticleDefinition}} = nothing
    position::G4ThreeVector = G4ThreeVector(0, 0, 0)
    direction::G4ThreeVector = G4ThreeVector(0, 0, 0)
end

function SSDGenerator(source::SolidStateDetectors.MonoenergeticSource;  kwargs...)

    @assert source.direction isa SolidStateDetectors.CartesianVector || source.direction == :isotropic

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
        if source.direction isa SolidStateDetectors.CartesianVector
          data.direction = G4ThreeVector(source.direction.x, source.direction.y, source.direction.z)
        end
        SetParticleMomentumDirection(gun, data.direction)
        SetParticleEnergy(gun, ustrip(u"MeV", source.energy))
        SetParticleDefinition(gun, particle)
    end

    function _gen(evt::G4Event, data::GeneratorData)::Nothing
      if source.direction == :isotropic
        ϕ = π/2 + rand()*π #rand()*2π 
        θ = acos(1 - 2*rand()) # because we want isotropic emission
        direction = G4ThreeVector(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ))
        data.direction = direction
        SetParticleMomentumDirection(data.gun, direction)
      end
      GeneratePrimaryVertex(data.gun, CxxPtr(evt))
    end

    G4JLPrimaryGenerator("SSDGenerator", data; init_method=_init, generate_method=_gen)
end

function SSDGenerator(source::SolidStateDetectors.IsotopeSource;  kwargs...)

    @assert source.direction isa SolidStateDetectors.CartesianVector || source.direction == :isotropic

    data = GeneratorData(;kwargs...)
    function _init(data::GeneratorData, ::Any)
        gun = data.gun = move!(G4ParticleGun())
        data.position = G4ThreeVector(
            source.position.x * Geant4.SystemOfUnits.meter,
            source.position.y * Geant4.SystemOfUnits.meter,
            source.position.z * Geant4.SystemOfUnits.meter
        )
        SetParticlePosition(gun, data.position)
        if source.direction isa SolidStateDetectors.CartesianVector
            data.direction = G4ThreeVector(source.direction.x, source.direction.y, source.direction.z)
        end
        SetParticleMomentumDirection(gun, data.direction)
        # SetParticleEnergy(gun, source.excitEnergy)
    end

    function _gen(evt::G4Event, data::GeneratorData)::Nothing
        if isnothing(data.particle)  # late initialize (after physics processes)
            data.particle = GetIon(source.Z, source.A, source.excitEnergy)
            SetParticleDefinition(data.gun, data.particle)
            SetParticleCharge(data.gun, source.ionCharge)
        end
        if source.direction == :isotropic
            ϕ = rand()*2π 
            θ = acos(1 - 2*rand())
            direction = G4ThreeVector(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ))
            data.direction = direction
        end
        SetParticleMomentumDirection(data.gun, data.direction)
        SetParticlePosition(data.gun, data.position) # needed ?
        GeneratePrimaryVertex(data.gun, CxxPtr(evt))
    end
    G4JLPrimaryGenerator("SSDGenerator", data; init_method=_init, generate_method=_gen)
end