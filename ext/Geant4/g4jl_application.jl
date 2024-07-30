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

function SSDGenerator(source::SolidStateDetectors.SSDSource;  kwargs...)

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
        SetParticleEnergy(gun, 2.615GeV)
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

#=
@with_kw mutable struct GeneratorC3aData <: G4JLGeneratorData
    gun::Union{Nothing, CxxPtr{G4ParticleGun}} = nothing
    ion::Union{Nothing, CxxPtr{G4ParticleDefinition}} = nothing
    Z::Int64 = 82
    A::Int64 = 212
    ionCharge::Float64 = 0eplus
    excitEnergy::Float64 = 0keV
    position::G4ThreeVector = G4ThreeVector(5cm, 0cm, 5cm)
    direction::G4ThreeVector = G4ThreeVector(-1, 0, 0)
end
function GeneratorC3a(;kwargs...)
    data = GeneratorC3aData(;kwargs...)
    function _init(data::GeneratorC3aData, ::Any)
        gun = data.gun = move!(G4ParticleGun())
        data.direction = direction = G4ThreeVector(-10,rand()*2 - 1,rand()*2 - 1)
        SetParticlePosition(gun, data.position)
        SetParticleMomentumDirection(gun, data.direction)
        SetParticleEnergy(gun, 0eV)
    end
    function _gen(evt::G4Event, data::GeneratorC3aData)::Nothing
        if isnothing(data.ion)  # late initialize (after physics processes)
            ion = data.ion = GetIon(data.Z, data.A, data.excitEnergy)
            SetParticleDefinition(data.gun, ion)
            SetParticleCharge(data.gun, data.ionCharge)
        end
        #position = data.position + G4ThreeVector((rand()-0.5)*1cm, (rand()-0.5)*1cm, (rand()-0.5)*1cm)
        #SetParticlePosition(data.gun, position)
        data.direction = direction = G4ThreeVector(-10,rand()*2 - 1,rand()*2 - 1)
        SetParticleMomentumDirection(data.gun, direction)
        #data.direction = direction
        SetParticlePosition(data.gun, data.position)
        GeneratePrimaryVertex(data.gun, CxxPtr(evt))
    end
    G4JLPrimaryGenerator("GeneratorC3a", data; init_method=_init, generate_method=_gen)
end
=#



