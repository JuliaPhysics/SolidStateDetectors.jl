using SolidStateDetectors
using SolidStateDetectors: Cylindrical, Cartesian, World, Simulation, SSDFloat
using SolidStateDetectors.ConstructiveSolidGeometry: CSGUnion, CSGDifference, CSGIntersection, 
    Box, Cone, Ellipsoid, Torus, RegularPrism, Polycone, AbstractVolumePrimitive, AbstractGeometry,
    AbstractConstructiveGeometry, CylindricalPoint, origin, rotation
using IntervalSets
using OrderedCollections: OrderedDict
using Unitful
using Rotations
using LightXML


# Returns position vector for Primitive relative to origin
@inline parse_origin(e::AbstractVolumePrimitive) = origin(e)

# Returns position vector for leftmost Primitive in geometry tree
@inline parse_origin(e::AbstractConstructiveGeometry) = parse_origin(has_volume(e.a) ? e.a : e.b)

# Returns rotation matrix for Primitive relative to standard basis
@inline parse_rotation(e::AbstractVolumePrimitive) = rotation(e)

# Returns rotation matrix for leftmost Primitive in geometry tree
@inline parse_rotation(e::AbstractConstructiveGeometry) = parse_rotation(has_volume(e.a) ? e.a : e.b)


# Add <position> to <define> section, referenced in the geometry definition (in <solids>) via the name
function create_position(e::AbstractConstructiveGeometry, x_define::XMLElement, name::String, v::Bool)
    # Calculate relative position of the volumes contained in the parent volume
    x, y, z = parse_origin(e.b) - parse_origin(e.a)
   
    xpos = new_child(x_define, "position")
    set_attributes(xpos, OrderedDict(
        "name" => name,
        "x" => x,
        "y" => y,
        "z" => z,
        "unit" => SolidStateDetectors.internal_length_unit
    ))
end


# Calculates position of world relative to [0, 0, 0], adds it to <define> section
# Cylindrical grid
function create_position(e::World{T, 3, Cylindrical}, x_define::XMLElement, name::String, v::Bool) where {T}
    x_pos = new_child(x_define, "position")
    set_attributes(x_pos, OrderedDict(
        "name" => name,
        # "x" => zero(T),
        # "y" => zero(T),
        "z" => -IntervalSets.mean(e.intervals[3]),
        "unit" => SolidStateDetectors.internal_length_unit,
    ))
end
# Cartesian grid
function create_position(e::World{T, 3, Cartesian}, x_define::XMLElement, name::String, v::Bool) where {T}
    x_pos = new_child(x_define, "position")
    set_attributes(x_pos, OrderedDict(
        "name" => name,
        "x" => -IntervalSets.mean(e.intervals[1]),
        "y" => -IntervalSets.mean(e.intervals[2]),
        "z" => -IntervalSets.mean(e.intervals[3]),
        "unit" => SolidStateDetectors.internal_length_unit
    ))
end
    

# Add <rotation> to <define> section, referenced in geometry definition (in <solids>) via the name
function create_rotation(e::AbstractConstructiveGeometry, x_define::XMLElement, name::String, v::Bool)
    R1 = RotMatrix{3}(parse_rotation(e.a))
    R2 = RotMatrix{3}(parse_rotation(e.b))

    R_rel = transpose(R1) * R2

    m = RotXYZ(R_rel)

    x_rot = new_child(x_define, "rotation")
    set_attributes(x_rot, OrderedDict(
        "name" => name,
        "x" => m.theta1,
        "y" => m.theta2,
        "z" => m.theta3,
        "unit" => "rad" 
    ))
end



# Multiple ways to define tubes/cones => needs parsing
@inline parse_tube_radius(e::Cone{T, <:Any, T}) where {T} = (zero(T), e.r)
@inline parse_tube_radius(e::Cone{T, <:Any, Tuple{T, T}}) where {T} = (e.r[1], e.r[2])
@inline parse_cone_radius(e::Cone{T, <:Any, Tuple{Tuple{T}, Tuple{T}}}) where {T} = (zero(T), e.r[1][1], zero(T), e.r[2][1])
@inline parse_cone_radius(e::Cone{T, <:Any, Tuple{Tuple{T, T}, Tuple{T, T}}}) where {T} = (e.r[1][1], e.r[1][2], e.r[2][1], e.r[2][2])


# Checks whether given geometry object has required dimension values (zero volume not supported in Geant4)
function has_volume(e::Box, v::Bool = false)
    if e.hX > 0 && e.hY > 0 && e.hZ > 0
        return true
    else
        v && @warn "Box: All dimensions must be positive"
        return false
    end
end

function has_volume(e::Cone{T,<:Any,TR}, v::Bool = false) where {T, TR <: Union{T, Tuple{T,T}}}
    rmin, rmax = parse_tube_radius(e)
    
    if rmin >= rmax
        v && @warn "Cone: Outer radius must be strictly bigger than inner radius"
        return false
    elseif e.hZ <= 0
        v && @warn "Cone: The height must be a positive number"
        return false
    end

    return true
end

function has_volume(e::Cone{T,<:Any,TR}, v::Bool = false) where {T, TR}    
    rmin1, rmax1, rmin2, rmax2 = parse_cone_radius(e)
    
    if rmin1 > rmax1 || rmin2 > rmax2 || (rmin1 == rmax1 && rmin2 == rmax2)
        v && @warn "Cone: The outer radii must be strictly bigger than the inner radii"
        return false
    elseif e.hZ <= 0
        v && @warn "Cone: The height must be a positive number"
        return false
    else 
        return true
    end
end

function has_volume(e::Ellipsoid{T,<:Any,T}, v::Bool = false) where {T}
    if e.r <= 0
        v && @warn "Ellipsoid: Radius must be a positive number"
        return false
    else
        return true
    end
end

function has_volume(e::Torus, v::Bool = false)
    if isa(e.r_tube, Tuple)
        if e.r_tube[1] >= e.r_tube[2]
            v && @warn "Torus: Outer tube radius must be bigger than the inner tube radius"
            return false
        elseif e.r_tube[1] < 0 || e.r_tube[2] < 0
            v && @warn "Torus: Inner tube radius must be non-negative, Outer tube radius must be positive"
            return false
        end
    elseif e.r_tube <= 0
        v && @warn "Torus: Tube radius must be a positive number"
        return false
    end
    if e.r_torus < 0
        v && @warn "Torus: Torus radius cannot be negative"
        return false
    end
    return true
end

function has_volume(e::RegularPrism, v::Bool = false)
    if e.r <= 0
        v && @warn "Prism: Radius of prism must be a positive number"
        return false
    end
    if e.hZ <= 0
        v && @warn "Prism: Height of prism must be a positive number"
        return false
    end
    return true
end

function has_volume(p::Polycone, v::Bool = false)
    if isapprox(PolygonOps.area(tuple.(p.r, p.z)), 0)
        v && @warn "Polycone: The points passed result in a zero-volume Polycone"
        return false
    end
    return true
end


# Check whether Boolean solid has any volume at all
function has_volume(e::AbstractConstructiveGeometry, v::Bool = false)
    return has_volume(e.a, v) || has_volume(e.b, v)
end


# Builds <solids> section by recursively iterating through geometry tree (DFS)

# Default function to call when parsing not yet implemented
@inline function parse_geometry(e, x_solids::XMLElement, x_define::XMLElement, id::Integer, pf::AbstractString, v::Bool)
    throw(ArgumentError("$(typeof(e).name.name) not implemented yet"))
end

@inline CSGtoGDML(::CSGUnion) = "union"
@inline CSGtoGDML(::CSGDifference) = "subtraction"
@inline CSGtoGDML(::CSGIntersection) = "intersection"
function parse_geometry(e::AbstractConstructiveGeometry, x_solids::XMLElement, x_define::XMLElement, id::Integer, pf::AbstractString, v::Bool)::Nothing
    # If composite geometry contains children with zero volume, omit its geometry definition
    # Omit the construction of a Boolean solid if some children have zero volume
    hva = has_volume(e.a, v)
    hvb = has_volume(e.b, v)
    if hva && hvb
        parse_geometry(e.a, x_solids, x_define, 2 * id, pf, v)
        parse_geometry(e.b, x_solids, x_define, 2 * id + 1, pf, v)
        
        y = new_child(x_solids, CSGtoGDML(e))
        set_attributes(y, OrderedDict("name" => pf * string(id)))
        
        z_1 = new_child(y, "first")
        z_2 = new_child(y, "second")
        set_attribute(z_1, "ref", pf * string(2 * id))
        set_attribute(z_2, "ref", pf * string(2 * id + 1))

        if !iszero(parse_origin(e))
            z_pos = new_child(y, "positionref")
            posname = pf * "pos_" * string(id)
            set_attribute(z_pos, "ref", posname)
            create_position(e, x_define, posname, v)
        end 

        if !(parse_rotation(e.a) ≈ parse_rotation(e.b))
            z_rot = new_child(y, "rotationref")
            rotname = pf * "rot_" * string(id)
            set_attribute(z_rot, "ref", rotname)
            create_rotation(e, x_define, rotname, v)
        end
    elseif hva
        parse_geometry(e.a, x_solids, x_define, id, pf, v)
    elseif hvb
        parse_geometry(e.b, x_solids, x_define, id, pf, v)
    end
    nothing
end

function parse_geometry(e::Box, x_solids::XMLElement, x_define::XMLElement, id::Integer, pf::AbstractString, v::Bool)::Nothing
    if has_volume(e, v)
        y = new_child(x_solids, "box")
        set_attributes(y, OrderedDict(
            "name" => pf * string(id),
            "x" => 2 * e.hX,
            "y" => 2 * e.hY,
            "z" => 2 * e.hZ,
            "lunit" => SolidStateDetectors.internal_length_unit
        ))
    end
    nothing
end

@inline parse_φ(::Type{T}, φ::Nothing) where {T} = (T(360), "deg")
@inline parse_φ(::Type{T}, φ::T) where {T} = (φ, string(SolidStateDetectors.internal_angle_unit))
@inline parse_φ(::Type, φ) = throw(ArgumentError("Cone: the type of φ is unexpected"))
function parse_geometry(e::Cone{T,<:Any, TR}, x_solids::XMLElement, x_define::XMLElement, id::Integer, pf::AbstractString, v::Bool)::Nothing where {T <: SSDFloat, TR <: Union{T, Tuple{T,T}}}
    if has_volume(e, v)
    
        y = new_child(x_solids, "tube")

        rmin, rmax = parse_tube_radius(e)
        phi, aunit = parse_φ(T, e.φ)
        
        set_attributes(y, OrderedDict(
            "name" => pf * string(id),
            "rmin" => rmin,
            "rmax" => rmax,
            "z" => 2 * e.hZ,
            "startphi" => zero(T),
            "deltaphi" => phi,
            "lunit" => SolidStateDetectors.internal_length_unit,
            "aunit" => aunit
        ))
    end
    nothing
end

function parse_geometry(e::Cone{T,<:Any, TR}, x_solids::XMLElement, x_define::XMLElement, id::Integer, pf::AbstractString, v::Bool)::Nothing where {T <: SSDFloat, TR <: Union{Tuple{Tuple{T}, Tuple{T}}, Tuple{Tuple{T, T}, Tuple{T, T}}}}
    if has_volume(e, v)
    
        y = new_child(x_solids, "cone")
        
        rmin1, rmax1, rmin2, rmax2 = parse_cone_radius(e)
        phi, aunit = parse_φ(T, e.φ)
        
        set_attributes(y, OrderedDict(
            "name" => pf * string(id),
            "rmin1" => rmin1,
            "rmax1" => rmax1,
            "rmin2" => rmin2,
            "rmax2" => rmax2,
            "z" => 2 * e.hZ,
            "startphi" => zero(T),
            "deltaphi" => phi,
            "lunit" => SolidStateDetectors.internal_length_unit,
            "aunit" => aunit
        ))
    end
    nothing
end

function parse_geometry(e::Ellipsoid{T,<:Any, TR, Nothing, Nothing}, x_solids::XMLElement, x_define::XMLElement, id::Integer, pf::AbstractString, v::Bool)::Nothing where {T, TR}
    if has_volume(e, v)
       
        y = new_child(x_solids, "orb")
        
        r = if TR == T
            e.r
        end
        isnothing(r) && throw(ArgumentError("Ellipsoid: the type of r is unexpected"))
        
        # if typeof(e.φ) == Nothing
        # end
        
        # if typeof(e.θ) == Nothing
        # end
        
        set_attributes(y, OrderedDict(
            "name" => pf * string(id),
            "r" => r,
            "lunit" => SolidStateDetectors.internal_length_unit
        ))
    end
    nothing
end

function parse_geometry(e::Torus{T,<:Any,<:Any,Nothing,Nothing}, x_solids::XMLElement, x_define::XMLElement, id::Integer, pf::AbstractString, v::Bool)::Nothing where {T}
    
    # if e.r_torus == 0, the torus is basically a sphere => convert it to a sphere
    if iszero(e.r_torus)
        parse_geometry(
            if isa(e.r_tube, Tuple)
                CSGDifference(
                    Ellipsoid{T}(r = e.r_tube[2], origin = e.origin, rotation = e.rotation), 
                    Ellipsoid{T}(r = e.r_tube[1], origin = e.origin, rotation = e.rotation)
                )
            else
                Ellipsoid{T}(r = e.r_tube, origin = e.origin, rotation = e.rotation)
            end,
            x_solids, x_define, id, pf, v
        )

    elseif has_volume(e, v)
        y = new_child(x_solids, "torus")
        
        rmin, rmax = isa(e.r_tube, Tuple) ? e.r_tube : (0, e.r_tube)

        set_attributes(y, OrderedDict(
            "name" => pf * string(id),
            "rmin" => rmin,
            "rmax" => rmax,
            "rtor" => e.r_torus,
            "startphi" => 0,
            "deltaphi" => 360,
            "lunit" => SolidStateDetectors.internal_length_unit,
            "aunit" => "deg"
        ))
    end
end

function _rectangle_point(x::T)::Tuple{T,T} where {T <: Real}
    return (
        clamp(tand(abs(mod(x-90,360)-180)-90),-1,1),
        clamp(tand(abs(mod(   x,360)-180)-90),-1,1)
    )
end

function parse_geometry(e::Torus{T,<:Any,<:Any,Nothing,Tuple{T,T}}, x_solids::XMLElement, x_define::XMLElement, id::Integer, pf::AbstractString, v::Bool)::Nothing where {T}
    if has_volume(e, v)
        
        theta1::T, theta2::T = rad2deg.(e.θ)
        rmax = isa(e.r_tube, Tuple) ? e.r_tube[2] : e.r_tube

        parse_geometry( 
            CSGIntersection(
                # create full-θ Torus
                Torus(r_torus = e.r_torus, r_tube = e.r_tube, φ = e.φ, θ = nothing, origin = e.origin, rotation = e.rotation),
                
                # create Polycone to intersect with
                begin
                    points = Tuple{T,T}[(zero(T), zero(T))]
                    push!(points, _rectangle_point(theta2))
                    tmp::T = theta2 - mod(theta2 - 45, 90)
                    while theta1 < tmp && tmp > 0
                        push!(points, _rectangle_point(tmp))
                        tmp -= 90
                    end
                    push!(points, _rectangle_point(theta1))
                    push!(points, (zero(T), zero(T)))
                    Polycone(
                        r = getindex.(points,1) .* rmax .+ e.r_torus,
                        z = getindex.(points,2) .* rmax,
                        origin = e.origin,
                        rotation = e.rotation
                    )
                end
            ),
            x_solids, x_define, id, pf, v
        )
    end
end


function parse_geometry(e::RegularPrism{N,T,<:Any}, x_solids::XMLElement, x_define::XMLElement, id::Integer, pf::AbstractString, v::Bool)::Nothing where {T, N}
    if has_volume(e, v)
        y = new_child(x_solids, "polyhedra")
        
        r = e.r
        hZ = e.hZ
        set_attributes(y, OrderedDict(
            "name" => pf * string(id),
            "startphi" => 0,
            "deltaphi" => 360,
            "numsides" => N,
            "aunit" => "deg",
            "lunit" => SolidStateDetectors.internal_length_unit
        ))
        
        z_top = new_child(y, "zplane")
        set_attributes(z_top, OrderedDict(
            "rmax" => r,
            "z" => hZ
        ))

        z_bottom = new_child(y, "zplane")
        set_attributes(z_bottom, OrderedDict(
            "rmax" => r,
            "z" => -hZ
        ))
    end
end

function parse_geometry(p::Polycone, x_solids::XMLElement, x_define::XMLElement, id::Integer, pf::AbstractString, v::Bool)::Nothing
    if has_volume(p, v)
        y = new_child(x_solids, "genericPolycone")

        for i in eachindex(p.r)
            rzpoint = new_child(y, "rzpoint")
            set_attributes(rzpoint, OrderedDict(
                "r" => p.r[i],
                "z" => p.z[i]
            ))
        end

        set_attributes(y, OrderedDict(
            "name" => pf * string(id),
            "startphi" => 0,
            "deltaphi" => 360,
            "lunit" => SolidStateDetectors.internal_length_unit,
            "aunit" => "deg"
        ))
    end
end


function parse_geometry(e::World{T, 3, Cylindrical}, x_solids::XMLElement, x_define::XMLElement, id::Integer, pf::AbstractString, v::Bool)::Nothing where {T}
    w = new_child(x_solids, "tube")
    rmin, rmax = endpoints(e.intervals[1])
    z = endpoints(e.intervals[3])
    
    !iszero(rmin) && throw(ArgumentError("The r-interval of the world must start at 0."))
    
    name = pf * string(id)
    set_attributes(w, OrderedDict(
        "name" => name,
        "rmin" => zero(T),
        "rmax" => rmax,
        "z" => 2*max(abs.(z)...),
        "startphi" => zero(T),
        "deltaphi" => T(360),
        "lunit" => SolidStateDetectors.internal_length_unit,
        "aunit" => "deg"
    ))
    
    create_position(e, x_define, name, v)
    nothing
end

function parse_geometry(e::World{T, 3, Cartesian}, x_solids::XMLElement, x_define::XMLElement, id::Integer, pf::AbstractString, v::Bool)::Nothing where {T}
    w = new_child(x_solids, "box")
    x, y, z = endpoints.(e.intervals)
    
    name = pf * string(id)
    set_attributes(w, OrderedDict(
        "name" => name,
        "x" => 2*max(abs.(x)...),
        "y" => 2*max(abs.(y)...),
        "z" => 2*max(abs.(z)...),
        "lunit" => SolidStateDetectors.internal_length_unit
    ))
    
    create_position(e, x_define, name, v)
    nothing
end


# Adds <structure> section for each geometry object from <solids> and references them 
function create_volume(x_structure::XMLElement, pf::String, material_ref::AbstractString)::XMLElement
    x_vol = new_child(x_structure, "volume")
    set_attribute(x_vol, "name", pf)

    x_mat = new_child(x_vol, "materialref")
    set_attribute(x_mat, "ref", material_ref)

    x_solid_ref = new_child(x_vol, "solidref")
    # The id of the top level parents is 1
    set_attribute(x_solid_ref, "ref", pf * "_1")
    
    return x_vol
end


# Places solid "name" inside world volume
function add_to_world(x_world::XMLElement, x_define::XMLElement, e::AbstractGeometry, name::AbstractString)::Nothing
    x_vol_ref = new_child(x_world, "volumeref")
    set_attribute(x_vol_ref, "ref", name)
    
    x_pos_ref = new_child(x_world, "positionref")
    set_attribute(x_pos_ref, "ref", "pos_" * name)
    
    x_pos_ref = new_child(x_define, "position")
    parsed_e = parse_origin(e)
    x, y, z = parsed_e.x, parsed_e.y, parsed_e.z
    set_attributes(x_pos_ref, OrderedDict(
        "name" => "pos_" * name,
        "x" => x,
        "y" => y,
        "z" => z,
        "unit" => SolidStateDetectors.internal_length_unit
    ))
    nothing
end


@inline function parse_material(material::String)
    if material == "High Purity Germanium" return "G4_Ge"
    elseif material == "Vacuum" return "G4_Vacuum"
    elseif material == "Silicon" return "G4_Si"
    elseif material == "Aluminium" return "G4_Al"
    elseif material == "liquid Argon" return "G4_lAr"
    elseif material == "Copper" return "G4_Cu"
    elseif material == "Polyethylene Naphthalate" return "G4_PEN"
    elseif material == "Polytetrafluorethylen" return "G4_PTFE"
    elseif material == "Cadmium zinc telluride" return "G4_CdZnTe"
    elseif material == "CsPbBr3" return "G4_CsPbBr3"
    elseif material == "Lead" return "G4_Pb"
    end
    throw(ArgumentError("Material characteristics for \"$(material)\" not defined yet"))
end
