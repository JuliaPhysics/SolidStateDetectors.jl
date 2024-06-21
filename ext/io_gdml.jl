using SolidStateDetectors
using SolidStateDetectors: Cylindrical, Cartesian, World, Simulation
using SolidStateDetectors.ConstructiveSolidGeometry: CSGUnion, CSGDifference, CSGIntersection, 
    Box, Cone, Ellipsoid, Torus, show_CSG_tree, AbstractVolumePrimitive, AbstractGeometry,
    AbstractConstructiveGeometry, CylindricalPoint, origin, rotation
using IntervalSets
using DataStructures: OrderedDict
using Unitful
using Rotations


# position relative to origin of a primitive is stored in its parameters
@inline parse_origin(e::AbstractVolumePrimitive) = origin(e)

# position relative to the origin of a boolean solid is defined as position of
# the leftmost child in the tree structure
@inline parse_origin(e::AbstractConstructiveGeometry) = parse_origin(e.a)

# Returns euler angles of the rotation matrix
function parse_rotation_matrix(e::AbstractVolumePrimitive, v::Bool)
    m = RotXYZ(RotMatrix{3}(rotation(e)))
    return m.theta1, m.theta2, m.theta3
end

function parse_rotation_matrix(e::AbstractConstructiveGeometry, v::Bool)
    v && @warn "Combining multiple rotations is an experimental feature"
    return parse_rotation_matrix(e.a, v) .- parse_rotation_matrix(e.b, v)
end


# Add <position> to <define> section, referenced in the geometry definition (in <solids>) via the name
function create_position(e::AbstractConstructiveGeometry, x_define::XMLElement, name::String, v::Bool)
    # Calculate relative position of the volumes contained in the parent volume
    x, y, z = parse_origin(e.b) .- parse_origin(e.a)
   
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
    # 
    x, y, z = parse_rotation_matrix(e, v)
    x_rot = new_child(x_define, "rotation")
    set_attributes(x_rot, OrderedDict(
        "name" => name,
        "x" => x,
        "y" => y,
        "z" => z,
        "unit" => "rad" 
    ))
end



# Multiple ways to define tubes/cones => needs parsing
@inline parse_tube_radius(e::Cone{T, <:Any, T}) where {T} = (zero(T), e.r)
@inline parse_tube_radius(e::Cone{T, <:Any, Tuple{T, T}}) where {T} = (e.r[1], e.r[2])
@inline parse_cone_radius(e::Cone{T, <:Any, Tuple{Tuple{T}, Tuple{T}}}) where {T} = (zero(T), e.r[1][1], zero(T), e.r[2][1])
@inline parse_cone_radius(e::Cone{T, <:Any, Tuple{Tuple{T, T}, Tuple{T, T}}}) where {T} = (e.r[1][1], e.r[1][2], e.r[2][1], e.r[2][2])


# Checks whether given geometry object has required dimension values (zero volume not supported in Geant4)
function has_volume(e::Box, ::XMLElement, ::XMLElement, id::Integer, pf::AbstractString, v::Bool; parse::Bool = true)    
    if e.hX > 0 && e.hY > 0 && e.hZ > 0
        return true
    else
        v && @warn "Box $(pf * string(id)): All dimensions must be positive"
        return false
    end
end

function has_volume(e::Cone{T,<:Any,TR}, ::XMLElement, ::XMLElement, id::Integer, pf::AbstractString, v::Bool; parse::Bool = true) where {T, TR <: Union{T, Tuple{T,T}}}
    rmin, rmax = parse_tube_radius(e)
    
    if rmin >= rmax
        v && @warn "Cone $(pf * string(id)): Outer radius must be strictly bigger than inner radius"
        return false
    else
        return true
    end
end

function has_volume(e::Cone{T,<:Any,TR}, ::XMLElement, ::XMLElement, id::Integer, pf::AbstractString, v::Bool; parse::Bool = true) where {T, TR}    
    rmin1, rmax1, rmin2, rmax2 = parse_cone_radius(e)
    
    if rmin1 >= rmax1 || rmin2 >= rmax2
        v && @warn "Cone $(pf * string(id)): The outer radii must be strictly bigger than the inner radii"
        return false
    else
        return true
    end
end


function has_volume(e::Ellipsoid{T,<:Any,T}, ::XMLElement, ::XMLElement, id::Integer, pf::AbstractString, v::Bool; parse::Bool = true) where {T}
    if e.r <= 0
        v && @warn "Ellipsoid $(pf * string(id)): Radius must be a positive number"
        return false
    else
        return true
    end
end

# If composite geometry contains children with zero volume, omit its geometry definition
# Omit the construction of a Boolean solid if some children have zero volume
function has_volume(e::AbstractConstructiveGeometry, x_solids::XMLElement, x_define::XMLElement, id::Integer, pf::AbstractString, v::Bool; parse::Bool = true)
    hva = has_volume(e.a, x_solids, x_define, 2 * id, pf, v, parse = false)
    hvb = has_volume(e.b, x_solids, x_define, 2 * id + 1, pf, v, parse = false)
    if hva && hvb
        return true
    elseif parse && hva 
        parse_geometry(e.a, x_solids, x_define, 2 * id, pf, v)
    elseif parse && hvb
        parse_geometry(e.b, x_solids, x_define, 2 * id + 1, pf, v)
    end
    return false
end



# Builds <solids> section by recursively iterating through geometry tree (DFS)

# Default function to call when parsing not yet implemented
function parse_geometry(e, x_solids::XMLElement, x_define::XMLElement, id::Integer, pf::AbstractString, v::Bool)
    throw(AssertionError("$(typeof(e).name.name) not implemented yet"))
end

@inline CSGtoGDML(::CSGUnion) = "union"
@inline CSGtoGDML(::CSGDifference) = "subtraction"
@inline CSGtoGDML(::CSGIntersection) = "intersection"
function parse_geometry(e::AbstractConstructiveGeometry, x_solids::XMLElement, x_define::XMLElement, id::Integer, pf::AbstractString, v::Bool)::Nothing
    if has_volume(e, x_solids, x_define, id, pf, v)
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

        if !all(iszero.(parse_rotation_matrix(e, v)))
            z_rot = new_child(y, "rotationref")
            rotname = pf * "rot_" * string(id)
            set_attribute(z_rot, "ref", rotname)
            create_rotation(e, x_define, rotname, v)
        end
    end
    nothing
end

function parse_geometry(e::Box, x_solids::XMLElement, x_define::XMLElement, id::Integer, pf::AbstractString, v::Bool)::Nothing
    if has_volume(e, x_solids, x_define, id, pf, v)
        y = new_child(x, "box")
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

parse_φ(::Type{T}, φ::Nothing) where {T} = (T(360), "deg")
parse_φ(::Type{T}, φ::T) where {T} = (φ, string(SolidStateDetectors.internal_angle_unit))
parse_φ(::Type, φ) = throw(AssertionError("Cone: the type of φ is unexpected"))
function parse_geometry(e::Cone{T,<:Any, TR}, x_solids::XMLElement, x_define::XMLElement, id::Integer, pf::AbstractString, v::Bool)::Nothing where {T, TR <: Union{T, Tuple{T,T}}}
    if has_volume(e, x_solids, x_define, id, pf, v)
    
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

function parse_geometry(e::Cone{T,<:Any,TR, <:Any}, x_solids::XMLElement, x_define::XMLElement, id::Integer, pf::AbstractString, v::Bool)::Nothing where {T, TR <: Union{Tuple{Tuple{T}, Tuple{T}}, Tuple{Tuple{T, T}, Tuple{T, T}}}}
    if has_volume(e, x_solids, x_define, id, pf, v)
    
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
    if has_volume(e, x_solids, x_define, id, pf, v)
       
        y = new_child(x_solids, "orb")
        
        r = if TR == T
            e.r
        end
        isnothing(r) && throw(AssertionError("Ellipsoid: the type of r is unexpected"))
        
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

function parse_geometry(e::World{T, 3, Cylindrical}, x_solids::XMLElement, x_define::XMLElement, id::Integer, pf::AbstractString, v::Bool)::Nothing where {T}

    w = new_child(x_solids, "tube")
    rmin, rmax = endpoints(e.intervals[1])
    z = width(e.intervals[3])
    
    !iszero(rmin) && throw(AssertionError("The r-interval of the world must start at 0."))
    
    name = pf * string(id)
    set_attributes(w, OrderedDict(
        "name" => name,
        "rmin" => zero(T),
        "rmax" => rmax,
        "z" => z,
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
    x, y, z = width.(e.intervals)
    
    name = pf * string(id)
    set_attributes(w, OrderedDict(
        "name" => name,
        "x" => x,
        "y" => y,
        "z" => z,
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
    x, y, z = parse_origin(e)
    set_attributes(x_pos_ref, OrderedDict(
        "name" => "pos_" * name,
        "x" => x,
        "y" => y,
        "z" => z,
        "unit" => SolidStateDetectors.internal_length_unit
    ))
    nothing
end
