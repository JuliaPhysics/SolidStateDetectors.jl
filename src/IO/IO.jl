"""
    write_to_hdf5(output, ep::ElectricPotential)

`output` should be and `HDF5.HDF5Group`.
"""
function write_to_hdf5(output, ep::ElectricPotential)
    output["cyclic"] = [ep.cyclic]
    output["r"] = ep.r
    output["phi"] = ep.φ
    output["z"] = ep.z
    output["potential"] = ep.potential
    nothing
end

"""
    read_from_hdf5(input, ::Type{ElectricPotential})::ElectricPotential

`input` should be and `HDF5.HDF5Group`.
"""
function read_from_hdf5(input, ::Type{ElectricPotential})::ElectricPotential
    T::Type = eltype(input["cyclic"][:])
    ep = ElectricPotential{T}(
        input["cyclic"][:][1],
        input["r"][:],
        input["phi"][:],
        input["z"][:],
        only_2d = length(input["phi"][:]) == 1
    )
    ep.potential = input["potential"][:,:,:]
    return ep
end


"""
    write_to_hdf5(output, ep::PointTypes)

`output` should be and `HDF5.HDF5Group`.
"""
function write_to_hdf5(output, pts::PointTypes)
    output["r"] = pts.r
    output["phi"] = pts.φ
    output["z"] = pts.z
    output["pointtypes"] = pts.pointtypes
    nothing
end

"""
    read_from_hdf5(input, ::Type{PointTypes})::PointTypes

`input` should be and `HDF5.HDF5Group`.
"""
function read_from_hdf5(input, ::Type{PointTypes})::PointTypes
    T::Type = eltype(input["r"][:])
    pts = PointTypes{T}(
        input["r"][:],
        input["phi"][:],
        input["z"][:],
        input["pointtypes"][:,:,:]
    )
    return pts
end
