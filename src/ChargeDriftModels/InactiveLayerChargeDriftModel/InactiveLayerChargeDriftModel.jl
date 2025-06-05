"""
    struct InactiveLayerChargeDriftModel{T <: SSDFloat} <: AbstractChargeDriftModel{T}
        
Charge drift model in which the electrons and holes drift along the electric field.
There factors are considered in the mobility calculation: ionized impurities, neutral impurities, and acoustic phonon.

## Fields
- `calculate_mobility::Function`: Mobility calculation function
"""
abstract type VolumeType end
abstract type BulkVolume <: VolumeType end 
abstract type SurfaceVolume <: VolumeType end


struct InactiveLayerChargeDriftModel{T <: SSDFloat} <: AbstractChargeDriftModel{T}
    calculate_mobility::Function
    pn_depth::T
    surface_imp_model::T
    pn_hole_mobility::T
    pn_electron_mobility::T
end

function calculate_mobility(
	pt::SolidStateDetectors.AbstractCoordinatePoint{T}, ::Type{Hole}) where {T}

	depth=calculate_depth2surface(pt)
	volume_type = depth>pn_depth ? BulkVolume : SurfaceVolume
	Temp = T_Deadlayer/K
	
	if (volume_type <: BulkVolume)
		u = hole_bulk_mobility

	elseif (volume_type <: SurfaceVolume)
		# pass the density to the following equation in unit cm^-3
		Nn = Neural_Impurity / cm^-3
		Ni = min(5e20*cm^-3, (p_type_density + the_li_concentration(pt))) / cm^-3 # in cm^-3 
		 
		uI = (2.35e17*Temp^1.5/Ni/log(9.13e13*Temp^2/Ni) + 1.51e18*Temp^1.5/Ni/log(5.82e14*Temp^2/Ni)) * (cm^2/V/s)
		uA = (7.77e7 * Temp^-1.5) * (cm^2/V/s)
		uN = (1/Nn * (2.31e18+2.36e20) * 0.82 * (0.228*Temp^0.5 + 0.976*Temp^-0.5)) * (cm^2/V/s)
		u = 1/(1/uI + 1/uA + 1/uN)
	end
	u
end

function calculate_mobility(
	pt::SolidStateDetectors.AbstractCoordinatePoint{T}, ::Type{Electron}) where {T}

	depth=calculate_depth2surface(pt)
	volume_type = depth>pn_depth ? BulkVolume : SurfaceVolume

	Temp = T_Deadlayer/K # pass to the following equation in unit K
	
	if (volume_type <: BulkVolume)
		# from dwh-sagepss: from ieee
		# u = 3.6e4 * (cm^2/s/V)

		# use the value at the pn point to substitute the value at the bulk (super close, save time)
		u = 34824.545  * (cm^2/s/V)

	elseif (volume_type <: SurfaceVolume)
		# from dwh-sagepss, 这个结果的范围是0.18-1.58 m2/s/V (从深度0mm-PN结)
		# A = 1.115e-7
		# B = 5.22e-18 *(cm^3)
		# Ni = max(p_type_density,the_li_concentration(pt))
		# u = 1 / (A*Temp^1.5+B*Ni*Temp^-1.5) * 1.5 * (cm^2/s/V) # FIX-ME 

		# from mdm-jinst-ntype # more ref to: #electron model from jinst-electron-MDM
		# pass to the following equation in unit cm^-3
		Nn = Neural_Impurity / cm^-3
		Ni = (p_type_density + the_li_concentration(pt)) / cm^-3 # in cm^-3 

		uI = 2.442e18*Temp^1.5/Ni/(log(2.496e14*Temp^2/Ni)) * (cm^2/V/s)
		uA = 9.32e7 * Temp^-1.5 * (cm^2/V/s)
		uN = 1.07e20/Nn * (0.28*Temp^0.5 + 0.54*Temp^-0.5) * (cm^2/V/s)
		u = 1/(1/uI + 1/uA + 1/uN)
	end
	u
end


function InactiveLayerChargeDriftModel{T}(config::AbstractDict) where {T <: SSDFloat}
    if !haskey(config, "mobilities")
        throw(ConfigFileError("InactiveLayerChargeDriftModel config file needs entry 'mobilities'."))
    end

    for axis in ("e", "h")
        if !haskey(config["mobilities"], axis)
            throw(ConfigFileError("InactiveLayerChargeDriftModel config file needs entry 'mobilities/$axis'."))
        end
    end
    
    InactiveLayerChargeDriftModel{T}(config["mobilities"]["e"], config["mobilities"]["h"])
end

@fastmath function getVe(fv::SVector{3, T}, cdm::InactiveLayerChargeDriftModel{T}, current_pos::CartesianPoint{T} = zero(CartesianPoint{T})) where {T <: SSDFloat}
    @inbounds begin 
        -cdm.calculate_mobility(current_pos, Electron)*fv
    end
end

@fastmath function getVh(fv::SVector{3, T}, cdm::InactiveLayerChargeDriftModel{T}, current_pos::CartesianPoint{T} = zero(CartesianPoint{T})) where {T <: SSDFloat}
    @inbounds begin 
        cdm.calculate_mobility(current_pos, Hole)*fv
    end
end