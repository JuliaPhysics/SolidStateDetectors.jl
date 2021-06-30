@recipe function f( ef::ElectricField{T, 3, Cylindrical};
                    r = missing,
                    φ = missing,
                    z = missing,
                    contours_equal_potential=false,
                    full_det = false ) where {T}

    g::Grid{T, 3, Cylindrical} = ef.grid
    ef_magn  = norm.(ef)

    seriescolor --> :inferno
    cross_section::Symbol, idx::Int, value::T = get_crosssection_idx_and_value(g, r, φ, z)
    cross_section != :r ? aspect_ratio --> 1 : nothing

    title --> "Electric Field (Magn.) @ $(cross_section) = $(round(value,sigdigits=2))"*(cross_section == :φ ? "°" : "m")
    colorbar_title --> "Electric Field Strength in V / m"

    ElectricPotential(ef_magn,g), cross_section, idx, value, contours_equal_potential, full_det

end

@recipe function f(ef::ElectricField{T, 3, Cartesian};
                    x = missing,
                    y = missing,
                    z = missing,
                    contours_equal_potential = false) where {T <: SSDFloat}

    g::Grid{T, 3, Cartesian} = ef.grid
    ef_magn  = norm.(ef)

    seriescolor --> :inferno
    cross_section::Symbol, idx::Int, value::T = get_crosssection_idx_and_value(g, x, y, z)
    aspect_ratio --> 1

    title --> "Electric Field (Magn.) @ $(cross_section) = $(round(value,sigdigits=2))m"
    colorbar_title --> "Electric Field Strength in V / m"

    ElectricPotential(ef_magn,g), cross_section, idx, value, contours_equal_potential
end



@userplot Plot_electric_fieldlines
@recipe function f(gdd::Plot_electric_fieldlines; φ = missing, r = missing, x = missing, y = missing, z = missing,
                    max_nsteps=3000,
                    sampling = 2u"mm", # Specifies in what density the Contacts are sampled to generate equally spaced surface charges. Also see spacing.
                    offset = 2u"mm", # should be at least as big as sampling. In doubt sampling can be reduced and the spacing keyword can be used to thin out the lines.
                    spacing = 1, # If, due to fine sampling, too many lines would clutter the plot, the spacing keyword allows to skip some fieldlines. Spacing = 2 means plot every second line. Spacing = 3 every third.
                    full_det = false,
                    skip_contact = 1) # Usually the "core" contact is skipped, and the other contacts are equally sampled for charges to drift
    sim = gdd.args[1]
    S = get_coordinate_system(sim.electric_field.grid)
    T = SolidStateDetectors.get_precision_type(sim.detector)

    dim_array = [φ, r, x, y, z]
    dim_symbols_array = [:φ, :r, :x, :y, :z]
    if isempty(skipmissing(dim_array))
        if S == Cartesian
            v::T = 0
            φ = 0
            dim_number = 2
            dim_symbol = :φ
        elseif S == Cylindrical
            dim_number = 1
            dim_symbol = :x
            v = sim.electric_field.grid[dim_number][ div(length(sim.electric_field.grid[dim_number]), 2) ]
        end
    elseif sum(ismissing.(dim_array)) == 4
        dim_idx = findfirst(x -> !ismissing(x), dim_array)
        dim_symbol = dim_symbols_array[dim_idx]
        v = dim_array[dim_idx]
        if dim_symbol == :φ v = deg2rad(v) end
    else
        throw(ArgumentError("Only one keyword for a certain dimension s allowed. One of 'φ', 'r', 'x', 'y', 'z'"))
    end

    show_full_det = full_det == true && dim_symbol == :φ ? true : false # The full_det keyword only makes sense for crossections in the xz plane in cylindrical grids

    (dim_symbol != :r && !(dim_symbol == :z && S == Cylindrical) ) ? aspect_ratio --> 1 : nothing
    title --> (show_full_det ? "Electric Field Lines @$(dim_symbol)=$(round(rad2deg(v),sigdigits = 3))°, (=$(round(rad2deg(T((v+π)%(2π))),sigdigits = 3))° on left side) " : "Electric Field Lines @$(dim_symbol)=$(round(dim_symbol == :φ ? rad2deg(v) : v, sigdigits=3))" * (dim_symbol == :φ ? "°" : "m"))
    xguide --> (S == Cylindrical ? (dim_symbol == :r ? "φ / rad" : "r / m") : (dim_symbol == :x ? "y / m" : "x / m"))
    yguide --> "z / m"
    (S == Cylindrical && dim_symbol == :z) ? xguide :=  "" : nothing
    (S == Cylindrical && dim_symbol == :z) ? yguide :=  "" : nothing

    contacts_to_spawn_charges_for = filter!(x -> x.id !=skip_contact, Contact{T}[c for c in sim.detector.contacts])
    spawn_positions = CartesianPoint{T}[]
    spawn_positions_mirror = CartesianPoint{T}[] # Will only be used if show_full_det == true
    grid = sim.electric_field.grid
    function get_closest_samples(pt::PT, sampled_points_pool) where PT # Looks for the closest neighbors of a sample point in pool of samplesthe three major directions and give pairs.
        pool1 = unique!(map(x->x[1], sampled_points_pool))
        pool2 = unique!(map(y->y[2], sampled_points_pool))
        pool3 = unique!(map(z->z[3], sampled_points_pool))
        idcs_closest = [searchsortednearest(pool1, pt[1]), searchsortednearest(pool2, pt[2]), searchsortednearest(pool3, pt[3])]

        pairs = [
          ( PT(pool1[maximum([idcs_closest[1]-1,1])], pool2[idcs_closest[2]], pool3[idcs_closest[3]] ), PT(pool1[minimum([idcs_closest[1]+1,length(pool1)])], pool2[idcs_closest[2]], pool3[idcs_closest[3]]) )
          ( PT(pool1[idcs_closest[1]], pool2[maximum([idcs_closest[2]-1,1])], pool3[idcs_closest[3]] ), PT(pool1[idcs_closest[1]], pool2[minimum([idcs_closest[2]+1,length(pool2)])], pool3[idcs_closest[3]]) )
          ( PT(pool1[idcs_closest[1]], pool2[idcs_closest[2]], pool3[maximum([idcs_closest[3]-1,1])] ), PT(pool1[idcs_closest[1]], pool2[idcs_closest[2]], pool3[minimum([idcs_closest[3]+1,length(pool3)])]) )
        ]
    end
    symbol2dimnumber_dict = Dict(:r => 1, :x => 1, :y => 2, :φ => 2, :z => 3)
    symbol2dimnumber(dim_symbol) = symbol2dimnumber_dict[dim_symbol] # give the corresponding dimension for each cross section symbol
    dim_number = symbol2dimnumber(dim_symbol)

    function get_vector_in_Xsec_dir(sample_points, sample_pool, dim_symbol) # get a vetor that points in the direction of the cross section
        local_samples_in_Xsec_dir = map(x->get_closest_samples(x,sample_pool)[symbol2dimnumber(dim_symbol)], sample_points)
        vec_in_Xsec_dir = map(x->CartesianPoint(x[2]) - CartesianPoint(x[1]), local_samples_in_Xsec_dir)
    end
    function get_offset_vector(sampled_points_Xsec_cart, vec_in_Xsec_dir) # Get the vector of two neighboring samples and take the cross product with the vector in direction of the xsection to generate a orthogonal offset vector.
        sample_forward_vec = CartesianPoint{T}[]
        for i in 1:length(sampled_points_Xsec_cart)
            i!=length(sampled_points_Xsec_cart) ? push!(sample_forward_vec, sampled_points_Xsec_cart[i+1] - sampled_points_Xsec_cart[i]) : push!(sample_forward_vec, sampled_points_Xsec_cart[i] - sampled_points_Xsec_cart[i-1])
        end
        offset_vector = CartesianVector{T}.(cross.(sample_forward_vec, vec_in_Xsec_dir))./norm.(cross.(sample_forward_vec, vec_in_Xsec_dir)) # calculate normaliezed cross product
    end

    for c in contacts_to_spawn_charges_for[:] # Main loop to generate the close-to-the-surface charges
        for positive_geometry in c.geometry_positive[:]
            sample_dummy = SolidStateDetectors.sample(positive_geometry, T[1000,1000,1000]) # used to extract the point type of the volume primitive
            PT = eltype(sample_dummy)
            sampling_vector_pool = T.(ustrip.([uconvert(u"m", sampling) for i in 1:3]))
            PT == CylindricalPoint{T} ? sampling_vector_pool[2] = sampling_vector_pool[2]/0.001 *2*π / 360 : nothing # rough translation of mm to radians; Might need some polish
            sample_pool = SolidStateDetectors.sample(positive_geometry, sampling_vector_pool)
            if length(sample_pool)<8
                @warn("The sampling step is to coarse ($(sampling)). Please use the 'sampling' keyword (e.g. 'sampling  = 0.1u\"mm\"') to specify smaller sampling steps. Also consider to reduce the 'offset' keyword accordingly. Attempting to automatically reduce the sampling size.")
                exponent = 1
                while length(unique!(map(x->x[dim_number],sample_pool))) < 8
                    exponent+=1
                    sampling_vector_pool ./=2^exponent
                    # sample_pool = SolidStateDetectors.sample(positive_geometry, sampling_vector_pool./(2^exponent))
                    sample_pool = SolidStateDetectors.sample(positive_geometry, sampling_vector_pool)
                end
                @info("Scaling down sampling steps by a factor of $(2^exponent). Now using sampling steps of $(sampling_vector_pool./2^exponent) m. Also scaling down offset by a factor $(2^exponent).")
                offset /= 2^exponent
            end
            # sample_pool = S == Cylindrical ? geom_round.(CylindricalPoint.(sample_pool)) : geom_round.(sample_pool)
            # sample_pool = S == Cartesian ? geom_round.(CartesianPoint.(sample_pool)) : geom_round.(sample_pool)
            sample_pool = S == Cylindrical ? CylindricalPoint.(sample_pool) : sample_pool
            sample_pool = S == Cartesian ? CartesianPoint.(sample_pool) : sample_pool
            sampled_planes = unique!(map(x->x[dim_number],sample_pool))
            v_Xsec_plane = sampled_planes[searchsortednearest(sampled_planes,v)]
            if abs(v-v_Xsec_plane) > sampling_vector_pool[dim_number] continue; end

            # v_Xsec_plane_mirror = show_full_det ? sampled_planes[searchsortednearest(sampled_planes,v+π)] : v_Xsec_plane
            sampled_points_Xsec = sample_pool[findall(x-> (abs(x[dim_number] - v_Xsec_plane)<sampling_vector_pool[dim_number]/3), sample_pool)]
            sampled_points_Xsec_cart = CartesianPoint.(sampled_points_Xsec)

            vec_in_Xsec_dir = get_vector_in_Xsec_dir(sampled_points_Xsec, sample_pool, dim_symbol)
            offset_vector::Vector{CartesianVector{T}} = get_offset_vector(sampled_points_Xsec_cart, vec_in_Xsec_dir) .* T.(ustrip.(uconvert(u"m", offset)))
            push!(spawn_positions, broadcast(+, sampled_points_Xsec_cart, offset_vector)...)
            push!(spawn_positions, broadcast(-, sampled_points_Xsec_cart, offset_vector)...)

            if show_full_det # generate a seperate pool of charges for the + 180 deg direction. This is necessary to later know which drift paths have to be mirrored.
                v_Xsec_plane_mirror::T = sampled_planes[searchsortednearest(sampled_planes,T((v+π)%(2π)))]
                if abs(T((v+π)%(2π)) - v_Xsec_plane_mirror) > sampling_vector_pool[dim_number] continue; end

                sampled_points_Xsec_mirror = sample_pool[findall(x-> (abs(x[dim_number] - v_Xsec_plane_mirror)<sampling_vector_pool[dim_number]/3), sample_pool)]
                sampled_points_Xsec_cart_mirror = CartesianPoint.(sampled_points_Xsec_mirror)

                vec_in_Xsec_dir_mirror = get_vector_in_Xsec_dir(sampled_points_Xsec_mirror, sample_pool, dim_symbol)
                offset_vector_mirror::Vector{CartesianVector{T}} = get_offset_vector(sampled_points_Xsec_cart_mirror, vec_in_Xsec_dir_mirror) .* T.(ustrip.(uconvert(u"m", offset)))

                push!(spawn_positions_mirror, broadcast(+, sampled_points_Xsec_cart_mirror, offset_vector_mirror)...)
                push!(spawn_positions_mirror, broadcast(-, sampled_points_Xsec_cart_mirror, offset_vector_mirror)...)
            end
        end
    end

    filter!(x -> x in sim.detector && !in(x, sim.detector.contacts), spawn_positions) # get rid of unneccessary spawnpositions
    show_full_det ? filter!(x ->x in sim.detector && !in(x, sim.detector.contacts), spawn_positions_mirror) : nothing
    show_full_det ? spawn_positions = vcat(spawn_positions,spawn_positions_mirror) : nothing
    #spawn_positions = unique!(geom_round.(spawn_positions))
    @info "$(round(Int,length(spawn_positions)/spacing)) drifts are now being simulated..."

    el_field_itp     = get_interpolated_drift_field(sim.electric_field.data      , sim.electric_field.grid) #Interpolate the Electric fields, in which the charges will drift. Is passed to the drift function.
    el_field_itp_inv = get_interpolated_drift_field(sim.electric_field.data .* -1, sim.electric_field.grid)

    @showprogress for (ipos, pos) in enumerate(spawn_positions) # Charge drift and plotting loop. Not optimized for speed, but it doesnt have to be. Uses low level drift function for more contorl.
        if ((spacing-1)+ipos)%spacing == 0
            path = CartesianPoint{T}[CartesianPoint{T}(0.0,0.0,0.0) for i in 1:max_nsteps]
            _drift_charge!(path, Vector{T}(undef, max_nsteps), sim.detector, sim.point_types, sim.electric_potential.grid, CartesianPoint(pos), T(1e-9), el_field_itp, verbose = false )
            filter!(x->x != CartesianPoint{T}(0.0,0.0,0.0), path)
            @series begin
                seriescolor --> :white
                if dim_symbol == :z && S == Cylindrical projection --> :polar end
                label --> ""
                x, y = if dim_symbol == :φ
                    map(x -> (pos in spawn_positions_mirror ? -1 : 1) * sqrt(x[1]^2+x[2]^2), path), # mirror the path for the charges that originate from v + 2 \pi
                    map(x -> x[3], path)
                elseif dim_symbol == :x
                    map(x -> x[2], path),
                    map(x -> x[3], path)
                elseif dim_symbol == :y
                    map(x -> x[1], path),
                    map(x -> x[3], path)
                elseif dim_symbol == :z
                    if S == Cylindrical
                        path = CylindricalPoint.(path)
                        map(x -> x[2], path),
                        map(x -> x[1], path)
                    else
                        map(x -> x[1], path),
                        map(x -> x[2], path)
                    end
                elseif dim_symbol == :r
                    path = CylindricalPoint.(path)
                    map(x -> x[2], path),
                    map(x -> x[3], path)
                end
                x, y
            end
            # Repeat for the other charge carrier type
            path = CartesianPoint{T}[CartesianPoint{T}(0.0,0.0,0.0) for i in 1:max_nsteps]
            _drift_charge!(path, Vector{T}(undef, max_nsteps), sim.detector, sim.point_types, sim.electric_potential.grid, CartesianPoint(pos), T(2e-9), el_field_itp_inv, verbose = false )
            filter!(x->x != CartesianPoint{T}(0.0,0.0,0.0), path)
            @series begin
                seriescolor --> :white
                if dim_symbol == :z && S == Cylindrical projection --> :polar end
                label --> ""
                x, y = if dim_symbol == :φ
                    map(x -> (pos in spawn_positions_mirror ? -1 : 1) * sqrt(x[1]^2+x[2]^2), path),
                    map(x -> x[3], path)
                elseif dim_symbol == :x
                    map(x -> x[2], path),
                    map(x -> x[3], path)
                elseif dim_symbol == :y
                    map(x -> x[1], path),
                    map(x -> x[3], path)
                elseif dim_symbol == :z
                    if S == Cylindrical
                        path = CylindricalPoint.(path)
                        ylims --> (0.0, sim.world.intervals[1].right)
                        map(x -> x[2], path),
                        map(x -> x[1], path)
                    else
                        map(x -> x[1], path),
                        map(x -> x[2], path)
                    end
                elseif dim_symbol == :r
                    ylims --> (sim.world.intervals[3].left, sim.world.intervals[3].right)
                    path = CylindricalPoint.(path)
                    map(x -> x[2], path),
                    map(x -> x[3], path)
                end
                x, y
            end

        end
    end
end
