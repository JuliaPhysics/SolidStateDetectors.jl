"""
    readsiggen(file_path::String; T::Type)
    
Reads a [SigGen](https://github.com/radforddc/icpc_siggen) configuration file (ending in `.config`) in 'file_path'
and returns a dictionary of all parameters. Non-existing parameteres are set to `0`.

## Arguments
* `file_path::String`: File path leading to the [SigGen](https://github.com/radforddc/icpc_siggen) configuration file.

## Keywords 
* `T::Type`: Type of the parameters in the output dictionary. Default is `Float64`.
"""
function readsiggen(file_path::String; dicttype::Type = Dict{String,Any}, T::Type=Float64)::dicttype
    
    detector_name = split(basename(file_path), ".config")[1]
    
    config = dicttype("name"       => detector_name,
        "verbosity_level"      => 0,
        "xtal_length"          => 0,
        "xtal_radius"          => 0,
        "top_bullet_radius"    => 0, # not implemented in later conversions
        "bottom_bullet_radius" => 0, # not implemented in later conversions
        "pc_length"            => 0,
        "pc_radius"            => 0,
        "bulletize_PC"         => 0,
        "wrap_around_radius"   => 0,
        "ditch_depth"          => 0,
        "ditch_thickness"      => 0,
        "bottom_taper_length"  => 0,
        "hole_length"          => 0,
        "hole_radius"          => 0,
        "outer_taper_length"   => 0,
        "inner_taper_length"   => 0,
        "outer_taper_width"    => 0,
        "inner_taper_width"    => 0,
        "taper_angle"          => 0,
        "taper_length"         => 0,
        "Li_thickness"         => 0,
        "energy"               => 0,
        "xtal_grid"            => 0,
        "impurity_z0"          => 0,
        "impurity_gradient"    => 0,
        "impurity_surface"     => 0,
        "xtal_HV"              => 0,
        "max_iterations"       => 0,
        "write_field"          => 0,
        "write_WP"             => 0,
        "drift_name"           => "",
        "field_name"           => "",
        "wp_name"              => "",
        "xtal_temp"            => 0,
        "preamp_tau"           => 0,
        "time_steps_calc"      => 0,
        "step_time_calc"       => 0,
        "step_time_out"        => 0,
        "charge_cloud_size"    => 0,
        "use_diffusion"        => 0,
        "surface_drift_vel_factor" => 0);

    file = open(file_path) do f
        temp          = readlines(f)
        config_file   = dicttype()
        for line in temp
            # parse lines that are neither empty, nor consist of spaces or comments only
            if isnothing(match(r"^\s*(#.*){0,1}$", line))
                line  = strip(split(line, "#")[1])
                let s = split(line, r"\s+")
                    if length(s) == 2
                        name, value = s
                        config_file[name] = value
                    else
                        @warn "Could not parse line '$line'"
                    end
                end
            end
        end
        return config_file
    end

    config_keys = keys(config)
    file_keys   = keys(file)


    for key in file_keys
        if key in config_keys
            if key == "wp_name" || key == "drift_name" || key == "field_name"
                config[key] = string(file[key])
            elseif key == "taper_angle" && parse(Float64, file["taper_angle"]) != 0
                config["outer_taper_width"] = tan(deg2rad(parse(T, file[key]))) * parse(T, file["outer_taper_length"])
                config["inner_taper_width"] = tan(deg2rad(parse(T, file[key]))) * parse(T, file["inner_taper_length"])
                config[key]                 = parse(T, file[key])
            else
                config[key] = parse(T, file[key])
            end
        end
    end
    return config
end



"""
    siggentodict(config::AbstractDict; dicttype::Type, units::AbstractDict)
    
Converts the dictionary containing the parameters from a [SigGen](https://github.com/radforddc/icpc_siggen) configuration file
to a dictionary that can be understood by SolidStateDetectors.jl. 

## Arguments
* `config::AbstractDict`: Dictionary containing [SigGen](https://github.com/radforddc/icpc_siggen) parameters (output of `readsiggen()``).

## Keywords
* `dicttype::Type`: Type of the dictionary returned.
* `units::AbstractDict`: Units used in [SigGen](https://github.com/radforddc/icpc_siggen) configuration file (set to `"mm"`, `"deg"`, `"V"` and `"K"`).
    The dictionary needs the fields `"length"`, `"angle"`, `"potential"` and `"temperature"`.
"""
function siggentodict(config::AbstractDict;
        dicttype::Type = Dict{String, Any},
        units::AbstractDict = Dict("length" => "mm","angle" => "deg","potential" => "V","temperature" => "K"))::dicttype

    #>----------------Settings START-------------------------------------------<

    #>----------------Setting the z and r limits-------------------------------<
    z_limit      = ceil(config["xtal_length"]; sigdigits=1);
    if z_limit - config["xtal_length"] <= 0
        z_limit += 10
    end
    r_limit      = ceil(config["xtal_radius"]; sigdigits=1);
    if r_limit - config["xtal_radius"] <= 0
        r_limit += 10
    end

    #>----------------Setting the Axes & Grid dictionaries---------------------<
    axes     = dicttype("r"           => dicttype("to"         => r_limit,
                                                  "boundaries" => "inf"),
                        "phi"         => dicttype("from"       => 0,
                                                  "to"         => 0,
                                                  "boundaries" => "periodic"),
                        "z"           => dicttype("from"       => -10,
                                                  "to"         => z_limit,
                                                  "boundaries" => dicttype("left" => "inf", "right" => "inf"))
                        )

    grid     = dicttype("coordinates" => "cylindrical",
                        "axes"        => axes);


    #>----------------Settings END---------------------------------------------<
    #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
    #>----------------Volume START---------------------------------------------<


    #>----------------Main Volume -> Initial Cylinder--------------------------<
    # r and L are given by the input file. Cylindrical symmetry is set to default
    main_vol = dicttype("translate"   => dicttype("z"    => config["xtal_length"]/2,
                        "tube"        => dicttype(
                        "name"        => "Initial Cylinder",
                        "r"           => dicttype("from" => 0.0, "to" => config["xtal_radius"]),
                        "phi"         => dicttype("from" => 0.0, "to" => 360.0),
                        "h"           => config["xtal_length"])));

    #>----------------Upper Cone / Taper on top--------------------------------<
    # If not existing, everything will be set to 0
    upper_cone = dicttype("translate" => dicttype("z"      => config["xtal_length"] -
                                            config["outer_taper_length"]/2,
                        "cone"        => dicttype(
                        "name"        => "Upper Cone",
                        "r"           => dicttype("bottom" => dicttype("from" => config["xtal_radius"],
                                                               "to"   => ceil(config["xtal_radius"] + 1.0)),
                                              "top"    => dicttype("from" => config["xtal_radius"] -
                                                                            config["outer_taper_width"],
                                                               "to"   => ceil(config["xtal_radius"] + 1.0))),
                        "phi"         => dicttype("from"   => 0.0, "to" => 360.0),
                        "h"           => config["outer_taper_length"])));

    #>----------------Lower Cone / Taper on bottom-----------------------------<
    # 45-degree taper at bottom of ORTEC-type crystal
    lower_cone = dicttype("translate" => dicttype("z"      => 0.4*config["taper_length"],
                        "cone"        => dicttype(
                        "name"        => "Lower Cone",
                        "r"           => dicttype("top"    => dicttype("from" => config["xtal_radius"],
                                                               "to"   => ceil(config["xtal_radius"] + 1.0)),
                                              "bottom" => dicttype("from" => config["xtal_radius"] -
                                                                            1.2 * config["taper_length"],
                                                               "to"   => ceil(config["xtal_radius"] + 1.0))),
                        "phi"         => dicttype("from"   => 0.0, "to" => 360.0),
                        "h"           => 1.2*config["taper_length"])));
                        

    #>----------------Ditch around the bottom contact--------------------------<
    # If not existing, everything will be set to 0
    ditch    = dicttype("translate"   => dicttype("z"    => config["ditch_depth"]/2,
                        "tube"        => dicttype(
                        "name"        => "Ditch",
                        "r"           => dicttype("from" => config["wrap_around_radius"] -
                                                            config["ditch_thickness"],
                                              "to"   => config["wrap_around_radius"]),
                        "phi"         => dicttype("from" => 0.0, "to" => 360.0),
                        "h"           => config["ditch_depth"])));

    #>----------------Borehole on bottom as contact----------------------------<
    # If not existing, everything will be set to 0
    borehole = if config["bulletize_PC"] == 1
        if config["pc_length"] < config["pc_radius"]
            dicttype( "union"   => [ 
            dicttype( "tube"    => dicttype(
                  "r"       => config["pc_radius"] - config["pc_length"],
                  "h"       => 1.2*config["pc_length"],
                  "origin"  => dicttype("z" => 0.4*config["pc_length"]))),
            dicttype( "torus"   => dicttype(
                  "r_torus" => config["pc_radius"] - config["pc_length"],
                  "r_tube"  => config["pc_length"]))])
        elseif config["pc_length"] > config["pc_radius"]
            dicttype( "union"  => [
            dicttype( "tube"   => dicttype(
                  "r"      => config["pc_radius"],
                  "h"      => 1.2*config["pc_length"] - config["pc_radius"],
                  "origin" => dicttype("z" => (0.8*config["pc_length"] - config["pc_radius"])/2))),
            dicttype( "sphere" => dicttype(
                  "r"      => config["pc_radius"],
                  "origin" => dicttype("z" => config["pc_length"] - config["pc_radius"])))])
        else
            dicttype( "sphere" => dicttype(
                       "r" => config["pc_radius"]))
        end
    else # not bulletized
            dicttype( "tube"   => dicttype(
                  "name"   => "Borehole",
                  "r"      => dicttype("from" => 0.0, "to" => config["pc_radius"]),
                  "h"      => config["pc_length"],
                  "origin" => dicttype("z" => config["pc_length"]/2)));
    end

    #>----------------Hole on top of the detector------------------------------<
    # If not existing, everything will be set to 0
    hole = if config["inner_taper_width"] != 0 && config["inner_taper_length"] != 0
        dicttype("translate"   => dicttype("z"      => config["xtal_length"] -
                               config["inner_taper_length"]/2,
             "cone"        => dicttype(
             "name"        => "Inner Cone",
             "r"           => dicttype("top"    => dicttype("from" => 0,
                                                   "to"   => config["hole_radius"] +
                                                                config["inner_taper_width"]),
                                  "bottom" => dicttype("from" => 0,
                                                   "to"   => config["hole_radius"])),
             "phi"         => dicttype("from"   => 0.0, "to" => 360.0),
             "h"           => config["inner_taper_length"])));
                            
    else
        dicttype("translate"   => dicttype("z"    => config["xtal_length"] -
                                               0.8*config["hole_length"]/2,
             "tube"        => dicttype(
             "name"        => "Top Hole",
             "r"           => dicttype("from" => 0.0, "to" => config["hole_radius"]),
             "phi"         => dicttype("from" => 0.0, "to" => 360.0),
             "h"           => 1.2*config["hole_length"])));
    end

    #>----------------Subtract volumes-----------------------------------------<
    geometry = main_vol;
    if config["outer_taper_length"] != 0 || config["outer_taper_width"] != 0
        geometry = dicttype("difference" => [geometry, upper_cone])
    end
    if config["inner_taper_length"] != 0 || config["inner_taper_width"] != 0 || config["hole_radius"] != 0 || config["hole_length"] != 0
        geometry = dicttype("difference" => [geometry, hole])
    end
    if config["taper_length"] != 0
        geometry = dicttype("difference" => [geometry, lower_cone])
    end
    if config["ditch_thickness"] != 0
        geometry = dicttype("difference" => [geometry, ditch])
    end
    if config["bulletize_PC"] == 1 || (config["pc_radius"] != 0 && config["pc_length"] >= 0.1) # below this, it is the thickness of the layer
        geometry = dicttype("difference" => [geometry, borehole])
    end


    #>----------------Volume END-----------------------------------------------<
    #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
    #>----------------Contacts START-------------------------------------------<
    

    #>----------------Area/volume of point contact-----------------------------<
    #>----------------Contact 1------------------------------------------------<
    # Depending on whether the contact is a borehole or not
    
    geometry_1 = if config["bulletize_PC"] == 1.0
        if config["pc_length"] < config["pc_radius"]
            dicttype( "union"  => [
            dicttype( "tube"   => dicttype(
                  "r"          => config["pc_radius"] - config["pc_length"],
                  "h"          => 0,
                  "origin"     => dicttype("z" => config["pc_length"]))),
            dicttype( "torus"  => dicttype(
                  "r_torus"    => config["pc_radius"] - config["pc_length"],
                  "r_tube"     => dicttype("from" => config["pc_length"], 
                                           "to"   => config["pc_length"]),
                  "theta"      => dicttype("from" => 0.0, "to" => 90.0)))])
        elseif config["pc_length"] > config["pc_radius"]
            dicttype( "union"  => [
            dicttype( "tube"   => dicttype(
                  "r"          => dicttype("from" => config["pc_radius"], "to" => config["pc_radius"]),
                  "h"          => config["pc_length"] - config["pc_radius"],
                  "origin"     => dicttype("z" => (config["pc_length"] - config["pc_radius"])/2))),
            dicttype( "difference" => [
            dicttype( "sphere"     => dicttype(
                  "r"          => dicttype("from" => config["pc_radius"], 
                                           "to"   => config["pc_radius"]),
                  "origin"     => dicttype("z" => config["pc_length"] - config["pc_radius"]))),
            dicttype( "tube"       => dicttype(
                  "r"          => 1.2*config["pc_radius"],
                  "h"          => 1.2*config["pc_radius"],
                  "origin"     => dicttype("z" => config["pc_length"] - 1.6*config["pc_radius"])))])])
        else
            dicttype( "difference" => [
            dicttype( "sphere"     => dicttype(
                  "r"          => dicttype("from" => config["pc_radius"], 
                                           "to"   => config["pc_radius"]))),
            dicttype( "tube"       => dicttype(
                  "r"          => 1.2*config["pc_radius"],
                  "h"          => 1.2*config["pc_radius"],
                  "origin"     => dicttype("z" => -0.6 * config["pc_radius"])))]) 
        end
    elseif config["pc_radius"] != 0 && config["pc_length"] >= 5 # below this, it is the thickness of the layer
            dicttype( "union"      => [
            # borehole_wall
            dicttype( "tube"       => dicttype(
                  "r"          => dicttype("from" => config["pc_radius"],
                                           "to"   => config["pc_radius"]),
                  "h"          => config["pc_length"],
                  "origin"     => dicttype("z" => config["pc_length"]/2))),
            # borehole_top
            dicttype( "tube"       => dicttype(
                  "r"          => dicttype("from" => 0.0,
                                           "to"   => config["pc_radius"]),
                  "h"          => 0.0,
                  "origin"     => dicttype("z" => config["pc_length"]))),
            # borehole_around
            dicttype( "tube"       => dicttype(
                  "r"          => dicttype("from" => config["pc_radius"],
                                           "to"   => config["wrap_around_radius"] - config["ditch_thickness"]),
                  "h"          => 0.0))])
    elseif config["pc_radius"] != 0 && config["pc_length"] < 5
            dicttype( "tube"       => dicttype(
                  "r"          => config["pc_radius"],
                  "h"          => config["pc_length"],
                  "origin"     => dicttype("z" => config["pc_length"]/2)))
    end

    #>----------------Define contact 1-----------------------------------------<

    contact_1 = dicttype("material"    => "HPGe",
                         "id"          => 1,
                         "potential"   => 0.0,
                         "geometry"    => geometry_1);

    #>----------------Contact 2------------------------------------------------<
    geometry_2 = dicttype[]

    if config["wrap_around_radius"] != 0.0
        bottom = dicttype("tube" => dicttype(
                      "r"    => dicttype("from" => config["wrap_around_radius"],
                                         "to" => config["xtal_radius"] - config["taper_length"]),
                      "phi"  => dicttype("from" => 0.0,
                                         "to" => 360.0),
                      "h"    => 0.0))
        push!(geometry_2, bottom)
    end

    wall   = dicttype("translate" => dicttype("z" => config["xtal_length"]/2 + config["taper_length"]/2 - config["outer_taper_length"]/2,
                  "tube" => dicttype(
                  "r"    => dicttype("from" => config["xtal_radius"],
                                     "to"   => config["xtal_radius"]),
                  "phi"  => dicttype("from" => 0.0,
                                     "to" => 360.0),
                  "h"    => config["xtal_length"] - config["taper_length"] - config["outer_taper_length"])))
    push!(geometry_2, wall)

    top    = dicttype("translate" => dicttype("z" => config["xtal_length"],
                  "tube" => dicttype(
                  "r"    => dicttype("from" => config["hole_radius"] + config["inner_taper_width"],
                                     "to" => config["xtal_radius"] - config["outer_taper_width"]),
                  "phi"  => dicttype("from" => 0.0,
                                     "to" => 360.0),
                  "h"    => 0.0)))
    push!(geometry_2, top)
    if config["hole_radius"] != 0 && config["hole_length"] != 0 && config["inner_taper_width"] == 0 && config["inner_taper_length"] == 0
        hole   = dicttype("translate" => dicttype("z" => config["xtal_length"] - config["hole_length"]/2,
                      "tube" => dicttype(
                      "r"    => dicttype("from" => config["hole_radius"],
                                         "to" => config["hole_radius"]),
                      "phi"  => dicttype("from" => 0.0,
                                         "to" => 360.0),
                      "h"    => config["hole_length"])))
        inner_hole_bottom =  dicttype("translate" => dicttype("z" => config["xtal_length"] - config["hole_length"],
                      "tube" => dicttype(
                      "r"    => dicttype("from" => 0.0,
                                         "to" => config["hole_radius"]),
                      "phi"  => dicttype("from" => 0.0,
                                         "to" => 360.0),
                      "h"    => 0.0)))
        push!(geometry_2, hole)
        push!(geometry_2, inner_hole_bottom)
    end
    if config["hole_radius"] != 0 && config["hole_length"] != 0 && config["inner_taper_width"] != 0 && config["inner_taper_length"] != 0
        inner_taper = dicttype("translate" => dicttype("z" => config["xtal_length"] - config["inner_taper_length"]/2,
                      "cone" => dicttype(
                      "r"    => dicttype("top" => dicttype("from" => config["hole_radius"] + config["inner_taper_width"],
                                                           "to"   => config["hole_radius"] + config["inner_taper_width"]),
                                      "bottom" => dicttype("from" => config["hole_radius"],
                                                           "to"   => config["hole_radius"])),
                      "phi"  => dicttype("from" => 0.0, "to" => 360.0),
                      "h"    => config["inner_taper_length"])))
        inner_hole_bottom = dicttype("translate" => dicttype("z" => config["xtal_length"] - config["hole_length"],
                      "tube" => dicttype(
                      "r"    => dicttype("from" => 0.0, "to" => config["hole_radius"]),
                      "phi"  => dicttype("from" => 0.0, "to" => 360.0),
                      "h"    => 0.0)))
        push!(geometry_2, inner_taper)
        push!(geometry_2, inner_hole_bottom)
    end
    if config["outer_taper_width"] != 0 && config["outer_taper_length"] != 0
        outer_taper = dicttype("translate" => dicttype("z" => config["xtal_length"] - config["outer_taper_length"]/2,
                      "cone" => dicttype(
                      "r"    => dicttype("top" => dicttype("from" => config["xtal_radius"] - config["outer_taper_width"],
                                                           "to"   => config["xtal_radius"] - config["outer_taper_width"]),
                                      "bottom" => dicttype("from" => config["xtal_radius"],
                                                           "to"   => config["xtal_radius"])),
                      "phi"  => dicttype("from"   => 0.0, "to" => 360.0),
                      "h"    => config["outer_taper_length"])))
        push!(geometry_2, outer_taper)
    end
    if config["taper_length"] != 0
        taper = dicttype("translate" => dicttype("z" => config["taper_length"]/2,
                      "cone" => dicttype(
                      "r"    => dicttype("top" => dicttype("from" => config["xtal_radius"],
                                                           "to"   => config["xtal_radius"]),
                                      "bottom" => dicttype("from" => config["xtal_radius"] - config["taper_length"],
                                                           "to"   => config["xtal_radius"] - config["taper_length"])),
                      "phi"   => dicttype("from"   => 0.0, "to" => 360.0),
                      "h"     => config["taper_length"])))
        push!(geometry_2, taper)
    end


    #>----------------Define contact 2-----------------------------------------<

    contact_2 = dicttype("material"    => "HPGe",
                     "id"          => 2,
                     "potential"   => config["xtal_HV"],
                     "geometry"    => dicttype("union" => geometry_2));

    #>----------------Contacts END---------------------------------------------<
    #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
    #>----------------Final configuration--------------------------------------<


    impurity_z0            = config["impurity_z0"]*1e10
    impurity_gradient      = config["impurity_gradient"]*1e10
    if units["length"]     == "mm"
        impurity_z0        /= 1e3
        impurity_gradient  /= 1e4
    elseif units["length"] == "m"
        impurity_z0        *= 1e6
        impurity_gradient  *= 1e8
    end
    
    objects = dicttype("semiconductor" => dicttype(
                    "material"    => "HPGe",
                    "temperature" => config["xtal_temp"],
                    "impurity_density" => dicttype(
                                "name" => "cylindrical",
                                "offset" => impurity_z0,
                                "gradient" => dicttype("z" => impurity_gradient)),
                    "geometry"    => geometry),
                    "contacts"    => [ contact_1, contact_2 ]);

    conf = dicttype("name"    => config["name"],
                "units"       => units,
                "grid"        => grid,
                "medium"      => "vacuum",
                "detectors"   => [objects]);

    return conf
end
