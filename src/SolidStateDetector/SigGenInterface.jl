"""
    readsiggen(file_path::String[, T::Type=Float64])
Read the '*.config' file in 'file_path' for SigGen and returns a dictionary of all parameters.
Non-existing parameteres are set to 0.
...
# Arguments
- `file_path::String`: file path for the SigGen config file.
- `T::Type=Float64`: type of the parameters in the output dictionary.
...
"""
function readsiggen(file_path::String; T::Type=Float64)
    #detector_name = split(basename(file_path), ".config")[1]
    #In order to use the plotting recipe:
    detector_name = "Public Inverted Coax"
    
    config = Dict("name"       => detector_name,
        "verbosity_level"      => 0,
        "xtal_length"          => 0,
        "xtal_radius"          => 0,
        "top_bullet_radius"    => 0, # not implemented in later conversions
        "bottom_bullet_radius" => 0, # not implemented in later conversions
        "pc_length"            => 0,
        "pc_radius"            => 0,
        "bulletize_PC"         => 0, # not implemented in later conversions
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
        config_file   = Dict()
        for line in temp
            if line  != "" && string(line[1]) != "#"
                line  = strip(split(line, "#")[1])
                name  = ""
                value = ""
                if length(split(line, "\t")) > 1
                    name  = split(line, "\t")[1]
                    value = strip(split(line, "\t")[2])
                else
                    name  = split(line, " ")[1]
                    value = strip(split(line, " ")[end])
                end

                config_file[name] = value

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
    siggentodict(config::Dict[, units::Dict, detector_name::String])
Converts the dictionary containing the parameters from a SigGen config file
to a SSD config dictionary. This dictionary can be saved as a JSON file using
the JSON package and 'JSON.print(file, config, 4)'.
The 'detector_name' is set to "Public Inverted Coax" by default to inherit
the colour scheme.
...
# Arguments
- `config::Dict`: dictionary containing SigGen parameters (output of readsiggen()).
- `units::Dict`: units used in SigGen file (set to 'mm', 'deg', 'V' and 'K').
- `detector_name::String`: name of the detector.
...
"""
function siggentodict(config::Dict;
        units::Dict = Dict("length" => "mm","angle" => "deg","potential" => "V","temperature" => "K"))

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
    axes         = Dict("r"           => Dict("to"         => r_limit,
                                              "boundaries" => "inf"),
                        "phi"         => Dict("from"       => 0,
                                              "to"         => 0,
                                              "boundaries" => "periodic"),
                        "z"           => Dict("from"       => -10,
                                              "to"         => z_limit,
                                              "boundaries" => Dict("left" => "inf", "right" => "inf"))
                        )

    grid         = Dict("coordinates" => "cylindrical",
                        "axes"        => axes);


    #>----------------Settings END---------------------------------------------<
    #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
    #>----------------Volume START---------------------------------------------<


    #>----------------Main Volume -> Initial Cylinder--------------------------<
    # r and L are given by the input file. Cylindrical symmetry is set to default
    main_vol     = Dict("translate"   => Dict("z"    => config["xtal_length"]/2,
                        "tube"        => Dict(
                        "name"        => "Initial Cylinder",
                        "r"           => Dict("from" => 0.0, "to" => config["xtal_radius"]),
                        "phi"         => Dict("from" => 0.0, "to" => 360.0),
                        "h"           => config["xtal_length"])));

    #>----------------Upper Cone / Taper on top--------------------------------<
    # If not existing, everything will be set to 0
    upper_cone   = Dict("translate"   => Dict("z"      => config["xtal_length"] -
                                            config["outer_taper_length"]/2,
                        "cone"        => Dict(
                        "name"        => "Upper Cone",
                        "r"           => Dict("bottom" => Dict("from" => config["xtal_radius"],
                                                               "to"   => ceil(config["xtal_radius"] + 1.0)),
                                              "top"    => Dict("from" => config["xtal_radius"] -
                                                                            config["outer_taper_width"],
                                                               "to"   => ceil(config["xtal_radius"] + 1.0))),
                        "phi"         => Dict("from"   => 0.0, "to" => 360.0),
                        "h"           => config["outer_taper_length"])));

    #>----------------Lower Cone / Taper on bottom-----------------------------<
    # 45-degree taper at bottom of ORTEC-type crystal
    lower_cone   = Dict("translate"   => Dict("z"      => config["taper_length"]/2,
                        "cone"        => Dict(
                        "name"        => "Lower Cone",
                        "r"           => Dict("top"    => Dict("from" => config["xtal_radius"],
                                                               "to"   => ceil(config["xtal_radius"] + 1.0)),
                                              "bottom" => Dict("from" => config["xtal_radius"] -
                                                                            config["taper_length"],
                                                               "to"   => ceil(config["xtal_radius"] + 1.0))),
                        "phi"         => Dict("from"   => 0.0, "to" => 360.0),
                        "h"           => config["taper_length"])));
                        

    #>----------------Ditch around the bottom contact--------------------------<
    # If not existing, everything will be set to 0
    ditch        = Dict("translate"   => Dict("z"    => config["ditch_depth"]/2,
                        "tube"        => Dict(
                        "name"        => "Ditch",
                        "r"           => Dict("from" => config["wrap_around_radius"] -
                                                            config["ditch_thickness"],
                                              "to"   => config["wrap_around_radius"]),
                        "phi"         => Dict("from" => 0.0, "to" => 360.0),
                        "h"           => config["ditch_depth"])));

    #>----------------Borehole on bottom as contact----------------------------<
    # If not existing, everything will be set to 0
    borehole     = Dict("translate"   => Dict("z"    => config["pc_length"]/2,
                        "tube"        => Dict(
                        "name"        => "Borehole",
                        "r"           => Dict("from" => 0.0, "to" => config["pc_radius"]),
                        "phi"         => Dict("from" => 0.0, "to" => 360.0),
                        "h"           => config["pc_length"])));

    #>----------------Hole on top of the detector------------------------------<
    # If not existing, everything will be set to 0
    hole = Dict()
    if config["inner_taper_width"] != 0 && config["inner_taper_length"] != 0
        hole         = Dict("translate"   => Dict("z"      => config["xtal_length"] -
                                               config["inner_taper_length"]/2,
                            "cone"        => Dict(
                            "name"        => "Inner Cone",
                            "r"           => Dict("top"    => Dict("from" => 0,
                                                                   "to"   => config["hole_radius"] +
                                                                                config["inner_taper_width"]),
                                                  "bottom" => Dict("from" => 0,
                                                                   "to"   => config["hole_radius"])),
                            "phi"         => Dict("from"   => 0.0, "to" => 360.0),
                            "h"           => config["inner_taper_length"])));
                            
    else
        hole         = Dict("translate"   => Dict("z"    => config["xtal_length"] -
                                               config["hole_length"]/2,
                            "tube"        => Dict(
                            "name"        => "Top Hole",
                            "r"           => Dict("from" => 0.0, "to" => config["hole_radius"]),
                            "phi"         => Dict("from" => 0.0, "to" => 360.0),
                            "h"           => config["hole_length"])));
    end

    #>----------------Subtract volumes-----------------------------------------<
    geometry = main_vol;
    if config["outer_taper_length"] != 0 || config["outer_taper_width"] != 0
        geometry = Dict("difference" => [geometry, upper_cone])
    end
    if config["inner_taper_length"] != 0 || config["inner_taper_width"] != 0 || config["hole_radius"] != 0 || config["hole_length"] != 0
        geometry = Dict("difference" => [geometry, hole])
    end
    if config["taper_length"] != 0
        geometry = Dict("difference" => [geometry, lower_cone])
    end
    if config["ditch_thickness"] != 0
        geometry = Dict("difference" => [geometry, ditch])
    end
    if config["pc_radius"] != 0 && config["pc_length"] >= 0.1 # below this, it is the thickness of the layer
        geometry = Dict("difference" => [geometry, borehole])
    end


    #>----------------Volume END-----------------------------------------------<
    #<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
    #>----------------Contacts START-------------------------------------------<


    geometry_1   = Any[] # gathers all parts of contact_1

    #>----------------Area/volume of point contact-----------------------------<
    #>----------------Contact 1------------------------------------------------<
    # Depending on whether the contact is a borehole or not
    if config["pc_radius"] != 0 && config["pc_length"] >= 5 # below this, it is the thickness of the layer
        borehole_wall = Dict("translate" => Dict("z" => config["pc_length"]/2,
                                 "tube" => Dict(
             "r"         => Dict("from" => config["pc_radius"],
                                   "to" => config["pc_radius"]),
             "phi"       => Dict("from" => 0.0,
                                   "to" => 360.0),
             "h"         => config["pc_length"])))
        borehole_top = Dict("translate" => Dict("z" => config["pc_length"],
                                 "tube" => Dict(
             "r"         => Dict("from" => 0.0,
                                   "to" => config["pc_radius"]),
             "phi"       => Dict("from" => 0.0,
                                   "to" => 360.0),
             "h"         => 0.0)))
        borehole_around = Dict(  "tube" => Dict(
             "r"         => Dict("from" => config["pc_radius"],
                                   "to" => config["wrap_around_radius"] - config["ditch_thickness"]),
             "phi"       => Dict("from" => 0.0,
                                   "to" => 360.0),
             "h"         => 0.0))
        push!(geometry_1, borehole_wall)
        push!(geometry_1, borehole_top)
        push!(geometry_1, borehole_around)

    elseif config["pc_radius"] != 0 && config["pc_length"] < 5
        pc_vol = Dict(      "translate" => Dict("z" => config["pc_length"]/2,
                                 "tube" => Dict(
             "r"         => Dict("from" => 0.0,
                                   "to" => config["pc_radius"]),
             "phi"       => Dict("from" => 0.0,
                                   "to" => 360.0),
             "h"         => config["pc_length"])))
        push!(geometry_1, pc_vol)
    end

    #>----------------Define contact 1-----------------------------------------<

    contact_1 = Dict("material"    => "HPGe",
                     "id"     => 1,
                     "potential"   => 0.0,
                     "geometry"    => Dict("union" => geometry_1));

    #>----------------Contact 2------------------------------------------------<
    geometry_2 = Any[]

    if config["wrap_around_radius"] != 0.0
        bottom = Dict("tube" => Dict(
                      "r"    => Dict("from" => config["wrap_around_radius"],
                                       "to" => config["xtal_radius"] - config["taper_length"]),
                      "phi"  => Dict("from" => 0.0,
                                       "to" => 360.0),
                      "h"    => 0.0))
        push!(geometry_2, bottom)
    end

    wall   = Dict("translate" => Dict("z" => config["xtal_length"]/2 + config["taper_length"]/2 - config["outer_taper_length"]/2,
                  "tube" => Dict(
                  "r"    => Dict("from" => config["xtal_radius"],
                                 "to"   => config["xtal_radius"]),
                  "phi"  => Dict("from" => 0.0,
                                   "to" => 360.0),
                  "h"    => config["xtal_length"] - config["taper_length"] - config["outer_taper_length"])))
    push!(geometry_2, wall)

    top    = Dict("translate" => Dict("z" => config["xtal_length"],
                  "tube" => Dict(
                  "r"    => Dict("from" => config["hole_radius"] + config["inner_taper_width"],
                                   "to" => config["xtal_radius"] - config["outer_taper_width"]),
                  "phi"  => Dict("from" => 0.0,
                                   "to" => 360.0),
                  "h"    => 0.0)))
    push!(geometry_2, top)
    if config["hole_radius"] != 0 && config["hole_length"] != 0 && config["inner_taper_width"] == 0 && config["inner_taper_length"] == 0
        hole   = Dict("translate" => Dict("z" => config["xtal_length"] - config["hole_length"]/2,
                      "tube" => Dict(
                      "r"    => Dict("from" => config["hole_radius"],
                                       "to" => config["hole_radius"]),
                      "phi"  => Dict("from" => 0.0,
                                       "to" => 360.0),
                      "h"    => config["hole_length"])))
        inner_hole_bottom =  Dict("translate" => Dict("z" => config["xtal_length"] - config["hole_length"],
                      "tube" => Dict(
                      "r"    => Dict("from" => 0.0,
                                       "to" => config["hole_radius"]),
                      "phi"  => Dict("from" => 0.0,
                                       "to" => 360.0),
                      "h"    => 0.0)))
        push!(geometry_2, hole)
        push!(geometry_2, inner_hole_bottom)
    end
    if config["hole_radius"] != 0 && config["hole_length"] != 0 && config["inner_taper_width"] != 0 && config["inner_taper_length"] != 0
        inner_taper = Dict("translate" => Dict("z" => config["xtal_length"] - config["inner_taper_length"]/2,
                      "cone" => Dict(
                      "r"    => Dict("top"    => Dict("from" => config["hole_radius"] + config["inner_taper_width"],
                                                      "to"   => config["hole_radius"] + config["inner_taper_width"]),
                                     "bottom" => Dict("from" => config["hole_radius"],
                                                      "to"   => config["hole_radius"])),
                      "phi"  => Dict("from" => 0.0, "to" => 360.0),
                      "h"    => config["inner_taper_length"])))
        inner_hole_bottom = Dict("translate" => Dict("z" => config["xtal_length"] - config["hole_length"],
                      "tube" => Dict(
                      "r"    => Dict("from" => 0.0, "to" => config["hole_radius"]),
                      "phi"  => Dict("from" => 0.0, "to" => 360.0),
                      "h"    => 0.0)))
        push!(geometry_2, inner_taper)
        push!(geometry_2, inner_hole_bottom)
    end
    if config["outer_taper_width"] != 0 && config["outer_taper_length"] != 0
        outer_taper = Dict("translate" => Dict("z" => config["xtal_length"] - config["outer_taper_length"]/2,
                      "cone" => Dict(
                      "r"    => Dict("top"    => Dict("from" => config["xtal_radius"] - config["outer_taper_width"],
                                                      "to"   => config["xtal_radius"] - config["outer_taper_width"]),
                                     "bottom" => Dict("from" => config["xtal_radius"],
                                                      "to"   => config["xtal_radius"])),
                      "phi"  => Dict("from"   => 0.0, "to" => 360.0),
                      "h"    => config["outer_taper_length"])))
        push!(geometry_2, outer_taper)
    end
    if config["taper_length"] != 0
        taper = Dict("translate" => Dict("z" => config["taper_length"]/2,
                      "cone" => Dict(
                      "r"    => Dict("top"    => Dict("from" => config["xtal_radius"],
                                                      "to"   => config["xtal_radius"]),
                                     "bottom" => Dict("from" => config["xtal_radius"] - config["taper_length"],
                                                      "to"   => config["xtal_radius"] - config["taper_length"])),
                      "phi"       => Dict("from"   => 0.0, "to" => 360.0),
                      "h"         => config["taper_length"])))
        push!(geometry_2, taper)
    end


    #>----------------Define contact 2-----------------------------------------<

    contact_2 = Dict("material"    => "HPGe",
                     "id"     => 2,
                     "potential"   => config["xtal_HV"],
                     "geometry"    => Dict("union" => geometry_2));

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
    
    objects = Dict("semiconductor" => Dict(
                    "material"    => "HPGe",
                    "temperature" => config["xtal_temp"],
                    "impurity_density" => Dict(
                                "name"     => "cylindrical",
                                "r"        => Dict("init"     => 0.0,
                                                   "gradient" => 0.0),
                                "z"        => Dict("init"     => impurity_z0,
                                                   "gradient" => impurity_gradient)),
                    "geometry"    => geometry),
                    "contacts" => [ contact_1, contact_2 ]);

    conf = Dict("name"    => config["name"],
                "units"   => units,
                "grid"    => grid,
                "medium"  => "vacuum",
                "detectors" => [objects]);

    return conf
end
