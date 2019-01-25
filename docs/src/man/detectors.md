# Detectors
Currently, three classes of detectors are supported: Coaxial, Inverted Coax, and BEGe type detectors. 

All detector properties are specified in a one-for-all .json file. 

## Example 1) Inverted Coax

Example minimum config file for an Inverted Coax detector (IVC) plus explanations. 
**Remember**, comments are not allowed in JSON files and have to be deleted if you want to use it.
```json
{
    "name":"ExampleInvertedCoax", // Arbitrary name of the detector
    "class":"InvertedCoax",  // either "Coax", "BEGe", "InvertedCoax"
    "type":"p", // either "p", "ptype", "p-type", "n", "ntype" or "n-type"
    "cyclic":0, // The periodicity of the detector in degree. 
                // `0` means complete symmetric in θ -> 2D simulation. 
                // The condition "360 / `cyclic` = n; n = 1, 2, 3, ..." must be fulfilled.
    "mirror_symmetry_θ": "true", // set to true, if a mirror symmetry exists within the periodicity specified via `cyclic`
    "materials":{
        "detector":"HPGe",  // Material of the detector 
        "environment":"Vacuum" // Material of the environment
    },
    "geometry":
    {
        "unit": "mm", // unit of the values which the user specifies in this config file
        "crystal":
        {
            "length": 80.0, // crystal goes from z=0mm to z=80mm
            "radius": 35.0  // crystal goes from r=0mm to r=35mm
        },
        "point_contact":
        {
            "endplate":"bot",
            "depth":3e-4,
            "radius":9.5
        },
        "groove":
        {
            "endplate":"bot",
            "rInner":13.5,
            "width":0.0,
            "depth":0.0
        },
        "borehole": // IVC: Borehole starts at the top of the crystal 
        {
          "length":55.0,
          "radius":5.0
        },
        "taper": 
        {
            "outer":
            {
                "rInner":24.42,
                "length":60
            },
            "inner":
            {
                "rOuter":0,
                "length":0
            }
        }
    },
    "charge_carrier_density":
    {
        "unit": "cm-3",
        "top": 0.94, // times 10^10
        "bot": 1.10  // just the absolute value. The sign is determined through the bulk material
    },
    "segmentation":
    {
        "n_contacts_total": 3, // n-contacts + p-contacts
        "core": // core = n-contact of the detector
        {
            "type": "Tubs", // "Box" volume in cylindrical coordinates
            "rStart": 0.0,
            "rStop": 3.0,
            "phiStart": 0,
            "phiStop": 360,
            "zStart": 0.0,
            "zStop": 0.0,
            "potential": 0.0 // Bias voltage of the core
        },        
        "n_individual_segments": 2, // now the p-contacts
        "S1":
        {
            "type": "Tubs",
            "rStart": 15.0,
            "rStop": 35.0,
            "phiStart": 0,
            "phiStop": 360,
            "zStart": 0,
            "zStop": 0,
            "potential": 3500.0, // Bias voltage of this segment
            "boundaryWidth":
            {
                "radial":0.0,
                "vertical": 0.0,
                "horizontal": 0.0
            },
            "repetitive":false,
            "repetitions":
            {
                "radial":0,
                "vertical": 0,
                "horizontal": 5
            }
        },
        "S2":
        {
            "type": "Tubs",
            "rStart": 35.0,
            "rStop": 35.0,
            "phiStart": 0,
            "phiStop": 360,
            "zStart": 0.0,
            "zStop": 20.0,
            "potential": 3500.0, // Bias voltage of this segment
            "boundaryWidth":
            {
                "radial":0.0,
                "vertical": 0.0,
                "horizontal": 0.0
            },
            "repetitive":false,
            "repetitions":
            {
                "radial": 0,
                "vertical": 0,
                "horizontal": 5
            }
        },
        "custom_grouping":false, // Multiple geometrical segments can build up one contact:
                                 // Specify which geometrical segments belong to which contact
                                 // To specify geometrical segments `3`, `4`, `5`, `6` use "3-4-5-6" (not only "3-6")
        "Chn1":"1",    // contact "Chn1" is build through the geometrical segment `1`. (This is the core in this case)
        "Chn2":"2-3"   // contact "Chn2" is build through the geometrical segments `2` & `3`. 
    }
}
```
