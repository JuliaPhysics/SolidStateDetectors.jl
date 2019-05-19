# Constructive Solid Geometry (CSG)

## Boolean operators

### Union

```json
"geometry": {
    "type": "union",
    "parts": [
        {
            "name":"Seg1 bottom",
            "type": "tube",
            "r": {
                "from": 13.5,
                "to": 39.5
            },
            "phi": {
                "from": 0.3582,
                "to": 59.6419
            },
            "h": 0
        },
        {
            "name": "Seg1 side",
            "type": "tube",
            "r": {
                "from": 39.5,
                "to": 39.5
            },
            "phi": {
                "from": 0.3582,
                "to": 59.6419
            },
            "h": 40
        },
    }
}
```

### Difference

```json
"geometry": {
    "type": "difference",
    "parts": [
        {
            "name": "Initial Cylinder",
            "type": "tube",
            "r": {
                "from": 0.0,
                "to":  35.0
            },
            "phi": {
                "from": 0.0,
                "to": 360.0
            },
            "z": {
                "from": 0,
                "to":  40
            }
        },
        {
            "name": "Borehole",
            "type": "tube",
            "r": {
                "from": 0.0,
                "to":   5.0
            },
            "phi": {
                "from": 0.0,
                "to": 360.0
            },
            "z": {
                "from": 0,
                "to":  40
            }
        }
    ]
}
```

### Intersection

ToDo...