# CSG Primitives

CSG = Constructive Solid Geometries

## Primitive: Box

```json
{
    "type": "box",
    "x": {
        "from": 0.0,
        "to":  10.0
    },
    "y":{
        "from": 0.0,
        "to":  10.0
    },
    "z": {
        "from": 0.0,
        "to":  10.0
    }
}
```

## Primitve: Tube

```json
{
    "type": "tube",
    "r": {
        "from": 3.0,
        "to":   7.0
    },
    "phi": {
        "from": 0.0,
        "to": 360.0
    },
    "z": {
        "from": 0.0,
        "to":  10.0
    }
}
```

## Primitive: Cone

```json
{
    "type": "cone",
    "r": {
        "bottom": {
            "from": 10.0,
            "to": 20.0
        },
        "top": {
            "from": 5.0,
            "to": 25.0
        }
    },
    "phi": {
        "from": 0.0,
        "to": 360.0
    },
    "z": {
        "from": 0.0,
        "to":  10.0
    }
}
```

## Primitive: Sphere

```json
{
    "type": "sphere",
    "r": 30,
    "translate": {
        "x": 0,
        "y": 0,
        "z": 0
    }
}
```