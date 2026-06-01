def energy_linear(*args, **kwargs): return {"name": "energy_linear"}
def distance_piecewise(*args, **kwargs): return {"name": "distance_piecewise"}
def gaussian(*args, **kwargs): return {"name": "gaussian"}

SCORING_TRANSFORMS = {
    "energy_linear": energy_linear,
    "distance_piecewise": distance_piecewise,
    "gaussian": gaussian,
}
