def distance(*args, **kwargs): return {"name": "distance"}
def angle(*args, **kwargs): return {"name": "angle"}
def angle_with_axis(*args, **kwargs): return {"name": "angle_with_axis"}

FEATURE_TYPES = {
    "distance": distance,
    "angle": angle,
    "angle_with_axis": angle_with_axis,
}
