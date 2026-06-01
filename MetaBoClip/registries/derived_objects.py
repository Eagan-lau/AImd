def vector(*args, **kwargs): return {"name": "vector"}
def plane(*args, **kwargs): return {"name": "plane"}
def normal_from_plane(*args, **kwargs): return {"name": "normal_from_plane"}
def oriented_normal_from_plane(*args, **kwargs): return {"name": "oriented_normal_from_plane"}

DERIVED_OBJECTS = {
    "vector": vector,
    "plane": plane,
    "normal_from_plane": normal_from_plane,
    "oriented_normal_from_plane": oriented_normal_from_plane,
}
