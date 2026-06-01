def collect_all(*args, **kwargs): return {"name": "collect_all"}
def cartesian_product(*args, **kwargs): return {"name": "cartesian_product"}
def pick_best_per_parent(*args, **kwargs): return {"name": "pick_best_per_parent"}

COLLECTION_OPS = {
    "collect_all": collect_all,
    "cartesian_product": cartesian_product,
    "pick_best_per_parent": pick_best_per_parent,
}
