def weighted_sum(*args, **kwargs): return {"name": "weighted_sum"}
def coupled_mean_positive(*args, **kwargs): return {"name": "coupled_mean_positive"}

SCORING_EXPRESSIONS = {
    "weighted_sum": weighted_sum,
    "coupled_mean_positive": coupled_mean_positive,
}
