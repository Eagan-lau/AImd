def hydroxyl(*args, **kwargs): return {"name": "hydroxyl"}
def h_abstractable_carbon(*args, **kwargs): return {"name": "h_abstractable_carbon"}
def protic_nucleophile(*args, **kwargs): return {"name": "protic_nucleophile"}
def acyl_carbonyl(*args, **kwargs): return {"name": "acyl_carbonyl"}

LIGAND_GROUPS = {
    "hydroxyl": hydroxyl,
    "h_abstractable_carbon": h_abstractable_carbon,
    "protic_nucleophile": protic_nucleophile,
    "acyl_carbonyl": acyl_carbonyl,
}
