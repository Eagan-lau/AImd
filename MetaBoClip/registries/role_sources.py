def selector(*args, **kwargs): return {"name": "selector"}
def selector_any(*args, **kwargs): return {"name": "selector_any"}
def selector_within(*args, **kwargs): return {"name": "selector_within"}
def template_anchor(*args, **kwargs): return {"name": "template_anchor"}
def nearby_from_role(*args, **kwargs): return {"name": "nearby_from_role"}

ROLE_SOURCES = {
    "selector": selector,
    "selector_any": selector_any,
    "selector_within": selector_within,
    "template_anchor": template_anchor,
    "nearby_from_role": nearby_from_role,
}
