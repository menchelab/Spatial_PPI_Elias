import json
from pathlib import Path

from constraints import (
    SphereConstraint,
    ShellConstraint,
    CylinderConstraint,
    EllipsoidConstraint,
    EllipsoidShellConstraint,
)

CONSTRAINTS = {
    "SphereConstraint": SphereConstraint,
    "ShellConstraint": ShellConstraint,
    "CylinderConstraint": CylinderConstraint,
    "EllipsoidConstraint": EllipsoidConstraint,
    "EllipsoidShellConstraint": EllipsoidShellConstraint,
}


class LocationSpec:
    def __init__(self, name, min_degree=0, center=(0, 0, 0), repulsion_strength=1.0,
                 constraint_type="SphereConstraint", constraint_params=None):
        self.name = str(name)
        self.min_degree = int(min_degree)
        self.center = tuple(float(x) for x in center)
        self.repulsion_strength = float(repulsion_strength)
        self.constraint_type = str(constraint_type)
        self.constraint_params = dict(constraint_params or {})


class AppConfig:
    def __init__(self, default_params=None, locations=None):
        self.default_params = dict(default_params or {})
        self.locations = list(locations or [])


def load_config(path):
    path = Path(path)
    data = json.loads(path.read_text(encoding="utf-8"))

    default_params = data.get("default", {}) or {}

    locations = []
    for loc in data.get("locations", []) or []:
        constraint = loc.get("constraint", {}) or {}
        ctype = constraint.get("type", "")

        if ctype not in CONSTRAINTS:
            raise ValueError(f"Unknown constraint type '{ctype}'. Allowed: {sorted(CONSTRAINTS)}")

        cparams = constraint.get("params", {}) or {}

        locations.append(
            LocationSpec(
                name=loc["name"],
                min_degree=loc.get("min_degree", 0),
                center=loc.get("center", (0, 0, 0)),
                repulsion_strength=loc.get("repulsion_strength", 1.0),
                constraint_type=ctype,
                constraint_params=cparams,
            )
        )

    if not locations:
        raise ValueError("Config contains no 'locations'.")

    return AppConfig(default_params=default_params, locations=locations)


def build_constraint(constraint_type, params):
    cls = CONSTRAINTS[constraint_type]
    return cls(**(params or {}))
