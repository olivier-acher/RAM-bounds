"""RAM-bounds: bounds on the thickness and magnetic mass of magneto-dielectric radar absorbers.

Companion code of the paper "Mass and thickness of magneto-dielectric radar absorbers:
fundamental constraints and practical guidelines" (O. Acher).

Two uses:
- Specification -> frontier:  Specification(F1, F2, level_dB), then .summary(material),
  .exact_frontier(D), .approximate_law(D), plots.plot_specification(spec)
- Reflection curve -> assessment:  assess_curve(f, gamma_dB, D, M, material), then .summary(),
  plots.plot_assessment(result)
"""
from .constants import C, K_DB, GAMMA_BAR, db_to_np, np_to_db
from .materials import Material, MATERIALS, get_material
from .specification import Specification
from .assessment import assess_curve, limit_curve, Assessment
from .io import read_curve

__version__ = "0.1.0"
