"""Checks against numbers given in the paper.

Run from the RAM-bounds folder with   python tests/test_paper_values.py
(or with pytest, if installed:        pytest tests)
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from rambounds import MATERIALS, Specification, read_curve, assess_curve
from rambounds.crossover import crossover_inductance_ideal


def close(a, b, rel):
    return abs(a - b) <= rel * abs(b)


def test_table_I():
    """Table I: I_mu,dB and K_mu,dB of the typical materials."""
    expected = {"MnZn ferrite": (37, 7.6), "NiZn ferrite": (24, 4.8), "carbonyl iron powder": (690, 88),
                "permalloy powder": (150, 17), "iron film": (2100, 260), "permalloy film": (450, 52),
                "CoNbZr film": (450, 51)}
    for name, (I_dB, K_dB) in expected.items():
        m = MATERIALS[name]
        assert close(m.I_mu_dB, I_dB, 0.03), (name, m.I_mu_dB)
        assert close(m.K_mu_dB, K_dB, 0.03), (name, m.K_mu_dB)


def test_flat_specification():
    """Sec. 2.2-2.3: -10 dB over 2-18 GHz with carbonyl iron powder."""
    spec = Specification(2.0, 18.0, 10.0)
    iron = MATERIALS["carbonyl iron powder"]
    assert close(spec.D_R, 7.8, 0.01)                       # Rozanov thickness, mm
    assert close(spec.M_A(iron), 1.8, 0.02)                 # thin magnetic absorber mass, kg/m^2
    D_b, J_b = spec.balanced_point()
    assert close(D_b, 1.94, 0.01) and close(iron.mass(J_b), 0.46, 0.01)


def test_ideal_absorber():
    """Crossover inductance of the ideal absorber W (-10 dB over 2-18 GHz), by the half-circle formula."""
    spec = Specification(2.0, 18.0, 10.0)
    assert close(crossover_inductance_ideal(6.0, spec.edges, spec.levels), 133.11, 1e-3)
    assert close(crossover_inductance_ideal(10.0, spec.edges, spec.levels), 209.24, 1e-3)


def test_stack_3():
    """Stack 3 (Fig. 10, Tables II and III): exact frontier at 3.89 mm, and efficiencies with the
    saturation factor, 0.88 for its own absorption over 0.3-18 GHz and 0.79 for the specification."""
    spec = Specification(1.0, 9.0, 10.0)
    J_ex, F_s = spec.exact_frontier(3.8857)
    assert close(J_ex, 31.6, 0.01) and close(F_s, 2.9, 0.02)
    f, g = read_curve(os.path.join(HERE, "..", "examples", "data", "stack3_attenuation.csv"))
    own = assess_curve(f, g, D=3.8857, J_mu=50.11, material="CoNbZr film", band=(0.3, 18.0))
    assert abs(own.eta_xi - 0.88) < 0.01, own.eta_xi
    tpl = assess_curve([1.0, 9.0], [10.0, 10.0], D=3.8857, J_mu=50.11, material="CoNbZr film", band=(1.0, 9.0))
    assert abs(tpl.eta_xi - 0.79) < 0.01, tpl.eta_xi


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok  ", name)
    print("all checks passed")
