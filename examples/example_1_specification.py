"""Example 1 - Specification -> frontier.

What is the minimal tradeoff between thickness and magnetic mass for an absorber that keeps
Gamma_dB >= 10 dB over 2-18 GHz, with carbonyl iron powder as magnetic material?
(Sec. 2.3 of the paper: D_R = 7.8 mm, M_A = 1.8 kg/m^2, balanced point D = 1.94 mm, M = 0.46 kg/m^2.)

Run from the RAM-bounds folder:   python examples/example_1_specification.py
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from rambounds import Specification, MATERIALS
from rambounds.plots import plot_attenuation, plot_specification

spec = Specification(2.0, 18.0, 10.0)            # F1, F2 in GHz, level in dB
material = MATERIALS["carbonyl iron powder"]     # or Material("my material", k_A, 4piMs_kG, rho)

print(spec.summary(material, thicknesses=[1.0, 2.0, 4.0]))

# the frontier as a table: thickness (mm) -> minimum magnetic mass (kg/m^2)
D = np.array([0.5, 1.0, 2.0, 3.0, 4.0, 6.0])
J, Fs = spec.exact_frontier(D)
print("\n  D (mm)   M_min (kg/m2)   F_s (GHz)")
for d, j, fs in zip(D, J, Fs):
    print("  %6.2f   %13.3f   %9.2f" % (d, material.mass(j), fs))

fig, (a, b) = plt.subplots(1, 2, figsize=(9.5, 3.8), gridspec_kw=dict(width_ratios=[1.2, 1]))
f = np.geomspace(1, 30, 400)
plot_attenuation(f, np.where((f >= 2) & (f <= 18), 10.0, np.nan), ax=a, spec=spec, color="none")
plot_specification(spec, material, ax=b)
fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "example_1_specification.png")
fig.savefig(out, dpi=150)
print("\nFigure saved:", out)
