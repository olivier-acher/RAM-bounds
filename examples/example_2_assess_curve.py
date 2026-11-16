"""Example 2 - Reflection curve alone -> assessment.

An absorber is known by its attenuation curve (no phase), its thickness and its magnetic mass.
How far is it from the limits? Here: stack 3 of the paper (Fig. 10, Table II), 3.89 mm thick,
with 0.56 kg/m^2 of CoNbZr films. Its efficiency with the saturation factor is 0.88 with respect
to its own absorption between 0.3 and 18 GHz, and 0.79 with respect to the specification
(-10 dB over 1-9 GHz), as in Table II of the paper.

To use your own curve: a text file with two columns, frequency (GHz) and attenuation (dB)
(or reflection in dB, negative), then change the file name, D, M and the material below.

Run from the RAM-bounds folder:   python examples/example_2_assess_curve.py
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from rambounds import read_curve, assess_curve, Specification
from rambounds.plots import plot_attenuation, plot_assessment

f, gamma_dB = read_curve(os.path.join(HERE, "data", "stack3_attenuation.csv"))
D = 3.89            # thickness, mm
M = 0.5625          # magnetic mass, kg/m^2
material = "CoNbZr film"

# 1) the absorber against its own reflection curve, credited between 0.3 and 18 GHz
own = assess_curve(f, gamma_dB, D=D, M=M, material=material, band=(0.3, 18.0))
print("=== Own reflection curve ===")
print(own.summary())

# 2) the same absorber against the specification it was designed for: -10 dB over 1-9 GHz
spec = Specification(1.0, 9.0, 10.0)
template = assess_curve([1.0, 9.0], [10.0, 10.0], D=D, M=M, material=material, band=(1.0, 9.0))
print("\n=== Specification: %s ===" % spec.describe())
print(template.summary())

fig, (a, b) = plt.subplots(1, 2, figsize=(9.5, 3.8), gridspec_kw=dict(width_ratios=[1.2, 1]))
plot_attenuation(f, gamma_dB, ax=a, spec=spec, label="stack 3")
a.set_xlim(0.3, 30)
plot_assessment(own, ax=b, label="stack 3")
fig.tight_layout()
out = os.path.join(HERE, "example_2_assess_curve.png")
fig.savefig(out, dpi=150)
print("\nFigure saved:", out)
