# RAM-bounds

Companion code of the paper

> O. Acher, *Mass and thickness of magneto-dielectric radar absorbers: fundamental constraints and practical guidelines*, preprint (2026), Zenodo, [doi:10.5281/zenodo.24345127](https://doi.org/10.5281/zenodo.24345127). Submitted to Physical Review Applied.

Equation, figure and section numbers below refer to this preprint.

The paper establishes a bound on the performance of any planar, passive, linear magneto-dielectric absorber on a metal plane at normal incidence: the absorption area of its reflection curve cannot exceed what its thickness and its amount of magnetic material can pay for. This code applies the bound in two ways.

| Use | You give | You get |
|---|---|---|
| **Specification → frontier** (Sec. 2.3) | a minimum attenuation over a band, e.g. −10 dB over 2–18 GHz, and a magnetic material | the Rozanov thickness `D_R`, the thin magnetic absorber mass `M_A`, the approximate law (32), the exact frontier (37): minimum magnetic mass at every thickness, the balanced point |
| **Reflection curve → assessment** (Secs. 3.1, 3.5) | an attenuation curve Γ_dB(f) (measured, computed or read on a data sheet, no phase needed), the thickness `D`, the magnetic mass `M` and the magnetic material | the position of the absorber in the reduced plane, its efficiency η (Eq. 51), the minimum magnetic mass allowed at its thickness by the simplified (20), local (21) and general (18) forms of the bound |

## Installation

You need Python 3.9 or later with `numpy` and `matplotlib` (both are included in Anaconda).

1. Download the code: green **Code** button at the top of this page → **Download ZIP**, then unzip it.
2. If needed, install the two packages: `pip install -r requirements.txt`
3. Check that everything works, from the `RAM-bounds` folder: `python tests/test_paper_values.py`

## Quick start

Run the two examples from the `RAM-bounds` folder:

```
python examples/example_1_specification.py
python examples/example_2_assess_curve.py
```

Each prints its results and saves a figure next to the script.

### Specification → frontier

```python
from rambounds import Specification, MATERIALS

spec = Specification(2.0, 18.0, 10.0)        # F1 (GHz), F2 (GHz), minimum attenuation (dB)
iron = MATERIALS["carbonyl iron powder"]
print(spec.summary(iron, thicknesses=[1.0, 2.0]))

J_ex, F_s = spec.exact_frontier(2.0)          # thickness in mm
print(iron.mass(J_ex), "kg/m2 at least, split frequency", F_s, "GHz")
```

Templates made of several flat steps are accepted: `Specification.from_steps([(1, 2, 10), (2, 18, 20)])` (GHz, GHz, dB).

### Reflection curve → assessment

```python
from rambounds import read_curve, assess_curve

f, gamma_dB = read_curve("my_absorber.csv")   # two columns: frequency (GHz), attenuation (dB)
result = assess_curve(f, gamma_dB, D=3.0, M=4.5, material="carbonyl iron powder", band=(2, 18))
print(result.summary())
```

The curve file is a text file with two columns, frequency in GHz and attenuation in dB, separated by commas, semicolons, tabs or spaces; header lines are skipped. A reflection coefficient in dB (negative values) is accepted as well. Only the attenuation within `band` is credited, which keeps all bounds valid.

### Figures

```python
from rambounds.plots import plot_attenuation, plot_specification, plot_assessment
```

They follow the conventions of the paper: attenuation increasing downward, frequency axis on top.

## Magnetic materials

`MATERIALS` contains the typical values of Table I of the paper:

| Name | k_A | 4πM_s (kG) | ρ (g/cm³) |
|---|---|---|---|
| MnZn ferrite | 1/3 | 5.0 | 4.9 |
| NiZn ferrite | 1/3 | 4.0 | 5.0 |
| carbonyl iron powder | 1/3 | 21.5 | 7.9 |
| permalloy powder | 1/3 | 10 | 8.7 |
| iron film | 1 | 21.5 | 7.9 |
| permalloy film | 1 | 10 | 8.7 |
| CoNbZr film | 1 | 10 | 8.8 |

Any other material: `Material("my ferrite", k_A=1/3, Ms_kG=4.5, rho=5.1)`. The geometric factor `k_A` is 1 for films or laminations magnetized in their plane and excited along their hard axis, 1/2 when the two transverse permeabilities are equal, 1/3 for isotropic materials such as sintered ferrites and random powders. The magnetic mass `M` is the mass of the magnetic material only, per unit area, not the total mass of the absorber.

## What the numbers mean

- **D_R** (Rozanov thickness): thickness of the thinnest dielectric absorber with this reflection curve.
- **M_A** (thin magnetic absorber mass): magnetic mass of the lightest magnetic absorber of negligible thickness with this reflection curve.
- **x = D/D_R, m = M/M_A**: reduced coordinates; the absorber must lie above the frontier.
- **η** (efficiency, Eq. 51): the largest fraction of the budget of the bound that the absorber uses; 1 on the frontier.
- **Forms of the bound** for a curve without phase: simplified (20), for quick estimates; local (21), with the effective attenuation 2 tanh(Γ/2), rigorous and simple; general (18), with the saturation factor ξ computed from Γ(f) alone, rigorous and the best of the three.
- **F_s** (split frequency): below F_s the absorption is paid with magnetic material, above F_s with thickness.

## Scope and limitations

The bound holds for passive, linear, time-invariant planar absorbers on a metal plane, at normal incidence. The exact frontier uses the ideal absorber of the specification, whose optimality is proven at the center of the band and near its ends (Sec. 2.4.5 of the paper).

## How to cite

Please cite the paper above and this code. The code is archived on Zenodo at each release, with its own DOI.

## License

MIT, see [LICENSE](LICENSE).
