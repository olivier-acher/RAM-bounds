"""Specification -> frontier (Sec. 2.3 of the paper).

A specification is a minimum attenuation, in dB, over a frequency band. The simplest one is
flat: Gamma_dB >= level over [F1, F2]. Templates made of several flat steps are also accepted.

For a specification, the module gives:
- the Rozanov thickness D_R (thinnest dielectric absorber, Eq. 27) and the thin magnetic
  absorber amount J_A and mass M_A (lightest magnetic absorber of negligible thickness, Eq. 28);
- the approximate law (32), and the exact frontier (37): the minimum amount of magnetic
  material J_ex(D), or magnetic mass M_min(D), of any absorber of thickness D that meets
  the specification, with the split frequency F_s at this thickness;
- the balanced point of the approximate law (Eq. 33).
"""
import numpy as np

from .constants import C, db_to_np, np_to_db
from .crossover import crossover_inductance_ideal, envelope
from .materials import get_material


class Specification:
    """Minimum attenuation over a band.

    Specification(F1, F2, level_dB)            flat: Gamma_dB >= level_dB over [F1, F2] (GHz)
    Specification.from_steps([(fa, fb, level_dB), ...])   piecewise-flat template
    """

    def __init__(self, F1, F2, level_dB):
        if not 0 < F1 < F2:
            raise ValueError("the band must satisfy 0 < F1 < F2")
        self._set([(F1, F2, level_dB)])

    @classmethod
    def from_steps(cls, steps):
        spec = cls.__new__(cls)
        spec._set(sorted(steps))
        return spec

    def _set(self, steps):
        edges = [steps[0][0]]
        levels = []
        for fa, fb, lev in steps:
            if fa < edges[-1] - 1e-12 or fb <= fa:
                raise ValueError("steps must be ordered, non-overlapping intervals (fa < fb)")
            if fa > edges[-1] + 1e-12:          # gap between steps: no attenuation required
                levels.append(0.0)
                edges.append(fa)
            levels.append(float(db_to_np(lev)))
            edges.append(fb)
        self.steps = list(steps)
        self.edges = np.array(edges, dtype=float)
        self.levels = np.array(levels, dtype=float)       # nepers
        self.F1, self.F2 = self.edges[0], self.edges[-1]
        self.flat = len(self.steps) == 1
        self._cache = {}

    # ---- the two pure solutions -------------------------------------------------
    @property
    def bandwidth_ratio(self):
        return self.F2 / self.F1

    @property
    def D_R(self):
        """Rozanov thickness, mm (Eq. 27)."""
        return C / (2 * np.pi**2) * np.sum(self.levels * (1 / self.edges[:-1] - 1 / self.edges[1:]))

    @property
    def J_A(self):
        """Amount of magnetic material of the thin magnetic absorber, mm.GHz^2 (Eq. 28)."""
        return C / (2 * np.pi**2) * np.sum(self.levels * (self.edges[1:] - self.edges[:-1]))

    def M_A(self, material):
        """Magnetic mass of the thin magnetic absorber, kg/m^2 (Eq. 28)."""
        return get_material(material).mass(self.J_A)

    # ---- crossover frequencies ----------------------------------------------------
    def crossover_grid(self, n=600):
        """Crossover frequencies used for the envelopes (GHz)."""
        Fc = np.geomspace(self.F1 / 50, self.F2 * 50, n)
        near_edges = (self.edges[None, :] * (1 + np.array([-1e-3, 0.0, 1e-3]))[:, None]).ravel()
        return np.unique(np.concatenate([Fc, near_edges]))

    def _lines(self):
        """Crossover inductances of the template (no saturation) and of the ideal absorber."""
        if "lines" not in self._cache:
            Fc = self.crossover_grid()
            # absorption area of the template, exact for piecewise-flat steps
            A = np.zeros_like(Fc)
            for g, fa, fb in zip(self.levels, self.edges[:-1], self.edges[1:]):
                lo = np.clip(Fc, fa, fb)
                A += g * ((lo - fa) + Fc**2 * (1 / lo - 1 / fb))
            L_lin = C / (2 * np.pi**2) * A
            # ideal absorber W; Z_in = 1 (no reflection) also meets any specification,
            # so the frontier never lies above the hyperbola (53)
            L_W = np.minimum(crossover_inductance_ideal(Fc, self.edges, self.levels), 2 * C * Fc / np.pi**2)
            self._cache["lines"] = (Fc, L_lin, L_W)
        return self._cache["lines"]

    # ---- frontiers -------------------------------------------------------------------
    def approximate_law(self, D):
        """Approximate frontier: minimum J (mm.GHz^2) at thickness D (mm), neglecting the
        saturation factor (Eqs. 30-32). For a flat specification it is the law (32),
        J = J_A (1 - x) / (1 + (r - 1) x), with x = D / D_R."""
        D = np.asarray(D, dtype=float)
        if self.flat:
            x = D / self.D_R
            r = self.bandwidth_ratio
            return self.J_A * np.clip((1 - x) / (1 + (r - 1) * x), 0, None)
        Fc, L_lin, _ = self._lines()
        return envelope(L_lin, Fc, D)[0]

    def exact_frontier(self, D):
        """Exact frontier (37): minimum amount of magnetic material J_ex(D) (mm.GHz^2) of an
        absorber of thickness D (mm) that meets the specification, and the split frequency
        F_s (GHz) at this thickness. Returns (J_ex, F_s)."""
        Fc, _, L_W = self._lines()
        return envelope(L_W, Fc, D)

    def minimum_mass(self, D, material):
        """Minimum magnetic mass M_min(D), kg/m^2, on the exact frontier."""
        return get_material(material).mass(self.exact_frontier(D)[0])

    def balanced_point(self):
        """Balanced point of the approximate law for a flat specification (Eq. 33):
        x_b = m_b = 1 / (1 + sqrt(r)), reached for F_s = sqrt(F1 F2). Returns (D_b, J_b)."""
        if not self.flat:
            raise ValueError("the balanced point is defined here for flat specifications")
        b = 1 / (1 + np.sqrt(self.bandwidth_ratio))
        return b * self.D_R, b * self.J_A

    # ---- report --------------------------------------------------------------------------
    def describe(self):
        if self.flat:
            return "Gamma_dB >= %g dB over %g-%g GHz" % (np_to_db(self.levels[0]), self.F1, self.F2)
        return "template: " + ", ".join("%g dB over %g-%g GHz" % (lev, fa, fb) for fa, fb, lev in self.steps)

    def summary(self, material, thicknesses=()):
        """Text summary for a magnetic material (a Material or a name from MATERIALS)."""
        mat = get_material(material)
        lines = ["Specification: " + self.describe(),
                 "Magnetic material: %s (k_A = %.3g, 4 pi Ms = %g kG, rho = %g g/cm3; K_mu,dB = %.3g dB.GHz per kg/m2)"
                 % (mat.name, mat.k_A, mat.Ms_kG, mat.rho, mat.K_mu_dB),
                 "Rozanov thickness            D_R = %.3g mm" % self.D_R,
                 "Thin magnetic absorber       J_A = %.4g mm.GHz2, M_A = %.3g kg/m2" % (self.J_A, mat.mass(self.J_A))]
        if self.flat:
            Db, Jb = self.balanced_point()
            Jx, Fs = self.exact_frontier(Db)
            lines.append("Balanced point (law 35)      D = %.3g mm, M = %.3g kg/m2  (exact frontier at this D: M_min = %.3g kg/m2, F_s = %.3g GHz)"
                         % (Db, mat.mass(Jb), mat.mass(Jx), Fs))
        for D in thicknesses:
            Jx, Fs = self.exact_frontier(D)
            lines.append("At D = %-6.3g mm               M_min = %.3g kg/m2 (exact frontier), %.3g kg/m2 (law), F_s = %.3g GHz"
                         % (D, mat.mass(Jx), mat.mass(self.approximate_law(D)), Fs))
        return "\n".join(lines)
