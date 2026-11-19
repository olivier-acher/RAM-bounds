"""Reflection curve alone -> assessment (Secs. 3.1 and 3.5 of the paper).

Input: an attenuation curve Gamma_dB(f), measured, computed or read on a data sheet, without
its phase; the thickness D of the absorber; its magnetic mass M (or its amount of magnetic
material J_mu) and its magnetic material.

Output:
- the Rozanov thickness D_R and the thin magnetic absorber mass M_A of the curve itself (Eq. 50),
  and the position (x, m) = (D/D_R, M/M_A) of the absorber in the reduced plane;
- its efficiency eta (Eq. 51), and the same with the saturation factor;
- the minimum magnetic mass that any absorber of thickness D needs to produce this curve,
  with the simplified form (20), the local form (21) and the general form (18).

Only the attenuation within the chosen band is credited (0 elsewhere), which keeps every
bound valid (Sec. 3.5).
"""
from dataclasses import dataclass, field

import numpy as np

from .constants import C, db_to_np
from .crossover import (trapezoid_weights, absorption_area, effective_attenuation,
                        saturation_factor, envelope)
from .materials import get_material


@dataclass
class Assessment:
    band: tuple
    D: float
    J_mu: float
    material: object
    D_R: float
    J_A: float
    M_A: float
    x: float
    m: float
    eta: float
    Fc_eta: float
    eta_xi: float
    M_min_simplified: float
    M_min_local: float
    M_min_rigorous: float
    Fs: float
    curves: dict = field(repr=False, default_factory=dict)

    @property
    def M(self):
        return self.material.mass(self.J_mu)

    def summary(self):
        mat = self.material
        return "\n".join([
            "Band credited: %g-%g GHz" % self.band,
            "Absorber: D = %.4g mm, M = %.4g kg/m2 of %s (J_mu = %.4g mm.GHz2)" % (self.D, self.M, mat.name, self.J_mu),
            "Own reflection curve: D_R = %.4g mm, M_A = %.4g kg/m2 (J_A = %.4g mm.GHz2)" % (self.D_R, self.M_A, self.J_A),
            "Reduced position: x = D/D_R = %.3g, m = M/M_A = %.3g" % (self.x, self.m),
            "Efficiency eta (Eq. 51) = %.3g, reached at Fc = %.3g GHz; with the saturation factor: %.3g" % (self.eta, self.Fc_eta, self.eta_xi),
            "Minimum magnetic mass at this thickness: %.4g kg/m2 (simplified form 22), %.4g kg/m2 (local form 23), "
            "%.4g kg/m2 (general form 18, rigorous), split frequency %.3g GHz" % (self.M_min_simplified, self.M_min_local, self.M_min_rigorous, self.Fs),
            "Ratio M / M_min (rigorous) = %s" % ("%.3g" % (self.M / self.M_min_rigorous) if self.M_min_rigorous > 0 else "infinite: no magnetic material is required")])


def assess_curve(f, gamma_dB, D, M=None, material="carbonyl iron powder", J_mu=None, band=None,
                 n_freq=3000, n_fc=120, nphi=96):
    """Assess an absorber from its attenuation curve.

    f, gamma_dB : frequencies (GHz) and attenuation (dB, positive = absorption)
    D           : thickness of the absorber, mm
    M           : magnetic mass, kg/m^2 (or give J_mu, mm.GHz^2, instead)
    material    : a Material, or a name from rambounds.MATERIALS
    band        : (fa, fb) in GHz, the band credited; default: the range of f
    """
    mat = get_material(material)
    if J_mu is None:
        if M is None:
            raise ValueError("give the magnetic mass M or the amount of magnetic material J_mu")
        J_mu = mat.amount(M)
    f = np.asarray(f, dtype=float)
    gamma_dB = np.asarray(gamma_dB, dtype=float)
    order = np.argsort(f)
    f, gamma_dB = f[order], gamma_dB[order]
    fa, fb = band if band is not None else (f[0], f[-1])
    fg = np.geomspace(fa, fb, n_freq)
    g = db_to_np(np.clip(np.interp(fg, f, gamma_dB), 0, None))      # nepers, negative values ignored
    w = trapezoid_weights(fg)

    D_R = C / (2 * np.pi**2) * np.sum(g / fg**2 * w)
    J_A = C / (2 * np.pi**2) * np.sum(g * w)
    Fc = np.geomspace(fa / 30, fb * 30, n_fc)
    A = absorption_area(fg, w, g, Fc)
    Ae = absorption_area(fg, w, effective_attenuation(g), Fc)
    xi = saturation_factor(fg, w, g, Fc, nphi=nphi)
    L22 = C / (2 * np.pi**2) * A
    L23 = C / (2 * np.pi**2) * Ae
    L18 = (1 - xi) * L22
    budget = Fc**2 * D + J_mu
    k = np.argmax(L22 / budget)
    J22, _ = envelope(L22, Fc, D)
    J23, _ = envelope(L23, Fc, D)
    J18, Fs = envelope(L18, Fc, D)
    return Assessment(band=(fa, fb), D=D, J_mu=J_mu, material=mat, D_R=D_R, J_A=J_A, M_A=mat.mass(J_A),
                      x=D / D_R, m=J_mu / J_A, eta=float(np.max(L22 / budget)), Fc_eta=float(Fc[k]),
                      eta_xi=float(np.max(L18 / budget)),
                      M_min_simplified=mat.mass(J22), M_min_local=mat.mass(J23), M_min_rigorous=mat.mass(J18), Fs=float(Fs),
                      curves=dict(f=fg, gamma_dB=gamma_dB, f_input=f, Fc=Fc, xi=xi, L22=L22, L23=L23, L18=L18))


def limit_curve(assessment, form="rigorous", n=200):
    """Reduced plane: lowest m = M/M_A as a function of x = D/D_R for the same reflection curve
    (the solid curves of Fig. 9b of the paper). form = 'simplified' (20), 'local' (21) or 'rigorous' (18)."""
    L = assessment.curves[{"simplified": "L22", "local": "L23", "rigorous": "L18"}[form]]
    x = np.linspace(0, 1.0, n)
    J, _ = envelope(L, assessment.curves["Fc"], x * assessment.D_R)
    return x, J / assessment.J_A
