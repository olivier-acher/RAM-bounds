"""Crossover quantities of the paper.

For a crossover frequency Fc, the crossover window is min(1, Fc^2/f^2): it counts
in frequency below Fc and in wavelength above Fc.

- absorption area      A(Fc)  = int Gamma(f) min(1, Fc^2/f^2) df          (Eq. 17)
- saturation factor    xi(Fc), computed from Gamma(f) alone               (Eq. B14)
- crossover inductance Lambda[Z](Fc) = (c/pi^2) int Re Z(f) min(1, Fc^2/f^2) df   (Eqs. 20, 45),
  computed here as an average of Re Z over a half-circle of complex frequencies (Eq. B4).

Attenuations Gamma are in nepers in this module.
"""
import numpy as np

from .constants import C

# ---------------------------------------------------------------------------
# integrals on a frequency grid
# ---------------------------------------------------------------------------


def trapezoid_weights(f):
    """Weights w such that sum(w * g) approximates the integral of g over the grid f."""
    f = np.asarray(f, dtype=float)
    w = np.empty_like(f)
    w[1:-1] = (f[2:] - f[:-2]) / 2
    w[0] = (f[1] - f[0]) / 2
    w[-1] = (f[-1] - f[-2]) / 2
    return w


def absorption_area(f, w, gamma, Fc):
    """Absorption area A(Fc), in Np.GHz (Eq. 17). Fc may be a scalar or an array."""
    Fc = np.atleast_1d(np.asarray(Fc, dtype=float))
    win = np.minimum(1.0, (Fc[:, None] / f[None, :]) ** 2)
    A = win @ (gamma * w)
    return A if A.size > 1 else A[0]


def effective_attenuation(gamma):
    """Effective attenuation Gamma_e = 2 tanh(Gamma/2) (Eq. 23)."""
    return 2.0 * np.tanh(np.asarray(gamma, dtype=float) / 2.0)


def saturation_factor(f, w, gamma, Fc, nphi=96):
    """Saturation factor xi(Fc) of a reflection curve Gamma(f) (Appendix B 5, Eq. B14).

    y(phi) is half the mean attenuation seen through a Lorentzian window of half-width
    Fc cos(phi) centered at Fc sin(phi); xi is the relative loss when y is replaced by tanh(y).
    No phase is needed: the bound (18) holds for any passive absorber with this attenuation.
    """
    Fc = np.atleast_1d(np.asarray(Fc, dtype=float))
    phi = (np.arange(nphi) + 0.5) / nphi * np.pi / 2
    gw = gamma * w
    out = np.empty(Fc.size)
    for i, F in enumerate(Fc):
        s = F * np.cos(phi)[:, None]
        a = F * np.sin(phi)[:, None]
        K = s / (s**2 + (f[None, :] - a) ** 2) + s / (s**2 + (f[None, :] + a) ** 2)
        y = (K @ gw) / (2 * np.pi)
        lin = np.sum(y * np.cos(phi))
        out[i] = 0.0 if lin <= 0 else 1.0 - np.sum(np.tanh(y) * np.cos(phi)) / lin
    return out if out.size > 1 else out[0]


# ---------------------------------------------------------------------------
# crossover inductance of an analytic impedance, by the half-circle formula
# ---------------------------------------------------------------------------

_xg, _wg = np.polynomial.legendre.leggauss(200)
_EDGES = np.concatenate([np.linspace(0, np.pi / 2 - 0.05, 20),
                         np.pi / 2 - np.logspace(np.log10(0.05), -9, 60)[1:], [np.pi / 2]])
_PHI = np.concatenate([(a + b) / 2 + (b - a) / 2 * _xg for a, b in zip(_EDGES[:-1], _EDGES[1:])])
_WPHI = np.concatenate([(b - a) / 2 * _wg for a, b in zip(_EDGES[:-1], _EDGES[1:])])


def crossover_inductance(Z_of_s, Fc):
    """Lambda[Z](Fc), in mm.GHz^2, for an impedance given as a function of the complex
    frequency s (rad/ns, right half-plane), by Eq. (B4):
        Lambda = (2 c Fc / pi^2) int_0^{pi/2} Re Z(2 pi Fc e^{j phi}) cos(phi) dphi.
    The integration panels are refined near phi = pi/2, where the half-circle meets the axis."""
    s = 2 * np.pi * Fc * np.exp(1j * _PHI)
    return 2 * C * Fc / np.pi**2 * np.sum(_WPHI * np.real(Z_of_s(s)) * np.cos(_PHI))


def beta_steps(s, edges, levels):
    """Analytic continuation of half the attenuation of a minimum-phase reflection coefficient
    whose attenuation is the piecewise-constant template `levels` (nepers) on the frequency
    intervals defined by `edges` (GHz), and 0 elsewhere:
        beta(s) = sum_k (Gamma_k / pi) [atan(2 pi f_{k+1} / s) - atan(2 pi f_k / s)].
    On the real-frequency axis, Re beta = Gamma/2, and R = -exp(-2 beta)."""
    beta = np.zeros_like(s, dtype=complex)
    for g, fa, fb in zip(levels, edges[:-1], edges[1:]):
        beta += g / np.pi * (np.arctan(2 * np.pi * fb / s) - np.arctan(2 * np.pi * fa / s))
    return beta


def crossover_inductance_ideal(Fc, edges, levels):
    """Crossover inductance of the ideal absorber W of a template (Eq. 37 for a flat template):
    its reflection coefficient has the modulus of the template within the band, 1 outside,
    and minimum phase, so that W = tanh(beta). Fc may be a scalar or an array."""
    Fc = np.atleast_1d(np.asarray(Fc, dtype=float))
    L = np.array([crossover_inductance(lambda s: np.tanh(beta_steps(s, edges, levels)), F) for F in Fc])
    return L if L.size > 1 else L[0]


def envelope(L, Fc, D):
    """Lowest amount of magnetic material compatible with the lines F_c^2 D + J >= L(F_c):
    J(D) = max over F_c of [L(F_c) - F_c^2 D], clipped at 0, and the F_c of the maximum.
    D may be a scalar or an array (mm)."""
    D = np.atleast_1d(np.asarray(D, dtype=float))
    G = L[None, :] - Fc[None, :] ** 2 * D[:, None]
    k = np.argmax(G, axis=1)
    J = np.maximum(G[np.arange(D.size), k], 0.0)
    Fs = Fc[k]
    return (J, Fs) if J.size > 1 else (J[0], Fs[0])
