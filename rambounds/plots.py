"""Figures, with the conventions of the paper:
attenuation increasing downward, frequency axis on top (at 0 dB), and no titles."""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from .constants import np_to_db
from .materials import get_material
from .assessment import limit_curve

INK, GRAY, GRID = "#0b0b0b", "#8a8a8a", "#d9d8d4"


def plot_attenuation(f, gamma_dB, ax=None, spec=None, label=None, logf=True, color="#b4541f"):
    """Attenuation Gamma_dB(f), increasing downward, frequency axis on top.
    If a Specification is given, its template is drawn, with the forbidden region in gray."""
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 3.2))
    ax.plot(f, gamma_dB, color=color, lw=1.4, label=label)
    ymax = max(np.nanmax(gamma_dB), 1.0)
    if spec is not None:
        for g, fa, fb in zip(spec.levels, spec.edges[:-1], spec.edges[1:]):
            lev = float(np_to_db(g))
            ax.fill_between([fa, fb], 0, lev, color=GRAY, alpha=0.25, lw=0)
            ax.plot([fa, fb], [lev, lev], color=INK, lw=1.0)
            ymax = max(ymax, lev)
    ax.set_ylim(1.15 * ymax, 0)
    if logf:
        ax.set_xscale("log")
        ax.xaxis.set_major_locator(mticker.LogLocator(subs=(1.0, 2.0, 5.0)))
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, p: "%g" % v))
        ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    ax.set_xlabel("frequency $f$ (GHz)")
    ax.set_ylabel("attenuation $\\Gamma_\\mathrm{dB}$ (dB)")
    ax.grid(color=GRID, lw=0.5)
    if label:
        ax.legend(frameon=False)
    return ax


def plot_specification(spec, material=None, ax=None, n=300, mark=()):
    """Plane (thickness, magnetic mass) for a Specification: approximate law (32), exact
    frontier (37), balanced point. In reduced units (x = D/D_R, m = M/M_A) if no material is
    given, in mm and kg/m^2 otherwise. `mark` = list of (D, M, label) points to add."""
    if ax is None:
        _, ax = plt.subplots(figsize=(4.6, 4))
    D = np.linspace(0, 1.05 * spec.D_R, n)
    Jlaw = spec.approximate_law(D)
    Jex, _ = spec.exact_frontier(D)
    if material is None:
        sx, sy, xl, yl = 1 / spec.D_R, 1 / spec.J_A, "reduced thickness $x=D/D_R$", "reduced magnetic mass $m=M/M_A$"
    else:
        mat = get_material(material)
        sx, sy, xl, yl = 1.0, 1 / mat.K_mu, "thickness $D$ (mm)", "magnetic mass $M$ (kg/m$^2$) of %s" % mat.name
    ax.fill_between(D * sx, 0, Jex * sy, color=GRAY, alpha=0.25, lw=0)
    ax.plot(D * sx, Jlaw * sy, "--", color=GRAY, lw=1.3, label="approximate law (32)")
    ax.plot(D * sx, Jex * sy, "-", color=INK, lw=1.5, label="exact frontier (37)")
    if spec.flat:
        Db, Jb = spec.balanced_point()
        ax.plot([Db * sx], [Jb * sy], "o", color="white", mec=INK, mew=1.2, ms=6, label="balanced point (law 35)")
    for Dm, Mm, lab in mark:
        Jm = Mm * get_material(material).K_mu if material is not None else Mm * spec.J_A
        ax.plot([Dm * sx], [Jm * sy], "*", color="#b4541f", ms=11)
        ax.annotate(lab, (Dm * sx, Jm * sy), (5, 4), textcoords="offset points")
    ax.set_xlim(0, 1.05 * spec.D_R * sx)
    ax.set_ylim(0, 1.05 * spec.J_A * sy)
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    ax.grid(color=GRID, lw=0.5)
    ax.legend(frameon=False)
    return ax


def plot_assessment(assessment, ax=None, label="absorber"):
    """Reduced plane for an assessed curve: lowest magnetic mass allowed for the same
    reflection curve (simplified form 22, dashed; general form 18, solid), and the absorber."""
    if ax is None:
        _, ax = plt.subplots(figsize=(4.6, 4))
    x, m22 = limit_curve(assessment, "simplified")
    _, m18 = limit_curve(assessment, "rigorous")
    ax.fill_between(x, 0, m18, color=GRAY, alpha=0.25, lw=0)
    ax.plot(x, m22, "--", color=GRAY, lw=1.3, label="limit, simplified form (20)")
    ax.plot(x, m18, "-", color=INK, lw=1.5, label="limit, general form (18)")
    ax.plot([assessment.x], [assessment.m], "*", color="#b4541f", ms=12, label=label)
    ax.set_xlim(0, max(1.05, 1.1 * assessment.x))
    ax.set_ylim(0, max(1.05, 1.15 * assessment.m))
    ax.set_xlabel("reduced thickness $x=D/D_R$")
    ax.set_ylabel("reduced magnetic mass $m=M/M_A$")
    ax.grid(color=GRID, lw=0.5)
    ax.legend(frameon=False)
    return ax
