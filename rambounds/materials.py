"""Magnetic materials: the constants I_mu and K_mu of the paper (Eqs. 9-13, Table I).

I_mu = k_A * (gamma_bar * 4 pi Ms)^2   (GHz^2), the generalized Snoek bound of the material,
K_mu = I_mu / rho                       (mm.GHz^2 per kg/m^2), its efficiency per unit mass.

A layer of thickness d (mm) containing a volume fraction p of the material brings an
amount of magnetic material J_mu = p d I_mu (mm.GHz^2) and a magnetic mass M = rho p d (kg/m^2),
so that J_mu = K_mu M.
"""
from dataclasses import dataclass

from .constants import C, K_DB, GAMMA_BAR


@dataclass(frozen=True)
class Material:
    """A magnetic material.

    name      : free text
    k_A       : geometric factor (1 for films excited along their hard axis,
                1/2 for equal transverse permeabilities, 1/3 for isotropic materials)
    Ms_kG     : saturation magnetization 4 pi Ms, in kG
    rho       : density of the magnetic material, in g/cm^3
    """
    name: str
    k_A: float
    Ms_kG: float
    rho: float

    @property
    def I_mu(self):
        """Generalized Snoek bound, GHz^2 (Eq. 9)."""
        return self.k_A * (GAMMA_BAR * self.Ms_kG) ** 2

    @property
    def K_mu(self):
        """Mass efficiency, mm.GHz^2 per kg/m^2 (Eq. 13)."""
        return self.I_mu / self.rho

    @property
    def I_mu_dB(self):
        """I_mu in dB.GHz/mm (Eq. 24, Table I)."""
        return self.I_mu * K_DB / C

    @property
    def K_mu_dB(self):
        """K_mu in dB.GHz per kg/m^2 (Eq. 24, Table I)."""
        return self.K_mu * K_DB / C

    def amount(self, M):
        """Amount of magnetic material J_mu (mm.GHz^2) for a magnetic mass M (kg/m^2)."""
        return self.K_mu * M

    def mass(self, J):
        """Magnetic mass M (kg/m^2) for an amount of magnetic material J (mm.GHz^2)."""
        return J / self.K_mu


# Typical values of Table I of the paper (4 pi Ms in kG, density in g/cm^3).
MATERIALS = {
    "MnZn ferrite": Material("MnZn ferrite", 1 / 3, 5.0, 4.9),
    "NiZn ferrite": Material("NiZn ferrite", 1 / 3, 4.0, 5.0),
    "carbonyl iron powder": Material("carbonyl iron powder", 1 / 3, 21.5, 7.9),
    "permalloy powder": Material("permalloy powder", 1 / 3, 10.0, 8.7),
    "iron film": Material("iron film", 1.0, 21.5, 7.9),
    "permalloy film": Material("permalloy film", 1.0, 10.0, 8.7),
    "CoNbZr film": Material("CoNbZr film", 1.0, 10.0, 8.8),
}


def get_material(material):
    """Accept a Material or the name of one of the MATERIALS."""
    if isinstance(material, Material):
        return material
    try:
        return MATERIALS[material]
    except KeyError:
        raise KeyError("unknown material %r; choose one of %s, or build a Material(name, k_A, Ms_kG, rho)"
                       % (material, ", ".join(MATERIALS)))
