"""The reduction proper: number density, absolute radius, refractivity and their companions.

SPEC_02 Step 3 (v0.6), the orchestration of `lib.reduction` on the loaded inputs and the frozen
anchor.

**Companions.** Every uncertainty here is the first order propagation of an uncertainty an input
file declares, added in quadrature with no covariance between inputs (SPEC_02 v0.4). Each
declared uncertainty is first converted to one standard deviation by its declared kind (SPEC_00
section 5 v0.13), so every companion is `1sigma`, and the conversions applied to the terms that
entered are listed with it. Unstated terms are left out and named. The partial derivatives are
kept, because the product records them (SPEC_02 Step 4).

**The anchor terms need marches.** The anchor radius and the label latitude reach every radius
through `r0`. Their derivatives are properties of the anchored march at the frozen pair:

* `dphi_c/dphi_g = 1 / (1 + dpsi/dphi_c)`, from differentiating the fixed point
  `phi_c = phi_g - psi(phi_c)`, with `dpsi/dphi_c` taken along the surface by a central
  difference on one march with `phi_c +- epsilon` as nodes;
* `dr0/dphi_c = r0 G_phi / g`, Eq. B3 itself at the solution;
* `dr0/dr_anchor`, the **total** derivative (v0.6): the central difference on the anchored
  march at the node `phi_c` with the latitude held, plus the coupling through the latitude,
  `dr0/dphi_c x dphi_c/dr_anchor`. `phi_c` is not an independent input; `psi` depends on radius,
  so moving the anchor moves the latitude too. `dphi_c/dr_anchor` comes from `dpsi/dr_anchor` at
  the node, from the same two marches, through the fixed point.

Since `phi_c + psi = phi_g` at the solution, `dpsi/d(.) = -dphi_c/d(.)` for the anchor radius
and `1 - dphi_c/dphi_g` for the label.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType

import numpy as np

from casspian.lib import reduction as red
from casspian.lib.constants import CODATA_RELEASE
from casspian.lib.control import ReductionInputs, ReductionManifest
from casspian.lib.gravity import G_phi_eff, g_eff_radial, g_eff_vector
from casspian.lib.schema import uncertainty_companion
from casspian.refrac.anchor import FrozenAnchor, geoid_setup

__all__ = ["Companion", "Reduction", "ReductionError", "anchor_partials", "composition_closure",
           "reduce_profile"]

#: Central difference steps of the anchor partials. The march converges to far below either
#: (anchoring residual about 1e-6 m), and Eq. B3 is nearly linear in both over these spans.
LATITUDE_STEP_RAD = 1.0e-4
ANCHOR_RADIUS_STEP_M = 1000.0

#: The closure rules of SPEC_00 section 6.2 v0.13, and the one this step propagates through.
CLOSURE_RULES = ("share_of_remainder", "normalized")
PROPAGATED_CLOSURE = "share_of_remainder"


class ReductionError(ValueError):
    """The inputs do not have a structure this step specifies the propagation for."""


@dataclass(frozen=True)
class Companion:
    """A `1sigma` uncertainty, the terms that entered and were left out, and the conversions."""

    value: np.ndarray
    included: tuple
    unstated: tuple
    conversions: tuple = ()


@dataclass(frozen=True)
class Reduction:
    """Everything Step 3 computes, per level and scalar. Angles in radians."""

    anchor: FrozenAnchor
    anchor_level: int
    h_ref_m: float
    pressure_Pa: np.ndarray
    temperature_K: np.ndarray
    height_m: np.ndarray
    radius_m: np.ndarray
    height_above_anchor_isobar_m: np.ndarray
    number_density_m3: np.ndarray
    mean_refractivity_m3: np.ndarray
    mean_molar_mass_kg_mol: np.ndarray
    refractivity: np.ndarray
    companions: MappingProxyType
    partials: MappingProxyType
    terms: MappingProxyType
    closure: MappingProxyType
    codata_release: str


def _declared(values, kind, label, where):
    """A declared uncertainty at one standard deviation, with its conversion text, if any.

    An uncertainty that is NaN everywhere enters nothing, so its kind does not matter. A finite
    one with no declared kind is refused, since there is then no rule to convert it by.
    """
    values = np.asarray(values, dtype="float64")
    if not np.any(np.isfinite(values)):
        return values, None
    if not kind:
        raise ReductionError(f"{where} is finite but declares no uncertainty_kind; SPEC_00 "
                             "section 5 v0.13 has no rule to convert it by.")
    converted, conversion = red.to_one_sigma(values, kind)
    return converted, (f"{label}: {conversion}" if conversion else None)


def _conversions(included, table):
    """The conversion texts of the terms that actually entered, in order, without repeats."""
    out = []
    for label in included:
        text = table.get(label)
        if text and text not in out:
            out.append(text)
    return tuple(out)


def composition_closure(composition, tolerance=1e-9):
    """Read the closure declaration of a kind C file and check it against the values.

    SPEC_00 section 6.2 v0.13: `closure_rule` and `closure_species` are declared by the tool
    that wrote the file. This step propagates through `share_of_remainder` only, with the share
    species named first and its partner second, every other species assigned first. The
    declaration is then checked against the values by the inference of the v0.5 build: the
    share species must carry a stated mole fraction uncertainty and its ratio to exactly the
    declared partner must be constant at every level. Any disagreement is refused.
    """
    root = composition.dataset
    rule = str(root.attrs.get("closure_rule", ""))
    declared = str(root.attrs.get("closure_species", "")).split()
    if rule not in CLOSURE_RULES:
        raise ReductionError(
            f"the composition file declares closure_rule = {rule!r}; SPEC_00 section 6.2 v0.13 "
            f"defines {list(CLOSURE_RULES)}, and a kind C file must declare one."
        )
    if rule != PROPAGATED_CLOSURE:
        raise ReductionError(
            f"the composition file declares closure_rule = {rule!r}; SPEC_02 Step 3 specifies "
            f"the constrained derivative for {PROPAGATED_CLOSURE!r} only."
        )
    names = [str(v) for v in composition["species"].dataset["species_name"].values]
    if len(declared) != 2 or any(d not in names for d in declared):
        raise ReductionError(
            f"closure_species = {' '.join(declared)!r} must name the share species then its "
            f"partner, both among {names}."
        )
    x = np.stack([np.asarray(root[f"x_{n}"].values, dtype="float64") for n in names])
    stated = [i for i, n in enumerate(names) if uncertainty_companion(f"x_{n}") in root.variables]

    found = []
    for i in stated:
        partners = []
        for j in range(len(names)):
            if j == i:
                continue
            remainder = x[i] + x[j]
            if np.any(remainder <= 0.0):
                continue
            s = x[i] / remainder
            if float(s.max() - s.min()) <= tolerance:
                partners.append(j)
        if len(partners) == 1:
            found.append((i, partners[0]))
    share, partner = names.index(declared[0]), names.index(declared[1])
    if found != [(share, partner)]:
        shown = [(names[i], names[j]) for i, j in found]
        raise ReductionError(
            f"the composition file declares a share of remainder of {declared[0]} with "
            f"{declared[1]}, but its values show {shown} among the species with a stated mole "
            f"fraction uncertainty {[names[i] for i in stated]}."
        )
    return names, x, share, partner


def anchor_partials(inputs: ReductionInputs, manifest: ReductionManifest, anchor: FrozenAnchor):
    """The derivatives of `phi_c`, `psi` and `r0` in the anchor radius and the label."""
    setup = geoid_setup(inputs, manifest)
    constants, u_of_phi = setup.constants, setup.u_of_phi
    phi_c, r0 = anchor.phi_c_rad, anchor.r0_m

    eps = LATITUDE_STEP_RAD
    nodes = np.array([phi_c - eps, phi_c, phi_c + eps])
    along = setup.march(nodes)
    psi_along = g_eff_vector(u_of_phi(nodes), along.radius, nodes, *constants)[3]
    dpsi_dphi_c = float((psi_along[2] - psi_along[0]) / (2.0 * eps))
    dphi_c_dphi_g = 1.0 / (1.0 + dpsi_dphi_c)
    dr0_dphi_c_difference = float((along.radius[2] - along.radius[0]) / (2.0 * eps))

    g = float(g_eff_radial(anchor.u_at_phi_c_ms, r0, phi_c, *constants))
    G_phi = float(G_phi_eff(anchor.u_at_phi_c_ms, r0, phi_c, *constants))
    dr0_dphi_c = r0 * G_phi / g

    delta = ANCHOR_RADIUS_STEP_M
    at = np.array([phi_c])
    plus = setup.march(at, r_anchor_m=setup.r_anchor_m + delta)
    minus = setup.march(at, r_anchor_m=setup.r_anchor_m - delta)
    dr0_dr_anchor_held = float((plus.radius[0] - minus.radius[0]) / (2.0 * delta))
    psi_plus = g_eff_vector(u_of_phi(at), plus.radius, at, *constants)[3]
    psi_minus = g_eff_vector(u_of_phi(at), minus.radius, at, *constants)[3]
    dpsi_dr_anchor_at_node = float((psi_plus[0] - psi_minus[0]) / (2.0 * delta))
    dphi_c_dr_anchor = -dpsi_dr_anchor_at_node * dphi_c_dphi_g
    coupling = dr0_dphi_c * dphi_c_dr_anchor

    return {
        "dr0_dr_anchor": dr0_dr_anchor_held + coupling,
        "dr0_dr_anchor_latitude_held": dr0_dr_anchor_held,
        "dr0_dr_anchor_latitude_coupling": coupling,
        "dr0_dphi_c_m_per_rad": dr0_dphi_c,
        "dr0_dphi_c_difference_m_per_rad": dr0_dphi_c_difference,
        "dpsi_dphi_c_along_surface": dpsi_dphi_c,
        "dphi_c_dphi_g": dphi_c_dphi_g,
        "dpsi_dphi_g": 1.0 - dphi_c_dphi_g,
        "dphi_c_dr_anchor_rad_per_m": dphi_c_dr_anchor,
        "dpsi_dr_anchor_rad_per_m": -dphi_c_dr_anchor,
        "latitude_step_rad": eps,
        "anchor_radius_step_m": delta,
    }


def reduce_profile(inputs: ReductionInputs, manifest: ReductionManifest,
                   anchor: FrozenAnchor) -> Reduction:
    """B3.1, B1 and B3.3 on every level, with the first order companions. SPEC_02 Step 3."""
    thermo, geodesy = inputs.thermo, inputs.geodesy
    p = np.asarray(thermo["pressure_Pa"].values, dtype="float64")
    T = np.asarray(thermo["temperature_K"].values, dtype="float64")
    h = np.asarray(thermo["height_m"].values, dtype="float64")
    level = int(np.flatnonzero(p == manifest.anchor_isobar_Pa)[0])
    h_ref = float(h[level])

    def thermo_uncertainty(name, label):
        variable = thermo[name]
        return _declared(variable.values, str(variable.attrs.get("uncertainty_kind", "")),
                         label, f"kind T {name}")

    dp, convert_p = thermo_uncertainty("pressure_uncertainty_Pa", "pressure")
    dT, convert_T = thermo_uncertainty("temperature_uncertainty_K", "temperature")
    dh, convert_h = thermo_uncertainty("height_uncertainty_m", "height")
    dh_ref = float(dh[level])

    # B3.1 and its relative companion.
    n = red.number_density(p, T)
    n_rel, n_included, n_unstated = red.quadrature({"pressure": dp / p, "temperature": dT / T})

    # B3.3 and the composition companion, with the declared closure applied.
    names, x, share, partner = composition_closure(inputs.composition)
    composition = inputs.composition
    root = composition.dataset
    species = composition["species"].dataset
    R_i = np.asarray(species["refractivity_per_molecule_m3"].values, dtype="float64")
    M_i = np.asarray(species["molar_mass_kg_mol"].values, dtype="float64")
    R_bar = red.mean_over_species(x, R_i)
    m_bar = red.mean_over_species(x, M_i)

    share_name = uncertainty_companion(f"x_{names[share]}")
    ds, convert_share = _declared(root[share_name].values,
                                  str(root[share_name].attrs.get("uncertainty_kind", "")),
                                  "share", f"kind C {share_name}")
    first_assigned = {}
    for k in range(len(names)):
        if k in (share, partner):
            continue
        name = uncertainty_companion(f"x_{names[k]}")
        if name in root.variables:
            first_assigned[k] = _declared(root[name].values,
                                          str(root[name].attrs.get("uncertainty_kind", "")),
                                          f"x_{names[k]}", f"kind C {name}")
        else:
            first_assigned[k] = (np.full(p.shape, np.nan), None)
    dR_i, convert_R_i = _declared(species["refractivity_uncertainty_m3"].values,
                                  str(species["refractivity_uncertainty_m3"].attrs.get(
                                      "uncertainty_kind", "")),
                                  "per_molecule", "kind C species refractivity_uncertainty_m3")

    def composition_term(q, dq_i):
        _, d_ds, d_dx_first = red.share_of_remainder_partials(x, q, share, partner)
        pieces = {"share": d_ds * ds}
        table = {"share": convert_share}
        for k, derivative in d_dx_first.items():
            dx_k, conversion = first_assigned[k]
            pieces[f"x_{names[k]}"] = derivative * dx_k
            table[f"x_{names[k]}"] = conversion
        if dq_i is not None:
            for i, name in enumerate(names):
                pieces[f"per_molecule_{name}"] = x[i] * dq_i[i]
                table[f"per_molecule_{name}"] = convert_R_i
        total, included, unstated = red.quadrature(pieces)
        return Companion(total, included, unstated, _conversions(included, table))

    R_bar_companion = composition_term(R_i, dR_i)
    m_bar_companion = composition_term(M_i, None)

    N = red.refractivity(n, R_bar)
    N_rel, N_included, N_unstated = red.quadrature({
        "pressure": dp / p, "temperature": dT / T,
        "composition": R_bar_companion.value / R_bar})
    N_conversions = _conversions(N_included, {"pressure": convert_p, "temperature": convert_T})
    if "composition" in N_included:
        N_conversions += tuple(c for c in R_bar_companion.conversions if c not in N_conversions)

    # B1 and the anchor terms.
    radius = red.absolute_radius(h, h_ref, anchor.r0_m)
    partials = anchor_partials(inputs, manifest, anchor)
    surfaces = np.asarray(geodesy["surface_pressure_Pa"].values, dtype="float64")
    row = int(np.flatnonzero(surfaces == manifest.anchor_surface_Pa)[0])
    anchor_name = uncertainty_companion(manifest.anchor_quantity)
    anchor_declared = float(geodesy[anchor_name].values[row])
    anchor_kind = str(geodesy[anchor_name].attrs.get("uncertainty_kind", ""))
    d_anchor, convert_anchor = _declared(anchor_declared, anchor_kind, "anchor_radius",
                                         f"kind D {anchor_name}")
    d_anchor = float(d_anchor)
    # The refrac boundary: the label uncertainty is stated in degrees.
    label_declared = float(np.radians(float(thermo.attrs.get(
        "latitude_planetographic_deg_uncertainty", np.nan))))
    label_kind = str(thermo.attrs.get("latitude_planetographic_deg_uncertainty_kind", ""))
    d_label, convert_label = _declared(label_declared, label_kind, "label_latitude",
                                       "kind T latitude_planetographic_deg_uncertainty")
    d_label = float(d_label)
    scalar_table = {"anchor_radius": convert_anchor, "label_latitude": convert_label}

    terms = {
        "r0_anchor_radius_m": abs(partials["dr0_dr_anchor"]) * d_anchor,
        "r0_label_latitude_m": abs(partials["dr0_dphi_c_m_per_rad"]
                                   * partials["dphi_c_dphi_g"]) * d_label,
        "phi_c_anchor_radius_rad": abs(partials["dphi_c_dr_anchor_rad_per_m"]) * d_anchor,
        "phi_c_label_latitude_rad": abs(partials["dphi_c_dphi_g"]) * d_label,
        "psi_anchor_radius_rad": abs(partials["dpsi_dr_anchor_rad_per_m"]) * d_anchor,
        "psi_label_latitude_rad": abs(partials["dpsi_dphi_g"]) * d_label,
    }
    scalar = {}
    for key in ("r0", "phi_c", "psi"):
        unit = "m" if key == "r0" else "rad"
        total, included, unstated = red.quadrature({
            "anchor_radius": terms[f"{key}_anchor_radius_{unit}"],
            "label_latitude": terms[f"{key}_label_latitude_{unit}"],
        })
        scalar[key] = Companion(np.float64(total), included, unstated,
                                _conversions(included, scalar_table))

    # At the anchor level h and h_ref are one tabulated value, so h - h_ref is exactly zero and
    # carries no height uncertainty; elsewhere the two are different table entries.
    at_anchor = np.arange(p.size) == level
    height_term = np.where(at_anchor & np.isfinite(dh), 0.0, dh)
    anchor_height_term = (np.where(at_anchor, 0.0, dh_ref) if np.isfinite(dh_ref)
                          else np.full(p.shape, np.nan))
    r_unc, r_included, r_unstated = red.quadrature({
        "anchor_radius": np.full(p.shape, terms["r0_anchor_radius_m"]),
        "label_latitude": np.full(p.shape, terms["r0_label_latitude_m"]),
        "height": height_term,
        "anchor_height": anchor_height_term,
    })
    radius_table = dict(scalar_table, height=convert_h, anchor_height=convert_h)

    companions = {
        "number_density_uncertainty_m3": Companion(
            n * n_rel, n_included, n_unstated,
            _conversions(n_included, {"pressure": convert_p, "temperature": convert_T})),
        "refractivity_uncertainty": Companion(N * N_rel, N_included, N_unstated, N_conversions),
        "radius_uncertainty_m": Companion(r_unc, r_included, r_unstated,
                                          _conversions(r_included, radius_table)),
        "mean_refractivity_uncertainty_m3": R_bar_companion,
        "mean_molar_mass_uncertainty_kg_mol": m_bar_companion,
        "anchor_isobar_radius_uncertainty_m": scalar["r0"],
        "latitude_planetocentric_uncertainty_rad": scalar["phi_c"],
        "psi_uncertainty_rad": scalar["psi"],
    }
    partials.update({
        "declared_anchor_radius_uncertainty_m": anchor_declared,
        "declared_anchor_radius_uncertainty_kind": anchor_kind,
        "anchor_radius_uncertainty_1sigma_m": d_anchor,
        "declared_label_latitude_uncertainty_rad": label_declared,
        "declared_label_latitude_uncertainty_kind": label_kind,
        "label_latitude_uncertainty_1sigma_rad": d_label,
    })

    return Reduction(
        anchor=anchor,
        anchor_level=level,
        h_ref_m=h_ref,
        pressure_Pa=p,
        temperature_K=T,
        height_m=h,
        radius_m=radius,
        height_above_anchor_isobar_m=h - h_ref,
        number_density_m3=n,
        mean_refractivity_m3=R_bar,
        mean_molar_mass_kg_mol=m_bar,
        refractivity=N,
        companions=MappingProxyType(companions),
        partials=MappingProxyType(partials),
        terms=MappingProxyType(terms),
        closure=MappingProxyType({"rule": str(root.attrs["closure_rule"]),
                                  "species": tuple(names), "share": names[share],
                                  "partner": names[partner],
                                  "assigned_first": tuple(n for i, n in enumerate(names)
                                                          if i not in (share, partner))}),
        codata_release=CODATA_RELEASE,
    )
