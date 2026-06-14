"""Pure-gemmi helpers for multi-domain NCS averaging with `dm`.

This module is deliberately CCP4-binary-free: it depends only on gemmi (which
is in the slim server bundle) and numpy. All it does is turn a model + domain
definitions into the operators and masks that `dm` consumes. The plugin
(dm_multidomain.py) owns the binary execution.

The operator-direction convention was validated empirically against the AHIR
driver case: `dm` wants operators that map the REFERENCE (masked) copy onto
each NCS copy (x' = R x + t, identity first). Inverting them corrupts the
density. `operator_ref_to_copy` is self-orienting, so it is robust regardless
of any library's internal fixed/moving convention.
"""
from __future__ import annotations

import gemmi
import numpy as np

# Identity 3x3 rotation as a dm ROTA MATRIX string (row-major) + zero TRAN.
IDENTITY_DM = ("ROTA MATRIX 1 0 0 0 1 0 0 0 1", "TRAN 0 0 0")


def parse_domain(spec):
    """'shell:340-485:average' -> dict(name, lo, hi, mode).

    mode is one of average | refine | exclude.
    """
    parts = spec.split(":")
    if len(parts) != 3:
        raise ValueError(
            f"domain spec {spec!r} must be name:lo-hi:mode")
    name, rng, mode = parts
    lo, hi = (int(x) for x in rng.split("-"))
    mode = mode.strip().lower()
    if mode not in ("average", "refine", "exclude"):
        raise ValueError(f"domain {name}: mode {mode!r} not in average/refine/exclude")
    return dict(name=name.strip(), lo=lo, hi=hi, mode=mode)


def _is_amino_acid(res):
    info = gemmi.find_tabulated_residue(res.name)
    return bool(info and info.is_amino_acid())


def protein_chains(model):
    """Ordered list of chain names that contain amino acids."""
    out = []
    for chain in model:
        if any(_is_amino_acid(r) for r in chain):
            out.append(chain.name)
    return out


def detect_reference_and_copies(model, reference=None, copies=None):
    """Resolve reference + copy chain ids, auto-detecting from the model.

    If reference/copies are not supplied, the first protein chain is the
    reference and the rest are the copies.
    """
    prot = protein_chains(model)
    if not prot:
        raise ValueError("no protein chains found in model")
    ref = reference or prot[0]
    if copies:
        cps = [c.strip() for c in copies if c.strip()]
    else:
        cps = [c for c in prot if c != ref]
    return ref, cps


def _chain(model, name):
    for c in model:
        if c.name == name:
            return c
    return None


def ca_positions(model, chain_name, lo, hi):
    """{seqid -> Position} for CA atoms of a chain within [lo, hi]."""
    out = {}
    chain = _chain(model, chain_name)
    if chain is None:
        return out
    for res in chain:
        n = res.seqid.num
        if lo <= n <= hi:
            a = res.find_atom("CA", "*")
            if a is not None:
                out[n] = a.pos
    return out


def all_atom_positions(model, chain_name, lo, hi):
    pts = []
    chain = _chain(model, chain_name)
    if chain is None:
        return pts
    for res in chain:
        if lo <= res.seqid.num <= hi:
            for a in res:
                pts.append(a.pos)
    return pts


def _rmsd(a, b):
    return float(np.sqrt(np.mean(
        [(p.x - q.x) ** 2 + (p.y - q.y) ** 2 + (p.z - q.z) ** 2
         for p, q in zip(a, b)])))


def operator_ref_to_copy(model, ref, copy, lo, hi):
    """(gemmi.Transform T, rmsd) with T mapping REFERENCE coords -> COPY coords.

    Self-orienting: compute a superposition, then keep whichever of T / T^-1
    actually maps ref->copy. Independent of gemmi's internal convention.
    Raises ValueError if fewer than 3 common CA atoms are available.
    """
    ref_ca = ca_positions(model, ref, lo, hi)
    cp_ca = ca_positions(model, copy, lo, hi)
    common = sorted(set(ref_ca) & set(cp_ca))
    if len(common) < 3:
        raise ValueError(
            f"<3 common CA atoms between {ref} and {copy} in {lo}-{hi}")
    ref_pts = [ref_ca[n] for n in common]
    cp_pts = [cp_ca[n] for n in common]

    sup = gemmi.superpose_positions(cp_pts, ref_pts)
    T = sup.transform
    Tinv = T.inverse()
    if _rmsd([Tinv.apply(p) for p in ref_pts], cp_pts) < \
       _rmsd([T.apply(p) for p in ref_pts], cp_pts):
        T = Tinv
    return T, _rmsd([T.apply(p) for p in ref_pts], cp_pts)


def transform_to_dm(T):
    """gemmi.Transform -> ('ROTA MATRIX ...', 'TRAN ...') lines (x' = R x + t)."""
    m = T.mat.tolist()
    v = T.vec
    rota = "ROTA MATRIX " + " ".join(
        f"{m[i][j]:.7f}" for i in range(3) for j in range(3))
    tran = f"TRAN {v.x:.5f} {v.y:.5f} {v.z:.5f}"
    return rota, tran


def domain_operators(model, ref, copies, lo, hi):
    """List of (rota, tran) dm lines for a domain: identity first, then ref->copy_n.

    Also returns the per-copy RMSDs for diagnostics/logging.
    """
    ops = [IDENTITY_DM]
    rmsds = {}
    for c in copies:
        T, rmsd = operator_ref_to_copy(model, ref, c, lo, hi)
        ops.append(transform_to_dm(T))
        rmsds[c] = rmsd
    return ops, rmsds


def write_domain_mask(model, ref, lo, hi, cell, spacegroup, path,
                      spacing=1.0, radius=2.5):
    """Write a CCP4 mode-0 mask covering the reference copy's domain.

    No symmetry expansion: the mask covers a single copy, as `dm` expects for
    an NCS averaging mask. Returns the number of grid points set.
    """
    g = gemmi.Int8Grid()
    g.unit_cell = cell
    g.spacegroup = spacegroup
    g.set_size_from_spacing(spacing, gemmi.GridSizeRounding.Up)
    for p in all_atom_positions(model, ref, lo, hi):
        g.set_points_around(p, radius=radius, value=1)
    ccp4 = gemmi.Ccp4Mask()
    ccp4.grid = g
    ccp4.update_ccp4_header()
    ccp4.write_ccp4_map(path)
    return int(np.count_nonzero(np.array(g, copy=False)))


def build_keyword_script(domains, operators_by_domain, solc, ncycle,
                         mode_solv=True, mode_hist=True,
                         labin="FP=FP SIGFP=SIGFP PHIO=PHIO FOMO=FOMO",
                         labout="FDM=FDM PHIDM=PHIDM FOMDM=FOMDM",
                         ncross=1):
    """Assemble the `dm` keyword (stdin) script.

    domains: ordered list of parsed-domain dicts (exclude ones are skipped).
    operators_by_domain: {name -> [(rota, tran), ...]}.
    The DOMAIN index emitted here MUST match the NCSIN<n> order on the command
    line (the caller is responsible for keeping them in step).
    """
    modes = ["AVER"]
    if mode_solv:
        modes.insert(0, "SOLV")
    if mode_hist:
        modes.insert(1 if mode_solv else 0, "HIST")
    # NB: dm computes free-R from the FREE column in LABIN automatically; it
    # does NOT accept "FREE <ncross>" as an NCYCLE sub-keyword (that is dmmulti
    # syntax) -- so ncross is not emitted here.
    lines = [f"SOLC {solc}", "MODE " + " ".join(modes), f"NCYCLE {ncycle}"]
    dnum = 0
    for d in domains:
        if d["mode"] == "exclude":
            continue
        dnum += 1
        refine = " REFINE" if d["mode"] == "refine" else ""
        for k, (rota, tran) in enumerate(operators_by_domain[d["name"]]):
            lines.append(f"AVERAGE DOMAIN {dnum}" + (refine if k == 0 else ""))
            lines.append(rota)
            lines.append(tran)
    lines.append("LABIN " + labin)
    lines.append("LABOUT " + labout)
    return lines


def estimate_solvent_fraction(structure, protein_density=1.23):
    """Rough Matthews-style solvent fraction from model mass and cell.

    Counts atoms in the ASU model, scales by the number of spacegroup
    operations to fill the cell, and compares occupied volume to cell volume.
    Returns a value in (0, 1) or None if it cannot be estimated.
    """
    cell = structure.cell
    sg = structure.find_spacegroup()
    if cell is None or sg is None or cell.volume <= 0:
        return None
    mass = 0.0
    for chain in structure[0]:
        for res in chain:
            for atom in res:
                el = atom.element
                mass += el.weight if el else 0.0
    nops = len([op for op in sg.operations()])
    occupied = protein_density * mass * nops  # A^3 (1.23 A^3 per Da)
    frac = 1.0 - occupied / cell.volume
    if not (0.0 < frac < 1.0):
        return None
    return round(frac, 3)
