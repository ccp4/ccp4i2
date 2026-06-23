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


# ---------------------------------------------------------------------------
# Role / instance / segment model
#
# A rigid body (dm "domain") is a set of SEGMENTS, each a (role, lo, hi) residue
# range. A role names an entity (e.g. "CDK", "cyclin"); in the homomer case
# there is a single implicit role. An NCS INSTANCE maps roles -> chain ids for
# one copy of the assembly (instances[0] is the reference). This is what lets a
# body cross chains (the CDK C-helix travelling with the cyclin N-lobe) and lets
# AmBn complexes work: a partial instance simply omits a role.
# ---------------------------------------------------------------------------
_IMPLICIT_ROLE = "_"


def parse_segments(spec):
    """'cyclin:10-95,CDK:45-60' -> [('cyclin',10,95),('CDK',45,60)].

    A token without a 'role:' prefix uses the single implicit role, so the
    homomer case stays terse: '340-485' -> [('_',340,485)].
    """
    out = []
    for tok in str(spec).split(","):
        tok = tok.strip()
        if not tok:
            continue
        if ":" in tok:
            role, rng = tok.split(":", 1)
            role = role.strip() or _IMPLICIT_ROLE
        else:
            role, rng = _IMPLICIT_ROLE, tok
        lo, hi = (int(x) for x in rng.split("-"))
        out.append((role, lo, hi))
    if not out:
        raise ValueError(f"no segments parsed from {spec!r}")
    return out


def parse_assembly_rows(rows):
    """List of 'role=chain' token rows -> list of instance dicts {role: chain}.

    Index 0 is the reference instance. A bare chain id (no '=') uses the
    implicit role, so a homomer assembly is just ['A','B','C',...]; a CDK/cyclin
    assembly is ['CDK=A cyclin=B', 'CDK=C cyclin=D', ...]. Blank/partial roles
    (AmBn) are simply omitted from a row.
    """
    instances = []
    for row in rows:
        inst = {}
        for tok in str(row).replace(",", " ").split():
            if "=" in tok:
                role, chain = tok.split("=", 1)
                role = role.strip() or _IMPLICIT_ROLE
            else:
                role, chain = _IMPLICIT_ROLE, tok
            chain = chain.strip()
            if chain:
                inst[role] = chain
        if inst:
            instances.append(inst)
    return instances


def instance_label(inst):
    """Stable human label for an instance, e.g. {'CDK':'A','cyclin':'B'}->'A+B'."""
    return "+".join(inst[r] for r in sorted(inst))


def body_segment_roles(segments):
    return {role for role, _, _ in segments}


def body_name(segments):
    """Filesystem/dm-safe name for a body, e.g. 'cyclin_10_95__CDK_45_60'
    (the implicit role is elided so the homomer case reads '340_485')."""
    return "__".join(
        (f"{role}_" if role != _IMPLICIT_ROLE else "") + f"{lo}_{hi}"
        for role, lo, hi in segments)


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


def _lattice_shift(positions, cell):
    """Integer lattice vector that brings the centroid of `positions` into the
    primary cell [0,1). Computed from CA-style positions so the mask and the
    operators use an IDENTICAL shift."""
    import math
    fr = [cell.fractionalize(p) for p in positions]
    n = len(fr)
    cen = (sum(f.x for f in fr) / n, sum(f.y for f in fr) / n,
           sum(f.z for f in fr) / n)
    return (-math.floor(cen[0]), -math.floor(cen[1]), -math.floor(cen[2]))


def _shift(p, s, cell):
    f = cell.fractionalize(p)
    return cell.orthogonalize(gemmi.Fractional(f.x + s[0], f.y + s[1],
                                               f.z + s[2]))


def body_ca_correspondence(model, ref_inst, copy_inst, segments):
    """Ordered (ref_pts, copy_pts) CA pairs for a (possibly multi-chain) rigid
    body, matched on (role, seqid) across all its segments.

    Matching on (role, seqid) -- not bare seqid -- is what keeps a cross-chain
    body honest: cyclin residue 50 and CDK residue 50 don't collide because they
    live under different roles. The two lists are element-aligned by
    construction, so they feed straight into superposition no matter how many
    chains the body spans. Segments whose role is absent in either instance are
    skipped (a partial complex contributes only the roles it has).
    """
    ref_pts, cp_pts = [], []
    for role, lo, hi in segments:
        ref_chain = ref_inst.get(role)
        cp_chain = copy_inst.get(role)
        if not ref_chain or not cp_chain:
            continue
        r = ca_positions(model, ref_chain, lo, hi)
        c = ca_positions(model, cp_chain, lo, hi)
        for n in sorted(set(r) & set(c)):
            ref_pts.append(r[n])
            cp_pts.append(c[n])
    return ref_pts, cp_pts


def _superpose_ref_to_copy(ref_pts, cp_pts):
    """gemmi.Transform mapping ref_pts -> cp_pts, self-orienting (keeps whichever
    of T / T^-1 actually maps ref onto copy, independent of gemmi's internal
    fixed/moving convention) + the achieved RMSD."""
    sup = gemmi.superpose_positions(cp_pts, ref_pts)
    T = sup.transform
    Tinv = T.inverse()
    if _rmsd([Tinv.apply(p) for p in ref_pts], cp_pts) < \
       _rmsd([T.apply(p) for p in ref_pts], cp_pts):
        T = Tinv
    return T, _rmsd([T.apply(p) for p in ref_pts], cp_pts)


def operator_ref_to_copy_body(model, ref_inst, copy_inst, segments, cell=None):
    """(gemmi.Transform T, rmsd) with T mapping the REFERENCE rigid body onto its
    copy. Generalises operator_ref_to_copy to a multi-segment, multi-chain body.

    Raises ValueError if fewer than 3 common CA atoms are available.

    T maps the reference at its REAL position (no lattice shift) onto the copy,
    so it shares the frame of the averaging mask -- which now carries a
    contiguous copy of the domain at its real place via a signed NSTART
    sub-volume (see write_partitioned_masks), not a PBC-wrapped image. Mask and
    operator move together in the crystal frame; dm averages along the operator
    from the masked region with no spurious lattice-vector displacement. `cell`
    is accepted for back-compat but no longer used.
    """
    ref_pts, cp_pts = body_ca_correspondence(model, ref_inst, copy_inst, segments)
    if len(ref_pts) < 3:
        raise ValueError(
            f"<3 common CA atoms for body {body_name(segments)} between "
            f"{instance_label(ref_inst)} and {instance_label(copy_inst)}")
    return _superpose_ref_to_copy(ref_pts, cp_pts)


def operator_ref_to_copy(model, ref, copy, lo, hi, cell=None):
    """Back-compat single-chain shim: one implicit-role segment over one chain.
    See operator_ref_to_copy_body for the general case.
    """
    return operator_ref_to_copy_body(
        model, {_IMPLICIT_ROLE: ref}, {_IMPLICIT_ROLE: copy},
        [(_IMPLICIT_ROLE, lo, hi)], cell=cell)


def transform_to_dm(T):
    """gemmi.Transform -> ('ROTA MATRIX ...', 'TRAN ...') lines (x' = R x + t)."""
    m = T.mat.tolist()
    v = T.vec
    rota = "ROTA MATRIX " + " ".join(
        f"{m[i][j]:.7f}" for i in range(3) for j in range(3))
    tran = f"TRAN {v.x:.5f} {v.y:.5f} {v.z:.5f}"
    return rota, tran


def body_operators(model, instances, segments, cell=None):
    """[(rota, tran), ...] dm lines for a rigid body: identity first, then
    ref->copy for each copy instance that contains EVERY role the body uses.

    instances[0] is the reference; the rest are copies. A copy missing one of
    the body's roles (a partial complex -- e.g. the spare A in A4B3 for a body
    that also needs B) is silently skipped for this body. Also returns the
    per-copy RMSDs (keyed by instance label) for diagnostics/logging.
    """
    ref = instances[0]
    roles = body_segment_roles(segments)
    ops = [IDENTITY_DM]
    rmsds = {}
    for copy in instances[1:]:
        if not all(copy.get(r) for r in roles):
            continue
        T, rmsd = operator_ref_to_copy_body(model, ref, copy, segments, cell=cell)
        ops.append(transform_to_dm(T))
        rmsds[instance_label(copy)] = rmsd
    return ops, rmsds


def domain_operators(model, ref, copies, lo, hi, cell=None):
    """Back-compat single-chain shim: one implicit-role segment averaged across
    single-chain copies. See body_operators for the general case.
    """
    instances = [{_IMPLICIT_ROLE: ref}] + [{_IMPLICIT_ROLE: c} for c in copies]
    return body_operators(model, instances, [(_IMPLICIT_ROLE, lo, hi)], cell=cell)


def body_ca_positions(model, instance, segments):
    """Reference CA positions for a body (union of its segments) -- used to
    compute the single lattice shift shared by mask and operators."""
    out = []
    for role, lo, hi in segments:
        chain = instance.get(role)
        if chain:
            out += list(ca_positions(model, chain, lo, hi).values())
    return out


def body_all_atom_positions(model, instance, segments):
    """All-atom positions for a body (union of its segments, across chains)."""
    pts = []
    for role, lo, hi in segments:
        chain = instance.get(role)
        if chain:
            pts += all_atom_positions(model, chain, lo, hi)
    return pts


def _parent_sampling(cell, spacegroup, spacing):
    """Symmetry-compatible full-cell sampling (MX, MY, MZ) -- the grid dm's
    working map uses. Sub-volume masks share this sampling (and the crystal cell)
    so they overlay that map voxel-for-voxel."""
    g = gemmi.Int8Grid()
    g.unit_cell = cell
    g.spacegroup = spacegroup
    g.set_size_from_spacing(spacing, gemmi.GridSizeRounding.Up)
    return (g.nu, g.nv, g.nw)


def _build_disjoint_box(model, ref_inst, segments_list, cell, spacegroup,
                        spacing, radius):
    """Build per-body masks on a SHARED padded sub-volume box and resolve their
    overlaps, returning (arrs, NS, M, sp, n_contested).

    The box is a contiguous parallelepiped in the parent grid (sampling M, signed
    origin NS in voxels) padded by enough voxels that the `radius` dilation never
    reaches a face -- so set_points_around's PBC wrap is a no-op and each body is
    placed at its REAL position (no lattice shift). All bodies share this one box
    so the nearest-atom overlap competition runs in a single index frame; arrs[b]
    is body b's disjoint int8 mask on that box.
    """
    M = _parent_sampling(cell, spacegroup, spacing)
    sp = (cell.a / M[0], cell.b / M[1], cell.c / M[2])
    body_atoms = [body_all_atom_positions(model, ref_inst, segs)
                  for segs in segments_list]

    # Box bounds in parent-voxel coords (frac * M), padded for the dilation shell
    # + a closure margin. min(sp) over-pads coarser axes, which is safe.
    pad = int(np.ceil(radius / min(sp))) + 2
    allv = np.array([[(lambda f: (f.x * M[0], f.y * M[1], f.z * M[2]))(
        cell.fractionalize(p))] for b in body_atoms for p in b]).reshape(-1, 3)
    NS = (np.floor(allv.min(0)) - pad).astype(int)
    NC = (np.ceil(allv.max(0)) + pad - NS + 1).astype(int)

    # Box-local cell preserves the parent per-voxel metric + angles, so a 2.5 A
    # set_points_around radius means the same physical distance as in the crystal.
    box_cell = gemmi.UnitCell(NC[0] * sp[0], NC[1] * sp[1], NC[2] * sp[2],
                              cell.alpha, cell.beta, cell.gamma)

    def to_box(pts):
        out = []
        for p in pts:
            f = cell.fractionalize(p)
            out.append(box_cell.orthogonalize(gemmi.Fractional(
                (f.x * M[0] - NS[0]) / NC[0], (f.y * M[1] - NS[1]) / NC[1],
                (f.z * M[2] - NS[2]) / NC[2])))
        return out

    grids, box_atoms = [], []
    for pts in body_atoms:
        g = gemmi.Int8Grid(int(NC[0]), int(NC[1]), int(NC[2]))
        g.set_unit_cell(box_cell)
        g.spacegroup = gemmi.SpaceGroup('P 1')
        bp = to_box(pts)
        for q in bp:
            g.set_points_around(q, radius=radius, value=1)
        grids.append(g)
        box_atoms.append(np.array([[q.x, q.y, q.z] for q in bp]) if bp
                         else np.empty((0, 3)))

    arrs = [np.array(g, copy=False) for g in grids]   # views: edits write through
    claimed = np.stack([a > 0 for a in arrs])         # (nbody, ni, nj, nk)
    overlap = claimed.sum(0) > 1
    n_contested = int(np.count_nonzero(overlap))
    if n_contested and len(grids) > 1:
        idx = np.argwhere(overlap)                    # (P, 3) box indices
        cellpos = np.array([list(grids[0].get_position(int(i), int(j), int(k)))
                            for i, j, k in idx])       # (P, 3) box-local orthogonal
        BIG = 1e30
        dists = np.full((len(grids), idx.shape[0]), BIG)
        contest = claimed[:, overlap]                  # (nbody, P)
        for b, ab in enumerate(box_atoms):
            if ab.shape[0] == 0 or not contest[b].any():
                continue
            d = np.sqrt(((cellpos[:, None, :] - ab[None, :, :]) ** 2)
                        .sum(-1)).min(1)               # (P,)
            dists[b] = np.where(contest[b], d, BIG)
        winner = dists.argmin(0)                       # (P,)
        for b in range(len(grids)):
            lose = (winner != b) & contest[b]
            for (i, j, k) in idx[lose]:
                arrs[b][int(i), int(j), int(k)] = 0
    return arrs, NS, M, sp, n_contested


def _write_subvolume_mask(arr, NS, M, sp, cell, path):
    """Extract the tight set-voxel sub-box of `arr` (+1 voxel closure margin) and
    write it as a mode-0, P1 CCP4 mask whose NSTART / MX-MY-MZ / cell place the
    contiguous domain copy at its REAL position in the parent grid.

    The header re-stamp happens AFTER update_ccp4_header (which resets NSTART to 0
    and MX/MY/MZ to the written block's own dims): we then set NSTART = the box
    origin in the parent grid, MX/MY/MZ = the parent sampling, and the crystal
    cell -- so gridding/origin match the map dm averages while NSTART + extent
    define the parallelepiped. Returns voxels set.
    """
    NSv = np.asarray(NS, dtype=int)
    sv = np.argwhere(arr > 0)
    if len(sv) == 0:
        tNS = NSv
        tarr = np.zeros((1, 1, 1), np.int8)
    else:
        dims = np.array(arr.shape)
        lo = np.maximum(sv.min(0) - 1, 0)             # 1-voxel zero margin for
        hi = np.minimum(sv.max(0) + 1, dims - 1)      # iso-surface closure
        tarr = np.ascontiguousarray(
            arr[lo[0]:hi[0] + 1, lo[1]:hi[1] + 1, lo[2]:hi[2] + 1]).astype(np.int8)
        tNS = NSv + lo
    tN = tarr.shape
    tcell = gemmi.UnitCell(tN[0] * sp[0], tN[1] * sp[1], tN[2] * sp[2],
                           cell.alpha, cell.beta, cell.gamma)
    tg = gemmi.Int8Grid(int(tN[0]), int(tN[1]), int(tN[2]))
    tg.set_unit_cell(tcell)
    tg.spacegroup = gemmi.SpaceGroup('P 1')
    np.array(tg, copy=False)[:] = tarr
    ccp4 = gemmi.Ccp4Mask()
    ccp4.grid = tg
    ccp4.update_ccp4_header()
    for w, v in [(5, int(tNS[0])), (6, int(tNS[1])), (7, int(tNS[2])),
                 (8, int(M[0])), (9, int(M[1])), (10, int(M[2]))]:
        ccp4.set_header_i32(w, v)
    for w, v in [(11, cell.a), (12, cell.b), (13, cell.c),
                 (14, cell.alpha), (15, cell.beta), (16, cell.gamma)]:
        ccp4.set_header_float(w, float(v))
    ccp4.write_ccp4_map(path)
    return int(np.count_nonzero(tarr))


def write_body_mask(model, ref_inst, segments, cell, spacegroup, path,
                    spacing=1.0, radius=2.5):
    """Write a CCP4 mode-0 mask covering the reference rigid body (union of its
    segments, possibly spanning several chains -- and possibly with a hole where
    a sub-range belongs to a different body, e.g. the C-helix excised from the
    CDK N-lobe).

    No symmetry expansion: the mask carries a single, CONTIGUOUS copy of the
    domain at its real position, as an NSTART sub-volume (see
    _write_subvolume_mask) -- the shape dm expects for an NCS averaging mask.
    Returns the number of grid points set.
    """
    arrs, NS, M, sp, _ = _build_disjoint_box(
        model, ref_inst, [segments], cell, spacegroup, spacing, radius)
    return _write_subvolume_mask(arrs[0], NS, M, sp, cell, path)


def write_domain_mask(model, ref, lo, hi, cell, spacegroup, path,
                      spacing=1.0, radius=2.5):
    """Back-compat single-chain shim: one implicit-role segment over one chain.
    See write_body_mask for the general case.
    """
    return write_body_mask(
        model, {_IMPLICIT_ROLE: ref}, [(_IMPLICIT_ROLE, lo, hi)],
        cell, spacegroup, path, spacing=spacing, radius=radius)


def write_partitioned_masks(model, ref_inst, segments_list, cell, spacegroup,
                            paths, spacing=1.0, radius=2.5):
    """Write DISJOINT averaging masks for several rigid bodies that share the
    reference instance.

    dm requires NCS averaging masks not to overlap: a grid point in two masks
    would be averaged under two operator sets. But masks built independently DO
    overlap -- the `radius` dilation bleeds across a shared boundary, so abutting
    bodies (the CDK N-lobe and the C-helix+cyclin body meeting at the helix, or
    AHIR's three stacked domains) contend for the boundary shell.

    Resolution is a nearest-atom competition: every contested grid point is kept
    only for the body whose reference atoms are closest, so the masks tile the
    molecule envelope without overlap. All bodies are masked on ONE shared,
    padded sub-volume box (so the competition runs in a single index frame, where
    contested cells -- being co-located by definition -- need no minimum-image
    correction); each body's disjoint result is then written as its own tight
    NSTART sub-volume placing a contiguous copy of the domain at its REAL
    position (no lattice shift, no PBC wrap). The masks and the operators
    (operator_ref_to_copy_body) thus share the crystal frame and move together.

    An NCS averaging mask is a shape that does not obey crystallographic symmetry
    (dm takes the crystal symmetry from the MTZ, not the mask), so each is written
    P1 -- only the NSTART/extent and parent sampling tie it back to the working
    map's grid.

    segments_list and paths are parallel (same order). Returns
    (nset_per_body, n_contested_points).
    """
    arrs, NS, M, sp, n_contested = _build_disjoint_box(
        model, ref_inst, segments_list, cell, spacegroup, spacing, radius)
    nsets = [_write_subvolume_mask(arr, NS, M, sp, cell, path)
             for arr, path in zip(arrs, paths)]
    return nsets, n_contested


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
