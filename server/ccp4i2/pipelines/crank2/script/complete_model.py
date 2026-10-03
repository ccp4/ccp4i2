"""A built model with its anomalous substructure, as one file.

Crank2 (and its SHELX route) hands on XYZOUT, the built model, and
XYZOUT_SUBSTR, the anomalous scatterers, separately; its own final REFMAC
refinement modelled both. Refining XYZOUT alone leaves the scatterers out of
the model, which matters as much as they scatter: for a mercury soak of a
10 kDa protein (HypF, two Hg sites) Servalcat went from R1-free 0.29 with the
sites to 0.44 without them, the protein distorting to explain the mercury.
Model building (ModelCraft) drops them the same way.

A site is not always a new atom. Where the scatterer is part of the model it
is already there: S-SAD sites sit on cysteine SG and methionine SD, and a
selenomethionine's Se sits where the model, built as methionine, has an SD.
So each site is
  - "on_model": left out, the model has an atom of the same element within
    SAME_ATOM of it (any symmetry copy);
  - "converted": for Se on a methionine SD, that residue becomes MSE, its
    SD replaced by SE;
  - "not_placed": an S or Se site that matches no model atom. These
    elements belong to residues, so a lone atom would be wrong: an unbuilt
    or misbuilt residue, a sulphate or chloride taken for S, or a
    misplaced site;
  - "clashes": left out, another model atom is within CLASH of it (no real
    scatterer sits that close to a protein atom);
  - "added": as a HETATM residue of its own in a new chain, with the
    occupancy and B factor the substructure gives.

Sites from the same job as the model share its origin and hand. Sites from
elsewhere (another job, another program) need not: with find_origin the
sites are first moved by whichever origin shift and hand of the model's
space group puts the most of them where scatterers can be (on a model atom
of their element, or 2-4.5 A from the model), if that is clearly better
than leaving them where they are.
"""
import math

import gemmi

SAME_ATOM = 2.0
CLASH = 1.5
NEAR = (2.0, 4.5)   # a bound scatterer sits this far from its nearest model atom
GRID = 12           # origin shifts tried on non-polar axes: multiples of 1/12
POLAR_STEP = 0.3    # A, the scan along a polar axis

# Monomer-library names for single-atom residues whose code is not the element
ION_RESIDUE = {"I": "IOD"}
# Elements that are part of residues: never added as lone atoms
RESIDUE_ELEMENTS = {"S", "Se"}


def _near(ns, model, cell, pos, radius):
    """(distance, cra) for every model atom within radius of pos, any image."""
    found = []
    for mark in ns.find_atoms(pos, "\0", radius=radius):
        cra = mark.to_cra(model)
        d = cell.find_nearest_image(pos, cra.atom.pos, gemmi.Asu.Any).dist()
        if d <= radius:
            found.append((d, cra))
    return sorted(found, key=lambda x: x[0])


def _where(cra):
    return f"{cra.chain.name}/{cra.residue.name}{cra.residue.seqid}/{cra.atom.name}"


def _classify(ns, model, cell, element, pos):
    """(kind, detail, cra or None, nearest distance or None) for one site."""
    near = _near(ns, model, cell, pos, max(SAME_ATOM, NEAR[1]))
    nearest = near[0][0] if near else None
    same = [(d, cra) for d, cra in near if d <= SAME_ATOM and cra.atom.element.name == element]
    if same:
        d, cra = same[0]
        return "on_model", f"{element} {_where(cra)} {d:.2f} A", cra, nearest
    met_sd = [(d, cra) for d, cra in near if d <= SAME_ATOM and element == "Se"
              and cra.residue.name == "MET" and cra.atom.name == "SD"]
    if met_sd:
        d, cra = met_sd[0]
        return "converted", f"{cra.chain.name}/MET{cra.residue.seqid}", cra, nearest
    clash = [(d, cra) for d, cra in near if d < CLASH]
    if clash:
        d, cra = clash[0]
        return "clashes", f"{element} {d:.2f} A from {_where(cra)}", cra, nearest
    if element in RESIDUE_ELEMENTS:
        where = f"nearest atom {nearest:.1f} A" if nearest is not None else "no model atom within 4.5 A"
        return "not_placed", f"{element} ({where})", None, nearest
    return "added", element, None, nearest


def _sites(sub):
    if len(sub) == 0:
        return []
    return [atom for chain in sub[0] for residue in chain for atom in residue]


def _score(ns, model, cell, sites, frac_positions):
    """(explained, on_model) for sites at these fractional positions."""
    explained = strong = 0
    for site, frac in zip(sites, frac_positions):
        kind, _, _, nearest = _classify(ns, model, cell, site.element.name, cell.orthogonalize(frac))
        if kind in ("on_model", "converted"):
            explained += 1
            strong += 1
        elif kind != "clashes" and nearest is not None and NEAR[0] <= nearest <= NEAR[1]:
            explained += 1
    return explained, strong


def _refine_polar(ns, model, cell, sites, fracs, hand, shift, axis):
    """The polar-axis shift, to 0.05 A, that best seats the sites on the
    model atoms of their element (the coarse scan only finds the band)."""
    length = cell.parameters[axis]
    best, best_error = shift[axis], None
    for k in range(-50, 51):
        trial = list(shift)
        trial[axis] = shift[axis] + k * 0.05 / length
        error = 0.0
        for site, f in zip(sites, fracs):
            pos = cell.orthogonalize(gemmi.Fractional(
                hand * f.x + trial[0], hand * f.y + trial[1], hand * f.z + trial[2]))
            near = [d for d, cra in _near(ns, model, cell, pos, SAME_ATOM)
                    if cra.atom.element.name == site.element.name
                    or (site.element.name == "Se" and cra.atom.name == "SD")]
            error += min(near + [SAME_ATOM]) ** 2
        if best_error is None or error < best_error:
            best, best_error = trial[axis], error
    return best % 1.0


def _ops(sg):
    """{(rotation, translation mod 1 in 24ths)} of every operation."""
    out = set()
    for op in sg.operations():
        out.add((tuple(tuple(r) for r in op.rot), tuple(t % 24 for t in op.tran)))
    return out


def _allowed(sg_sites, sg_model, hand, t24):
    """Does x -> hand*x + t map the sites' group onto the model's?"""
    target = _ops(sg_model)
    for op in sg_sites.operations():
        rot = [[v // 24 for v in row] for row in op.rot]
        shift = [t24[i] - sum(rot[i][j] * t24[j] for j in range(3)) for i in range(3)]
        tran = tuple((hand * op.tran[i] + shift[i]) % 24 for i in range(3))
        if (tuple(tuple(r) for r in op.rot), tran) not in target:
            return False
    return True


def _polar_axes(sg):
    axes = []
    for i in range(3):
        if all(op.rot[j][i] == (24 if j == i else 0) for op in sg.operations() for j in range(3)):
            axes.append(i)
    return axes


def fit_origin(st, sub):
    """Move the sites of ``sub`` to the model's origin and hand, in place.

    Returns what was found, for the report.
    """
    model = st[0]
    cell = st.cell
    sites = _sites(sub)
    sg_model = st.find_spacegroup()
    sg_sites = sub.find_spacegroup() or sg_model
    result = {"searched": False, "moved": False, "sites": len(sites)}
    if sg_model is None:
        result["note"] = "the model has no space group"
        return result
    ns = gemmi.NeighborSearch(model, cell, 5).populate()
    fracs = [cell.fractionalize(s.pos) for s in sites]
    before = _score(ns, model, cell, sites, fracs)
    result["explained_before"] = before[0]
    if len(sites) < 3:
        result["note"] = "fewer than 3 sites: too few to fix an origin, left as given"
        return result
    if before[0] >= math.ceil(0.6 * len(sites)):
        result["note"] = "the sites already fit the model"
        return result

    polar = _polar_axes(sg_model)
    centric = sg_model.is_centrosymmetric()
    hands = (1,) if centric else (1, -1)
    fixed = [i for i in range(3) if i not in polar]
    candidates = []
    steps = range(0, 24, 24 // GRID)
    for hand in hands:
        for a in (steps if 0 in fixed else (0,)):
            for b in (steps if 1 in fixed else (0,)):
                for c in (steps if 2 in fixed else (0,)):
                    t24 = (a, b, c)
                    if _allowed(sg_sites, sg_model, hand, t24):
                        candidates.append((hand, [v / 24 for v in t24]))
    scan = [0.0]
    if len(polar) == 1:
        length = cell.parameters[polar[0]]
        scan = [k * POLAR_STEP / length for k in range(int(length / POLAR_STEP))]
    elif len(polar) > 1:
        result["note"] = "more than one polar axis: the origin along them was not searched"

    # Every placement scored; a site on a model atom of its element counts
    # three times one merely near the model (2-4.5 A), which a wrong origin
    # often manages by chance (Xe in a cavity, half the sites near protein).
    scored = []
    for hand, t in candidates:
        for p in scan:
            shift = list(t)
            if len(polar) == 1:
                shift[polar[0]] = p
            moved = [gemmi.Fractional(hand * f.x + shift[0], hand * f.y + shift[1], hand * f.z + shift[2])
                     for f in fracs]
            explained, strong = _score(ns, model, cell, sites, moved)
            scored.append((explained + 2 * strong, explained, hand, shift))
    result.update(searched=True, tried=len(scored))
    scored.sort(key=lambda s: -s[0])
    key, explained, hand, shift = scored[0]
    result["explained_after"] = explained

    def same_placement(other):
        # the scan along a polar axis finds one placement over a band of
        # steps: a site counts as on its atom up to SAME_ATOM either side
        _, _, h, s = other
        if h != hand:
            return False
        for i in range(3):
            d = abs(s[i] - shift[i]) % 1.0
            d = min(d, 1.0 - d) * cell.parameters[i]
            if d > (2 * SAME_ATOM + 0.5 if i in polar else 0.01):
                return False
        return True

    runner_up = next((s for s in scored[1:] if not same_placement(s)), None)
    margin = key - runner_up[0] if runner_up else key
    identity = (hand, shift) == (1, [0.0, 0.0, 0.0])
    if explained < max(3, math.ceil(0.6 * len(sites))):
        result["note"] = ("no origin or hand of the model's space group puts the sites "
                          "where scatterers can be: are they from this structure?")
    elif identity:
        result["note"] = "the sites fit best where they are"
    elif explained < before[0] + 2:
        result["note"] = ("no origin or hand of the model's space group fits the sites "
                          "clearly better: are they from this structure?")
    elif margin < 2:
        result["note"] = ("several origins or hands fit the sites about equally well, so "
                          "they were left as given; sites that sit on model atoms (S, Se) "
                          "or more of them would decide it")
    else:
        if len(polar) == 1:
            shift[polar[0]] = _refine_polar(ns, model, cell, sites, fracs, hand, shift, polar[0])
        for site, f in zip(sites, fracs):
            site.pos = cell.orthogonalize(gemmi.Fractional(
                hand * f.x + shift[0], hand * f.y + shift[1], hand * f.z + shift[2]))
        result["moved"] = True
        sign = "" if hand == 1 else "-"
        result["transform"] = ", ".join(f"{sign}{axis}{'+' if s >= 0 else ''}{s:.4f}"
                                        for axis, s in zip("xyz", shift))
        result["note"] = ("sites moved to the model's origin" +
                          (" and hand (inverted)" if hand == -1 else ""))
    return result


def _remove_clashing(model, cell, sites, report):
    """Remove the residues a heavy scatterer clashes with (see complete_model)."""
    ns = gemmi.NeighborSearch(model, cell, 5).populate()
    doomed = {}
    for site in sites:
        element = site.element.name
        if element in RESIDUE_ELEMENTS:
            continue
        near = [(d, cra) for d, cra in _near(ns, model, cell, site.pos, CLASH) if d < CLASH]
        if not near or any(cra.atom.element.name in RESIDUE_ELEMENTS or
                           cra.atom.element.name == element for _, cra in near):
            continue
        for _, cra in near:
            doomed[(cra.chain.name, str(cra.residue.seqid))] = (
                f"{cra.chain.name}/{cra.residue.name}{cra.residue.seqid} (under {element})")
    for chain in model:
        for i in reversed(range(len(chain))):
            key = (chain.name, str(chain[i].seqid))
            if key in doomed:
                report["removed"].append(doomed[key])
                del chain[i]
    for i in reversed(range(len(model))):
        if len(model[i]) == 0:
            del model[i]


def _free_chain_name(model):
    used = {ch.name for ch in model}
    for name in "WXYZHIJKLMNOPQRSTUVABCDEFG" + "abcdefghijklmnopqrstuvwxyz0123456789":
        if name not in used:
            return name
    raise ValueError("no free single-character chain name")


def _to_mse(residue):
    residue.name = "MSE"
    residue.het_flag = "H"
    for atom in residue:
        if atom.name == "SD":
            atom.name = "SE"
            atom.element = gemmi.Element("Se")


def complete_model(model_path, substr_path, out_path, all_met_to_mse=False, find_origin=False,
                   remove_clashing=False):
    """Write the model plus its substructure to out_path; return what was done.

    The result is a dict of lists ("added", "on_model", "converted",
    "not_placed", "clashes", "removed"), "all_mse" (methionines made MSE by
    all_met_to_mse) and, with find_origin, "origin" (fit_origin's account).

    remove_clashing: a heavy scatterer (not S or Se) with model atoms within
    CLASH of it, none of them S or Se, is taken to be right and the model
    wrong there: model building traced chain into the heavy-atom peak (HypF:
    ModelCraft put Cys67 N 0.31 A and an UNK atom 0.13 A from the two Hg;
    removing those residues and adding the Hg took R-free 0.409 to 0.297).
    Those residues are removed ("removed") and the site added. A site on an
    S or Se atom is that atom taken for the soaked element, so it is left.
    """
    st = gemmi.read_structure(str(model_path))
    sub = gemmi.read_structure(str(substr_path))
    if len(st) == 0 or len(st[0]) == 0:
        raise ValueError(f"{model_path}: no model to complete")
    if not st.cell.is_crystal():
        st.cell = sub.cell
    if not st.spacegroup_hm:
        st.spacegroup_hm = sub.spacegroup_hm
    st.setup_entities()
    model = st[0]
    cell = st.cell

    report = {"added": [], "on_model": [], "converted": [], "not_placed": [], "clashes": [],
              "removed": [], "all_mse": 0}
    if find_origin:
        report["origin"] = fit_origin(st, sub)
    origin = report.get("origin") or {}
    # Searched and could not put the sites in the model's frame: adding them
    # (or removing residues for them) would put heavy atoms where none are
    unplaced = bool(origin.get("searched") and not origin.get("moved")
                    and origin.get("note") != "the sites fit best where they are")
    if remove_clashing and not unplaced:
        _remove_clashing(model, cell, _sites(sub), report)
    ns = gemmi.NeighborSearch(model, cell, 5).populate()
    additions, to_mse = [], []
    for site in _sites(sub):
        kind, detail, cra, _ = _classify(ns, model, cell, site.element.name, site.pos)
        if unplaced and kind in ("added", "converted"):
            report["not_placed"].append(f"{site.element.name} (frame not established)")
            continue
        if kind == "added":
            additions.append(site)
        elif kind == "converted":
            to_mse.append((cra.chain.name, str(cra.residue.seqid)))
        else:
            report[kind].append(detail)

    for chain_name, seqid in dict.fromkeys(to_mse):
        for residue in model[chain_name]:
            if str(residue.seqid) == seqid and residue.name == "MET":
                _to_mse(residue)
                report["converted"].append(f"{chain_name}/MET{seqid}")
    if all_met_to_mse:
        for chain in model:
            for residue in chain:
                if residue.name == "MET":
                    _to_mse(residue)
                    report["all_mse"] += 1

    if additions:
        chain = gemmi.Chain(_free_chain_name(model))
        for number, site in enumerate(additions, start=1):
            residue = gemmi.Residue()
            residue.name = ION_RESIDUE.get(site.element.name.upper(), site.element.name.upper())
            residue.seqid = gemmi.SeqId(number, " ")
            residue.het_flag = "H"
            atom = site.clone()
            atom.name = site.element.name.upper()
            residue.add_atom(atom)
            chain.add_residue(residue)
            report["added"].append(f"{residue.name} occupancy {site.occ:.2f}")
        model.add_chain(chain)

    st.setup_entities()
    if str(out_path).lower().endswith((".cif", ".mmcif")):
        st.make_mmcif_document().write_file(str(out_path))
    else:
        st.write_pdb(str(out_path))
    return report


KINDS = ("added", "on_model", "converted", "not_placed", "clashes", "removed")


def report_element(report):
    """complete_model's account as a CompleteModel element, for program.xml:
    a count per kind as attributes, each item as a child, the origin search
    as an Origin child."""
    from lxml import etree

    element = etree.Element("CompleteModel")
    for kind in KINDS:
        element.set(kind, str(len(report[kind])))
        for text in report[kind]:
            etree.SubElement(element, kind).text = text
    element.set("all_mse", str(report.get("all_mse", 0)))
    origin = report.get("origin")
    if origin:
        node = etree.SubElement(element, "Origin")
        for key, value in origin.items():
            node.set(key, str(value))
    return element
