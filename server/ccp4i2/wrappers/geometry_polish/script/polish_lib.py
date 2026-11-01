"""NCS-aware geometry polish: the library behind the geometry_polish task.

Three tools, three roles:

* MolProbity (molprobity.rotalyze / ramalyze / clashscore) is the judge. It is
  what the wwPDB validation report runs, so a residue it no longer flags is a
  residue the report will not flag either.
* gemmi reads the maps and the model, scores density at atoms, finds contacts
  and transplants side chains or backbone windows between NCS copies.
* Coot's headless API ("chapi") does the building: auto-fit rotamer, library
  rotamers, peptide flips and torsion- or Ramachandran-restrained real-space
  refinement.

A candidate change moves a small window of residues (Coot's TRIPLE or
QUINTUPLE refinement moves the neighbours too). It is accepted only when,
over everything that moved, MolProbity stops flagging the target without
flagging anything new, no new contact or clash appears, and the 2mFo-DFc
density at both side chains and backbone is kept (see :func:`density_ok`).
What cannot be fixed that way is listed for a human rather than forced.

Residues are identified by (chain, number, insertion code) throughout, and
windows are taken by position in the chain, so insertion codes and numbering
gaps cannot make two residues collide.

Only the standard library is imported at module level: gemmi, numpy and
coot_headless_api are execution-time dependencies and every task plugin must
import with none of them present.
"""
import csv
import os
import subprocess
import time

BB_ONLY = ("N", "CA", "C", "O")
# Atoms a side-chain transplant leaves alone: the backbone and the CB it is
# superposed on. For density the CB counts as side chain, so an alanine has
# something to score.
TRANSPLANT_KEEP = {"N", "CA", "C", "O", "CB"}
# Side-chain atoms for the clash stage: anything not on this list, so HA and
# the amide H count as backbone and a GLY never has a side chain to move.
BACKBONE_H = {"N", "CA", "C", "O", "H", "HA", "HA2", "HA3"}
STANDARD_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
    "LYS", "MET", "MSE", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HOH",
}
MOLPROBITY_TOOLS = ("molprobity.rotalyze", "molprobity.ramalyze", "molprobity.clashscore")

ROTAMER_LIBRARY_TRIALS = 8          # library rotamers tried for a strained side chain
NCS_DONORS = 3                      # best NCS copies tried before anything else
ABSOLUTE_SLACK = 0.3                # sigma a change may lose outright (above MIN_DENSITY)
CONTACT_DISTANCE = 2.4              # A; a new heavy-atom contact closer than this vetoes a change
TORSION_WEIGHTS = (5.0, 10.0, 20.0)
RAMA_WEIGHT = 10.0
RAMA_MARGIN = 2                     # Ramachandran is scored this many residues beyond what moved
NEIGHBOURHOOD_RADIUS = 6.0          # A; residues this close to a moved one are re-scored for clashes


class PolishError(RuntimeError):
    """A tool the polish depends on is missing or failed outright."""


# --------------------------------------------------------------------------
# Residue keys
# --------------------------------------------------------------------------

def rkey(chain, num, icode=""):
    """The (chain, number, insertion code) key of a residue; blank icode is ''."""
    return (str(chain), int(num), (icode or "").strip())


def residue_key(chain_name, residue):
    return rkey(chain_name, residue.seqid.num, residue.seqid.icode)


def key_label(key, resname=""):
    """'SER A 100A' style label for logs and reports."""
    chain, num, icode = key
    return "{} {} {}{}".format(resname, chain, num, icode).strip()


def parse_residue_field(field):
    """MolProbity's fixed-width residue field: `` A1005A LYS``, ``AB  11  SER``.

    Columns: 0-1 chain (right-justified), 2-5 residue number, 6 insertion
    code, 8-10 residue name. Whitespace splitting breaks on 4-digit numbers
    and on 2-character chains, which is why this is positional.
    """
    if len(field) < 11:
        return None
    try:
        num = int(field[2:6])
    except ValueError:
        return None
    return rkey(field[0:2].strip(), num, field[6:7]), field[8:11].strip()


# --------------------------------------------------------------------------
# MolProbity
# --------------------------------------------------------------------------

def program(name):
    """Resolve a MolProbity tool through the program-location preferences."""
    try:
        from ccp4i2.config.program_discovery import resolve_program
        found = resolve_program(name)
    except Exception:  # pragma: no cover - discovery is optional outside i2
        found = None
    return found or name


def run_tool(name, pdb_path):
    """Run one MolProbity tool on a PDB file and return its stdout."""
    try:
        done = subprocess.run([program(name), pdb_path], capture_output=True, text=True)
    except OSError as exc:
        raise PolishError("cannot run {}: {}".format(name, exc)) from exc
    if done.returncode != 0 and not done.stdout.strip():
        raise PolishError("{} failed on {}: {}".format(name, os.path.basename(pdb_path),
                                                      done.stderr.strip()[-400:]))
    return done.stdout


def parse_rotalyze(text):
    """{key: {score, evaluation, rotamer, resname}} from rotalyze output.

    Lines are ``residue:occupancy:score%:chi1:chi2:chi3:chi4:evaluation:rotamer``.
    """
    result = {}
    for line in text.splitlines():
        if ":" not in line or line.startswith(("residue", "SUMMARY")):
            continue
        parts = line.split(":")
        if len(parts) < 9:
            continue
        parsed = parse_residue_field(parts[0])
        if parsed is None:
            continue
        try:
            score = float(parts[2])
        except ValueError:
            continue
        key, resname = parsed
        result[key] = dict(score=score, evaluation=parts[7].strip(), rotamer=parts[8].strip(),
                           resname=resname)
    return result


def parse_ramalyze(text):
    """{key: {score, evaluation, type, resname}} from ramalyze output.

    Lines are ``residue:score%:phi:psi:evaluation:type``.
    """
    result = {}
    for line in text.splitlines():
        if ":" not in line or line.startswith(("residue", "SUMMARY")):
            continue
        parts = line.split(":")
        if len(parts) < 6:
            continue
        parsed = parse_residue_field(parts[0])
        if parsed is None:
            continue
        try:
            score = float(parts[1])
        except ValueError:
            continue
        key, resname = parsed
        result[key] = dict(score=score, evaluation=parts[4].strip(), type=parts[5].strip(),
                           resname=resname)
    return result


def parse_clash_line(line):
    """One clashscore contact line, positionally.

    `` C1106  PHE  HB3  C1127A SER  HB3 :1.068`` and
    ``CX 116  PHE  HB3 CX 137A SER  HB3 :1.068``: residue field (11), space,
    atom (4), space, residue field (11), space, atom (4), `` :overlap``.
    """
    if len(line) < 35 or ":" not in line:
        return None
    first = parse_residue_field(line[0:11])
    second = parse_residue_field(line[17:28])
    if first is None or second is None:
        return None
    try:
        overlap = float(line.split(":")[-1])
    except ValueError:
        return None
    (key1, resname1), (key2, resname2) = first, second
    return dict(key1=key1, resname1=resname1, atom1=line[11:16].strip(),
                key2=key2, resname2=resname2, atom2=line[28:33].strip(), overlap=overlap)


def parse_clashscore(text):
    """(clashscore or None, [clash dicts]) from clashscore output."""
    score, clashes = None, []
    for line in text.splitlines():
        if line.startswith("clashscore"):
            try:
                score = float(line.split("=")[1])
            except (IndexError, ValueError):
                pass
        elif ":" in line and not line.startswith(("Bad", "Using", "Adding")):
            clash = parse_clash_line(line)
            if clash:
                clashes.append(clash)
    return score, clashes


def clash_pair(clash):
    """The unordered residue pair a clash belongs to."""
    return tuple(sorted([clash["key1"], clash["key2"]]))


def clash_id(clash):
    """The unordered atom pair of a clash, for diffing two clash lists."""
    return tuple(sorted([(clash["key1"], clash["atom1"]), (clash["key2"], clash["atom2"])]))


def new_clashes(before, after, threshold):
    """Clashes at or above ``threshold`` in ``after`` that ``before`` did not have."""
    known = {clash_id(c) for c in before}
    return [c for c in after if c["overlap"] >= threshold and clash_id(c) not in known]


def rotalyze(pdb_path):
    return parse_rotalyze(run_tool("molprobity.rotalyze", pdb_path))


def ramalyze(pdb_path):
    return parse_ramalyze(run_tool("molprobity.ramalyze", pdb_path))


def clashscore(pdb_path):
    return parse_clashscore(run_tool("molprobity.clashscore", pdb_path))


def summarize(pdb_path):
    """The headline MolProbity numbers for a model, as the validation report counts them."""
    rot = rotalyze(pdb_path)
    ram = ramalyze(pdb_path)
    n_rot = sum(1 for v in rot.values() if v["evaluation"] == "OUTLIER")
    n_ram = sum(1 for v in ram.values() if v["evaluation"] == "OUTLIER")
    n_fav = sum(1 for v in ram.values() if v["evaluation"] == "Favored")
    score, clashes = clashscore(pdb_path)
    return dict(rot_outliers=n_rot, rot_pct=round(100.0 * n_rot / max(1, len(rot)), 2),
                rama_outliers=n_ram, rama_pct=round(100.0 * n_ram / max(1, len(ram)), 2),
                rama_favored_pct=round(100.0 * n_fav / max(1, len(ram)), 2),
                clashscore=score, n_clashes=len(clashes))


# --------------------------------------------------------------------------
# Decisions
# --------------------------------------------------------------------------

def density_ok(before, after, guard_fraction, min_density, slack=ABSOLUTE_SLACK):
    """Does a change keep enough of the density it started with?

    The result must stay above ``min_density`` and either within ``slack``
    sigma of the original or above ``guard_fraction`` of it. The fraction is
    what let strained side chains in strong density (4 sigma) move to a
    rotamer at 2.5 sigma; the floor stops a side chain at 1 sigma dropping
    below 0.8, and stops the slack excusing a fall from 0.5 to 0.2.
    """
    if after != after or before != before:  # NaN: no atoms to score
        return True
    if after < min_density:
        return False
    return after >= before - slack or after >= guard_fraction * before


# --------------------------------------------------------------------------
# gemmi: maps, density, contacts, transplants, windows
# --------------------------------------------------------------------------

_map_cache = {}


def map_labels(labels):
    """(amplitude, phase) column labels of a map-coefficient MTZ, or Nones."""
    f_label = "F" if "F" in labels else next((l for l in labels if l.startswith("F")), None)
    phi_label = "PHI" if "PHI" in labels else next((l for l in labels if l.startswith("PH")), None)
    return f_label, phi_label


def mtz_labels(mtz_path):
    import gemmi
    labels = [c.label for c in gemmi.read_mtz_file(mtz_path).columns]
    f_label, phi_label = map_labels(labels)
    if f_label is None or phi_label is None:
        raise PolishError("{} has no F,PHI columns (found {})".format(mtz_path, ",".join(labels)))
    return f_label, phi_label


def load_map(mtz_path, sample_rate=3.0):
    """Sigma-scaled map from an F,PHI mini-MTZ (cached per path)."""
    import gemmi
    if mtz_path in _map_cache:
        return _map_cache[mtz_path]
    f_label, phi_label = mtz_labels(mtz_path)
    grid = gemmi.read_mtz_file(mtz_path).transform_f_phi_to_map(f_label, phi_label, sample_rate=sample_rate)
    grid.normalize()
    _map_cache[mtz_path] = grid
    return grid


def read_structure(path):
    import gemmi
    st = gemmi.read_structure(path)
    st.setup_entities()
    return st


def write_structure(st, path):
    """Write PDB or mmCIF according to the extension."""
    if path.lower().endswith((".cif", ".mmcif")):
        st.make_mmcif_document().write_file(path)
    else:
        st.write_pdb(path)


def chain_residues(st, chain_id):
    chain = st[0].find_chain(chain_id)
    return list(chain) if chain is not None else None


def locate(st, key):
    """Index of a residue in its chain, or None."""
    residues = chain_residues(st, key[0])
    if residues is None:
        return None
    for index, res in enumerate(residues):
        if residue_key(key[0], res) == key:
            return index
    return None


def residue_at(st, chain_id, index):
    residues = chain_residues(st, chain_id)
    if residues is None or index < 0 or index >= len(residues):
        return None
    return residues[index]


def find_residue(st, chain_id, resno, icode=""):
    index = locate(st, rkey(chain_id, resno, icode))
    return None if index is None else residue_at(st, chain_id, index)


def find_by_key(st, key):
    return find_residue(st, key[0], key[1], key[2])


def window(st, chain_id, centre, before, after):
    """Chain positions centre-before..centre+after, clamped to the chain."""
    residues = chain_residues(st, chain_id)
    if residues is None:
        return []
    return list(range(max(0, centre - before), min(len(residues), centre + after + 1)))


def keys_at(st, chain_id, indices):
    return [residue_key(chain_id, residue_at(st, chain_id, i)) for i in indices]


def polymer_chains(st):
    return [ch.name for ch in st[0] if ch.get_polymer()]


def is_amino_acid(resname):
    import gemmi
    info = gemmi.find_tabulated_residue(resname)
    return bool(info and info.is_amino_acid())


def has_altlocs(residue):
    return any(atom.altloc != "\0" for atom in residue)


def residue_density(grid, residue, keep=lambda atom: True):
    """Mean sigma-scaled density at the centres of the atoms ``keep`` selects."""
    values = [grid.interpolate_value(a.pos) for a in residue if keep(a)]
    if not values:
        return float("nan")
    return float(sum(values) / len(values))


def sidechain_density(grid, residue):
    return residue_density(grid, residue, lambda a: a.name not in BB_ONLY and a.element.name != "H")


def backbone_density(grid, residue):
    return residue_density(grid, residue, lambda a: a.name in BB_ONLY)


def clash_contacts(st, key, dist=CONTACT_DISTANCE):
    """Heavy-atom contacts closer than ``dist`` involving one residue.

    Symmetry aware; adjacent residues ignored. Returns a set of
    (other key, atom, other atom) so two models can be compared.
    """
    import gemmi
    ns = gemmi.NeighborSearch(st[0], st.cell, 5.0).populate(include_h=False)
    cs = gemmi.ContactSearch(dist)
    cs.ignore = gemmi.ContactSearch.Ignore.AdjacentResidues
    found = set()
    for contact in cs.find_contacts(ns):
        for a, b in ((contact.partner1, contact.partner2), (contact.partner2, contact.partner1)):
            if residue_key(a.chain.name, a.residue) == key:
                found.add((residue_key(b.chain.name, b.residue), a.atom.name, b.atom.name))
    return found


def copy_positions(src, dst):
    """Move every atom of ``dst`` onto ``src``'s position, adding atoms ``dst``
    lacks (with ``dst``'s own occupancy and B)."""
    import gemmi
    for atom in src:
        target = dst.find_atom(atom.name, "*")
        if target is not None:
            target.pos = atom.pos
        elif atom.element.name != "H":
            new_atom = gemmi.Atom()
            new_atom.name = atom.name
            new_atom.element = atom.element
            new_atom.pos = atom.pos
            anchor = dst.find_atom("CB", "*") or dst.find_atom("CA", "*")
            new_atom.occ = anchor.occ if anchor is not None else atom.occ
            new_atom.b_iso = anchor.b_iso if anchor is not None else atom.b_iso
            dst.add_atom(new_atom)


def copy_window(src_st, dst_st, chain_id, indices):
    """Bring every atom of the given chain positions (matched by key) from one model into another."""
    for index in indices:
        dst = residue_at(dst_st, chain_id, index)
        if dst is None:
            continue
        src = find_by_key(src_st, residue_key(chain_id, dst))
        if src is not None and src.name == dst.name:
            copy_positions(src, dst)


def _paired(st, chain_id, index, donor_chain):
    """(target residue, donor residue with the same number, icode and name) or None."""
    target = residue_at(st, chain_id, index)
    if target is None:
        return None
    key = residue_key(chain_id, target)
    donor = find_residue(st, donor_chain, key[1], key[2])
    if donor is None or donor.name != target.name:
        return None
    return target, donor


def transplant_sidechain(st, chain_id, resno, donor_chain, icode="", flank=1):
    """Copy the donor copy's side chain onto the target after superposing the
    local backbone (N, CA, C, O, CB of the residue and ``flank`` neighbours).

    Keeps the target's occupancies. Refuses residues with alternate
    conformations on either side, since a transplant would mix them.
    Returns the superposition rmsd, or None when the copies cannot be paired.
    """
    import gemmi
    index = locate(st, rkey(chain_id, resno, icode))
    if index is None:
        return None
    pair = _paired(st, chain_id, index, donor_chain)
    if pair is None:
        return None
    target, donor = pair
    if has_altlocs(target) or has_altlocs(donor):
        return None
    target_pos, donor_pos = [], []
    for i in window(st, chain_id, index, flank, flank):
        neighbours = _paired(st, chain_id, i, donor_chain)
        if neighbours is None:
            continue
        for name in TRANSPLANT_KEEP:
            ta, da = neighbours[0].find_atom(name, "*"), neighbours[1].find_atom(name, "*")
            if ta is not None and da is not None:
                target_pos.append(ta.pos)
                donor_pos.append(da.pos)
    if len(target_pos) < 4:
        return None
    sup = gemmi.superpose_positions(target_pos, donor_pos)
    anchor = target.find_atom("CB", "*") or target.find_atom("CA", "*")
    for da in donor:
        if da.name in TRANSPLANT_KEEP or da.element.name == "H":
            continue
        new_pos = gemmi.Position(sup.transform.apply(da.pos))
        ta = target.find_atom(da.name, "*")
        if ta is not None:
            ta.pos = new_pos
        else:
            atom = gemmi.Atom()
            atom.name, atom.element, atom.pos = da.name, da.element, new_pos
            atom.occ = anchor.occ if anchor is not None else 1.0
            atom.b_iso = anchor.b_iso if anchor is not None else da.b_iso
            target.add_atom(atom)
    return sup.rmsd


def transplant_backbone_indices(st, chain_id, lo, hi, donor_chain, flank=2):
    """Replace all atoms of chain positions lo..hi with the donor copy's,
    superposed on the backbone of ``flank`` residues either side.
    Returns rmsd or None."""
    import gemmi
    residues = chain_residues(st, chain_id)
    if residues is None:
        return None
    target_pos, donor_pos = [], []
    flanks = [i for i in range(lo - flank, lo) if i >= 0] + \
             [i for i in range(hi + 1, hi + flank + 1) if i < len(residues)]
    for i in flanks:
        pair = _paired(st, chain_id, i, donor_chain)
        if pair is None:
            continue
        for name in BB_ONLY:
            ta, da = pair[0].find_atom(name, "*"), pair[1].find_atom(name, "*")
            if ta is not None and da is not None:
                target_pos.append(ta.pos)
                donor_pos.append(da.pos)
    if len(target_pos) < 6:
        return None
    pairs = []
    for i in range(lo, hi + 1):
        pair = _paired(st, chain_id, i, donor_chain)
        if pair is None or has_altlocs(pair[0]) or has_altlocs(pair[1]):
            return None
        pairs.append(pair)
    sup = gemmi.superpose_positions(target_pos, donor_pos)
    for target, donor in pairs:
        for da in donor:
            if da.element.name == "H":
                continue
            ta = target.find_atom(da.name, "*")
            if ta is not None:
                ta.pos = gemmi.Position(sup.transform.apply(da.pos))
    return sup.rmsd


def transplant_backbone_window(st, chain_id, start, end, donor_chain, flank=2):
    """Number-based wrapper of :func:`transplant_backbone_indices`."""
    lo, hi = locate(st, rkey(chain_id, start)), locate(st, rkey(chain_id, end))
    if lo is None or hi is None:
        return None
    return transplant_backbone_indices(st, chain_id, lo, hi, donor_chain, flank)


def _mini_structure(st):
    import gemmi
    mini = gemmi.Structure()
    mini.cell = st.cell
    mini.spacegroup_hm = st.spacegroup_hm
    return mini


def write_index_window(st, chain_id, lo, hi, path):
    """A PDB holding chain positions lo..hi of one chain: enough for rotalyze/ramalyze."""
    import gemmi
    mini = _mini_structure(st)
    model = gemmi.Model("1")
    chain = gemmi.Chain(chain_id)
    residues = chain_residues(st, chain_id) or []
    for index in range(max(0, lo), min(len(residues), hi + 1)):
        chain.add_residue(residues[index])
    model.add_chain(chain)
    mini.add_model(model)
    mini.setup_entities()
    mini.write_pdb(path)
    return path


def write_window(st, chain_id, lo, hi, path):
    """A PDB holding residues numbered lo..hi of one chain (every insertion code)."""
    import gemmi
    mini = _mini_structure(st)
    model = gemmi.Model("1")
    chain = gemmi.Chain(chain_id)
    for res in chain_residues(st, chain_id) or []:
        if lo <= res.seqid.num <= hi:
            chain.add_residue(res)
    model.add_chain(chain)
    mini.add_model(model)
    mini.setup_entities()
    mini.write_pdb(path)
    return path


def neighbourhood_keys(st, keys, radius=NEIGHBOURHOOD_RADIUS):
    """Keys of every residue with an atom within ``radius`` of the given residues."""
    import gemmi
    ns = gemmi.NeighborSearch(st[0], st.cell, max(5.0, radius)).populate(include_h=False)
    keep = set()
    for key in keys:
        res = find_by_key(st, key)
        if res is None:
            continue
        keep.add(key)
        for atom in res:
            for mark in ns.find_atoms(atom.pos, "\0", radius=radius):
                cra = mark.to_cra(st[0])
                keep.add(residue_key(cra.chain.name, cra.residue))
    return keep


def write_neighbourhood(st, keys, path, radius=NEIGHBOURHOOD_RADIUS, selection=None):
    """A PDB holding the residues around ``keys`` (or exactly ``selection`` when
    given, so two models can be scored on the same residue set)."""
    import gemmi
    keep = selection if selection is not None else neighbourhood_keys(st, keys, radius)
    mini = _mini_structure(st)
    model = gemmi.Model("1")
    for source in st[0]:
        chain = gemmi.Chain(source.name)
        for res in source:
            if residue_key(source.name, res) in keep:
                chain.add_residue(res)
        if len(chain) > 0:
            model.add_chain(chain)
    mini.add_model(model)
    mini.setup_entities()
    mini.write_pdb(path)
    return path


def nonstandard_residues(st):
    return sorted({res.name for ch in st[0] for res in ch if res.name not in STANDARD_RESIDUES})


# --------------------------------------------------------------------------
# Coot headless API
# --------------------------------------------------------------------------

def coot_environment():
    """Point Coot at CCP4's data if the environment does not already."""
    ccp4 = os.environ.get("CCP4")
    if not ccp4:
        return
    share = os.path.join(ccp4, "share", "coot")
    os.environ.setdefault("COOT_DATA_DIR", share)
    os.environ.setdefault("SYMINFO", os.path.join(share, "data", "syminfo.lib"))
    os.environ.setdefault("COOT_REF_STRUCTS", os.path.join(share, "reference-structures"))
    os.environ.setdefault("COOT_STANDARD_RESIDUES", os.path.join(share, "standard-residues.pdb"))


def monomer_dictionary(resname):
    """CCP4 monomer-library file for a residue name, if the library is on hand."""
    lib = os.environ.get("CLIBD_MON")
    if not lib or not resname:
        return None
    path = os.path.join(lib, resname[0].lower(), resname + ".cif")
    return path if os.path.exists(path) else None


def coot_cid(key):
    """MMDB selection for one residue: ``//A/100`` or ``//A/100.A`` with an icode."""
    chain, num, icode = key
    return "//{}/{}{}".format(chain, num, "." + icode if icode else "")


class CootSession:
    """One Coot molecules container holding the model and the refinement map.

    A fresh session per trial is deliberate: Coot keeps state (restraints
    weights, backups) and the trial must start from the current model.
    """

    def __init__(self, pdb_path, map_mtz, dictionaries=()):
        coot_environment()
        from ccp4i2.lib.coot_api import molecules_container
        mc = molecules_container(False)
        mc.set_make_backups(False)
        mc.set_use_gemmi(False)
        for path in dictionaries:
            mc.import_cif_dictionary(path, mc.get_imol_enc_any())
        reader = getattr(mc, "read_coordinates", None) or mc.read_pdb
        self.mc = mc
        self.imol = reader(pdb_path)
        if self.imol < 0:
            raise PolishError("Coot could not read {}".format(pdb_path))
        f_label, phi_label = mtz_labels(map_mtz)
        self.imap = mc.read_mtz(map_mtz, f_label, phi_label, "", False, False)
        if self.imap < 0:
            raise PolishError("Coot could not read map coefficients {} as {},{}".format(
                map_mtz, f_label, phi_label))
        mc.set_imol_refinement_map(self.imap)

    def auto_fit(self, key):
        return self.mc.auto_fit_rotamer(self.imol, key[0], key[1], key[2], "", self.imap)

    def set_rotamer(self, key, number):
        return bool(self.mc.set_residue_to_rotamer_number(self.imol, coot_cid(key), "", number))

    def use_torsion_restraints(self, weight):
        self.mc.set_use_torsion_restraints(True)
        self.mc.set_torsion_restraints_weight(weight)

    def use_rama_restraints(self, weight):
        self.mc.set_use_rama_plot_restraints(True)
        self.mc.set_rama_plot_restraints_weight(weight)

    def refine(self, key, mode, cycles):
        """Real-space refine around a residue: mode TRIPLE (n-1..n+1) or QUINTUPLE (n-2..n+2)."""
        return self.mc.refine_residues(self.imol, key[0], key[1], key[2], "", mode, cycles)

    def flip_peptide(self, key):
        return bool(self.mc.flip_peptide_using_cid(self.imol, coot_cid(key), ""))

    def write(self, path):
        self.mc.write_coordinates(self.imol, path)
        return path


# --------------------------------------------------------------------------
# The polish itself
# --------------------------------------------------------------------------

class Options:
    """Knobs of the polish, with the defaults that worked."""

    def __init__(self, **kw):
        self.do_rotamers = True
        self.do_rama = True
        self.do_flips = True
        self.do_clashes = True
        self.use_ncs = True
        self.guard = 0.6
        self.min_density = 0.8
        self.cycles = 200
        self.max_clash_rotamers = 6
        self.min_overlap = 0.8
        self.chains = None          # None = every polymer chain
        for key, value in kw.items():
            if not hasattr(self, key):
                raise TypeError("unknown option {}".format(key))
            setattr(self, key, value)


class Verdict:
    """What the judge found for one candidate."""

    def __init__(self):
        self.ok = True
        self.reasons = []
        self.favoured = 0
        self.mean_sidechain = float("nan")
        self.target_rot = None
        self.target_rama = None
        self.clashes = []

    def veto(self, reason):
        self.ok = False
        self.reasons.append(reason)

    def __str__(self):
        return "ok" if self.ok else "no: " + "; ".join(self.reasons)


def judge(base, trial, target_key, guard, min_density, min_overlap):
    """Score a candidate against the baseline over everything that moved.

    ``base`` and ``trial`` are the dicts :meth:`Polisher._observe` builds:
    per-key (side-chain density, backbone density, contacts), rotalyze and
    ramalyze entries over the window, and the neighbourhood clash list.
    Pure, so it can be tested without Coot or MolProbity.
    """
    verdict = Verdict()
    verdict.target_rot = trial["rot"].get(target_key)
    verdict.target_rama = trial["ram"].get(target_key)
    verdict.clashes = trial["clashes"]
    side_values = []
    for key in base["keys"]:
        if key not in trial["density"]:
            verdict.veto("{} missing".format(key_label(key)))
            continue
        side0, back0, contacts0 = base["density"][key][0], base["density"][key][1], base["contacts"][key]
        side1, back1, contacts1 = trial["density"][key][0], trial["density"][key][1], trial["contacts"][key]
        if side1 == side1:
            side_values.append(side1)
        if not density_ok(side0, side1, guard, min_density):
            verdict.veto("side-chain density {} {:.2f} -> {:.2f}".format(key_label(key), side0, side1))
        if not density_ok(back0, back1, guard, min_density):
            verdict.veto("backbone density {} {:.2f} -> {:.2f}".format(key_label(key), back0, back1))
        if contacts1 - contacts0:
            verdict.veto("new contact at {}".format(key_label(key)))
        entry = trial["rot"].get(key)
        if key != target_key and entry and entry["evaluation"] == "OUTLIER" \
                and base["rot"].get(key, {}).get("evaluation") != "OUTLIER":
            verdict.veto("new rotamer outlier {}".format(key_label(key)))
    for key, entry in trial["ram"].items():
        if entry["evaluation"] == "OUTLIER" and key != target_key \
                and base["ram"].get(key, {}).get("evaluation") != "OUTLIER":
            verdict.veto("new Ramachandran outlier {}".format(key_label(key)))
    verdict.favoured = sum(1 for key in base["keys"]
                           if trial["ram"].get(key, {}).get("evaluation") == "Favored")
    verdict.mean_sidechain = sum(side_values) / len(side_values) if side_values else float("nan")
    added = new_clashes(base["clashes"], trial["clashes"], min_overlap)
    if added:
        worst = max(added, key=lambda c: c["overlap"])
        verdict.veto("new clash {} {} - {} {} ({:.2f})".format(
            key_label(worst["key1"], worst["resname1"]), worst["atom1"],
            key_label(worst["key2"], worst["resname2"]), worst["atom2"], worst["overlap"]))
    return verdict


class Polisher:
    """Run the stages on one model, keeping a change log and a review list.

    ``progress`` is called with a short ASCII line after every decision so a
    running report has something to show; ``snapshot`` is called with the
    current state dict at the same moments.
    """

    def __init__(self, xyzin, map_mtz, workdir, options=None, diff_mtz=None,
                 dict_path=None, progress=print, snapshot=None):
        self.xyzin = xyzin
        self.map_mtz = map_mtz
        self.diff_mtz = diff_mtz
        self.workdir = workdir
        self.opt = options or Options()
        self.progress = progress
        self.snapshot = snapshot
        self.changes = []           # accepted changes
        self.review = []            # what a human should look at
        self.stages = []            # per-stage counts
        self.before = None
        self.after = None
        self.dictionaries = []
        if dict_path:
            self.dictionaries.append(dict_path)
        self.grid = None
        self.diff_grid = None
        self.st = None
        self.cur = None             # path of the current model (PDB, for Coot)

    # ---- bookkeeping ------------------------------------------------------

    def say(self, text):
        self.progress(text)

    def _snapshot(self):
        if self.snapshot is not None:
            self.snapshot(self.state())

    def state(self):
        return dict(before=self.before, after=self.after, stages=list(self.stages),
                    changes=list(self.changes), review=list(self.review))

    def _path(self, name):
        return os.path.join(self.workdir, name)

    def _commit(self, note):
        """Write the current structure as the model every later trial starts from."""
        self.cur = self._path("polish_current.pdb")
        self.st.write_pdb(self.cur)
        self.say(note)
        self._snapshot()

    def _record(self, stage, key, resname, action, before, after, detail=""):
        self.changes.append(dict(stage=stage, chain=key[0], resno=key[1], icode=key[2], resname=resname,
                                 action=action, density_before=_r(before), density_after=_r(after),
                                 detail=detail))

    def _review(self, kind, key, resname, detail=""):
        self.review.append(dict(kind=kind, chain=key[0], resno=key[1], icode=key[2], resname=resname,
                                detail=detail))

    def _selected(self, chain_id):
        return self.opt.chains is None or chain_id in self.opt.chains

    def _coot(self, pdb_path):
        return CootSession(pdb_path, self.map_mtz, self.dictionaries)

    # ---- setup --------------------------------------------------------------

    def prepare(self):
        self.grid = load_map(self.map_mtz)
        if self.diff_mtz:
            try:
                self.diff_grid = load_map(self.diff_mtz)
            except PolishError as exc:
                self.say("WARNING: difference map not used: {}".format(exc))
        self.st = read_structure(self.xyzin)
        for name in nonstandard_residues(self.st):
            path = monomer_dictionary(name)
            if path and path not in self.dictionaries:
                self.dictionaries.append(path)
        self.cur = self._path("polish_input.pdb")
        self.st.write_pdb(self.cur)
        self.say("model {} chains {}; {} dictionaries".format(
            os.path.basename(self.xyzin), ",".join(polymer_chains(self.st)), len(self.dictionaries)))

    def run(self, xyzout):
        started = time.time()
        self.prepare()
        self.say("MolProbity before")
        self.before = summarize(self.cur)
        self.say("before: {}".format(_fmt(self.before)))
        self._snapshot()
        if self.opt.do_rotamers:
            self._stage("rotamers", self.stage_rotamers)
        if self.opt.do_rama:
            self._stage("ramachandran", self.stage_rama)
        if self.opt.do_clashes:
            self._stage("clashes", self.stage_clashes)
        write_structure(self.st, xyzout)
        self.say("MolProbity after")
        self.after = summarize(self.cur)
        self.say("after: {}".format(_fmt(self.after)))
        self.finish_review()
        self._snapshot()
        self.say("total {:.0f} s".format(time.time() - started))
        return xyzout

    def _stage(self, name, method):
        started = time.time()
        self.say("stage {}".format(name))
        try:
            targets, fixed = method()
        except PolishError:
            raise
        except Exception as exc:  # one stage failing must not lose the others
            self.say("ERROR: stage {} failed: {}: {}".format(name, type(exc).__name__, exc))
            targets, fixed = 0, 0
        self.stages.append(dict(name=name, targets=targets, fixed=fixed,
                                left=max(0, targets - fixed), seconds=int(time.time() - started)))
        self.say("stage {}: {} targets, {} fixed, {:.0f} s".format(
            name, targets, fixed, time.time() - started))
        self._snapshot()

    # ---- observing a model where a candidate moves atoms ------------------------

    def _observe(self, st, chain_id, moved, selection=None):
        """Everything the judge compares, for one model over chain positions ``moved``."""
        lo, hi = moved[0], moved[-1]
        keys = keys_at(st, chain_id, moved)
        rot = rotalyze(write_index_window(st, chain_id, lo - 1, hi + 1, self._path("polish_window.pdb")))
        ram = ramalyze(write_index_window(st, chain_id, lo - RAMA_MARGIN - 1, hi + RAMA_MARGIN + 1,
                                          self._path("polish_window.pdb")))
        density, contacts = {}, {}
        for key in keys:
            res = find_by_key(st, key)
            if res is None:
                continue
            density[key] = (sidechain_density(self.grid, res), backbone_density(self.grid, res))
            contacts[key] = clash_contacts(st, key)
        if selection is None:
            selection = neighbourhood_keys(st, keys)
        _, clashes = clashscore(write_neighbourhood(st, keys, self._path("polish_neighbourhood.pdb"),
                                                    selection=selection))
        return dict(keys=keys, lo=lo, hi=hi, rot=rot, ram=ram, density=density, contacts=contacts,
                    clashes=clashes, selection=selection)

    def _judge(self, trial_st, chain_id, base, target_key):
        trial = self._observe(trial_st, chain_id, list(range(base["lo"], base["hi"] + 1)), base["selection"])
        return judge(base, trial, target_key, self.opt.guard, self.opt.min_density, self.opt.min_overlap)

    def _accept(self, trial_st, chain_id, base):
        copy_window(trial_st, self.st, chain_id, range(base["lo"], base["hi"] + 1))

    # ---- stage 1: rotamers --------------------------------------------------

    def _still_outlier(self, key):
        index = locate(self.st, key)
        entry = rotalyze(write_index_window(self.st, key[0], index, index, self._path("polish_window.pdb"))).get(key)
        return entry is not None and entry["evaluation"] == "OUTLIER"

    def stage_rotamers(self):
        rot0 = rotalyze(self.cur)
        targets = sorted(k for k, v in rot0.items()
                         if v["evaluation"] == "OUTLIER" and self._selected(k[0]))
        self.say("{} rotamer outliers".format(len(targets)))
        fixed = 0
        for position, key in enumerate(targets, 1):
            index = locate(self.st, key)
            if index is None:
                continue
            res0 = residue_at(self.st, key[0], index)
            label = key_label(key, res0.name)
            if position > 1 and not self._still_outlier(key):
                fixed += 1
                self._record("rotamers", key, res0.name, "resolved by a neighbour's change",
                             float("nan"), float("nan"))
                self.say("[{}/{}] {}: no longer an outlier".format(position, len(targets), label))
                continue
            d0 = sidechain_density(self.grid, res0)
            self.say("[{}/{}] {} density {:.2f}".format(position, len(targets), label, d0))
            base = self._observe(self.st, key[0], window(self.st, key[0], index, 1, 1))
            trials = []
            if self.opt.use_ncs and not has_altlocs(res0):
                for donor in self._donors_for(key, res0.name, rot0):
                    trials.append(("ncs_" + donor, dict(donor=donor, weight=TORSION_WEIGHTS[0])))
            trials.append(("autofit", dict(autofit=True)))
            for weight in TORSION_WEIGHTS:
                trials.append(("torsion{:g}".format(weight), dict(weight=weight)))
            for k in range(ROTAMER_LIBRARY_TRIALS):
                trials.append(("rotamer{}".format(k), dict(rotamer=k, weight=TORSION_WEIGHTS[0])))
            tried = []
            accepted = None
            for name, how in trials:
                try:
                    trial_st = self._sidechain_trial(key, how)
                except PolishError:
                    raise
                except Exception as exc:
                    self.say("   {}: failed ({}: {})".format(name, type(exc).__name__, exc))
                    continue
                if trial_st is None:
                    if "rotamer" in how:
                        break       # no such library rotamer: nothing further to try
                    continue
                verdict = self._judge(trial_st, key[0], base, key)
                entry = verdict.target_rot
                if entry is None or entry["evaluation"] == "OUTLIER":
                    verdict.veto("still an outlier")
                d1 = sidechain_density(self.grid, find_by_key(trial_st, key))
                tried.append("{}:{:.2f}:{}".format(name, d1, entry["evaluation"] if entry else "?"))
                self.say("   {}: {}".format(name, verdict))
                if verdict.ok:
                    self._accept(trial_st, key[0], base)
                    accepted = (name, d1, entry["rotamer"])
                    break
            if accepted:
                fixed += 1
                self._record("rotamers", key, res0.name, accepted[0], d0, accepted[1],
                             "now {}".format(accepted[2]))
                self._commit("   fixed {} via {} ({:.2f} -> {:.2f})".format(label, accepted[0], d0, accepted[1]))
            else:
                self._review("rotamer outlier kept", key, res0.name,
                             "density {:.2f}; tried {}".format(d0, " ".join(tried) or "nothing"))
                self.say("   left {}: {}".format(label, " ".join(tried)))
                self._snapshot()
        return len(targets), fixed

    def _donors_for(self, key, resname, scores, favoured_only=False):
        """Other chains' copies of this residue, best MolProbity score first."""
        rank = {"Favored": 0, "Allowed": 1}
        donors = []
        for other in polymer_chains(self.st):
            if other == key[0]:
                continue
            other_key = rkey(other, key[1], key[2])
            entry = scores.get(other_key)
            res = find_by_key(self.st, other_key)
            if entry is None or res is None or res.name != resname or has_altlocs(res):
                continue
            if entry["evaluation"] not in rank:
                continue
            if favoured_only and entry["evaluation"] != "Favored":
                continue
            donors.append((rank[entry["evaluation"]], -entry["score"], other))
        donors.sort()
        return [d[2] for d in donors[:NCS_DONORS]]

    def _sidechain_trial(self, key, how):
        """Build one candidate for a side chain: the whole refined triple.

        Returns the trial structure, or None when the requested library rotamer
        does not exist or the NCS transplant could not be made.
        """
        start = self.cur
        if "donor" in how:
            st = read_structure(self.cur)
            if transplant_sidechain(st, key[0], key[1], how["donor"], icode=key[2]) is None:
                return None
            start = self._path("polish_trial_in.pdb")
            st.write_pdb(start)
        coot = self._coot(start)
        if how.get("autofit"):
            coot.auto_fit(key)
        if "rotamer" in how:
            if not coot.set_rotamer(key, how["rotamer"]):
                return None
        if "weight" in how:
            coot.use_torsion_restraints(how["weight"])
        coot.refine(key, "TRIPLE", self.opt.cycles)
        return read_structure(coot.write(self._path("polish_trial.pdb")))

    # ---- stage 2: Ramachandran ------------------------------------------------

    def stage_rama(self):
        ram0 = ramalyze(self.cur)
        targets = sorted(k for k, v in ram0.items()
                         if v["evaluation"] == "OUTLIER" and self._selected(k[0]))
        self.say("{} Ramachandran outliers".format(len(targets)))
        fixed = 0
        for position, key in enumerate(targets, 1):
            index = locate(self.st, key)
            if index is None:
                continue
            res0 = residue_at(self.st, key[0], index)
            label = key_label(key, res0.name)
            base = self._observe(self.st, key[0], window(self.st, key[0], index, 2, 2))
            if position > 1 and base["ram"].get(key, {}).get("evaluation") != "OUTLIER":
                fixed += 1
                self._record("ramachandran", key, res0.name, "resolved by a neighbour's change",
                             float("nan"), float("nan"))
                self.say("[{}/{}] {}: no longer an outlier".format(position, len(targets), label))
                continue
            d0 = base["density"][key][1]
            self.say("[{}/{}] {} backbone density {:.2f}".format(position, len(targets), label, d0))
            candidates = []
            if self.opt.use_ncs and not has_altlocs(res0):
                for donor in self._donors_for(key, res0.name, ram0, favoured_only=True):
                    candidates.append(("ncs_" + donor, dict(donor=donor)))
            candidates.append(("rama_refine", dict(rama=True)))
            if self.opt.do_flips:
                before = residue_at(self.st, key[0], index - 1) if index > 0 else None
                if before is not None:
                    candidates.append(("flip_before", dict(flip=residue_key(key[0], before))))
                candidates.append(("flip_after", dict(flip=key)))
            scored = []
            for name, how in candidates:
                try:
                    trial_st = self._backbone_trial(key, index, how)
                except PolishError:
                    raise
                except Exception as exc:
                    self.say("   {}: failed ({}: {})".format(name, type(exc).__name__, exc))
                    continue
                if trial_st is None:
                    continue
                verdict = self._judge(trial_st, key[0], base, key)
                entry = verdict.target_rama
                if entry is None or entry["evaluation"] == "OUTLIER":
                    verdict.veto("still an outlier")
                if verdict.target_rot is not None and verdict.target_rot["evaluation"] == "OUTLIER" \
                        and base["rot"].get(key, {}).get("evaluation") != "OUTLIER":
                    verdict.veto("new rotamer outlier at the target")
                self.say("   {}: {}".format(name, verdict))
                if verdict.ok:
                    scored.append(((verdict.favoured, _nz(verdict.mean_sidechain)), name, trial_st))
            if scored:
                scored.sort(key=lambda item: item[0], reverse=True)
                _, name, trial_st = scored[0]
                self._accept(trial_st, key[0], base)
                d1 = backbone_density(self.grid, find_by_key(self.st, key))
                fixed += 1
                self._record("ramachandran", key, res0.name, name, d0, d1)
                self._commit("   fixed {} via {} ({:.2f} -> {:.2f})".format(label, name, d0, d1))
            else:
                self._review("Ramachandran outlier kept", key, res0.name,
                             "backbone density {:.2f}; no candidate removed it without cost".format(d0))
                self._snapshot()
        return len(targets), fixed

    def _backbone_trial(self, key, index, how):
        start = self.cur
        if "donor" in how:
            st = read_structure(self.cur)
            if transplant_backbone_indices(st, key[0], index - 1, index + 1, how["donor"]) is None:
                return None
            start = self._path("polish_trial_in.pdb")
            st.write_pdb(start)
        coot = self._coot(start)
        if how.get("rama"):
            coot.use_rama_restraints(RAMA_WEIGHT)
        if "flip" in how:
            if not coot.flip_peptide(how["flip"]):
                return None
        coot.refine(key, "QUINTUPLE", max(50, self.opt.cycles // 2))
        return read_structure(coot.write(self._path("polish_trial.pdb")))

    # ---- stage 3: clashes -------------------------------------------------------

    def stage_clashes(self):
        score, clashes = clashscore(self.cur)
        seen, targets = set(), []
        for clash in clashes:
            if clash["overlap"] < self.opt.min_overlap:
                continue
            pair = clash_pair(clash)
            if pair in seen or not (self._selected(pair[0][0]) or self._selected(pair[1][0])):
                continue
            seen.add(pair)
            partners = []
            for side in ("1", "2"):
                key, resname, atom = clash["key" + side], clash["resname" + side], clash["atom" + side]
                if atom in BACKBONE_H or resname == "GLY" or not self._selected(key[0]):
                    continue
                if not is_amino_acid(resname):
                    continue        # ligands and ions are not rotamer problems
                res = find_by_key(self.st, key)
                if res is None or has_altlocs(res):
                    continue
                partners.append((key, resname))
            targets.append((pair, partners, clash))
        self.say("clashscore {}; {} pairs with overlap >= {}".format(score, len(targets), self.opt.min_overlap))
        fixed = 0
        for position, (pair, partners, clash) in enumerate(targets, 1):
            label = "{} {} - {} {}".format(key_label(clash["key1"], clash["resname1"]), clash["atom1"],
                                           key_label(clash["key2"], clash["resname2"]), clash["atom2"])
            if not partners:
                self._review("clash left for refinement", clash["key1"], clash["resname1"],
                             "overlap {:.2f} with {} {}; no movable side chain".format(
                                 clash["overlap"], key_label(clash["key2"], clash["resname2"]), clash["atom2"]))
                continue
            self.say("[{}/{}] {} overlap {:.2f}".format(position, len(targets), label, clash["overlap"]))
            done = False
            for key, resname in partners:
                index = locate(self.st, key)
                res0 = residue_at(self.st, key[0], index)
                d0 = sidechain_density(self.grid, res0)
                base = self._observe(self.st, key[0], window(self.st, key[0], index, 1, 1))
                if self._pair_overlap(base["clashes"], pair) < self.opt.min_overlap:
                    done = True     # already resolved by an earlier change
                    self.say("   no longer clashing")
                    break
                trials = [("autofit", dict(autofit=True))]
                trials += [("rotamer{}".format(k), dict(rotamer=k, weight=TORSION_WEIGHTS[0]))
                           for k in range(self.opt.max_clash_rotamers)]
                for name, how in trials:
                    try:
                        trial_st = self._sidechain_trial(key, how)
                    except PolishError:
                        raise
                    except Exception as exc:
                        self.say("   {}: failed ({}: {})".format(name, type(exc).__name__, exc))
                        continue
                    if trial_st is None:
                        break
                    verdict = self._judge(trial_st, key[0], base, key)
                    entry = verdict.target_rot
                    if entry is None or entry["evaluation"] == "OUTLIER":
                        verdict.veto("rotamer outlier")
                    worst = self._pair_overlap(verdict.clashes, pair)
                    if worst >= self.opt.min_overlap:
                        verdict.veto("still overlapping {:.2f}".format(worst))
                    self.say("   {} {}: {}".format(key_label(key, resname), name, verdict))
                    if not verdict.ok:
                        continue
                    self._accept(trial_st, key[0], base)
                    d1 = sidechain_density(self.grid, find_by_key(self.st, key))
                    fixed += 1
                    other = clash["key2"] if clash["key1"] == key else clash["key1"]
                    other_name = clash["resname2"] if clash["key1"] == key else clash["resname1"]
                    self._record("clashes", key, resname, name, d0, d1,
                                 "overlap {:.2f} with {}".format(clash["overlap"], key_label(other, other_name)))
                    self._commit("   fixed {} via {} {} ({:.2f} -> {:.2f})".format(
                        label, key_label(key, resname), name, d0, d1))
                    done = True
                    break
                if done:
                    break
            if not done:
                self._review("clash kept", clash["key1"], clash["resname1"],
                             "overlap {:.2f} with {} {}; no side-chain fix kept the density".format(
                                 clash["overlap"], key_label(clash["key2"], clash["resname2"]), clash["atom2"]))
                self.say("   left {}".format(label))
                self._snapshot()
        return len(targets), fixed

    @staticmethod
    def _pair_overlap(clashes, pair):
        worst = 0.0
        for clash in clashes:
            if clash_pair(clash) == pair:
                worst = max(worst, clash["overlap"])
        return worst

    # ---- wrap-up --------------------------------------------------------------

    def finish_review(self):
        rot, ram = rotalyze(self.cur), ramalyze(self.cur)
        listed = {(item["chain"], item["resno"], item["icode"], item["kind"]) for item in self.review}
        for key, v in sorted(rot.items()):
            if v["evaluation"] == "OUTLIER" and key + ("rotamer outlier kept",) not in listed:
                self._review("rotamer outlier", key, v["resname"], "still an outlier after polish")
        for key, v in sorted(ram.items()):
            if v["evaluation"] == "OUTLIER" and key + ("Ramachandran outlier kept",) not in listed:
                self._review("Ramachandran outlier", key, v["resname"], "still an outlier after polish")

    def write_changes_csv(self, path, original_path):
        """Every residue whose atoms moved, with density and MolProbity status before and after."""
        st0, st1 = read_structure(original_path), self.st
        rot0, rot1 = rotalyze(self._path("polish_input.pdb")), rotalyze(self.cur)
        ram0, ram1 = ramalyze(self._path("polish_input.pdb")), ramalyze(self.cur)
        by_key = {}
        for change in self.changes:
            by_key.setdefault((change["chain"], change["resno"], change["icode"]), []).append(
                "{}:{}".format(change["stage"], change["action"]))
        rows = changed_residue_rows(st0, st1, self.grid, self.diff_grid, rot0, rot1, ram0, ram1, by_key)
        write_changes(path, rows)
        return rows


def changed_residue_rows(st0, st1, grid, diff_grid, rot0, rot1, ram0, ram1, actions):
    rows = []
    for chain in st0[0]:
        for r0 in chain:
            key = residue_key(chain.name, r0)
            r1 = find_by_key(st1, key)
            if r1 is None or r0.name != r1.name:
                continue
            shift = 0.0
            for atom in r0:
                other = r1.find_atom(atom.name, "*")
                if other is not None:
                    shift = max(shift, atom.pos.dist(other.pos))
            if shift < 0.05:
                continue
            row = dict(chain=chain.name, resno=r0.seqid.num, icode=key[2], resname=r0.name,
                       max_shift=round(shift, 2),
                       action=" ".join(actions.get(key, [])) or "refined neighbour",
                       density_before=_r(sidechain_density(grid, r0)),
                       density_after=_r(sidechain_density(grid, r1)),
                       rotamer_before=rot0.get(key, {}).get("rotamer", ""),
                       rotamer_after=rot1.get(key, {}).get("rotamer", ""),
                       rama_before=ram0.get(key, {}).get("evaluation", ""),
                       rama_after=ram1.get(key, {}).get("evaluation", ""))
            if diff_grid is not None:
                row["diff_density_before"] = _r(sidechain_density(diff_grid, r0))
                row["diff_density_after"] = _r(sidechain_density(diff_grid, r1))
            rows.append(row)
    return rows


CHANGES_COLUMNS = ("chain", "resno", "icode", "resname", "max_shift", "action", "density_before",
                   "density_after", "diff_density_before", "diff_density_after", "rotamer_before",
                   "rotamer_after", "rama_before", "rama_after")


def write_changes(path, rows):
    """The per-residue change log as CSV with a fixed column order."""
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(CHANGES_COLUMNS), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in CHANGES_COLUMNS})
    return path


def _r(value):
    if value is None or value != value:
        return ""
    return round(float(value), 2)


def _nz(value):
    return 0.0 if value != value else value


def _fmt(summary):
    return ("rotamer outliers {rot_outliers} ({rot_pct}%), Ramachandran outliers {rama_outliers} "
            "({rama_pct}%), favoured {rama_favored_pct}%, clashscore {clashscore}").format(**summary)
