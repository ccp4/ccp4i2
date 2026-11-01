"""NCS-aware geometry polish: the library behind the geometry_polish task.

Three tools, three roles:

* MolProbity (molprobity.rotalyze / ramalyze / clashscore) is the judge. It is
  what the wwPDB validation report runs, so a residue it no longer flags is a
  residue the report will not flag either.
* gemmi reads the maps and the model, scores density at atoms, finds clashes
  and transplants side chains or backbone windows between NCS copies.
* Coot's headless API ("chapi") does the building: auto-fit rotamer, library
  rotamers, peptide flips and torsion- or Ramachandran-restrained real-space
  refinement.

Every candidate change is accepted only if MolProbity stops flagging the
residue, no new heavy-atom contact appears, and the 2mFo-DFc density at the
moved atoms is kept (see :func:`density_ok`). What cannot be fixed that way is
listed for a human rather than forced.

Only the standard library is imported at module level: gemmi, numpy and
coot_headless_api are execution-time dependencies and every task plugin must
import with none of them present.
"""
import csv
import os
import re
import subprocess
import tempfile
import time

BACKBONE = {"N", "CA", "C", "O", "CB"}
BB_ONLY = ("N", "CA", "C", "O")
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
ABSOLUTE_SLACK = 0.3                # sigma a change may lose outright
CONTACT_DISTANCE = 2.4              # A; a new heavy-atom contact closer than this vetoes a change
TORSION_WEIGHTS = (5.0, 10.0, 20.0)
RAMA_WEIGHT = 10.0


class PolishError(RuntimeError):
    """A tool the polish depends on is missing or failed outright."""


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


_RESIDUE_KEY = re.compile(r"\s*(\S+)\s+(-?\d+)\S*\s+(\S+)")


def _residue_key(text):
    """(chain, resno, resname) from the residue field MolProbity prints."""
    found = _RESIDUE_KEY.match(text)
    if not found:
        return None
    return found.group(1), int(found.group(2)), found.group(3)


def parse_rotalyze(text):
    """{(chain, resno): {score, evaluation, rotamer, resname}} from rotalyze output.

    Lines are ``residue:occupancy:score%:chi1:chi2:chi3:chi4:evaluation:rotamer``.
    """
    result = {}
    for line in text.splitlines():
        if ":" not in line or line.startswith(("residue", "SUMMARY")):
            continue
        parts = line.split(":")
        if len(parts) < 9:
            continue
        key = _residue_key(parts[0])
        if key is None:
            continue
        try:
            score = float(parts[2])
        except ValueError:
            continue
        result[(key[0], key[1])] = dict(score=score, evaluation=parts[7].strip(),
                                        rotamer=parts[8].strip(), resname=key[2])
    return result


def parse_ramalyze(text):
    """{(chain, resno): {score, evaluation, type, resname}} from ramalyze output.

    Lines are ``residue:score%:phi:psi:evaluation:type``.
    """
    result = {}
    for line in text.splitlines():
        if ":" not in line or line.startswith(("residue", "SUMMARY")):
            continue
        parts = line.split(":")
        if len(parts) < 6:
            continue
        key = _residue_key(parts[0])
        if key is None:
            continue
        try:
            score = float(parts[1])
        except ValueError:
            continue
        result[(key[0], key[1])] = dict(score=score, evaluation=parts[4].strip(),
                                        type=parts[5].strip(), resname=key[2])
    return result


_CLASH_LINE = re.compile(
    r"\s*(\S+)\s+(-?\d+)\S?\s+(\S+)\s+(\S+)\s+(\S+)\s+(-?\d+)\S?\s+(\S+)\s+(\S+)\s*:([\d.]+)")


def parse_clash_line(line):
    """One clashscore contact line, e.g. `` C 117  PHE  HB3  C 138  SER  HB3 :1.153``."""
    found = _CLASH_LINE.match(line)
    if not found:
        return None
    return dict(chain1=found.group(1), resno1=int(found.group(2)), resname1=found.group(3),
                atom1=found.group(4), chain2=found.group(5), resno2=int(found.group(6)),
                resname2=found.group(7), atom2=found.group(8), overlap=float(found.group(9)))


def parse_clashscore(text):
    """(clashscore or None, [clash dicts]) from clashscore output."""
    score, clashes = None, []
    for line in text.splitlines():
        if line.startswith("clashscore"):
            try:
                score = float(line.split("=")[1])
            except (IndexError, ValueError):
                pass
        elif ":" in line and line[:1] == " ":
            clash = parse_clash_line(line)
            if clash:
                clashes.append(clash)
    return score, clashes


def clash_pair(clash):
    """The unordered residue pair a clash belongs to."""
    return tuple(sorted([(clash["chain1"], clash["resno1"]), (clash["chain2"], clash["resno2"])]))


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

    Accepted when the density is within ``slack`` sigma of the original, or
    when it stays above both ``min_density`` and ``guard_fraction`` of the
    original. The second clause is what let strained side chains in strong
    density (4 sigma) move to a rotamer at 2.5 sigma, while a side chain at
    1 sigma may not drop below 0.8.
    """
    if after != after:  # NaN: nothing to score (no side-chain atoms)
        return True
    if before != before:
        return True
    return after >= before - slack or after >= max(min_density, guard_fraction * before)


# --------------------------------------------------------------------------
# gemmi: maps, density, contacts, transplants, windows
# --------------------------------------------------------------------------

_map_cache = {}


def load_map(mtz_path, sample_rate=3.0):
    """Sigma-scaled map from an F,PHI mini-MTZ (cached per path)."""
    import gemmi
    if mtz_path in _map_cache:
        return _map_cache[mtz_path]
    mtz = gemmi.read_mtz_file(mtz_path)
    labels = [c.label for c in mtz.columns]
    f_label = "F" if "F" in labels else next((l for l in labels if l.startswith("F")), None)
    phi_label = "PHI" if "PHI" in labels else next((l for l in labels if l.startswith("PH")), None)
    if f_label is None or phi_label is None:
        raise PolishError("{} has no F,PHI columns (found {})".format(mtz_path, ",".join(labels)))
    grid = mtz.transform_f_phi_to_map(f_label, phi_label, sample_rate=sample_rate)
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


def find_residue(st, chain_id, resno):
    chain = st[0].find_chain(chain_id)
    if chain is None:
        return None
    for res in chain:
        if res.seqid.num == resno:
            return res
    return None


def polymer_chains(st):
    return [ch.name for ch in st[0] if ch.get_polymer()]


def residue_density(grid, residue, keep=lambda atom: True):
    """Mean sigma-scaled density at the centres of the atoms ``keep`` selects."""
    values = [grid.interpolate_value(a.pos) for a in residue if keep(a)]
    if not values:
        return float("nan")
    return float(sum(values) / len(values))


def sidechain_density(grid, residue):
    return residue_density(grid, residue, lambda a: a.name not in BACKBONE and a.element.name != "H")


def backbone_density(grid, residue):
    return residue_density(grid, residue, lambda a: a.name in BB_ONLY)


def clash_contacts(st, chain_id, resno, dist=CONTACT_DISTANCE):
    """Heavy-atom contacts closer than ``dist`` involving one residue.

    Symmetry aware; adjacent residues ignored. Returns a set of
    (other chain, other resno, atom, other atom) so two models can be compared.
    """
    import gemmi
    ns = gemmi.NeighborSearch(st[0], st.cell, 5.0).populate(include_h=False)
    cs = gemmi.ContactSearch(dist)
    cs.ignore = gemmi.ContactSearch.Ignore.AdjacentResidues
    found = set()
    for contact in cs.find_contacts(ns):
        for a, b in ((contact.partner1, contact.partner2), (contact.partner2, contact.partner1)):
            if a.chain.name == chain_id and a.residue.seqid.num == resno:
                found.add((b.chain.name, b.residue.seqid.num, a.atom.name, b.atom.name))
    return found


def copy_positions(src, dst, include_backbone=True):
    """Move ``dst``'s atoms onto ``src``'s positions, adding atoms ``dst`` lacks."""
    import gemmi
    for atom in src:
        if not include_backbone and atom.name in BACKBONE:
            continue
        target = dst.find_atom(atom.name, "*")
        if target is not None:
            target.pos = atom.pos
        elif atom.element.name != "H":
            new_atom = gemmi.Atom()
            new_atom.name = atom.name
            new_atom.element = atom.element
            new_atom.pos = atom.pos
            new_atom.occ = atom.occ
            anchor = dst.find_atom("CB", "*") or dst.find_atom("CA", "*")
            new_atom.b_iso = anchor.b_iso if anchor is not None else atom.b_iso
            dst.add_atom(new_atom)


def transplant_sidechain(st, chain_id, resno, donor_chain, window=1):
    """Copy the donor copy's side chain onto the target after superposing the
    local backbone (N, CA, C, O, CB of resno-window..resno+window).

    Returns the superposition rmsd, or None when the copies cannot be paired.
    """
    import gemmi
    target_pos, donor_pos = [], []
    for n in range(resno - window, resno + window + 1):
        tr, dr = find_residue(st, chain_id, n), find_residue(st, donor_chain, n)
        if tr is None or dr is None:
            if n == resno:
                return None
            continue
        for name in ("N", "CA", "C", "O", "CB"):
            ta, da = tr.find_atom(name, "*"), dr.find_atom(name, "*")
            if ta is not None and da is not None:
                target_pos.append(ta.pos)
                donor_pos.append(da.pos)
    if len(target_pos) < 4:
        return None
    tr, dr = find_residue(st, chain_id, resno), find_residue(st, donor_chain, resno)
    if tr.name != dr.name:
        return None
    sup = gemmi.superpose_positions(target_pos, donor_pos)
    for da in dr:
        if da.name in BACKBONE or da.element.name == "H":
            continue
        new_pos = gemmi.Position(sup.transform.apply(da.pos))
        ta = tr.find_atom(da.name, "*")
        if ta is not None:
            ta.pos = new_pos
        else:
            atom = gemmi.Atom()
            atom.name, atom.element, atom.pos, atom.occ = da.name, da.element, new_pos, 1.0
            anchor = tr.find_atom("CB", "*")
            atom.b_iso = anchor.b_iso if anchor is not None else da.b_iso
            tr.add_atom(atom)
    return sup.rmsd


def transplant_backbone_window(st, chain_id, start, end, donor_chain, flank=2):
    """Replace residues start..end of the target with the donor copy's atoms after
    superposing on the flanking residues' backbone. Returns rmsd or None."""
    import gemmi
    target_pos, donor_pos = [], []
    for n in list(range(start - flank, start)) + list(range(end + 1, end + flank + 1)):
        tr, dr = find_residue(st, chain_id, n), find_residue(st, donor_chain, n)
        if tr is None or dr is None:
            continue
        for name in BB_ONLY:
            ta, da = tr.find_atom(name, "*"), dr.find_atom(name, "*")
            if ta is not None and da is not None:
                target_pos.append(ta.pos)
                donor_pos.append(da.pos)
    if len(target_pos) < 6:
        return None
    pairs = []
    for n in range(start, end + 1):
        tr, dr = find_residue(st, chain_id, n), find_residue(st, donor_chain, n)
        if tr is None or dr is None or tr.name != dr.name:
            return None
        pairs.append((tr, dr))
    sup = gemmi.superpose_positions(target_pos, donor_pos)
    for tr, dr in pairs:
        for da in dr:
            if da.element.name == "H":
                continue
            ta = tr.find_atom(da.name, "*")
            if ta is not None:
                ta.pos = gemmi.Position(sup.transform.apply(da.pos))
    return sup.rmsd


def write_window(st, chain_id, lo, hi, path):
    """A PDB holding residues lo..hi of one chain: enough for rotalyze/ramalyze."""
    import gemmi
    mini = gemmi.Structure()
    mini.cell = st.cell
    mini.spacegroup_hm = st.spacegroup_hm
    model = gemmi.Model("1")
    chain = gemmi.Chain(chain_id)
    source = st[0].find_chain(chain_id)
    if source is not None:
        for res in source:
            if lo <= res.seqid.num <= hi:
                chain.add_residue(res)
    model.add_chain(chain)
    mini.add_model(model)
    mini.setup_entities()
    mini.write_pdb(path)
    return path


def write_neighbourhood(st, residues, path, radius=6.0):
    """A PDB holding every residue with an atom within ``radius`` of the given
    (chain, resno) residues: enough for clashscore to re-score those contacts."""
    import gemmi
    ns = gemmi.NeighborSearch(st[0], st.cell, max(5.0, radius)).populate(include_h=False)
    keep = set()
    for chain_id, resno in residues:
        res = find_residue(st, chain_id, resno)
        if res is None:
            continue
        keep.add((chain_id, resno))
        for atom in res:
            for mark in ns.find_atoms(atom.pos, "\0", radius=radius):
                cra = mark.to_cra(st[0])
                keep.add((cra.chain.name, cra.residue.seqid.num))
    mini = gemmi.Structure()
    mini.cell = st.cell
    mini.spacegroup_hm = st.spacegroup_hm
    model = gemmi.Model("1")
    for source in st[0]:
        chain = gemmi.Chain(source.name)
        for res in source:
            if (source.name, res.seqid.num) in keep:
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
        self.imap = mc.read_mtz(map_mtz, "F", "PHI", "", False, False)
        mc.set_imol_refinement_map(self.imap)

    def auto_fit(self, chain_id, resno):
        return self.mc.auto_fit_rotamer(self.imol, chain_id, resno, "", "", self.imap)

    def set_rotamer(self, chain_id, resno, number):
        return bool(self.mc.set_residue_to_rotamer_number(
            self.imol, "//{}/{}".format(chain_id, resno), "", number))

    def use_torsion_restraints(self, weight):
        self.mc.set_use_torsion_restraints(True)
        self.mc.set_torsion_restraints_weight(weight)

    def use_rama_restraints(self, weight):
        self.mc.set_use_rama_plot_restraints(True)
        self.mc.set_rama_plot_restraints_weight(weight)

    def refine_triple(self, chain_id, resno, cycles):
        return self.mc.refine_residues(self.imol, chain_id, resno, "", "", "TRIPLE", cycles)

    def refine_range(self, chain_id, lo, hi, cycles):
        return self.mc.refine_residue_range(self.imol, chain_id, lo, hi, cycles)

    def flip_peptide(self, chain_id, resno):
        return bool(self.mc.flip_peptide_using_cid(self.imol, "//{}/{}".format(chain_id, resno), ""))

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

    def _record(self, stage, chain_id, resno, resname, action, before, after, detail=""):
        self.changes.append(dict(stage=stage, chain=chain_id, resno=resno, resname=resname,
                                 action=action, density_before=_r(before), density_after=_r(after),
                                 detail=detail))

    def _review(self, kind, chain_id, resno, resname, detail=""):
        self.review.append(dict(kind=kind, chain=chain_id, resno=resno, resname=resname, detail=detail))

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

    # ---- evaluation helpers ----------------------------------------------------

    def _rot_window(self, st, chain_id, resno):
        path = write_window(st, chain_id, resno - 1, resno + 1, self._path("polish_window.pdb"))
        return rotalyze(path).get((chain_id, resno))

    def _rama_window(self, st, chain_id, resno):
        path = write_window(st, chain_id, resno - 3, resno + 3, self._path("polish_window.pdb"))
        return ramalyze(path), rotalyze(path)

    def _donors_for(self, chain_id, resno, resname, scores, favoured_only=False):
        """Other chains' copies of this residue, best MolProbity score first."""
        rank = {"Favored": 0, "Allowed": 1}
        donors = []
        for other in polymer_chains(self.st):
            if other == chain_id:
                continue
            entry = scores.get((other, resno))
            res = find_residue(self.st, other, resno)
            if entry is None or res is None or res.name != resname:
                continue
            if entry["evaluation"] not in rank:
                continue
            if favoured_only and entry["evaluation"] != "Favored":
                continue
            donors.append((rank[entry["evaluation"]], -entry["score"], other))
        donors.sort()
        return [d[2] for d in donors[:NCS_DONORS]]

    # ---- stage 1: rotamers --------------------------------------------------

    def stage_rotamers(self):
        rot0 = rotalyze(self.cur)
        targets = sorted(k for k, v in rot0.items()
                         if v["evaluation"] == "OUTLIER" and self._selected(k[0]))
        self.say("{} rotamer outliers".format(len(targets)))
        fixed = 0
        for index, (chain_id, resno) in enumerate(targets, 1):
            res0 = find_residue(self.st, chain_id, resno)
            if res0 is None:
                continue
            resname = res0.name
            d0 = sidechain_density(self.grid, res0)
            self.say("[{}/{}] {} {}{} density {:.2f}".format(index, len(targets), resname, chain_id, resno, d0))
            trials = []
            if self.opt.use_ncs:
                for donor in self._donors_for(chain_id, resno, resname, rot0):
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
                    outcome = self._rotamer_trial(chain_id, resno, how)
                except PolishError:
                    raise
                except Exception as exc:
                    self.say("   {}: failed ({}: {})".format(name, type(exc).__name__, exc))
                    continue
                if outcome is None:
                    break           # no such library rotamer: nothing further to try
                trial_st, d1, new_contacts, evaluation = outcome
                tried.append("{}:{:.2f}:{}".format(name, d1, evaluation["evaluation"] if evaluation else "?"))
                ok = (evaluation is not None and evaluation["evaluation"] != "OUTLIER"
                      and not new_contacts
                      and density_ok(d0, d1, self.opt.guard, self.opt.min_density))
                if ok:
                    copy_positions(find_residue(trial_st, chain_id, resno), res0)
                    accepted = (name, d1, evaluation["rotamer"])
                    break
            if accepted:
                fixed += 1
                self._record("rotamers", chain_id, resno, resname, accepted[0], d0, accepted[1],
                             "now {}".format(accepted[2]))
                self._commit("   fixed {}{} {} via {} ({:.2f} -> {:.2f})".format(
                    chain_id, resno, resname, accepted[0], d0, accepted[1]))
            else:
                self._review("rotamer outlier kept", chain_id, resno, resname,
                             "density {:.2f}; tried {}".format(d0, " ".join(tried) or "nothing"))
                self.say("   left {}{} {}: {}".format(chain_id, resno, resname, " ".join(tried)))
                self._snapshot()
        return len(targets), fixed

    def _rotamer_trial(self, chain_id, resno, how):
        """Build one candidate for a side chain and score it.

        Returns (structure, density, new contacts, rotalyze entry), or None when
        the requested library rotamer does not exist.
        """
        start = self.cur
        if "donor" in how:
            st = read_structure(self.cur)
            if transplant_sidechain(st, chain_id, resno, how["donor"]) is None:
                return (st, float("nan"), set(), None)
            start = self._path("polish_trial_in.pdb")
            st.write_pdb(start)
        coot = self._coot(start)
        if how.get("autofit"):
            coot.auto_fit(chain_id, resno)
        if "rotamer" in how:
            if not coot.set_rotamer(chain_id, resno, how["rotamer"]):
                return None
        if "weight" in how:
            coot.use_torsion_restraints(how["weight"])
        coot.refine_triple(chain_id, resno, self.opt.cycles)
        path = coot.write(self._path("polish_trial.pdb"))
        st = read_structure(path)
        res = find_residue(st, chain_id, resno)
        d1 = sidechain_density(self.grid, res)
        new_contacts = clash_contacts(st, chain_id, resno) - clash_contacts(self.st, chain_id, resno)
        evaluation = self._rot_window(st, chain_id, resno)
        return st, d1, new_contacts, evaluation

    # ---- stage 2: Ramachandran ------------------------------------------------

    def stage_rama(self):
        ram0 = ramalyze(self.cur)
        rot_now = rotalyze(self.cur)
        targets = sorted(k for k, v in ram0.items()
                         if v["evaluation"] == "OUTLIER" and self._selected(k[0]))
        self.say("{} Ramachandran outliers".format(len(targets)))
        fixed = 0
        for index, (chain_id, resno) in enumerate(targets, 1):
            res0 = find_residue(self.st, chain_id, resno)
            if res0 is None:
                continue
            resname = res0.name
            window = range(resno - 2, resno + 3)
            d0 = self._window_density(self.st, chain_id, window)
            base = self._window_contacts(self.st, chain_id, window)
            self.say("[{}/{}] {} {}{} backbone density {:.2f}".format(
                index, len(targets), resname, chain_id, resno, d0))
            candidates = []
            if self.opt.use_ncs:
                for donor in self._donors_for(chain_id, resno, resname, ram0, favoured_only=True):
                    candidates.append(("ncs_" + donor, dict(donor=donor)))
            candidates.append(("rama_refine", dict(rama=True)))
            if self.opt.do_flips:
                candidates.append(("flip_before", dict(flip=resno - 1)))
                candidates.append(("flip_after", dict(flip=resno)))
            scored = []
            for name, how in candidates:
                try:
                    outcome = self._rama_trial(chain_id, resno, how)
                except PolishError:
                    raise
                except Exception as exc:
                    self.say("   {}: failed ({}: {})".format(name, type(exc).__name__, exc))
                    continue
                if outcome is None:
                    continue
                st, ram, rot = outcome
                outliers = [m for m in window if ram.get((chain_id, m), {}).get("evaluation") == "OUTLIER"]
                favoured = sum(1 for m in window if ram.get((chain_id, m), {}).get("evaluation") == "Favored")
                new_rot = [m for m in window
                           if rot.get((chain_id, m), {}).get("evaluation") == "OUTLIER"
                           and rot_now.get((chain_id, m), {}).get("evaluation") != "OUTLIER"]
                d1 = self._window_density(st, chain_id, window)
                new_contacts = self._window_contacts(st, chain_id, window) - base
                ok = not outliers and not new_rot and not new_contacts and d1 >= d0 - ABSOLUTE_SLACK
                self.say("   {}: outliers {} favoured {} new rotamer outliers {} density {:.2f} contacts {} -> {}".format(
                    name, len(outliers), favoured, len(new_rot), d1, len(new_contacts), "ok" if ok else "no"))
                if ok:
                    scored.append(((favoured, d1), name, st, d1))
            if scored:
                scored.sort(key=lambda item: item[0], reverse=True)
                _, name, st, d1 = scored[0]
                for m in window:
                    src, dst = find_residue(st, chain_id, m), find_residue(self.st, chain_id, m)
                    if src is not None and dst is not None:
                        copy_positions(src, dst)
                fixed += 1
                self._record("ramachandran", chain_id, resno, resname, name, d0, d1)
                self._commit("   fixed {}{} {} via {} ({:.2f} -> {:.2f})".format(
                    chain_id, resno, resname, name, d0, d1))
            else:
                self._review("Ramachandran outlier kept", chain_id, resno, resname,
                             "backbone density {:.2f}; no candidate removed it without cost".format(d0))
                self._snapshot()
        return len(targets), fixed

    def _window_density(self, st, chain_id, window):
        values = []
        for m in window:
            res = find_residue(st, chain_id, m)
            if res is not None:
                values.append(backbone_density(self.grid, res))
        values = [v for v in values if v == v]
        return sum(values) / len(values) if values else float("nan")

    def _window_contacts(self, st, chain_id, window):
        found = set()
        for m in window:
            for contact in clash_contacts(st, chain_id, m):
                found.add((m,) + contact)
        return found

    def _rama_trial(self, chain_id, resno, how):
        start = self.cur
        if "donor" in how:
            st = read_structure(self.cur)
            if transplant_backbone_window(st, chain_id, resno - 1, resno + 1, how["donor"]) is None:
                return None
            start = self._path("polish_trial_in.pdb")
            st.write_pdb(start)
        coot = self._coot(start)
        if how.get("rama"):
            coot.use_rama_restraints(RAMA_WEIGHT)
        if "flip" in how:
            if not coot.flip_peptide(chain_id, how["flip"]):
                return None
        coot.refine_range(chain_id, resno - 2, resno + 2, max(50, self.opt.cycles // 2))
        path = coot.write(self._path("polish_trial.pdb"))
        st = read_structure(path)
        ram, rot = self._rama_window(st, chain_id, resno)
        return st, ram, rot

    # ---- stage 3: clashes -------------------------------------------------------

    def stage_clashes(self):
        score, clashes = clashscore(self.cur)
        rot_now = rotalyze(self.cur)
        ram_now = ramalyze(self.cur)
        seen, targets = set(), []
        for clash in clashes:
            if clash["overlap"] < self.opt.min_overlap:
                continue
            pair = clash_pair(clash)
            if pair in seen or not (self._selected(pair[0][0]) or self._selected(pair[1][0])):
                continue
            seen.add(pair)
            partners = []
            if clash["atom1"] not in BACKBONE_H and clash["resname1"] != "GLY" and self._selected(clash["chain1"]):
                partners.append((clash["chain1"], clash["resno1"]))
            if clash["atom2"] not in BACKBONE_H and clash["resname2"] != "GLY" and self._selected(clash["chain2"]):
                partners.append((clash["chain2"], clash["resno2"]))
            targets.append((pair, partners, clash))
        self.say("clashscore {}; {} pairs with overlap >= {}".format(score, len(targets), self.opt.min_overlap))
        fixed = 0
        for index, (pair, partners, clash) in enumerate(targets, 1):
            label = "{}{}-{}{}".format(pair[0][0], pair[0][1], pair[1][0], pair[1][1])
            if not partners:
                self._review("clash left for refinement", pair[0][0], pair[0][1], clash["resname1"],
                             "backbone-backbone overlap {:.2f} with {}{}".format(
                                 clash["overlap"], pair[1][0], pair[1][1]))
                continue
            self.say("[{}/{}] {} overlap {:.2f}".format(index, len(targets), label, clash["overlap"]))
            done = False
            for chain_id, resno in partners:
                res0 = find_residue(self.st, chain_id, resno)
                if res0 is None:
                    continue
                d0 = sidechain_density(self.grid, res0)
                trials = [("autofit", dict(autofit=True))]
                trials += [("rotamer{}".format(k), dict(rotamer=k, weight=TORSION_WEIGHTS[0]))
                           for k in range(self.opt.max_clash_rotamers)]
                for name, how in trials:
                    try:
                        outcome = self._rotamer_trial(chain_id, resno, how)
                    except PolishError:
                        raise
                    except Exception as exc:
                        self.say("   {}: failed ({}: {})".format(name, type(exc).__name__, exc))
                        continue
                    if outcome is None:
                        break
                    st, d1, new_contacts, evaluation = outcome
                    if not density_ok(d0, d1, self.opt.guard, self.opt.min_density) or new_contacts:
                        continue
                    if evaluation is None or evaluation["evaluation"] == "OUTLIER":
                        continue
                    if self._pair_overlap(st, pair) >= self.opt.min_overlap:
                        continue
                    if self._new_outliers(st, chain_id, resno, rot_now, ram_now):
                        continue
                    copy_positions(find_residue(st, chain_id, resno), res0)
                    fixed += 1
                    self._record("clashes", chain_id, resno, res0.name, name, d0, d1,
                                 "overlap {:.2f} with {}{}".format(clash["overlap"], *(
                                     pair[1] if pair[0] == (chain_id, resno) else pair[0])))
                    self._commit("   fixed {} via {}{} {} ({:.2f} -> {:.2f})".format(
                        label, chain_id, resno, name, d0, d1))
                    done = True
                    break
                if done:
                    break
            if not done:
                self._review("clash kept", pair[0][0], pair[0][1], clash["resname1"],
                             "overlap {:.2f} with {}{}; no side-chain fix kept the density".format(
                                 clash["overlap"], pair[1][0], pair[1][1]))
                self.say("   left {}".format(label))
                self._snapshot()
        return len(targets), fixed

    def _pair_overlap(self, st, pair):
        path = write_neighbourhood(st, list(pair), self._path("polish_neighbourhood.pdb"))
        _, clashes = clashscore(path)
        worst = 0.0
        for clash in clashes:
            if clash_pair(clash) == pair:
                worst = max(worst, clash["overlap"])
        return worst

    def _new_outliers(self, st, chain_id, resno, rot_now, ram_now):
        ram, rot = self._rama_window(st, chain_id, resno)
        for m in range(resno - 1, resno + 2):
            if (rot.get((chain_id, m), {}).get("evaluation") == "OUTLIER"
                    and rot_now.get((chain_id, m), {}).get("evaluation") != "OUTLIER"):
                return True
            if (ram.get((chain_id, m), {}).get("evaluation") == "OUTLIER"
                    and ram_now.get((chain_id, m), {}).get("evaluation") != "OUTLIER"):
                return True
        return False

    # ---- wrap-up --------------------------------------------------------------

    def finish_review(self):
        rot, ram = rotalyze(self.cur), ramalyze(self.cur)
        listed = {(item["chain"], item["resno"], item["kind"]) for item in self.review}
        for (chain_id, resno), v in sorted(rot.items()):
            if v["evaluation"] == "OUTLIER" and (chain_id, resno, "rotamer outlier kept") not in listed:
                self._review("rotamer outlier", chain_id, resno, v["resname"], "still an outlier after polish")
        for (chain_id, resno), v in sorted(ram.items()):
            if v["evaluation"] == "OUTLIER" and (chain_id, resno, "Ramachandran outlier kept") not in listed:
                self._review("Ramachandran outlier", chain_id, resno, v["resname"], "still an outlier after polish")

    def write_changes_csv(self, path, original_path):
        """Every residue whose atoms moved, with density and MolProbity status before and after."""
        st0, st1 = read_structure(original_path), self.st
        rot0, rot1 = rotalyze(self._path("polish_input.pdb")), rotalyze(self.cur)
        ram0, ram1 = ramalyze(self._path("polish_input.pdb")), ramalyze(self.cur)
        by_key = {}
        for change in self.changes:
            by_key.setdefault((change["chain"], change["resno"]), []).append(
                "{}:{}".format(change["stage"], change["action"]))
        rows = changed_residue_rows(st0, st1, self.grid, self.diff_grid, rot0, rot1, ram0, ram1, by_key)
        write_changes(path, rows)
        return rows


def changed_residue_rows(st0, st1, grid, diff_grid, rot0, rot1, ram0, ram1, actions):
    rows = []
    for chain in st0[0]:
        for r0 in chain:
            r1 = find_residue(st1, chain.name, r0.seqid.num)
            if r1 is None or r0.name != r1.name:
                continue
            shift = 0.0
            for atom in r0:
                other = r1.find_atom(atom.name, "*")
                if other is not None:
                    shift = max(shift, atom.pos.dist(other.pos))
            if shift < 0.05:
                continue
            key = (chain.name, r0.seqid.num)
            row = dict(chain=chain.name, resno=r0.seqid.num, resname=r0.name, max_shift=round(shift, 2),
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


CHANGES_COLUMNS = ("chain", "resno", "resname", "max_shift", "action", "density_before", "density_after",
                   "diff_density_before", "diff_density_after", "rotamer_before", "rotamer_after",
                   "rama_before", "rama_after")


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


def _fmt(summary):
    return ("rotamer outliers {rot_outliers} ({rot_pct}%), Ramachandran outliers {rama_outliers} "
            "({rama_pct}%), favoured {rama_favored_pct}%, clashscore {clashscore}").format(**summary)
