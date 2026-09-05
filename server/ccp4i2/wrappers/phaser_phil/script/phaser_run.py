"""Running Phaser from a working phil, and recording what it says.

Phaser's driver turns the working phil into the phaser.Input* object; we run
the mode function ourselves so that our own callback object hears the run.
The recorder writes what it hears into program.xml as it happens, and after
the run the Result object -- not the log -- supplies the solutions.

Nothing here is inferred from prose. The module timeline comes from the
MODE callbacks, progress from the progress-bar callbacks, the verdict from
the `result` callback, solutions from the mr_solution object, and the
search-strategy narrative from the fixed control sentences Phaser prints in
its AUTOMATED MOLECULAR REPLACEMENT summary blocks -- with a count of the
blocks that matched none of them, so a change of wording shows rather than
misleads.

Phaser and libtbx are imported inside functions: this module is imported by
the task registry on a CCP4-free server.
"""
import os
import re
import time

from lxml import etree

AUTO = "AUTOMATED MOLECULAR REPLACEMENT"
ROTATION = "MOLECULAR REPLACEMENT ROTATION FUNCTION"
TRANSLATION = "MOLECULAR REPLACEMENT TRANSLATION FUNCTION"


class PhaserRecorder:
    """Phaser's callback object, writing program.xml as the run proceeds.

    `flush` is called with the XML root whenever something worth showing
    changed; progress increments are rate-limited to `min_interval` seconds
    because a rotation function increments hundreds of times a second.
    """

    def __init__(self, xmlroot, flush=None, min_interval=1.0, clock=time.monotonic):
        self.xmlroot = xmlroot
        self.flush = flush
        self.min_interval = min_interval
        self.clock = clock
        self.t0 = clock()
        self._last_flush = None
        self.modules = etree.SubElement(xmlroot, "Modules")
        self.attempts = etree.SubElement(xmlroot, "Attempts")
        self._attempt = None
        self._placed = 0

    # -- what Phaser calls -------------------------------------------------
    def startProgressBar(self, label, size):
        node = self._activity()
        node.find("label").text = str(label)
        node.find("max").text = str(size)
        node.find("value").text = "0"
        self._notify(force=True)

    def incrementProgressBar(self):
        node = self._activity()
        value = node.find("value")
        value.text = str(int(value.text or 0) + 1)
        self._notify(force=False)

    def endProgressBar(self):
        pass

    def warn(self, message):
        message = str(message).strip()
        if not message:
            return
        warnings = self._child("PhaserWarnings")
        etree.SubElement(warnings, "Warning").text = message
        if self._attempt is not None:
            etree.SubElement(self._attempt, "Warning").text = message
        self._notify(force=True)

    def loggraph(self, title, data):
        try:
            from ccp4i2.pimple.logtable import CCP4LogToEtree
            self.xmlroot.append(CCP4LogToEtree(data))
        except Exception:
            return
        self._notify(force=True)

    def call_back(self, message, data):
        if message == "MODE":
            self._module(str(data).strip())
        elif message == "current best solution":
            self._best(str(data))
        elif message == "result":
            self._child("Verdict").text = str(data).strip()
        elif message in ("phaser_start", "phaser_end"):
            self._child("Run").set(message.split("_")[1], f"{self.clock() - self.t0:.1f}")
        else:
            return
        self._notify(force=True)

    # -- what it records ---------------------------------------------------
    def _module(self, name):
        node = etree.SubElement(self.modules, "Module")
        node.set("name", name)
        node.set("t", f"{self.clock() - self.t0:.1f}")
        if name == ROTATION:
            # A rotation function opens a search attempt for the next
            # component; a placement (below) closes it
            self._attempt = etree.SubElement(self.attempts, "Attempt")
            self._attempt.set("component", str(self._placed + 1))
            self._attempt.set("number", str(1 + sum(
                1 for a in self.attempts if a.get("component") == str(self._placed + 1)) - 1))
            self._attempt.set("outcome", "running")
        self._child("CurrentModule").text = name

    def _best(self, text):
        for old in self.xmlroot.findall("PhaserCurrentBestSolution"):
            self.xmlroot.remove(old)
        node = etree.SubElement(self.xmlroot, "PhaserCurrentBestSolution")
        node.text = text
        if self._attempt is not None:
            self._attempt.set("outcome", "placed")
            self._attempt = None
        self._placed += 1

    def _activity(self):
        node = self.xmlroot.find("CurrentActivity")
        if node is None:
            node = etree.SubElement(self.xmlroot, "CurrentActivity")
            for tag in ("label", "max", "value"):
                etree.SubElement(node, tag).text = "0"
        return node

    def _child(self, tag):
        node = self.xmlroot.find(tag)
        if node is None:
            node = etree.SubElement(self.xmlroot, tag)
        return node

    def _notify(self, force):
        if self.flush is None:
            return
        now = self.clock()
        if not force and self._last_flush is not None and now - self._last_flush < self.min_interval:
            return
        self._last_flush = now
        self.flush(self.xmlroot)

    def finish(self):
        """The run is over: no attempt is still running.

        What became of an attempt the callbacks never closed -- a model that
        was already at the origin is placed without a `current best solution`
        message -- is not known here; the post-run strategy narrative says.
        """
        for attempt in self.attempts:
            if attempt.get("outcome") == "running":
                attempt.set("outcome", "ended")
        self._attempt = None


# ---------------------------------------------------------------------------
# After the run
# ---------------------------------------------------------------------------

_MODULE_LINE = re.compile(r"Phaser Module:\s*(.*?)\s+\d+\.\d+(?:\.\d+)?\s*\**\s*$")


def summary_blocks(text):
    """r.summary() as [(module name, block text)] in order."""
    blocks = []
    name, lines = None, []
    for line in text.splitlines():
        m = _MODULE_LINE.search(line)
        if m:
            if name is not None:
                blocks.append((name, "\n".join(lines)))
            name, lines = m.group(1).strip(), []
        elif name is not None:
            lines.append(line)
    if name is not None:
        blocks.append((name, "\n".join(lines)))
    return blocks


# The control sentences an AUTO block may carry. Anything else in an AUTO
# block is counted as unparsed, never guessed at.
_SEARCH_ORDER = "Search Order (next search *)"
_NEW_BEST = re.compile(r"New Best LLG\s*=\s*([-\d.]+)\s*\(resolution\s*=\s*([\d.]+)\)")
_LOWERED = re.compile(r"High resolution limit lowered by expected LLG\s*=\s*([\d.]+)")
_UNALTERED = "High resolution limit unaltered by expected LLG"
_NO_LOWER = "No solutions found at lower resolution"
_STORED = "Solution stored for possible reversion if next search fails"
_NSOL = re.compile(r"Number of solutions\s*=\s*(\d+)")
_AT_ORIGIN = "Molecular replacement is not necessary and is aborted"
_GIVEN_UP = "End of possible ensembles for this component"
# Blocks that inform rather than decide. Recognising them keeps the unparsed
# count honest: it should rise only for wording this code has never seen.
_INFORMATIONAL = (
    "Steps:", "Input Search Order", "Number of search ensembles", "Z-score test",
    "Composition Table", "Solutions will be refined", "Resolution for refinement",
    "Refinement may have introduced clashes", "$TEXT:MR Result", "Search: Next component",
    "Current is Best Solution", "No more components", "End with solutions",
    "There was only 1 possible ensemble", "Increase resolution",
)


def strategy_attempts(blocks):
    """The search attempts Phaser narrates, from its own control sentences.

    Returns (attempts, unparsed): attempts as dicts with component,
    resolution, outcome and llg; unparsed the number of AUTO blocks that
    said nothing this function recognises.
    """
    attempts, unparsed = [], 0
    current = None
    for name, text in blocks:
        if name == AUTO:
            recognised = False
            if _SEARCH_ORDER in text:
                recognised = True
                current = {"component": _next_component(text), "resolution": None,
                           "outcome": "searching", "llg": None}
                attempts.append(current)
            m = _NEW_BEST.search(text)
            if m and current is not None:
                recognised = True
                current["outcome"] = "placed"
                current["llg"] = float(m.group(1))
                current = None
            if _NO_LOWER in text and current is not None:
                recognised = True
                current["outcome"] = "no solution at lower resolution"
                current = None
            if _AT_ORIGIN in text and current is not None:
                # The input model already sits at the origin: placed without a search
                recognised = True
                current["outcome"] = "placed (model already at the origin, no search needed)"
                current = None
            if _GIVEN_UP in text and current is not None:
                recognised = True
                current["outcome"] = current["outcome"] if current["outcome"] != "searching" \
                    else "no definite solution"
                current = None
            if _STORED in text or _NSOL.search(text) or any(w in text for w in _INFORMATIONAL):
                recognised = True
            if not recognised and text.strip():
                unparsed += 1
        elif name == ROTATION and current is not None:
            m = _LOWERED.search(text)
            if m:
                current["resolution"] = float(m.group(1))
            elif _UNALTERED in text:
                current["resolution"] = "full"
        elif name == TRANSLATION and current is not None:
            if "No Signal in Translation Function" in text:
                current["outcome"] = "no signal in translation function"
            elif "NOT over Z-score cutoff" in text:
                current["outcome"] = "no definite solution"
    return attempts, unparsed


def _next_component(text):
    """The '*' line under 'Search Order (next search *)': '#2  model *'."""
    for line in text.splitlines():
        if "*" in line and "#" in line and "Search Order" not in line:
            words = [w for w in line.replace("*", " ").split() if not w.startswith("#")]
            return words[0] if words else None
    return None


# The annotation grammar, as documented at
# phenix-online.org/documentation/reference/phaser.html
_TOKEN_MEANING = {
    "RF*0": "rotation 000 by the input model's low R-factor",
    "TF*0": "translation 000 by the input model's low R-factor",
    "RF++": "rotation from a previous strong solution reused",
    "TFZ=*": "first molecule in P1, no translation function",
    "+TNCS": "added by a tNCS relation",
}


def annotation_tokens(annotation):
    """The per-component records in a solution annotation.

    Returns (components, unknown): each component a dict with the values
    it carried (RFZ, TFZ, TFZeq, PAK, LLG -- LLG last of possibly several,
    the documented low-then-high resolution refinement), plus notes for
    the flag tokens and amalgamation groups; unknown lists tokens the
    documented grammar does not cover.
    """
    components, unknown = [], []
    current = None
    depth = 0

    def new():
        nonlocal current
        current = {"notes": [], "amalgamated": []}
        components.append(current)

    for tok in annotation.replace("(", " ( ").replace(")", " ) ").split():
        if tok == "(":
            depth += 1
            continue
        if tok == ")":
            depth = max(0, depth - 1)
            continue
        if depth:
            if current is not None:
                current["amalgamated"].append(tok)
            continue
        if tok.startswith("RFZ=") or tok in ("RF*0", "RF++"):
            new()
            if tok.startswith("RFZ="):
                current["RFZ"] = _num(tok[4:])
            else:
                current["notes"].append(_TOKEN_MEANING[tok])
        elif tok.startswith("TFZ=="):
            if current is None:
                new()
            current["TFZeq"] = _num(tok[5:])
        elif tok == "TFZ=*" or tok.startswith("TFZ="):
            if current is None:
                new()
            if tok == "TFZ=*":
                current["notes"].append(_TOKEN_MEANING[tok])
            else:
                current["TFZ"] = _num(tok[4:])
        elif tok == "TF*0":
            if current is None:
                new()
            current["notes"].append(_TOKEN_MEANING[tok])
        elif tok.startswith("PAK="):
            if current is None:
                new()
            current["PAK"] = _num(tok[4:])
        elif tok.startswith("LLG+="):
            if current is not None:
                current["notes"].append("LLG during amalgamation")
        elif tok.startswith("LLG="):
            if current is None:
                new()
            current["LLG"] = _num(tok[4:])
        elif tok == "+TNCS":
            if current is not None:
                current["notes"].append(_TOKEN_MEANING[tok])
        elif tok.startswith("*T="):
            if current is not None:
                current["notes"].append(f"matches template solution {tok[3:]}")
        else:
            unknown.append(tok)
    return components, unknown


def _num(text):
    try:
        return float(text)
    except ValueError:
        return text


def solutions_xml(result, parent):
    """The solutions from the mr_solution object, as <Solutions><Solution/>.

    Typed fields are read as fields; only the annotation is tokenised.
    """
    solutions = etree.SubElement(parent, "Solutions")
    for i, sol_set in enumerate(result.getDotSol()):
        node = etree.SubElement(solutions, "Solution")
        etree.SubElement(node, "Number").text = str(i + 1)
        etree.SubElement(node, "Annotation").text = str(sol_set.ANNOTATION).strip()
        for field in ("LLG", "TFZ", "TFZeq", "R", "PAK", "NUM", "ORIG_LLG", "ORIG_R"):
            value = getattr(sol_set, field, None)
            if value is not None:
                etree.SubElement(node, field).text = f"{value:.2f}" if isinstance(value, float) else str(value)
        etree.SubElement(node, "spaceGroup").text = str(sol_set.getSpaceGroupName()).strip()
        for line in sol_set.unparse().splitlines():
            if line.startswith("SOLU HISTORY"):
                etree.SubElement(node, "History").text = line[len("SOLU HISTORY"):].strip()
        components = etree.SubElement(node, "Components")
        for comp in sol_set.KNOWN:
            c = etree.SubElement(components, "Component")
            etree.SubElement(c, "Name").text = str(comp.MODLID)
            euler = comp.getEuler()
            etree.SubElement(c, "Euler").text = " ".join(f"{x:.3f}" for x in euler)
            etree.SubElement(c, "Frac").text = " ".join(f"{x:.4f}" for x in comp.TRA)
            etree.SubElement(c, "Bfac").text = f"{comp.BFAC:.3f}"
            etree.SubElement(c, "Mult").text = str(comp.MULT)
        tokens, unknown = annotation_tokens(str(sol_set.ANNOTATION))
        placements = etree.SubElement(node, "Placements")
        for t in tokens:
            p = etree.SubElement(placements, "Placement")
            for key in ("RFZ", "TFZ", "TFZeq", "PAK", "LLG"):
                if key in t:
                    etree.SubElement(p, key).text = str(t[key])
            for note in t["notes"]:
                etree.SubElement(p, "Note").text = note
            if t["amalgamated"]:
                etree.SubElement(p, "Amalgamated").text = " ".join(t["amalgamated"])
        if unknown:
            etree.SubElement(node, "UnknownTokens").text = " ".join(unknown)
    return solutions


def run_mode(master_phil, working_phil_path, mode, work_directory, recorder, log_path):
    """Run phaser.run<MODE> on the Input object Phaser's driver builds from
    the working phil, with `recorder` as the callback. Returns the Result."""
    import iotbx.phil
    import phaser
    from libtbx.phil import interface as phil_interface
    from phaser.phenix_interface import driver

    working = master_phil.fetch(sources=[iotbx.phil.parse(file_name=working_phil_path)])
    index = phil_interface.index(master_phil, working, parse=iotbx.phil.parse)
    with open(log_path, "a") as log:
        interpreter = driver.phaser_parameter_interpreter(index, work_directory, out=log)
        output = phaser.Output()
        output.setPackagePhenix(log)
        output.setPhenixCallback(recorder)
        run = getattr(phaser, f"run{mode}")
        result = run(interpreter.input, output)
    recorder.finish()
    return result


# ---------------------------------------------------------------------------
# Experimental phasing
# ---------------------------------------------------------------------------

SAD = "SINGLE-WAVELENGTH ANOMALOUS DISPERSION"
_EP_LLG = re.compile(r"^\s*(?:Final )?Log-Likelihood\s*=\s*([-\d.]+)")
_EP_NEW_ATOMS = re.compile(r"Number of new atoms identified this cycle\s*=\s*(\d+)")
_EP_DELETED = re.compile(r"Deleted Atom #'s:\s*(.*)$")
_EP_FOM_ALL = re.compile(r"^\s*ALL\s+[\d.]+-\s*[\d.]+(?:\s+\d+\s+[\d.]+){3}\s+(\d+)\s+([\d.]+)")
_EP_COMPLETION = re.compile(r"SUB-STRUCTURE COMPLETION #(\d+)")
_EP_ANALYSIS = re.compile(r"SUB-STRUCTURE ANALYSIS #(\d+)")


def ep_cycles(blocks):
    """The substructure-completion cycles Phaser narrates in its SAD block.

    Each cycle: the LLG after refinement, the atoms found and deleted, the
    overall FOM, and whether completion converged. The block's vocabulary
    is small and fixed; anything else is left alone.
    """
    cycles = []
    final = {"llg": None, "fom": None, "converged": None}
    for name, text in blocks:
        if name != SAD:
            continue
        current = None
        for line in text.splitlines():
            m = _EP_ANALYSIS.search(line)
            if m:
                # "ANALYSIS #n: LLG MAPS" then "ANALYSIS #n: WEAK SITES" --
                # one cycle, two headings
                number = int(m.group(1))
                if current is None or current["cycle"] != number:
                    current = {"cycle": number, "llg": None, "added": 0,
                               "deleted": [], "fom": None, "converged": None}
                    cycles.append(current)
                continue
            m = _EP_NEW_ATOMS.search(line)
            if m and current is not None:
                current["added"] = int(m.group(1))
                continue
            m = _EP_DELETED.search(line)
            if m and current is not None:
                current["deleted"] = [int(x) for x in m.group(1).split() if x.isdigit()]
                continue
            if "Substructure completion has NOT converged" in line and current is not None:
                current["converged"] = False
            elif "Substructure completion has converged" in line and current is not None:
                current["converged"] = True
            m = _EP_LLG.match(line)
            if m:
                value = float(m.group(1))
                if line.strip().startswith("Final"):
                    final["llg"] = value
                elif current is not None:
                    current["llg"] = value
                continue
            m = _EP_FOM_ALL.match(line)
            if m:
                fom = float(m.group(2))
                if current is not None:
                    current["fom"] = fom
                final["fom"] = fom
    final["converged"] = cycles[-1]["converged"] if cycles else None
    return cycles, final


def ep_results_xml(result, parent):
    """The hands, sites and figures of merit, from the ResultEP object."""
    hands = etree.SubElement(parent, "Hands")
    count = 2 if getattr(result, "second_hand", False) else 1
    for i in range(count):
        hand = result.getHand(i)
        node = etree.SubElement(hands, "Hand")
        etree.SubElement(node, "Number").text = str(i + 1)
        etree.SubElement(node, "PDB").text = str(hand.PDBfile)
        etree.SubElement(node, "MTZ").text = str(hand.MTZfile)
        etree.SubElement(node, "LLG").text = f"{hand.getLogLikelihood():.1f}"
        stats = hand.allStats
        try:
            fom = list(stats.FOM_bin)
            num = list(stats.NUM_bin)
            hires = list(stats.HiRes_bin)
            lores = list(stats.LoRes_bin)
            total_num = sum(num)
            etree.SubElement(node, "FOM").text = f"{sum(fom) / total_num:.3f}" if total_num else ""
            bins = etree.SubElement(node, "FOMByResolution")
            for f, n, hi, lo in zip(fom, num, hires, lores):
                b = etree.SubElement(bins, "Bin")
                etree.SubElement(b, "LoRes").text = f"{lo:.2f}"
                etree.SubElement(b, "HiRes").text = f"{hi:.2f}"
                etree.SubElement(b, "Number").text = str(int(n))
                etree.SubElement(b, "FOM").text = f"{f / n:.3f}" if n else ""
        except Exception:
            pass
    overall = etree.SubElement(parent, "Overall")
    for name in ("stats_fom", "stats_acentric_fom", "stats_centric_fom", "stats_num",
                 "stats_hires", "stats_lores"):
        fn = getattr(result, name, None)
        if fn is None:
            continue
        try:
            value = fn()
            etree.SubElement(overall, name.replace("stats_", "")).text = (
                f"{value:.3f}" if isinstance(value, float) else str(value))
        except Exception:
            pass
    sites = etree.SubElement(parent, "Sites")
    try:
        atoms = list(result.getTopAtoms())
    except Exception:
        atoms = []
    for i, atom in enumerate(atoms):
        s = etree.SubElement(sites, "Site")
        etree.SubElement(s, "Number").text = str(i + 1)
        etree.SubElement(s, "Element").text = str(atom.scattering_type)
        x, y, z = atom.site
        etree.SubElement(s, "Frac").text = f"{x:.3f} {y:.3f} {z:.3f}"
        etree.SubElement(s, "Occupancy").text = f"{atom.occupancy:.2f}"
        b_iso = atom.b_iso() if callable(getattr(atom, "b_iso", None)) else getattr(atom, "b_iso", None)
        if b_iso is not None and atom.u_iso >= 0:
            etree.SubElement(s, "B").text = f"{b_iso:.1f}"
        else:
            etree.SubElement(s, "B").text = "aniso"
        etree.SubElement(s, "History").text = str(atom.label).strip()
    return hands


def ep_cycles_xml(cycles, final, parent):
    node = etree.SubElement(parent, "Completion")
    for key, value in final.items():
        if value is not None:
            node.set(key, str(value))
    for cycle in cycles:
        c = etree.SubElement(node, "Cycle")
        for key, value in cycle.items():
            if value is None:
                continue
            c.set(key, " ".join(str(v) for v in value) if isinstance(value, list) else str(value))
    return node
