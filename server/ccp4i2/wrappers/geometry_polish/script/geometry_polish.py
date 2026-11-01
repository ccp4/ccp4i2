"""Polish geometry for deposition: an NCS-aware MolProbity repair.

Takes a refined model and its maps, works through the MolProbity rotamer
outliers, Ramachandran outliers and side-chain clashes, and accepts only the
changes the density supports. Everything else goes on a review list. The
science lives in polish_lib; this plugin wires it to the job, writes the
program XML the report reads, and registers the outputs.
"""
import os
import pathlib
import sys

from lxml import etree

from ccp4i2.core import CCP4ErrorHandling, CCP4Utils
from ccp4i2.core.CCP4ModelData import CPdbDataFile
from ccp4i2.core.CCP4PluginScript import CPluginScript


def _value(param, default, cast):
    """A control parameter as a plain value, or ``default`` when it is unset.

    ``param.get()`` returns the coerced primitive (the idiom the other
    wrappers use for CBoolean); ``bool(param)`` would answer "set and true".
    """
    try:
        raw = param.get()
    except Exception:
        raw = None
    if raw is None or raw == "":
        return default
    try:
        return cast(raw)
    except (TypeError, ValueError):
        return default


def _chains(text):
    """'A, B,C' -> ['A', 'B', 'C']; blank -> None (every chain)."""
    if not text:
        return None
    chains = [c.strip() for c in str(text).replace(";", ",").split(",") if c.strip()]
    return chains or None


class geometry_polish(CPluginScript):
    TASKNAME = "geometry_polish"
    # Driven in-process through their CLI; declared so a missing MolProbity is
    # reported before the job starts rather than minutes into it.
    AUXILIARY_PROGRAMS = ("molprobity.rotalyze", "molprobity.ramalyze", "molprobity.clashscore")
    PROGRAM_LOCATED_BY = "the CCP4 environment (MolProbity ships with CCP4)"
    WHATNEXT = ["prosmart_refmac", "servalcat_pipe", "coot1"]
    ERROR_CODES = {
        201: {"description": "The polish could not run"},
        202: {"description": "The polished model was not written"},
    }

    def startProcess(self):
        # Execution-time imports: gemmi, numpy and Coot's API live on the worker only.
        from . import polish_lib

        inp = self.container.inputData
        out = self.container.outputData
        par = self.container.controlParameters

        xyzin = str(inp.XYZIN.fullPath)
        if inp.XYZIN.isMMCIF():
            path = pathlib.Path(str(out.XYZOUT.fullPath))
            out.XYZOUT.setFullPath(str(path.with_suffix(".cif")))
            out.XYZOUT.contentFlag.set(CPdbDataFile.CONTENT_FLAG_MMCIF)
        xyzout = str(out.XYZOUT.fullPath)
        workdir = self.getWorkDirectory()
        changes_csv = os.path.join(workdir, "geometry_polish_changes.csv")
        out.CHANGES_CSV.setFullPath(changes_csv)

        options = polish_lib.Options(
            do_rotamers=_value(par.DO_ROTAMERS, True, bool),
            do_rama=_value(par.DO_RAMA, True, bool),
            do_flips=_value(par.DO_PEPTIDE_FLIPS, True, bool),
            do_clashes=_value(par.DO_CLASHES, True, bool),
            use_ncs=_value(par.USE_NCS, True, bool),
            guard=_value(par.DENSITY_GUARD, 0.6, float),
            min_density=_value(par.MIN_DENSITY, 0.8, float),
            cycles=_value(par.REFINE_CYCLES, 200, int),
            max_clash_rotamers=_value(par.MAX_CLASH_ROTAMERS, 6, int),
            min_overlap=_value(par.MIN_CLASH_OVERLAP, 0.8, float),
            chains=_chains(_value(par.CHAIN_SELECTION, "", str)),
        )
        diff_mtz = str(inp.DIFFPHIOUT.fullPath) if inp.DIFFPHIOUT.isSet() else None
        dict_path = str(inp.DICT.fullPath) if inp.DICT.isSet() else None

        polisher = polish_lib.Polisher(
            xyzin, str(inp.FPHIOUT.fullPath), workdir, options=options, diff_mtz=diff_mtz,
            dict_path=dict_path, progress=self._say, snapshot=self._write_xml)

        # Coot writes its scratch files in the working directory.
        original_cwd = os.getcwd()
        os.chdir(workdir)
        try:
            polisher.run(xyzout)
            polisher.write_changes_csv(changes_csv, xyzin)
        except polish_lib.PolishError as exc:
            self._write_xml(polisher.state())
            return self._fail(201, "Geometry polish could not run: {}".format(exc))
        except Exception as exc:
            self._write_xml(polisher.state())
            return self._fail(201, "Geometry polish failed ({}: {})".format(type(exc).__name__, exc))
        finally:
            os.chdir(original_cwd)

        self._write_xml(polisher.state())
        if not os.path.exists(xyzout):
            return self._fail(202, "The polished model {} was not written".format(xyzout))
        out.XYZOUT.annotation.set("Polished model: {} changes, {} items for review".format(
            len(polisher.changes), len(polisher.review)))
        out.CHANGES_CSV.annotation.set("Per-residue change log of the geometry polish")
        return CPluginScript.SUCCEEDED

    # ---- plumbing -----------------------------------------------------------

    def _say(self, text):
        print(text)
        sys.stdout.flush()

    def _fail(self, code, message):
        print("ERROR:", message)
        sys.stdout.flush()
        self.appendErrorReport(code, message, severity=CCP4ErrorHandling.SEVERITY_ERROR)
        return CPluginScript.FAILED

    def _write_xml(self, state):
        """Program XML for the (running) report: summaries, stages, changes, review."""
        root = etree.Element("geometry_polish")
        summary = etree.SubElement(root, "summary")
        for tag in ("before", "after"):
            numbers = state.get(tag)
            if numbers:
                node = etree.SubElement(summary, tag)
                for key, value in numbers.items():
                    node.set(key, "" if value is None else str(value))
        stages = etree.SubElement(root, "stages")
        for stage in state.get("stages", []):
            node = etree.SubElement(stages, "stage")
            for key, value in stage.items():
                node.set(key, str(value))
        changes = etree.SubElement(root, "changes")
        for change in state.get("changes", []):
            node = etree.SubElement(changes, "change")
            for key, value in change.items():
                node.set(key, str(value))
        review = etree.SubElement(root, "review")
        for item in state.get("review", []):
            node = etree.SubElement(review, "item")
            for key, value in item.items():
                node.set(key, str(value))
        with open(self.makeFileName("PROGRAMXML"), "w") as handle:
            CCP4Utils.writeXML(handle, etree.tostring(root, pretty_print=True))
