"""Automated molecular replacement with Phaser, driven from its PHIL.

The parameters are Phaser's own (phenix_interface/__init__.params), read at
run time and filtered to MR_AUTO. The typed inputs -- reflections, search
models, composition -- come from CCP4i2's data model through shims, and
Phaser's driver builds the phaser.Input object from the resulting working
phil, so the translation from parameters to Phaser's API is Phaser's, not
ours. Phaser runs in this process, as the classic wrapper ran it, with a
recorder as its callback so that program.xml is written as it goes.
"""
import os
import pickle

from lxml import etree

from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.utils.phil_shims import FixedPhilShim
from ccp4i2.wrappers.phaser_phil.script.phaser_phil import phaser_phil
from ccp4i2.wrappers.phaser_phil.script.phaser_shims import (
    CompositionShim, EnsembleListShim, ObsDataShim, SolutionHook)
from ccp4i2.wrappers.phaser_phil.script import phaser_run


class phaser_mr_auto_phil(phaser_phil):

    TASKNAME = "phaser_mr_auto_phil"
    TASKCOMMAND = None          # Phaser runs in-process
    PHIL_MODE = "MR_AUTO"
    #: Whether the mode searches for the ensembles' copies; a mode that works
    #: on placed solutions (RNP) does not, and need not be asked for copies
    SEARCHES_ENSEMBLES = True
    PHIL_MODE_PATH = "phaser.mode"
    WHATNEXT = ["prosmart_refmac", "modelcraft", "coot_rebuild"]

    #: Beyond what the shims and the mode take out: the driver's own switches
    #: and the "simple" shortcut inputs the ensemble list supersedes.
    PHIL_EXCLUDE_SCOPES = [
        "phaser.sad_mode", "phaser.run_control", "phaser.output_dir",
        # The driver's own switches sit at the top level of its PHIL
        "dry_run", "show_script", "show_defaults", "test_mode", "verbose",
        "phaser.model", "phaser.seq_file", "phaser.mol_weight",
        "phaser.model_rmsd", "phaser.model_identity", "phaser.model_idhi",
        "phaser.model_idlo", "phaser.component_copies", "phaser.search_copies",
        "phaser.chain_type", "phaser.solution",
    ]

    ERROR_CODES = {
        201: {"description": "Expected output file missing"},
        202: {"description": "Phaser reported an error"},
        203: {"description": "Failed to prepare a search model"},
        204: {"description": "Failed to prepare the reflection data"},
        115: {"description": "An ensemble holds more than one model"},
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.xmlroot = etree.Element("PhaserMrResults")
        self.resultObject = None
        #: Callables given the XML root whenever program.xml is written; a
        #: pipeline embeds the live record in its own this way
        self.xml_responders = []
        self._input_space_group = None

    @property
    def _obs_shim(self):
        # Made on first use: the base constructor asks for the shims (to keep
        # their targets out of the parameter tree) before __init__ here runs
        if getattr(self, "_obs_shim_instance", None) is None:
            self._obs_shim_instance = ObsDataShim(self, "F_SIGF", "phaser.hklin", "phaser.labin")
        return self._obs_shim_instance

    @property
    def _model_paths(self):
        if getattr(self, "_model_paths_map", None) is None:
            self._model_paths_map = {}
        return self._model_paths_map

    # -- validation ----------------------------------------------------------
    def validity(self):
        error = super().validity()
        inp = self.container.inputData
        name = f"{self.TASKNAME}.container.inputData"
        ensembles = getattr(inp, "ENSEMBLES", None)
        if (self.SEARCHES_ENSEMBLES and ensembles is not None and len(ensembles) > 0
                and ensembles.copiesToPlace() == 0):
            error.append(
                klass=self.TASKNAME, code=111,
                details=(f"None of the {len(ensembles)} search model(s) is in use and asks "
                         "for copies to be placed, so Phaser has nothing to search for. "
                         "Tick 'use' and set the number of copies on at least one."),
                name=f"{name}.ENSEMBLES", severity=CCP4ErrorHandling.SEVERITY_ERROR)
        for i, ensemble in enumerate(ensembles or []):
            models = [item for item in ensemble.pdbItemList if item.structure.isSet()]
            if len(models) > 1:
                # Phaser checks the models agree in scattering within 20% and
                # refuses otherwise ("Molecular scattering of beta.pdb deviates
                # more than 20% from the mean"); say so before it does
                error.append(
                    klass=self.TASKNAME, code=115,
                    details=(f"Search model {i + 1} holds {len(models)} coordinate files. "
                             "An ensemble is a set of alternative models of the same "
                             "molecule, superposed; Phaser rejects it if their scattering "
                             "differs by more than 20%. Different molecules go in "
                             "separate search models."),
                    name=f"{name}.ENSEMBLES", severity=CCP4ErrorHandling.SEVERITY_WARNING)
            for item in ensemble.pdbItemList:
                if item.identity_to_target.isSet() and item.rms_to_target.isSet():
                    error.append(
                        klass=self.TASKNAME, code=112,
                        details=(f"Search model {i + 1}: give either the sequence identity "
                                 "or the RMS to the target, not both -- Phaser accepts one."),
                        name=f"{name}.ENSEMBLES", severity=CCP4ErrorHandling.SEVERITY_ERROR)
        comp_by = str(inp.COMP_BY)
        if comp_by == "ASU" and not inp.ASUFILE.isSet():
            error.append(klass=self.TASKNAME, code=113,
                         details="Composition by AsuContent file: choose the file.",
                         name=f"{name}.ASUFILE", severity=CCP4ErrorHandling.SEVERITY_ERROR)
        if comp_by == "MW" and not (inp.ASU_PROTEIN_MW.isSet() or inp.ASU_NUCLEICACID_MW.isSet()):
            error.append(klass=self.TASKNAME, code=113,
                         details="Composition by molecular weight: give the protein and/or "
                                 "nucleic-acid weight in the asymmetric unit.",
                         name=f"{name}.ASU_PROTEIN_MW", severity=CCP4ErrorHandling.SEVERITY_ERROR)
        if comp_by == "SEQUENCES" and not any(
                item.seqFile.isSet() for item in inp.SEQUENCES):
            error.append(klass=self.TASKNAME, code=113,
                         details="Composition by sequence files: add at least one.",
                         name=f"{name}.SEQUENCES", severity=CCP4ErrorHandling.SEVERITY_ERROR)
        return error

    # -- shims ----------------------------------------------------------------
    @classmethod
    def phil_shim_targets(cls):
        """The PHIL paths this task's shims write: a pipeline hosting the same
        PHIL leaves them out too, so the two trees have the same shape."""
        return [t for shim in (
            ObsDataShim(None, "F_SIGF", "phaser.hklin", "phaser.labin"),
            EnsembleListShim("ENSEMBLES", "phaser.ensemble", "phaser.search"),
            CompositionShim(),
            FixedPhilShim({"phaser.keywords.general.root": ""}),
        ) for t in shim.phil_targets()]

    def get_shim_definitions(self):
        root = os.path.join(str(self.getWorkDirectory()), "PHASER")
        return [
            self._obs_shim,
            EnsembleListShim("ENSEMBLES", "phaser.ensemble", "phaser.search",
                             path_map=self._model_paths),
            CompositionShim(),
            FixedPhilShim({"phaser.keywords.general.root": root}),
        ]

    # -- the run --------------------------------------------------------------
    def processInputFiles(self):
        error = self._obs_shim.prepare()
        if error.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
            self.appendErrorReport(204, str(error), severity=CCP4ErrorHandling.SEVERITY_ERROR)
            return CPluginScript.FAILED
        try:
            import gemmi
            self._input_space_group = gemmi.read_mtz_file(self._obs_shim.hklin).spacegroup.hm
        except Exception:
            self._input_space_group = None
        # Phaser reads PDB; a model given as mmCIF is converted alongside
        for ensemble in self.container.inputData.ENSEMBLES:
            for item in ensemble.pdbItemList:
                if not item.structure.isSet():
                    continue
                src = str(item.structure.getFullPath())
                if src.lower().endswith((".cif", ".mmcif")):
                    try:
                        import gemmi
                        dst = os.path.join(str(self.getWorkDirectory()),
                                           os.path.splitext(os.path.basename(src))[0] + ".pdb")
                        gemmi.read_structure(src).write_pdb(dst)
                        self._model_paths[src] = dst
                    except Exception as err:
                        self.appendErrorReport(203, f"{src}: {err}", severity=CCP4ErrorHandling.SEVERITY_ERROR)
                        return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self):
        self._phil_path = self.build_working_phil()
        return CPluginScript.SUCCEEDED

    def startProcess(self):
        recorder = phaser_run.PhaserRecorder(self.xmlroot, flush=self.flushXML)
        try:
            self.resultObject = self._run(recorder)
        except phaser_run.PhaserInputError as err:
            # Phaser refused the input; its sentence is the whole story
            self.appendErrorReport(202, str(err), severity=CCP4ErrorHandling.SEVERITY_ERROR)
            self.flushXML(self.xmlroot)
            return CPluginScript.FAILED
        if not self.resultObject.Success():
            self.appendErrorReport(
                202, f"{self.resultObject.ErrorName()}: {self.resultObject.ErrorMessage()}",
                severity=CCP4ErrorHandling.SEVERITY_ERROR)
            self.flushXML(self.xmlroot)
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def _run(self, recorder):
        return phaser_run.run_mode(
            self.get_master_phil(), self._phil_path, self.PHIL_MODE,
            str(self.getWorkDirectory()), recorder, self.makeFileName("LOG"),
            input_hooks=[SolutionHook(self.container)])

    # -- what came out ---------------------------------------------------------
    def processOutputFiles(self):
        result = self.resultObject
        out = self.container.outputData
        work_dir = str(self.getWorkDirectory())
        for i in range(1, len(result.getPdbFiles()) + 1):
            xyz = os.path.join(work_dir, f"PHASER.{i}.pdb")
            hkl = os.path.join(work_dir, f"PHASER.{i}.mtz")
            for path in (xyz, hkl):
                if not os.path.exists(path):
                    self.appendErrorReport(201, path, severity=CCP4ErrorHandling.SEVERITY_ERROR)
                    return CPluginScript.FAILED
            out.XYZOUT.append(out.XYZOUT.makeItem())
            out.XYZOUT[-1].setFullPath(xyz)
            out.XYZOUT[-1].annotation.set(f"Positioned coordinates for solution {i}")
            out.HKLOUT.append(out.HKLOUT.makeItem())
            out.HKLOUT[-1].setFullPath(hkl)
        if len(out.HKLOUT) > 0:
            from ccp4i2.core import CCP4XtalData
            self.splitHkloutList(
                miniMtzsOut=["MAPOUT", "DIFMAPOUT", "PHASEOUT"],
                programColumnNames=["FWT,PHWT", "DELFWT,PHDELWT", "PHIC,FOM"],
                outputBaseName=["MAPOUT", "DIFMAPOUT", "PHASEOUT"],
                outputContentFlags=[1, 1, CCP4XtalData.CPhsDataFile.CONTENT_FLAG_PHIFOM],
                infileList=out.HKLOUT)
            for i in range(len(out.MAPOUT)):
                out.MAPOUT[i].annotation.set(f"Map for solution {i + 1}")
                out.MAPOUT[i].contentFlag.set(1)
                out.MAPOUT[i].subType.set(1)
                out.DIFMAPOUT[i].annotation.set(f"Difference map for solution {i + 1}")
                out.DIFMAPOUT[i].contentFlag.set(1)
                out.DIFMAPOUT[i].subType.set(2)
                out.PHASEOUT[i].annotation.set(f"Calculated phases for solution {i + 1}")
        solutions = result.getDotSol()
        if len(solutions) > 0:
            with open(str(out.SOLOUT.fullPath), "wb") as handle:
                pickle.dump(solutions, handle)
            out.SOLOUT.annotation.set("Solutions from Phaser")
            # Phaser may have solved in another setting of the space group:
            # a pipeline then reindexes the data to match
            solved = str(solutions[0].getSpaceGroupName()).strip()
            given = (self._input_space_group or "").strip()
            reindexed = bool(given) and _same_symbol(solved) != _same_symbol(given)
            out.dataReindexed.set(reindexed)
            if reindexed:
                warnings = self.xmlroot.find("PhaserWarnings")
                if warnings is None:
                    warnings = etree.SubElement(self.xmlroot, "PhaserWarnings")
                etree.SubElement(warnings, "Warning").text = (
                    f"Space group of the best solution ({solved}) differs from the "
                    f"input data ({given})")

        # The record of the run, from the Result object and Phaser's summary
        phaser_run.solutions_xml(result, self.xmlroot)
        blocks = phaser_run.summary_blocks(result.summary())
        summaries = etree.SubElement(self.xmlroot, "Summaries")
        for name, text in blocks:
            node = etree.SubElement(summaries, "Summary")
            node.set("module", name)
            node.text = text
        attempts, unparsed = phaser_run.strategy_attempts(blocks)
        strategy = etree.SubElement(self.xmlroot, "Strategy")
        strategy.set("unparsed", str(unparsed))
        for attempt in attempts:
            node = etree.SubElement(strategy, "Attempt")
            for key, value in attempt.items():
                if value is not None:
                    node.set(key, str(value))
        warnings = self.xmlroot.find("PhaserWarnings")
        if warnings is None:
            warnings = etree.SubElement(self.xmlroot, "PhaserWarnings")
        seen = {w.text for w in warnings}
        for text in result.warnings():
            if text not in seen:
                etree.SubElement(warnings, "Warning").text = text
        self.flushXML(self.xmlroot)
        return CPluginScript.SUCCEEDED

    def flushXML(self, xmlroot):
        """program.xml, written whole then renamed, so a reader never sees half."""
        target = self.makeFileName("PROGRAMXML")
        tmp = target + "_tmp"
        etree.ElementTree(xmlroot).write(tmp, pretty_print=True, xml_declaration=True,
                                         encoding="utf-8")
        os.replace(tmp, target)
        for responder in getattr(self, "xml_responders", ()):
            responder(xmlroot)


def _same_symbol(symbol):
    """Space-group symbols compared without spaces or case."""
    return "".join(str(symbol).split()).upper()
