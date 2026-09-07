"""Shims from CCP4i2's typed inputs to Phaser's PHIL.

Phaser's driver (phaser.phenix_interface.driver) builds the phaser.Input*
object from the working phil, so what these write is the whole of the
translation: a reflection file and its labels, the ensembles and what to
search for, and the composition. Each shim owns the paths it writes, which
PhilPluginScript then keeps out of the generic parameter tree.
"""
import logging
import os

from ccp4i2.utils.phil_shims import PhilShim

logger = logging.getLogger(f"ccp4i2:{__name__}")


class ObsDataShim(PhilShim):
    """F_SIGF -> phaser.hklin and phaser.labin.

    The mini-MTZ is rewritten as mean intensities when it holds intensities
    and mean amplitudes otherwise -- Phaser's own preference, and the choice
    the driver then makes for itself from the file (setREFL_I_SIGI or
    setREFL_F_SIGF). The rewrite is the plugin's makeHklin(), run from
    processInputFiles() through prepare(); convert() only names the result.
    """

    def __init__(self, plugin, input_name, hklin_path, labin_path):
        self.plugin = plugin
        self.input_name = input_name
        self.phil_hklin_path = hklin_path
        self.phil_labin_path = labin_path
        self.hklin = None
        self.labin = None

    def prepare(self):
        """Write the reflection file Phaser will read; returns a CErrorReport."""
        from ccp4i2.core.CCP4XtalData import CObsDataFile
        obs = getattr(self.plugin.container.inputData, self.input_name)
        flag = int(obs.contentFlag) if obs.contentFlag.isSet() else 0
        intensities = flag in (CObsDataFile.CONTENT_FLAG_IPAIR,
                               CObsDataFile.CONTENT_FLAG_IMEAN)
        target = (CObsDataFile.CONTENT_FLAG_IMEAN if intensities
                  else CObsDataFile.CONTENT_FLAG_FMEAN)
        hklin, error = self.plugin.makeHklin([[self.input_name, target]], hklin="hklin")
        self.hklin = hklin
        self.labin = "I,SIGI" if intensities else "F,SIGF"
        return error

    def convert(self, container, work_directory):
        if not self.hklin:
            return []
        return [(self.phil_hklin_path, self.hklin), (self.phil_labin_path, self.labin)]


class EnsembleListShim(PhilShim):
    """ENSEMBLES (CEnsembleList) -> one phaser.ensemble block per ensemble,
    and one phaser.search block per ensemble asking for copies.

    A coordinate set carries identity OR rmsd, as Phaser requires; one with
    neither reads its variance from the PDB remarks, which is what the
    classic wrapper's "CARD ON" meant. `path_map` rewrites a model path to
    a prepared copy (a CIF converted to PDB, say).
    """

    def __init__(self, input_name, ensemble_path, search_path, path_map=None,
                 fixed_name="FIXENSEMBLES"):
        self.input_name = input_name
        self.phil_ensemble_path = ensemble_path
        self.phil_search_path = search_path
        self.path_map = path_map if path_map is not None else {}
        #: An optional list of ensemble labels already placed, at the origin
        #: of their own coordinates: Phaser's solution_at_origin
        self.fixed_name = fixed_name

    def convert(self, container, work_directory):
        entries = []
        ensembles = getattr(container.inputData, self.input_name, None)
        if ensembles is None:
            return entries
        fixed = getattr(container.inputData, self.fixed_name, None)
        fixed = {str(x) for x in fixed} if fixed is not None else set()
        for i, ensemble in enumerate(ensembles):
            label = str(ensemble.label) if ensemble.label.isSet() else f"ensemble_{i + 1}"
            fields = [("model_id", label)]
            if label in fixed:
                fields.append(("solution_at_origin", True))
            for item in ensemble.pdbItemList:
                if not item.structure.isSet():
                    continue
                path = str(item.structure.getFullPath())
                path = self.path_map.get(path, path)
                coords = [("pdb", path)]
                if item.identity_to_target.isSet():
                    coords.append(("identity", float(item.identity_to_target)))
                elif item.rms_to_target.isSet():
                    coords.append(("rmsd", float(item.rms_to_target)))
                else:
                    coords.append(("read_variance_from_pdb_remarks", True))
                fields.append(("coordinates", coords))
            entries.append((self.phil_ensemble_path, fields))
            copies = int(ensemble.copiesToPlace())
            if copies > 0:
                entries.append((self.phil_search_path,
                                [("ensembles", label), ("copies", copies)]))
        return entries


class CompositionShim(PhilShim):
    """The composition, from whichever source COMP_BY names.

    DEFAULT   average solvent content, 50% -- Phaser's COMP AVERAGE, and it
              raises the composition itself if the search asks for more
    ASU       the ASU file: one composition.chain per selected sequence,
              with its copies (the file's per-sequence tick boxes are honoured)
    SEQUENCES sequence files with copies, one composition.chain each
    MW        protein and/or nucleic-acid molecular weight in the ASU
    SOLVENT   a solvent fraction the user gives

    Phaser reads a sequence from a file, so ASU sequences are written out
    one FASTA each into the job directory.
    """

    def __init__(self, comp_by="COMP_BY", asu="ASUFILE", sequences="SEQUENCES",
                 solvent="SOLVENT_FRACTION", protein_mw="ASU_PROTEIN_MW",
                 nucleic_mw="ASU_NUCLEICACID_MW",
                 chain_path="phaser.composition.chain",
                 solvent_path="phaser.composition.solvent"):
        self.comp_by = comp_by
        self.asu = asu
        self.sequences = sequences
        self.solvent = solvent
        self.protein_mw = protein_mw
        self.nucleic_mw = nucleic_mw
        self.phil_chain_path = chain_path
        self.phil_solvent_path = solvent_path

    def convert(self, container, work_directory):
        inp = container.inputData
        mode = str(getattr(inp, self.comp_by)) if hasattr(inp, self.comp_by) else "DEFAULT"
        if mode == "SOLVENT":
            return [(self.phil_solvent_path, float(getattr(inp, self.solvent)))]
        if mode == "MW":
            entries = []
            for name, chain_type in ((self.protein_mw, None), (self.nucleic_mw, "na")):
                mw = getattr(inp, name, None)
                if mw is not None and mw.isSet() and float(mw) > 0:
                    fields = [("mw", float(mw)), ("num", 1)]
                    if chain_type:
                        fields.append(("chain_type", chain_type))
                    entries.append((self.phil_chain_path, fields))
            return entries
        if mode == "ASU":
            return self._chains_from_asu(getattr(inp, self.asu, None), work_directory)
        if mode == "SEQUENCES":
            return self._chains_from_sequences(getattr(inp, self.sequences, None))
        return [(self.phil_solvent_path, 0.5)]

    def _chains_from_sequences(self, seq_list):
        entries = []
        for item in seq_list or []:
            if not item.seqFile.isSet():
                continue
            copies = int(item.nCopies) if item.nCopies.isSet() else 1
            entries.append((self.phil_chain_path, [
                ("sequence_file", str(item.seqFile.getFullPath())),
                ("num", max(1, copies))]))
        return entries

    def _chains_from_asu(self, asu, work_directory):
        entries = []
        if asu is None or not asu.isSet():
            return entries
        try:
            asu.loadFile()
        except Exception as err:
            logger.warning("Could not load ASU file %s: %s", self.asu, err)
            return entries
        content = getattr(asu, "fileContent", None)
        seq_list = getattr(content, "seqList", None) if content is not None else None
        for i, seq in enumerate(seq_list or []):
            if hasattr(asu, "isSelected") and not asu.isSelected(seq):
                continue
            sequence = _text(seq.sequence).replace(" ", "").replace("\n", "")
            if not sequence:
                continue
            name = _text(getattr(seq, "name", None)) or f"sequence_{i + 1}"
            path = os.path.join(work_directory, f"composition_{i + 1}.fasta")
            with open(path, "w", encoding="utf-8") as fasta:
                fasta.write(f">{name}\n")
                for j in range(0, len(sequence), 60):
                    fasta.write(sequence[j:j + 60] + "\n")
            try:
                copies = max(1, int(seq.nCopies))
            except (TypeError, ValueError):
                copies = 1
            polymer = _text(getattr(seq, "polymerType", None)).upper()
            fields = [("sequence_file", path), ("num", copies)]
            if polymer in ("DNA", "RNA"):
                fields.append(("chain_type", "na"))
            entries.append((self.phil_chain_path, fields))
        return entries


def _text(obj):
    if obj is None:
        return ""
    value = getattr(obj, "value", obj)
    return "" if value is None else str(value)


class EpCrystalShim(PhilShim):
    """The one crystal EP_AUTO takes: anomalous data, wavelength, and the
    substructure -> one phaser.crystal block.

    The mini-MTZ is rewritten as anomalous amplitude pairs (FPAIR) by the
    plugin's makeHklin() from prepare(), as the classic wrapper did, and the
    labels are read back from that file with iotbx: the driver matches its
    own label_string() exactly ("Fplus,SIGFplus,Fminus,SIGFminus,merged"),
    which is not something to compose by hand. The substructure goes to
    crystal.pdb_file, which the driver passes to setATOM_PDB -- ha_file is
    Phaser's own heavy-atom format, not a PDB.
    """

    def __init__(self, plugin, obs_name="F_SIGF", sites_name="XYZIN_HA",
                 wavelength_name="WAVELENGTH", crystal_path="phaser.crystal"):
        self.plugin = plugin
        self.obs_name = obs_name
        self.sites_name = sites_name
        self.wavelength_name = wavelength_name
        self.phil_crystal_path = crystal_path
        self.hklin = None
        self.labin = None

    def prepare(self):
        from ccp4i2.core.CCP4XtalData import CObsDataFile
        hklin, error = self.plugin.makeHklin(
            [[self.obs_name, CObsDataFile.CONTENT_FLAG_FPAIR]], hklin="hklin")
        if hklin:
            from iotbx import file_reader
            arrays = file_reader.any_file(hklin, force_type="hkl").file_server.miller_arrays
            anomalous = [a.info().label_string() for a in arrays if a.anomalous_flag()]
            self.hklin = hklin
            self.labin = anomalous[0] if anomalous else None
        return error

    def convert(self, container, work_directory):
        if not self.hklin or not self.labin:
            return []
        inp = container.inputData
        fields = [("xtal_id", "crystal1")]
        sites = getattr(inp, self.sites_name, None)
        if sites is not None and sites.isSet():
            fields.append(("pdb_file", str(sites.getFullPath())))
        dataset = [("wave_id", "wave1"), ("hklin", self.hklin), ("labin", self.labin)]
        wavelength = getattr(inp, self.wavelength_name, None)
        if wavelength is not None and wavelength.isSet():
            dataset.append(("wavelength", float(wavelength)))
        fields.append(("dataset", dataset))
        return [(self.phil_crystal_path, fields)]


class PartialModelShim(PhilShim):
    """A partial model to phase from, when there is one."""

    def __init__(self, partial_by="PARTIAL_BY", model_name="XYZIN_PARTIAL",
                 path="phaser.keywords.partial"):
        self.partial_by = partial_by
        self.model_name = model_name
        self.phil_partial_path = path

    def convert(self, container, work_directory):
        inp = container.inputData
        if str(getattr(inp, self.partial_by, "NONE")) != "MODEL":
            return []
        model = getattr(inp, self.model_name, None)
        if model is None or not model.isSet():
            return []
        return [(self.phil_partial_path, [("mode", "model"), ("pdb", str(model.getFullPath()))])]


class LlgCompletionShim(PhilShim):
    """Log-likelihood-gradient completion of the substructure: which
    elements to look for, how many cycles, and whether they are pure
    anomalous scatterers (Phaser's IMAGINARY method)."""

    def __init__(self, elements_name="ELEMENTS", cycles_name="LLGC_CYCLES",
                 pure_name="PURE_ANOMALOUS", path="phaser.keywords.llgcompletion"):
        self.elements_name = elements_name
        self.cycles_name = cycles_name
        self.pure_name = pure_name
        self.phil_llgc_path = path

    def convert(self, container, work_directory):
        inp = container.inputData
        elements = [str(e).strip() for e in getattr(inp, self.elements_name, []) or []]
        elements = [e for e in elements if e]
        if not elements:
            return []
        fields = [("complete", True), ("scatterer", " ".join(elements))]
        cycles = getattr(inp, self.cycles_name, None)
        if cycles is not None and cycles.isSet():
            fields.append(("ncycle", int(cycles)))
        pure = getattr(inp, self.pure_name, None)
        if pure is not None and pure.isSet() and bool(pure):
            fields.append(("method", "imaginary"))
        return [(self.phil_llgc_path, fields)]


def read_solutions(path):
    """The mr_solution object pickled by an earlier job: a set of solutions,
    or a rotation list (is_rlist()) from a rotation-function run."""
    import pickle
    with open(path, "rb") as handle:
        return pickle.load(handle)


def is_rotation_list(solutions):
    """Whether these are a rotation list rather than solutions with
    translations. A plain sequence of sets (a test double, an old pickle)
    is taken to be the latter."""
    check = getattr(solutions, "is_rlist", None)
    return bool(check()) if callable(check) else False


def solution_model_ids(solutions):
    """The ensemble labels the solutions name, of either kind."""
    if is_rotation_list(solutions):
        return {str(r.MODLID) for s in solutions for r in s.RLIST}
    return {str(k.MODLID) for s in solutions for k in s.KNOWN}


class SolutionHook:
    """SOLIN, a pickled mr_solution from an earlier job, handed to Phaser as
    the solutions to start from -- or, for the translation function, the
    rotation list to try. Not a shim: the PHIL's own `solution` keyword
    wants a pickled Result, so this goes in through setSOLU on the built
    phaser.Input (run_mode's input_hooks)."""

    def __init__(self, container, input_name="SOLIN"):
        self.container = container
        self.input_name = input_name

    def __call__(self, phaser_input):
        solin = getattr(self.container.inputData, self.input_name, None)
        if solin is None or not solin.isSet():
            return
        phaser_input.setSOLU(read_solutions(str(solin.getFullPath())))
