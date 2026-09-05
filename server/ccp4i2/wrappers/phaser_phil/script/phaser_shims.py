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

    def __init__(self, input_name, ensemble_path, search_path, path_map=None):
        self.input_name = input_name
        self.phil_ensemble_path = ensemble_path
        self.phil_search_path = search_path
        self.path_map = path_map if path_map is not None else {}

    def convert(self, container, work_directory):
        entries = []
        ensembles = getattr(container.inputData, self.input_name, None)
        if ensembles is None:
            return entries
        for i, ensemble in enumerate(ensembles):
            label = str(ensemble.label) if ensemble.label.isSet() else f"ensemble_{i + 1}"
            fields = [("model_id", label)]
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
    ASU       the ASU file: one composition.chain per sequence, with its copies
    SEQUENCES sequence files with copies, one composition.chain each
    SOLVENT   a solvent fraction the user gives

    Phaser reads a sequence from a file, so ASU sequences are written out
    one FASTA each into the job directory.
    """

    def __init__(self, comp_by="COMP_BY", asu="ASUFILE", sequences="SEQUENCES",
                 solvent="SOLVENT_FRACTION",
                 chain_path="phaser.composition.chain",
                 solvent_path="phaser.composition.solvent"):
        self.comp_by = comp_by
        self.asu = asu
        self.sequences = sequences
        self.solvent = solvent
        self.phil_chain_path = chain_path
        self.phil_solvent_path = solvent_path

    def convert(self, container, work_directory):
        inp = container.inputData
        mode = str(getattr(inp, self.comp_by)) if hasattr(inp, self.comp_by) else "DEFAULT"
        if mode == "SOLVENT":
            return [(self.phil_solvent_path, float(getattr(inp, self.solvent)))]
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
