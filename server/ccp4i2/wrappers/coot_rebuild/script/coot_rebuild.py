"""Interactive Coot 0.9 session, database-connected via the cootbridge.

The 0.9 twin of the coot1 wrapper: launch Coot with a Python-2 stub
script and an environment carrying only connection details and the job
identity. The stub loads the shared data layer (api_client) and the
0.9 adapter (coot09_loader) directly by file path - Coot 0.9's
embedded interpreter is Python 2.7 and cannot import the ccp4i2
package - hands them the flat-namespace Coot functions they need, and
the adapter fetches and loads the job's input data.

Retained legacy behaviour: the classic COOT_FILE_DROP/output<N> save
contract for harvest, keybindings, COOTSTATEFILE seeding, and best-
effort inlining of COOTSCRIPTFILE scripts (with a no-op shim for the
retired ccp4i2Interface object those scripts referenced, so they
degrade instead of crashing). The legacy in-Coot GTK2 menu system
(ccp4i2CootInterface) is gone by design; a GTK2 browser over the
shared browse model is a later tier.
"""

import os
from pathlib import Path

from ccp4i2.core import CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript


class coot_rebuild(CPluginScript):
    TASKNAME = "coot_rebuild"
    TASKCOMMAND = "coot"
    ASYNCHRONOUS = True

    ERROR_CODES = {
        200: {"description": "Coot exited with error status"},
        201: {"description": "Failed in harvest operation"},
        202: {"description": "Failed in processOutputFiles"},
    }

    def makeCommandAndScript(self):
        work_dir = Path(self.getWorkDirectory())
        self.dropDir = str(work_dir / "COOT_FILE_DROP")
        try:
            os.makedirs(self.dropDir)
        except OSError:
            if not os.path.isdir(self.dropDir):
                self.dropDir = str(work_dir)

        from ccp4i2.cootbridge import COOT09_STARTUP_STUB
        from ccp4i2.cootbridge.handshake import export_handshake

        export_handshake(self, work_dir, Path(self.dropDir))

        # Per-job extras the static stub reads from the environment, rather
        # than being interpolated into generated source.
        if self.container.inputData.USEKEYBINDINGS.isSet() and \
                self.container.inputData.USEKEYBINDINGS:
            os.environ["CCP4I2_COOT_KEYBINDINGS"] = "1"
        else:
            os.environ.pop("CCP4I2_COOT_KEYBINDINGS", None)
        if self.container.inputData.COOTSCRIPTFILE.isSet():
            os.environ["CCP4I2_COOTSCRIPTFILE"] = \
                self.container.inputData.COOTSCRIPTFILE.fullPath.__str__()
        else:
            os.environ.pop("CCP4I2_COOTSCRIPTFILE", None)

        if self.container.inputData.COOTSTATEFILE.isSet():
            self.copyStateFile()

        cl_args = ["--no-state-script", "--python"]
        if self.container.inputData.DICT.isSet():
            cl_args += ["--dictionary",
                        self.container.inputData.DICT.fullPath.__str__()]
        cl_args += ["--script", str(COOT09_STARTUP_STUB)]
        self.appendCommandLine(cl_args)
        return CPluginScript.SUCCEEDED

    # -- state file seeding (legacy behaviour, bug fixed) -------------------

    def copyStateFile(self):
        """Seed the session from a saved state file, re-pointing its
        molecule load at XYZIN_LIST[0]."""
        text = CCP4Utils.readFile(
            self.container.inputData.COOTSTATEFILE.fullPath.__str__())
        new_text = ""
        for line in text.split("\n"):
            if "handle-read-draw-molecule" in line and \
                    self.container.inputData.XYZIN_LIST.isSet() and \
                    len(self.container.inputData.XYZIN_LIST) > 0:
                new_text += ('(handle-read-draw-molecule "' +
                             self.container.inputData.XYZIN_LIST[0].__str__() +
                             '" 1)\n')
            else:
                new_text += line + "\n"
        # The legacy version computed new_text and then saved the
        # unpatched original - documented as defect #5 in
        # docs/interrupt-and-resume.md. Save the patched text.
        CCP4Utils.saveFile(
            os.path.join(self.dropDir, "0-coot-history.scm"), new_text)

    # -- harvesting ---------------------------------------------------------

    def numberOfOutputFiles(self):
        from ccp4i2.cootbridge import api_client

        outputs = api_client.harvestable_outputs(self.dropDir)
        return outputs[-1][0] if outputs else 0

    def processOutputFiles(self):
        from lxml import etree

        from ccp4i2.cootbridge import api_client
        from ccp4i2.core.CCP4ModelData import CPdbDataFile

        work_dir = Path(self.getWorkDirectory())
        self.xmlroot = etree.Element("coot_rebuild")
        n_models = 0
        n_dicts = 0
        try:
            # Models saved through the drop-dir contract.
            xyzout = self.container.outputData.XYZOUT
            for number, path in api_client.harvestable_outputs(self.dropDir):
                source = Path(path)
                target = work_dir / f"XYZOUT_{n_models}{source.suffix}"
                while target.exists():
                    target = work_dir / \
                        f"XYZOUT_{n_models}_{target.stem}{source.suffix}"
                os.replace(source, target)
                while n_models >= len(xyzout):
                    xyzout.append(xyzout.makeItem())
                xyzout[n_models].setFullPath(str(target))
                xyzout[n_models].annotation.set(
                    f"Coot output file number {number}")
                xyzout[n_models].subType.set(CPdbDataFile.SUBTYPE_MODEL)
                xyzout[n_models].contentFlag.set(
                    CPdbDataFile.CONTENT_FLAG_MMCIF if source.suffix == ".cif"
                    else CPdbDataFile.CONTENT_FLAG_PDB)
                n_models += 1
            # Truncate in place; XYZOUT.set(slice) would deep-copy the
            # items through CDataFile.get()/set() and drop annotation and
            # subType (see coot1.py).
            while len(xyzout) > n_models:
                xyzout.pop()

            # Dictionaries left behind by Coot's ligand tools. Two passes:
            # the classic name patterns (acedrg/pyrogen/prodrg), then a
            # content sniff of any other loose CIF in the work/drop dirs -
            # so a builder dict whose name matches no pattern is still
            # caught. output<N>.cif saved coordinates classify as models
            # and are skipped.
            import glob as _glob

            from ccp4i2.cootbridge.harvest import cif_is_restraint_dictionary

            cif_list = []
            seen = set()

            def _add(candidate):
                real = os.path.normpath(candidate)
                if real not in seen and os.path.isfile(real):
                    seen.add(real)
                    cif_list.append(real)

            for pattern in (
                os.path.join(self.dropDir, "coot-ccp4", "prodrg-*.cif"),
                os.path.join(str(work_dir), "coot-ccp4", "prodrg-*.cif"),
                os.path.join(str(work_dir), "*pyrogen.cif"),
                os.path.join(str(work_dir), "acedrg-*.cif"),
            ):
                for hit in _glob.glob(os.path.normpath(pattern)):
                    _add(hit)
            for extra in (
                _glob.glob(os.path.join(str(work_dir), "*.cif"))
                + _glob.glob(os.path.join(self.dropDir, "*.cif"))
                + _glob.glob(os.path.join(self.dropDir, "coot-ccp4", "*.cif"))
            ):
                if cif_is_restraint_dictionary(extra):
                    _add(extra)

            dictout = self.container.outputData.DICTOUT
            for output_cif in cif_list:
                name = os.path.basename(output_cif)
                target = work_dir / f"DICTOUT_{n_dicts}.cif"
                os.replace(output_cif, target)
                while n_dicts >= len(dictout):
                    dictout.append(dictout.makeItem())
                dictout[n_dicts].setFullPath(str(target))
                if "acedrg" in name:
                    producer = "Acedrg"
                elif "pyrogen" in name:
                    producer = "Pyrogen"
                elif "prodrg" in name:
                    producer = "Prodrg"
                else:
                    producer = "Coot"
                dictout[n_dicts].annotation.set(
                    f"Coot/{producer} created geometry for ligand")
                n_dicts += 1
            while len(dictout) > n_dicts:
                dictout.pop()

            etree.SubElement(self.xmlroot, "number_output_files").text = \
                str(n_models)
            etree.SubElement(self.xmlroot, "number_output_dicts").text = \
                str(n_dicts)

            # Merge new dictionaries into the project library; failures
            # warn rather than fail the job ("no sad face of doom").
            for dict_file in dictout[:n_dicts]:
                try:
                    self.mergeDictToProjectLib(fileName=dict_file.__str__())
                except Exception:
                    self.addReportWarning(
                        "mergeDictToProjectLib raised exception: does not "
                        "compromise output dictionary")
                try:
                    lig_nodes = self.xmlroot.xpath("//LIGANDS")
                    lig_node = (lig_nodes[0] if lig_nodes else
                                etree.SubElement(self.xmlroot, "LIGANDS"))
                    for item in dict_file.fileContent.monomerList:
                        etree.SubElement(lig_node, "ligand").text = \
                            str(item.three_letter_code)
                except Exception:
                    self.addReportWarning(
                        "fileContent.monomerList raised exception: does not "
                        "compromise output dictionary")
        except Exception:
            self.appendErrorReport(202, "Data harvesting failed")

        CCP4Utils.saveEtreeToFile(self.xmlroot,
                                  self.makeFileName("PROGRAMXML"))
        if n_models + n_dicts > 0:
            return CPluginScript.SUCCEEDED
        # Nothing was saved: the job self-deletes rather than litter the
        # project (classic coot_rebuild behaviour).
        return CPluginScript.MARK_TO_DELETE

    def addReportWarning(self, text):
        from lxml import etree

        nodes = self.xmlroot.xpath("//Warnings")
        parent = nodes[0] if nodes else etree.SubElement(self.xmlroot,
                                                         "Warnings")
        etree.SubElement(parent, "Warning").text = text
