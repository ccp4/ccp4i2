"""Interactive Coot 1.x session, database-connected via the cootbridge.

The wrapper's job is deliberately small: it launches Coot with a stub
script and an environment carrying ONLY connection details and the job
identity (the handshake). Everything else - deciding what to load from
the job's input parameters, the in-Coot CCP4i2 menu, the project
browser, saving - happens inside Coot in ccp4i2.cootbridge (see that
package). Harvesting picks up the COOT_FILE_DROP contract the bridge's
save action writes to, plus Coot's own --show-ccp4i2-save-button
directory and loose saves in the work directory.
"""

import os
import sys
from pathlib import Path

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4ModelData import CPdbDataFile
from ccp4i2.cootbridge.harvest import cif_is_restraint_dictionary


class coot1(CPluginScript):
    TASKNAME = "coot1"
    TASKCOMMAND = "coot-1"
    ASYNCHRONOUS = True
    WHATNEXT = ["prosmart_refmac", "coot_rebuild", "coot1", "modelcraft"]

    ERROR_CODES = {}

    def makeCommandAndScript(self):
        work_dir = Path(self.getWorkDirectory())
        drop_dir = work_dir / "COOT_FILE_DROP"
        drop_dir.mkdir(parents=True, exist_ok=True)

        from ccp4i2.cootbridge import COOT1_STARTUP_STUB
        from ccp4i2.cootbridge.handshake import export_handshake

        export_handshake(self, work_dir, drop_dir)

        # A real static stub, not generated source: it reads the
        # environment handshake and hands over to the bridge.
        self.appendCommandLine(
            ["--no-state-script", "--script", str(COOT1_STARTUP_STUB)])
        return CPluginScript.SUCCEEDED

    # -- harvesting ---------------------------------------------------------

    def processOutputFiles(self):
        work_dir = Path(self.getWorkDirectory())
        candidates = []  # (sort_key, path)

        # 1. The bridge's save contract: COOT_FILE_DROP/output<N>.pdb|cif,
        #    in save order.
        from ccp4i2.cootbridge import api_client

        for number, path in api_client.harvestable_outputs(
            str(work_dir / "COOT_FILE_DROP")
        ):
            candidates.append(((0, number), Path(path)))

        # 2. Coot 1's own --show-ccp4i2-save-button writes into
        #    ./coot-ccp4i2/ under the CWD (the work directory).
        button_dir = work_dir / "coot-ccp4i2"
        if button_dir.is_dir():
            for path in sorted(
                list(button_dir.glob("*.pdb")) + list(button_dir.glob("*.cif")),
                key=lambda p: p.stat().st_mtime,
            ):
                candidates.append(((1, path.stat().st_mtime), path))

        # 3. Loose files in the work directory: saved coordinates AND
        #    ligand-builder restraint dictionaries (acedrg/get-monomer
        #    write a *.cif here). They are told apart by content below,
        #    not by extension -- a dictionary CIF is not a model.
        for path in sorted(
            list(work_dir.glob("*.pdb")) + list(work_dir.glob("*.cif"))
        ):
            candidates.append(((2, path.name), path))

        # Split coordinates from restraint dictionaries.
        model_paths = []
        dict_paths = []
        for _key, path in candidates:
            if path.suffix == ".cif" and cif_is_restraint_dictionary(path):
                dict_paths.append(path)
            else:
                model_paths.append(path)

        n_models = self._file_list_into(
            self.container.outputData.XYZOUT, model_paths, work_dir, "XYZOUT",
            self._annotate_model)
        n_dicts = self._file_list_into(
            self.container.outputData.DICTOUT, dict_paths, work_dir, "DICTOUT",
            self._annotate_dict)

        # Merge harvested dictionaries into the project monomer library so
        # downstream tasks see the ligand geometry. Best-effort.
        for dict_file in self.container.outputData.DICTOUT[:n_dicts]:
            try:
                self.mergeDictToProjectLib(fileName=dict_file.__str__())
            except Exception:
                pass
        return CPluginScript.SUCCEEDED

    def _annotate_model(self, item, path):
        item.annotation.set(f"Coot output: {path.name}")
        item.subType.set(CPdbDataFile.SUBTYPE_MODEL)
        item.contentFlag.set(
            CPdbDataFile.CONTENT_FLAG_MMCIF if path.suffix == ".cif"
            else CPdbDataFile.CONTENT_FLAG_PDB)

    def _annotate_dict(self, item, path):
        item.annotation.set(f"Coot ligand dictionary: {path.name}")

    def _file_list_into(self, out_list, paths, work_dir, stem, annotate):
        """File ``paths`` into the ``out_list`` COutputFileList, moving
        files from outside the work dir to canonical names first, and
        setting metadata via ``annotate(item, path)``. Truncates spare
        slots in place with pop() -- NOT out_list.set(slice), which
        deep-copies items through CDataFile.get()/set() and drops the
        annotation/subType just set (the gleaner then falls back to the
        bare param name). Returns the number filed."""
        index = 0
        for path in paths:
            if path.parent != work_dir:
                target = work_dir / f"{stem}_{index}{path.suffix}"
                while target.exists():
                    target = work_dir / \
                        f"{stem}_{index}_{target.stem}{path.suffix}"
                os.replace(path, target)
                path = target
            while index >= len(out_list):
                out_list.append(out_list.makeItem())
            out_list[index].setFullPath(str(path))
            annotate(out_list[index], path)
            index += 1
        while len(out_list) > index:
            out_list.pop()
        return index
