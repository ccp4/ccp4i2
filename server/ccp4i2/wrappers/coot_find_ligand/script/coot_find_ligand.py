import os
import pathlib
import shutil
import xml.etree.ElementTree as ET

import coot_headless_api

from ccp4i2.core.CCP4ModelData import CPdbDataFile
from ccp4i2.core.CCP4PluginScript import CPluginScript


class coot_find_ligand(CPluginScript):
    TASKMODULE = "model_building"
    TASKTITLE = "Find ligand with Coot API"
    TASKNAME = "coot_find_ligand"
    TASKVERSION = 0.0
    WHATNEXT = ["prosmart_refmac"]
    TIMEOUT_PERIOD = 9999999.9
    RUNEXTERNALPROCESS = False

    def startProcess(self):
        inputData = self.container.inputData
        outputData = self.container.outputData
        params = self.container.controlParameters

        if inputData.XYZIN.isMMCIF():
            path = pathlib.Path(str(outputData.XYZOUT.fullPath))
            outputData.XYZOUT.setFullPath(str(path.with_suffix(".cif")))
            outputData.XYZOUT.contentFlag.set(CPdbDataFile.CONTENT_FLAG_MMCIF)

        xyzin = str(inputData.XYZIN.fullPath)
        mtzin = str(inputData.FPHI.fullPath)
        cifin = str(inputData.DICT.fullPath)
        xyzout = str(outputData.XYZOUT.fullPath)

        mc = coot_headless_api.molecules_container_py(True)
        mc.set_make_backups(False)
        mc.set_use_gemmi(False)
        imol = mc.read_pdb(xyzin)
        imap = mc.read_mtz(mtzin, "F", "PHI", "", False, False)
        mc.import_cif_dictionary(cifin, mc.get_imol_enc_any())
        ilig = mc.get_monomer(str(inputData.COMPID))

        result = mc.fit_ligand(
            imol_protein=imol,
            imol_map=imap,
            imol_ligand=ilig,
            n_rmsd=float(params.THRESHOLD),
            use_conformers=bool(params.FLEXIBLE),
            n_conformers=int(params.CONFORMERS),
        )

        xmlRoot = ET.Element("coot_find_ligand")
        clusterMols = {}
        for info in result:
            if info.cluster_idx not in clusterMols:
                clusterMols[info.cluster_idx] = str(info.imol)
                element = ET.SubElement(xmlRoot, "AddedLigand")
                element.set("clusterVolume", str(info.get_cluster_volume()))
                element.set("fittingScore", str(info.get_fitting_score()))
        mc.merge_molecules(imol, ":".join(clusterMols.values()))
        mc.write_coordinates(imol, xyzout)

        ET.ElementTree(xmlRoot).write(self.makeFileName("PROGRAMXML"))

        shutil.rmtree("coot-backup", ignore_errors=True)

        status = CPluginScript.FAILED
        if os.path.exists(xyzout):
            status = CPluginScript.SUCCEEDED
        self.reportStatus(status)
        return status
