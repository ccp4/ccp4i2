import os
import pathlib
import shutil

import coot_headless_api

from core.CCP4ModelData import CPdbDataFile
from core.CCP4PluginScript import CPluginScript


class coot_rsr_morph(CPluginScript):
    TASKMODULE = "refinement"
    TASKTITLE = "Real space refinement morphing with Coot API"
    TASKNAME = "coot_rsr_morph"
    TASKVERSION = 202110261437
    WHATNEXT = ["prosmart_refmac"]
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 9999999.9
    MAINTAINER = "stuart.mcnicholas@york.ac.uk"
    RUNEXTERNALPROCESS = False

    def startProcess(self, command=None, handler=None, **kw):
        outFormat = "cif" if self.container.inputData.XYZIN.isMMCIF() else "pdb"
        oldFullPath = pathlib.Path(str(self.container.outputData.XYZOUT.fullPath))
        if outFormat == "cif":
            self.container.outputData.XYZOUT.setFullPath(str(oldFullPath.with_suffix('.cif')))
            self.container.outputData.XYZOUT.contentFlag.set(CPdbDataFile.CONTENT_FLAG_MMCIF)

        xyzin = str(self.container.inputData.XYZIN.fullPath)
        mtzin = str(self.container.inputData.FPHIIN.fullPath)
        xyzout = os.path.normpath(str(self.container.outputData.XYZOUT))
        local_radius = self.container.controlParameters.LOCAL_RADIUS
        gm_alpha = self.container.controlParameters.GM_ALPHA
        blur_b_factor = self.container.controlParameters.BLUR_B_FACTOR

        mc = coot_headless_api.molecules_container_py(True)
        mc.set_make_backups(False)
        mc.set_use_gemmi(True)
        imol = mc.read_pdb(xyzin)
        imap = mc.read_mtz(mtzin, "F", "PHI", "", False, False)
        imap_blurred = mc.sharpen_blur_map(imap, blur_b_factor, False)
        mc.set_imol_refinement_map(imap_blurred)
        mc.generate_self_restraints(imol, local_radius)
        mc.set_refinement_geman_mcclure_alpha(gm_alpha)
        success = mc.refine_residues_using_atom_cid(imol, "//", "ALL", 4000)
        mc.write_coordinates(imol, xyzout)

        shutil.rmtree("coot-backup", ignore_errors=True)

        status = CPluginScript.FAILED
        if success and os.path.exists(str(self.container.outputData.XYZOUT)):
            status = CPluginScript.SUCCEEDED
        self.reportStatus(status)
        return status
