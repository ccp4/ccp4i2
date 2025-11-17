import os
import pathlib

import coot_headless_api
import gemmi

from core.CCP4ModelData import CPdbDataFile
from core.CCP4PluginScript import CPluginScript


class coot_rsr_morph(CPluginScript):
    TASKMODULE = "refinement"
    TASKTITLE = "Real space refinement morphing with Coot API"
    TASKNAME = "coot_rsr_morph"
    TASKVERSION = 202511171635
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

        pdb_path = str(self.container.inputData.XYZIN.fullPath)
        mtz_path = str(self.container.inputData.FPHIIN.fullPath)
        out_path = os.path.normpath(str(self.container.outputData.XYZOUT))
        local_radius = self.container.controlParameters.LOCAL_RADIUS
        gm_alpha = self.container.controlParameters.GM_ALPHA
        blur_b_factor = self.container.controlParameters.BLUR_B_FACTOR

        mc = coot_headless_api.molecules_container_py(True)
        mc.set_use_gemmi(True)
        imol = mc.read_pdb(pdb_path)
        gemmiStructure = gemmi.read_structure(pdb_path)
        imap = mc.read_mtz(mtz_path,"F","PHI","",False,False)
        mc.generate_self_restraints(imol, local_radius)

        mc.set_refinement_geman_mcclure_alpha(gm_alpha)

        imap_blurred = mc.sharpen_blur_map(imap, blur_b_factor, False)
        mc.set_imol_refinement_map(imap_blurred)

        for model in gemmiStructure:
            for chain in model:
                firstResidue, lastResidue = chain[0].seqid.num,chain[-1].seqid.num
                mc.refine_residue_range(imol,chain.name,firstResidue, lastResidue,4000)

        mc.write_coordinates(imol,out_path)

        status = CPluginScript.FAILED
        if os.path.exists(str(self.container.outputData.XYZOUT)):
            status = CPluginScript.SUCCEEDED
        self.reportStatus(status)
        return status
