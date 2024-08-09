from __future__ import print_function

import pathlib
import sys,os
import textwrap
import lxml

from core import CCP4File
from core.CCP4ModelData import CPdbDataFile
from core.CCP4PluginScript import CPluginScript


class coot_rsr_morph(CPluginScript):

    TASKMODULE = "refinement"
    TASKTITLE = "Real space refinement morphing with coot"
    TASKNAME = "coot_rsr_morph"
    TASKCOMMAND = "coot"
    TASKVERSION = 202110261437
    WHATNEXT = ["prosmart_refmac"]
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 9999999.9
    MAINTAINER = "stuart.mcnicholas@york.ac.uk"

    def makeCommandAndScript(self):
        outFormat = "cif" if self.container.inputData.XYZIN.isMMCIF() else "pdb"
        oldFullPath = pathlib.Path(str(self.container.outputData.XYZOUT.fullPath))
        if outFormat == "cif":
            self.container.outputData.XYZOUT.setFullPath(str(oldFullPath.with_suffix('.cif')))
            self.container.outputData.XYZOUT.contentFlag.set(CPdbDataFile.CONTENT_FLAG_MMCIF)

        pdb_path = str(self.container.inputData.XYZIN.fullPath)
        mtz_path = str(self.container.inputData.FPHIIN.fullPath)
        out_path = os.path.normpath(str(self.container.outputData.XYZOUT))
        if sys.platform.startswith("win"):
            out_path = out_path.replace("\\","\\\\")
        local_radius = self.container.controlParameters.LOCAL_RADIUS
        gm_alpha = self.container.controlParameters.GM_ALPHA
        blur_b_factor = self.container.controlParameters.BLUR_B_FACTOR
        script = textwrap.dedent(
            f"""\
            try:
                imol = read_pdb("{pdb_path}")
                imap = make_and_draw_map("{mtz_path}", "F", "PHI", "", 0, 0)
                generate_self_restraints(imol, {local_radius})
                set_show_extra_restraints(imol, 0)
                set_refinement_geman_mcclure_alpha({gm_alpha})
                rmsd = map_sigma_py(imap)
                if rmsd is not None:
                    imap_blurred = sharpen_blur_map(imap, {blur_b_factor})
                    set_imol_refinement_map(imap_blurred)
                    set_matrix(15.0 / rmsd)
                    residues = fit_protein_make_specs(imol, "all-chains")
                    with AutoAccept():
                        refine_residues_py(imol, residues)
                write_{outFormat}_file(imol, "{out_path}")
            except Exception:
                import traceback
                print(traceback.format_exc())
            coot_real_exit(0)
            """
        )
        script_path = os.path.join(self.workDirectory, "script.py")
        with open(script_path, "w") as stream:
            stream.write(script)
        self.appendCommandLine(["--no-graphics", "--script", script_path])
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        status = CPluginScript.FAILED
        if os.path.exists(self.container.outputData.XYZOUT.__str__()):
            status = CPluginScript.SUCCEEDED
        print(
            "coot_rsr_morph.handleFinish",
            self.container.outputData.XYZOUT,
            os.path.exists(self.container.outputData.XYZOUT.__str__()),
        )
        # Create a trivial xml output file
        root = lxml.etree.Element("coot_rsr_morph")
        self.container.outputData.XYZOUT.subType = 1
        xml_file = CCP4File.CXmlDataFile(fullPath=self.makeFileName("PROGRAMXML"))
        xml_file.saveFile(root)
        return status
