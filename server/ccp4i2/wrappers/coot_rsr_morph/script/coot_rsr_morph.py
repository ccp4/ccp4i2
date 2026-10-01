import math
import os
import shutil
import pathlib

from lxml import etree

from ccp4i2.core.CCP4ModelData import CPdbDataFile
from ccp4i2.core.CCP4PluginScript import CPluginScript


def shifts(before_path, after_path, top=5):
    """How far morphing moved the model, as an XML element: the number of
    atoms matched by chain, residue and name, their RMS, mean and largest
    shift, and the residues that moved most."""
    import gemmi
    def atoms(path):
        structure = gemmi.read_structure(str(path))
        return {(c.name, str(r.seqid.num) + r.seqid.icode.strip(), a.name):
                (a.pos, f"{c.name}/{r.name} {r.seqid.num}")
                for c in structure[0] for r in c for a in r}
    before, after = atoms(before_path), atoms(after_path)
    common = [k for k in before if k in after]
    root = etree.Element("Shifts")
    if not common:
        return root
    d = {k: before[k][0].dist(after[k][0]) for k in common}
    values = list(d.values())
    root.set("atoms", str(len(values)))
    root.set("rms", "%.2f" % math.sqrt(sum(v * v for v in values) / len(values)))
    root.set("mean", "%.2f" % (sum(values) / len(values)))
    root.set("max", "%.2f" % max(values))
    by_residue = {}
    for k, v in d.items():
        label = before[k][1]
        by_residue[label] = max(v, by_residue.get(label, 0.0))
    for label, v in sorted(by_residue.items(), key=lambda kv: -kv[1])[:top]:
        etree.SubElement(root, "Residue", name=label, max="%.2f" % v)
    return root


class coot_rsr_morph(CPluginScript):
    TASKNAME = "coot_rsr_morph"
    WHATNEXT = ["prosmart_refmac"]
    ASYNCHRONOUS = True

    def startProcess(self):
        # lazy: the helper imports the external Coot API at execution (worker) only
        from ccp4i2.lib.coot_api import molecules_container
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

        mc = molecules_container(True)
        mc.set_make_backups(False)
        mc.set_use_gemmi(False)
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
            # Say how far the model moved (the report said only "finished",
            # the output was "XYZOUT.pdb"): a few tenths of an Angstrom is
            # the tidy local correction morphing is for.
            moved = shifts(xyzin, str(self.container.outputData.XYZOUT))
            root = etree.Element("coot_rsr_morph")
            root.append(moved)
            with open(self.makeFileName("PROGRAMXML"), "wb") as f:
                f.write(etree.tostring(root, pretty_print=True))
            if moved.get("atoms"):
                self.container.outputData.XYZOUT.annotation = (
                    "RSR morph: %s atoms moved, RMS %s A, largest %s A" % (
                        moved.get("atoms"), moved.get("rms"), moved.get("max")))
        self.reportStatus(status)
        return status
