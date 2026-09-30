"""Build the project the refinement route's pages are illustrated from: data
reduction (Aimless), a ligand dictionary (Acedrg), refinement (ProSMART and
Refmac), and the density analysis that follows it (EDSTATS).

MDM2 with Nutlin-3a, from the demo data that ships with CCP4i2
(demo_data/mdm2, the Ligand Tutorial): the unmerged data of 4hg7, reduced
with the indexing matched to the deposited model; that model refined against
them, restrained by ProSMART to the independent MDM2 structure 4qo4, with
waters rebuilt from scratch.

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_refine.py

Makes the project MDM2. Each task's last job gets an unrun clone, for the
input figures.
"""
import xml.etree.ElementTree as ET
from pathlib import Path

import gemmi

from scenario_common import clone_last, i2run as _i2run, scratch_home

PROJECT = "MDM2"
DEMO = Path(__file__).resolve().parents[3] / "server/ccp4i2/demo_data/mdm2"
# Nutlin-3a, with its stereochemistry (demo_data/mdm2/README.html).
NUTLIN_3A = ("COc1ccc(c(OC(C)C)c1)C2=N[C@H]([C@H](N2C(=O)N3CCNC(=O)C3)"
             "c4ccc(Cl)cc4)c5ccc(Cl)cc5")


def i2run(*args):
    _i2run(PROJECT, *args)


def job_dir(number):
    return scratch_home() / f"projects/mdm2/CCP4_JOBS/job_{number}"


def program_xml(number):
    return ET.parse(job_dir(number) / "program.xml").getroot()


def main():
    # 1. Data reduction, indexed to match the model to be refined. The
    #    images were integrated to 1.25 A, well past where the data end
    #    (outer-shell completeness under 5%), so the automatic cutoff is
    #    the point: a first Aimless run estimates the limit, a second
    #    applies it.
    i2run("aimless_pipe",
          "--UNMERGEDFILES", "crystalName=hg7", "dataset=DS1",
          f"file={DEMO / 'mdm2_unmerged.mtz'}",
          "--XYZIN_REF", f"fullPath={DEMO / '4hg7.pdb'}",
          "--MODE", "MATCH", "--REFERENCE_DATASET", "XYZ",
          "--AUTOCUTOFF", "True")
    out = gemmi.read_mtz_file(str(job_dir(1) / "HKLOUT_0-observed_data.mtz"))
    names = [(d.project_name, d.crystal_name, d.dataset_name)
             for d in out.datasets if d.dataset_name != "HKL_base"]
    assert names == [("MDM2", "hg7", "DS1")], f"dataset names: {names}"

    # 2. A dictionary for Nutlin-3a from its SMILES string.
    i2run("LidiaAcedrgNew",
          "--MOLSMILESORSKETCH", "SMILES", "--TLC", "NUT",
          "--SMILESIN", NUTLIN_3A)

    # 3. Refinement: the deposited model without its waters, ProSMART
    #    restraints to 4qo4, waters added back.
    i2run("prosmart_refmac",
          "--F_SIGF", "fileOut=aimless_pipe[-1].HKLOUT[0]",
          "--FREERFLAG", "fileOut=aimless_pipe[-1].FREEROUT",
          "--XYZIN", f"fullPath={DEMO / '4hg7.pdb'}",
          "selection/text=not (HOH)",
          "--prosmartProtein.TOGGLE", "True",
          "--prosmartProtein.REFERENCE_MODELS", f"fullPath={DEMO / '4qo4.cif'}",
          "--ADD_WATERS", "True")

    # 4. How well the refined model fits its density. (Not PISA: qtpisa
    #    opens the QtPISA window on the desktop, so it cannot be scripted.)
    i2run("edstats",
          "--XYZIN", "fileOut=prosmart_refmac[-1].XYZOUT",
          "--FPHIIN1", "fileOut=prosmart_refmac[-1].FPHIOUT",
          "--FPHIIN2", "fileOut=prosmart_refmac[-1].DIFFPHIOUT")

    for task in ("aimless_pipe", "LidiaAcedrgNew", "prosmart_refmac",
                 "edstats"):
        clone_last(PROJECT, task)


if __name__ == "__main__":
    main()
