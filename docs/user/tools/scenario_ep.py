"""Build the project the experimental-phasing route's pages are illustrated
from: SAD phasing with Phaser, alone and as a pipeline.

Gamma-adaptin ear domain (PDB 1gyu) soaked in xenon, from the demo data
that ships with CCP4i2 (demo_data/gamma: home-source Cu K-alpha data,
1.54 A, to 1.8 A). The four xenon sites are given, as a substructure search
would have found them (SHELX, licensed separately, is not needed).

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_ep.py

Makes the project GammaXe: import, Phaser EP (a single run), the Phaser EP
pipeline (both hands, parrot, ModelCraft), and an unrun clone of each
Phaser task for the input figures.
"""
import xml.etree.ElementTree as ET
from pathlib import Path

from scenario_common import clone_last, i2run as _i2run, scratch_home

PROJECT = "GammaXe"
DEMO = Path(__file__).resolve().parents[3] / "server/ccp4i2/demo_data/gamma"


def i2run(*args):
    _i2run(PROJECT, *args)


def program_xml(number):
    path = scratch_home() / f"projects/gammaxe/CCP4_JOBS/job_{number}/program.xml"
    return ET.parse(path).getroot()


def main():
    # 1. Import the anomalous intensities and the crystal's Free R set.
    i2run("import_merged",
          "--HKLIN", f"fullPath={DEMO / 'merged_intensities_Xe.mtz'}",
          "--FREERFLAG", f"fullPath={DEMO / 'freeR.mtz'}")

    common = ["--F_SIGF", "fileOut=import_merged[-1].OBSOUT",
              "--XYZIN_HA", f"fullPath={DEMO / 'heavy_atoms.pdb'}",
              "--ELEMENTS", "Xe",
              "--WAVELENGTH", "1.542",
              "--LLGC_CYCLES", "20",
              "--COMP_BY", "ASU",
              "--ASUFILE", f"fullPath={DEMO / 'gamma.asu.xml'}"]

    # 2. One Phaser run: phase from the sites and complete the substructure.
    i2run("phaser_ep_auto_phil", *common)
    fom = float(program_xml(2).findtext(".//Overall/fom") or 0)
    assert fom > 0.3, f"Phaser EP: overall FOM {fom}"

    # 3. The pipeline: both hands, parrot on each, then a short ModelCraft.
    i2run("phaser_ep_phil", *common,
          "--FREERFLAG", "fileOut=import_merged[-1].FREEOUT",
          "--RUNPARROT", "True",
          "--RUNMODELCRAFT", "True",
          "--MODELCRAFT_ITERATIONS", "3")

    for task in ("phaser_ep_auto_phil", "phaser_ep_phil"):
        clone_last(PROJECT, task)


if __name__ == "__main__":
    main()
