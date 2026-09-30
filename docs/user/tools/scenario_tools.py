"""Add to the MDM2 project (scenario_refine.py) the runs the small data and
coordinate tools' pages are illustrated from: AU contents and the Matthews
analysis, a free R set, intensities to amplitudes, a structural superposition
and a coordinate selection.

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_tools.py

Run scenario_refine.py first. Each task's last job gets an unrun clone.
"""
from pathlib import Path

from scenario_common import clone_last, i2run as _i2run
from scenario_refine import DEMO, PROJECT


def i2run(*args):
    _i2run(PROJECT, *args)


def sequence():
    lines = (DEMO / "4hg7.seq").read_text().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


def main():
    data = "fileOut=aimless_pipe[-1].HKLOUT[0]"

    # What is in the asymmetric unit, and how much solvent that implies.
    i2run("ProvideAsuContents",
          "--ASU_CONTENT", f"sequence={sequence()}", "nCopies=1",
          "name=MDM2", "description=MDM2 N-terminal domain",
          "polymerType=PROTEIN", "source/baseName=4hg7.seq",
          f"source/relPath={DEMO}")
    i2run("matthews", "--HKLIN", data, "--MODE", "asu_components",
          "--ASUIN", "fileOut=ProvideAsuContents[-1].ASUCONTENTFILE")

    # A free R set of its own for the reduced data; intensities to amplitudes.
    i2run("freerflag", "--F_SIGF", data, "--FRAC", "0.05")
    i2run("ctruncate", "--OBSIN", data)

    # 4qo4, an independent MDM2 structure, superposed on the refined model;
    # and the ligand alone, taken out of the refined model.
    i2run("gesamt",
          "--XYZIN_QUERY", f"fullPath={DEMO / '4qo4.cif'}",
          "--XYZIN_TARGET", "fileOut=prosmart_refmac[-1].XYZOUT")
    i2run("coordinate_selector",
          "--XYZIN", "fileOut=prosmart_refmac[-1].XYZOUT", "selection/text=(NUT)")

    for task in ("ProvideAsuContents", "matthews", "freerflag", "ctruncate",
                 "gesamt", "coordinate_selector"):
        clone_last(PROJECT, task)


if __name__ == "__main__":
    main()
