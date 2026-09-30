"""Build the project the Molrep, DIMPLE and Csymmatch pages are illustrated
from: the native data of the gamma demo (P212121, 1.8 A), solved by
molecular replacement with the model supplied with it.

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_mr2.py

Makes the project Gamma: import, Molrep, DIMPLE, then Csymmatch moving the
Molrep solution onto the supplied model's origin. Each task's last job gets
an unrun clone.
"""
from pathlib import Path

from scenario_common import clone_last, i2run as _i2run

PROJECT = "Gamma"
DEMO = Path(__file__).resolve().parents[3] / "server/ccp4i2/demo_data/gamma"


def i2run(*args):
    _i2run(PROJECT, *args)


def main():
    i2run("import_merged", "--HKLIN",
          f"fullPath={DEMO / 'merged_intensities_native.mtz'}")
    # The free set the model was refined against, not a new one: reflections
    # the model has already been fitted to are not free, and with a new set
    # DIMPLE's R-free rose as its R fell.
    data = ["--F_SIGF", "fileOut=import_merged[-1].OBSOUT",
            "--FREERFLAG", f"fullPath={DEMO / 'freeR.mtz'}"]
    model = f"fullPath={DEMO / 'gamma_model.pdb'}"

    i2run("molrep_pipe", *data, "--XYZIN", model)
    i2run("i2Dimple", *data, "--XYZIN", model)
    i2run("csymmatch", "--XYZIN_QUERY", "fileOut=molrep_pipe[-1].XYZOUT",
          "--XYZIN_TARGET", model)

    for task in ("molrep_pipe", "i2Dimple", "csymmatch"):
        clone_last(PROJECT, task)


if __name__ == "__main__":
    main()
