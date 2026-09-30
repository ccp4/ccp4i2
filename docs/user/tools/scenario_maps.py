"""Add to the Gamma project (scenario_mr2.py) the runs the reflection, phase
and map tools' pages are illustrated from, using the rest of the gamma demo:
a xenon derivative, and the phases supplied with it (initial_phases.mtz;
the demo does not say how they were made).

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_maps.py

Adds: the xenon data; Scaleit (native against derivative); a native
Patterson; HL coefficients to phase and FOM; a short refinement of the
Molrep solution, for model phases; Cphasematch (how good the experimental
phases are against the model); an anomalous difference map (where the xenons are); and a map
calculated from the model. Each documented task's last job gets an unrun
clone.
"""
from scenario_common import clone_last, i2run as _i2run
from scenario_mr2 import DEMO, PROJECT


def i2run(*args):
    _i2run(PROJECT, *args)


def main():
    i2run("import_merged", "--HKLIN", f"fullPath={DEMO / 'merged_intensities_Xe.mtz'}")
    native = "fileOut=[1].OBSOUT"
    xenon = "fileOut=import_merged[-1].OBSOUT"
    phases = f"fullPath={DEMO / 'initial_phases.mtz'}"

    i2run("scaleit", "--MERGEDFILES", f"fullPath={DEMO / 'merged_intensities_native.mtz'}",
          "--MERGEDFILES", f"fullPath={DEMO / 'merged_intensities_Xe.mtz'}")
    i2run("cpatterson", "--F_SIGF", native)
    i2run("chltofom", "--HKLIN", phases)
    i2run("prosmart_refmac", "--F_SIGF", native,
          "--FREERFLAG", f"fullPath={DEMO / 'freeR.mtz'}",
          "--XYZIN", "fileOut=molrep_pipe[-1].XYZOUT", "--NCYCLES", "5",
          "--VALIDATE_MOLPROBITY", "False")
    i2run("cphasematch", "--F_SIGF", native, "--ABCD1", phases,
          "--ABCD2", "fileOut=prosmart_refmac[-1].ABCDOUT")
    i2run("cmapcoeff", "--MAPTYPE", "anom", "--F_SIGF1", xenon, "--ABCD1", phases)
    i2run("density_calculator", "--XYZIN", "fileOut=prosmart_refmac[-1].XYZOUT",
          "--D_MIN", "2.0")

    for task in ("scaleit", "cpatterson", "chltofom", "cphasematch", "cmapcoeff",
                 "density_calculator"):
        clone_last(PROJECT, task)


if __name__ == "__main__":
    main()
