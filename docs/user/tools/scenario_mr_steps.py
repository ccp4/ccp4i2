"""Add to the MDM2 project the runs the Phaser step-by-step page is
illustrated from: molecular replacement taken apart into its steps, each a
job of its own, with the automated search beside them for comparison.

The search model is the Chainsaw model of MDMX (scenario_prep.py), 57%
identical to MDM2, which the automated search places with LLG 104 and TFZ
12.7. The steps: the rotation function (a list of orientations), the
translation function on that list (placed solutions), the packing test
(solutions that do not clash with their symmetry mates), and rigid-body
refinement of what is left.

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_mr_steps.py

Run scenario_refine.py, scenario_tools.py and scenario_prep.py first. Each
task's last job gets an unrun clone.
"""
from scenario_common import clone_last, i2run as _i2run, output_file_id
from scenario_refine import PROJECT


def i2run(*args):
    _i2run(PROJECT, *args)


def common():
    # One copy in the asymmetric unit. Only the rotation function and the
    # automated search use the number; the steps after work on what they
    # are given.
    model = output_file_id(PROJECT, "chainsaw", "XYZOUT")
    return ["--F_SIGF", "fileOut=[1].HKLOUT[0]",
            "--ENSEMBLES", "label=MDMX_chainsaw", "use=True", "number=1",
            "pdbItemList/identity_to_target=0.57", f"pdbItemList/structure/dbFileId={model}",
            "--COMP_BY", "ASU",
            "--ASUFILE", "fileOut=ProvideAsuContents[-1].ASUCONTENTFILE"]


# How far down the rotation function's peaks the search goes. The step's
# default keeps those within 75% of the top (12 orientations here) and the
# right one is not among them: run so, the steps find nothing (best LLG 12,
# in the wrong space group). The automated search keeps those within 60%
# and searches 15% further down ("deep search": 33 orientations), and finds
# the solution from the 23rd. These are its settings.
AS_DEEP_AS_AUTO = ["--phaser__keywords__peaks_rotation__cutoff", "60",
                   "--phaser__keywords__peaks_rotation__down", "15"]


def steps(*rotation_options):
    i2run("phaser_mr_frf_phil", *common(), *rotation_options)
    i2run("phaser_mr_ftf_phil", *common(),
          "--RFILEIN", "fileOut=phaser_mr_frf_phil[-1].RFILEOUT")
    i2run("phaser_mr_pak_phil", *common(),
          "--SOLIN", "fileOut=phaser_mr_ftf_phil[-1].SOLOUT")
    i2run("phaser_mr_rnp_phil", *common(),
          "--SOLIN", "fileOut=phaser_mr_pak_phil[-1].SOLOUT")


def main():
    steps()                    # the defaults: the solution is missed
    i2run("phaser_mr_auto_phil", *common())   # all of it in one job
    steps(*AS_DEEP_AS_AUTO)    # the steps again, as deep as the automated search

    for task in ("phaser_mr_frf_phil", "phaser_mr_ftf_phil", "phaser_mr_pak_phil",
                 "phaser_mr_rnp_phil", "phaser_mr_auto_phil"):
        clone_last(PROJECT, task)


if __name__ == "__main__":
    main()
