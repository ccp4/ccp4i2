"""Build the project the servalcat_pipe help page is illustrated from.

1h1s, phospho-CDK2/cyclin A with the inhibitor NU6102 (4SP): the deposited
model re-refined against the PDB-REDO reflections, the everyday use of the
task. Validation is left on, so the report shows all it can.

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_servalcat_pipe.py

Makes the project CDK2_NU6102: one refinement, and an unrun clone of it.
"""
from scenario_common import clone_last, fetch, i2run

PROJECT = "CDK2_NU6102"


def main():
    model = fetch("https://www.ebi.ac.uk/pdbe/entry-files/download/pdb1h1s.ent")
    mtz = fetch("https://pdb-redo.eu/db/1h1s/1h1s_final.mtz")
    i2run(PROJECT, "servalcat_pipe",
          "--XYZIN", f"fullPath={model}",
          "--HKLIN", f"fullPath={mtz}", "columnLabels=/*/*/[FP,SIGFP]",
          "--FREERFLAG", f"fullPath={mtz}", "columnLabels=/*/*/[FREE]",
          "--F_SIGF_OR_I_SIGI", "F_SIGF")
    clone_last(PROJECT, "servalcat_pipe")


if __name__ == "__main__":
    main()
