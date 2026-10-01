"""Build the project the Privateer page is illustrated from: the glyco demo
(demo_data/glyco), 4iid, a fungal glycoprotein with many extended N-glycans,
some of whose sugars sit in higher-energy ring conformations. That is what
Privateer is for: it validates each sugar's ring conformation, anomer and
linkage against what the chemistry allows, and draws the glycans.

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_privateer.py

Makes the project Glyco. The job gets an unrun clone, for the input figure.
"""
from pathlib import Path

from scenario_common import clone_last, i2run as _i2run

PROJECT = "Glyco"
DEMO = Path(__file__).resolve().parents[3] / "server/ccp4i2/demo_data/glyco"


def main():
    _i2run(PROJECT, "privateer",
           "--XYZIN", f"fullPath={DEMO / '4iid.pdb'}",
           "--F_SIGF", f"fullPath={DEMO / '4iid.mtz'}",
           "columnLabels=/*/*/[FP,SIGFP]")
    clone_last(PROJECT, "privateer")


if __name__ == "__main__":
    main()
