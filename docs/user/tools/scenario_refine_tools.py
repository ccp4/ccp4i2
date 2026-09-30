"""Add to the MDM2 project (scenario_refine.py, scenario_tools.py,
scenario_prep.py) the runs the refinement and coordinate tools' pages are
illustrated from, each asking the question the tool exists to answer:

- Sheetbend: how far can shift-field refinement move the molecular
  replacement solution from MDMX (scenario_prep.py) towards the data?
- Zanuda: is P6(5)22 right, or is the model fitting a pseudosymmetry?
- TLS groups for the refined model.
- AREAIMOL: which protein atoms does Nutlin-3a bury? The refined model
  compared with the same model without its ligand.
- The model against the AU contents; fractional coordinates.
- PAIREF: the data were integrated to 1.25 A and Aimless's automatic cutoff
  stopped at 1.35 A. Do the shells beyond it improve the model? Paired
  refinement needs the uncut data, so Aimless is run again without the
  cutoff, completing the free set it already made rather than drawing a new
  one.

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_refine_tools.py

Run the three scenarios above first. Each task's last job gets an unrun clone.
"""
from scenario_common import clone_last, i2run as _i2run
from scenario_refine import DEMO, PROJECT


def i2run(*args):
    _i2run(PROJECT, *args)


def main():
    data = ["--F_SIGF", "fileOut=aimless_pipe[-1].HKLOUT[0]",
            "--FREERFLAG", "fileOut=aimless_pipe[-1].FREEROUT"]
    refined = "fileOut=prosmart_refmac[-1].XYZOUT"

    i2run("sheetbend", *data,
          "--XYZIN", "fileOut=phaser_simple_phil[-1].XYZOUT[0]")
    i2run("zanuda", *data, "--XYZIN", refined)
    i2run("ProvideTLS", "--XYZIN", refined,
          "--TLSGROUPS", "groupId=1", "chainId=A", "firstRes=18",
          "lastRes=108", "selection=ALL")

    # The ligand's footprint: the same model with and without Nutlin-3a.
    i2run("coordinate_selector", "--XYZIN", refined,
          "selection/text=not (NUT)")
    i2run("areaimol", "--XYZIN", refined,
          "--XYZIN2", "fileOut=coordinate_selector[-1].XYZOUT",
          "--DIFFMODE", "COMPARE")

    i2run("modelASUCheck", "--XYZIN", refined,
          "--ASUIN", "fileOut=ProvideAsuContents[-1].ASUCONTENTFILE")
    i2run("add_fractional_coords", "--XYZIN", refined)

    # Last, because it makes aimless_pipe[-1] the uncut data.
    i2run("aimless_pipe",
          "--UNMERGEDFILES", "crystalName=hg7", "dataset=DS1",
          f"file={DEMO / 'mdm2_unmerged.mtz'}",
          "--XYZIN_REF", f"fullPath={DEMO / '4hg7.pdb'}",
          "--MODE", "MATCH", "--REFERENCE_DATASET", "XYZ",
          "--FREERFLAG", "fileOut=aimless_pipe[-1].FREEROUT",
          "--COMPLETE", "True")
    # No dictionary: the model was refined with the monomer library's NUT,
    # whose atom names it carries. The Acedrg dictionary made from SMILES
    # (scenario_refine.py) names the same atoms differently, and Refmac
    # then finds no restraints for any of them.
    i2run("pairef", *data, "--XYZIN", refined,
          "--INIRES", "1.35", "--NSHELL", "2", "--WSHELL", "0.05")

    for task in ("sheetbend", "zanuda", "ProvideTLS", "areaimol",
                 "modelASUCheck", "add_fractional_coords", "pairef"):
        clone_last(PROJECT, task)


if __name__ == "__main__":
    main()
