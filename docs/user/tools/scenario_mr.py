"""Build the project the molecular-replacement route's pages are illustrated
from: Import merged, Define AU contents, Phaser (basic and expert),
ModelCraft and Validation.

The beta-lactamase / BLIP complex from the demo data that ships with CCP4i2
(demo_data/beta_blip, 3.0 A, P3221): a complex where searching with one
component is not enough. Basic Phaser places beta-lactamase alone; expert
Phaser places both; ModelCraft rebuilds from that solution.

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_mr.py

Makes the project BetaBlip. Each task's last job gets an unrun clone, for
the input figures.
"""
import xml.etree.ElementTree as ET
from pathlib import Path

from scenario_common import clone_last, i2run as _i2run, scratch_home

PROJECT = "BetaBlip"
DEMO = Path(__file__).resolve().parents[3] / "server/ccp4i2/demo_data/beta_blip"


def i2run(*args):
    _i2run(PROJECT, *args)


def sequence(name):
    """The one-letter sequence from a demo .seq (FASTA) file."""
    lines = (DEMO / name).read_text().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


def job_dir(number):
    return scratch_home() / f"projects/betablip/CCP4_JOBS/job_{number}"


def verdict(number):
    xml = ET.parse(job_dir(number) / "program.xml").getroot()
    return xml.findtext(".//Verdict") or ""


def main():
    mtz = DEMO / "beta_blip_P3221.mtz"

    # 1. Import: the merged data; a Free R set is generated.
    i2run("import_merged", "--HKLIN", f"fullPath={mtz}")

    # 2. AU contents, with the data for the Matthews analysis.
    i2run("ProvideAsuContents",
          "--ASU_CONTENT", f"sequence={sequence('beta.seq')}", "nCopies=1",
          "name=BETA", "description=Beta-lactamase TEM-1",
          "polymerType=PROTEIN", "source/baseName=beta.seq",
          f"source/relPath={DEMO}",
          "--ASU_CONTENT", f"sequence={sequence('blip.seq')}", "nCopies=1",
          "name=BLIP", "description=Beta-lactamase inhibitory protein",
          "polymerType=PROTEIN", "source/baseName=blip.seq",
          f"source/relPath={DEMO}",
          "--HKLIN", f"fullPath={mtz}")

    data = ["--F_SIGF", "fileUse=import_merged[-1].OBSOUT",
            "--FREERFLAG", "fileUse=import_merged[-1].FREEOUT",
            "--COMP_BY", "ASU",
            "--ASUFILE", "fileUse=ProvideAsuContents[-1].ASUCONTENTFILE"]

    # 3. Basic Phaser, searching with beta-lactamase alone.
    i2run("phaser_simple_phil", *data,
          "--XYZIN", f"fullPath={DEMO / 'beta.pdb'}",
          "--SEARCHSEQUENCEIDENTITY", "0.9")
    assert "Single" in verdict(3), f"Basic Phaser: {verdict(3)!r}"

    # 4. Expert Phaser, both components, then sheetbend and refmac.
    i2run("phaser_pipeline_phil", *data,
          "--ENSEMBLES", "label=beta", "use=True", "number=1",
          "pdbItemList/identity_to_target=0.9",
          f"pdbItemList/structure={DEMO / 'beta.pdb'}",
          "--ENSEMBLES", "label=blip", "use=True", "number=1",
          "pdbItemList/identity_to_target=0.9",
          f"pdbItemList/structure={DEMO / 'blip.pdb'}")
    assert "Single" in verdict(4), f"Expert Phaser: {verdict(4)!r}"

    # 5. ModelCraft from the expert solution. Fewer cycles than the default
    # (25, stopping when R-free stops improving) to keep the scenario short.
    i2run("modelcraft",
          "--F_SIGF", "fileUse=import_merged[-1].OBSOUT",
          "--FREERFLAG", "fileUse=import_merged[-1].FREEOUT",
          "--ASUIN", "fileUse=ProvideAsuContents[-1].ASUCONTENTFILE",
          "--XYZIN", "fileUse=phaser_pipeline_phil[-1].XYZOUT_REFMAC",
          "--CYCLES", "5")
    assert (job_dir(5) / "XYZOUT.cif").exists(), "ModelCraft built nothing"

    # 6. Validation of the built model, against the data.
    i2run("validate_protein",
          "--XYZIN_1", "fileUse=modelcraft[-1].XYZOUT",
          "--F_SIGF_1", "fileUse=import_merged[-1].OBSOUT")

    for task in ("import_merged", "ProvideAsuContents", "phaser_simple_phil",
                 "phaser_pipeline_phil", "modelcraft", "validate_protein"):
        clone_last(PROJECT, task)


if __name__ == "__main__":
    main()
