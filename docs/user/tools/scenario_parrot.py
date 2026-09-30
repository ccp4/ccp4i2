"""Build the project the parrot help page is illustrated from.

1h1s, phospho-CDK2/cyclin A: two complexes in the asymmetric unit. The
starting phases come from refining a partial model (the two CDK2 chains
only), so parrot has something to recover (the cyclin density) and an NCS
two-fold to find from the model.

Run from server/ with ccp4-python, against a scratch home, never a live one:

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_parrot.py

It makes the project CDK2_CyclinA with four jobs: ASU contents,
refinement of the partial model, parrot, and an unrun clone of parrot.
"""
import os
import subprocess
import sys
import urllib.request
from pathlib import Path

import gemmi

PROJECT = "CDK2_CyclinA"
CODE = "1h1s"
URLS = {
    "model": f"https://www.ebi.ac.uk/pdbe/entry-files/download/pdb{CODE}.ent",
    "mtz": f"https://pdb-redo.eu/db/{CODE}/{CODE}_final.mtz",
}


def scratch_home() -> Path:
    home = os.environ.get("CCP4I2_HOME")
    if not home:
        sys.exit("Set CCP4I2_HOME to a scratch directory: this script "
                 "creates a project and must not touch a live database.")
    for live in (".ccp4i2", ".ccp4i2-django", ".ccp4i2x"):
        if Path(home).resolve() == (Path.home() / live).resolve():
            sys.exit(f"CCP4I2_HOME is the live home {home}; use a scratch one.")
    return Path(home)


def fetch(kind: str, work: Path) -> Path:
    path = work / Path(URLS[kind]).name
    if not path.exists():
        print(f"Fetching {URLS[kind]}")
        urllib.request.urlretrieve(URLS[kind], path)
    return path


def one_letter(entity) -> str:
    """Entity sequence, with modified residues (TPO160) as their parent."""
    letters = []
    for name in entity.full_sequence:
        info = gemmi.find_tabulated_residue(gemmi.Entity.first_mon(name))
        code = info.one_letter_code.upper() if info else "X"
        letters.append(code if code.isalpha() and code != " " else "X")
    return "".join(letters)


def prepare(work: Path):
    """Sequences of the two entities, and the CDK2-only partial model."""
    structure = gemmi.read_structure(str(fetch("model", work)))
    structure.setup_entities()
    polymers = [e for e in structure.entities
                if e.entity_type == gemmi.EntityType.Polymer]
    names = {"A": "CDK2", "B": "CyclinA"}
    sequences = {names[e.name]: one_letter(e) for e in polymers}
    for name, seq in sequences.items():
        (work / f"{name}.seq").write_text(seq + "\n", encoding="ascii")
    structure.remove_ligands_and_waters()
    for model in structure:
        for chain in ("B", "D"):  # the cyclins
            model.remove_chain(chain)
    structure.remove_empty_chains()
    partial = work / "cdk2_only.cif"
    structure.make_mmcif_document().write_file(str(partial))
    return sequences, partial


def i2run(*args: str):
    cmd = [sys.executable, "-m", "ccp4i2.cli.i2run", *args,
           "--project_name", PROJECT]
    print("i2run", args[0])
    subprocess.run(cmd, check=True, env={
        **os.environ, "DJANGO_SETTINGS_MODULE": "ccp4i2.config.settings"})


def main():
    home = scratch_home()
    work = home / "scenario_inputs"
    work.mkdir(parents=True, exist_ok=True)
    sequences, partial = prepare(work)
    mtz = fetch("mtz", work)

    asu = []
    for name, description in (("CDK2", "Cyclin-dependent kinase 2"),
                              ("CyclinA", "Cyclin A2 (171-432)")):
        asu += ["--ASU_CONTENT", f"sequence={sequences[name]}", "nCopies=2",
                f"name={name}", f"description={description}",
                "polymerType=PROTEIN", f"source/baseName={name}.seq",
                f"source/relPath={work}"]
    i2run("ProvideAsuContents", *asu)

    i2run("prosmart_refmac",
          "--XYZIN", f"fullPath={partial}",
          "--F_SIGF", f"fullPath={mtz}", "columnLabels=/*/*/[FP,SIGFP]",
          "--FREERFLAG", f"fullPath={mtz}", "columnLabels=/*/*/[FREE]",
          "--NCYCLES", "5", "--ADD_WATERS", "False",
          "--VALIDATE_MOLPROBITY", "False")

    i2run("parrot",
          "--F_SIGF", "fileUse=prosmart_refmac[-1].F_SIGF",
          "--FREERFLAG", "fileUse=prosmart_refmac[-1].FREERFLAG",
          "--ABCD", "fileUse=prosmart_refmac[-1].ABCDOUT",
          "--ASUIN", "fileUse=ProvideAsuContents[-1].ASUCONTENTFILE",
          "--XYZIN_MODE", "mr",
          "--XYZIN_MR", "fileUse=prosmart_refmac[-1].XYZOUT",
          "--CYCLES", "10")

    # A clone of the parrot job: the input screenshots show a job being set
    # up, not one that has run.
    subprocess.run([sys.executable, "manage.py", "shell", "-c", (
        "from ccp4i2.db.models import Job\n"
        "from ccp4i2.lib.utils.jobs.clone import clone_job\n"
        f"job = Job.objects.filter(project__name='{PROJECT}', "
        "task_name='parrot').order_by('-id').first()\n"
        "clone_job(str(job.uuid))\n")],
        check=True, env={**os.environ,
                         "DJANGO_SETTINGS_MODULE": "ccp4i2.config.settings"})


if __name__ == "__main__":
    main()
