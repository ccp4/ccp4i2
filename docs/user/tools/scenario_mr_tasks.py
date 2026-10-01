"""Build the runs the remaining molecular-replacement pages are illustrated from.

- comit (MDM2): an omit map for job 25, Phaser's placement of a pruned MDMX
  model, to see the density without the model's bias.
- mrbump_basic (MDM2): MrBUMP's whole pipeline (alignment, model editing,
  Phaser, refinement) on a search model it is given, MDMX chain A from PDB
  3dab, with every online search off: nothing is sent anywhere.
- SIMBAD (Gamma): sequence-free molecular replacement by lattice search:
  CCP4's database of known unit cells against the gamma data's cell.
- molrep_map (CAK, a new project): place the CDK-activating kinase model
  (PDB 7b5o) in its cryo-EM map (EMD-12042), both hands, as the test does.

AMPLE is not run: its helical-ensemble mode fails in CCP4 9 (see test_ample).
MoRDa is not installed here; its page is written from the interface.

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_mr_tasks.py

Run scenario_refine.py, scenario_prep.py and scenario_mr_steps.py (MDM2) and
scenario_import.py (Gamma) first. Each task's last job gets an unrun clone.
"""
import gzip
import shutil

from scenario_common import clone_last, fetch, i2run, last_job_dir

EMDB = "https://ftp.ebi.ac.uk/pub/databases/emdb/structures"


def gunzipped(path):
    out = path.with_suffix("")
    if not out.exists():
        with gzip.open(path, "rb") as fin, open(out, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    return out


def main():
    # 1. Omit map for the MR solution's refinement (job 25).
    i2run("MDM2", "comit",
          "--F_SIGF", "fileOut=[1].HKLOUT[0]",
          "--F_PHI_IN", "fileOut=[25].MAPOUT_REFMAC")
    assert list(last_job_dir("MDM2", "comit").glob("*.mtz")), "comit wrote no map"

    # 2. MrBUMP with one local search model and no online search.
    mdmx = fetch("https://www.ebi.ac.uk/pdbe/entry-files/download/pdb3dab.ent", "3dab.pdb")
    i2run("MDM2", "mrbump_basic",
          "--F_SIGF", "fileOut=[1].HKLOUT[0]", "--FREERFLAG", "fileOut=[1].FREEROUT",
          "--ASUIN", "fileOut=[9].ASUCONTENTFILE",
          "--XYZIN_LIST", f"fullPath={mdmx}", "selection/text=A/",
          "--LOCAL", "True", "--LOCALONLY", "True",
          "--SEARCH_PDB", "False", "--SEARCH_AFDB", "False",
          "--MRMAX", "5", "--PJOBS", "2", "--NCYC", "10")
    assert list(last_job_dir("MDM2", "mrbump_basic").glob("*.pdb")) or \
        list(last_job_dir("MDM2", "mrbump_basic").glob("*.cif")), "MrBUMP produced no model"

    # 3. SIMBAD lattice search on the gamma data.
    i2run("Gamma", "SIMBAD", "--F_SIGF", "fileOut=[1].OBSOUT",
          "--SIMBAD_SEARCH_LEVEL", "Lattice", "--SIMBAD_NPROC", "4")

    # 4. A model into its cryo-EM map.
    mapin = gunzipped(fetch(f"{EMDB}/EMD-12042/map/emd_12042.map.gz"))
    model = fetch("https://files.rcsb.org/download/7b5o.pdb", "7b5o.pdb")
    i2run("CAK", "molrep_map", "--MAPIN", f"fullPath={mapin}", "--XYZIN", f"fullPath={model}",
          "--SEARCH_RESOLUTION", "4.0", "--DOWNSAMPLE", "2")
    job = last_job_dir("CAK", "molrep_map")
    for name in ("ORIGINALMODEL.pdb", "INVERTEDMODEL.pdb"):
        assert (job / name).exists(), f"molrep_map wrote no {name}"

    for project, task in (("MDM2", "comit"), ("MDM2", "mrbump_basic"),
                          ("Gamma", "SIMBAD"), ("CAK", "molrep_map")):
        clone_last(project, task)


if __name__ == "__main__":
    main()
