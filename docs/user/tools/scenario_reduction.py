"""Build the runs the remaining data-reduction and pipeline pages are
illustrated from.

Thaumatin (a new project), from 300 images (45 degrees, 0.15 degree each) of
the xia2 test sweep th_8_2 (Zenodo record 10271; the i2run test uses 20,
which give 14% of the data):
- xia2_dials: the 300 images, integrated, scaled and merged by xia2 with DIALS;
- xia2_dials twice more, on images 1-150 and 151-300, as two wedges for
- xia2_multiplex, which scales and merges them together;
- AlternativeImportXIA2: the full xia2 run's directory, imported.
- xia2_ssx_reduce and import_serial_pipe: set up, not run (no serial data).

- dr_mr_modelbuild_pipeline: from the full xia2 run's unmerged reflections
  to a built model: scaling, MR with thaumatin (PDB 1rqw), then ModelCraft
  (5 cycles); the AU contents from 1rqw's sequence.

MDM2: findmyseq and arp_warp_classic, set up and not run (neither is
installed here).

Gamma:
- mrparse: search models for the gamma sequence, from CCP4's local PDB
  sequence search (USEAPI False). It downloads the public PDB entries it
  finds; it sends nothing of the project.

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_reduction.py

Needs scenario_refine.py / scenario_mr_steps.py (MDM2) and scenario_import.py
(Gamma) first. Name tasks on the command line to run only those.
"""
import sys
import tarfile
from pathlib import Path

import ccp4i2
from scenario_common import clone_last, fetch, i2run, inputs_dir, last_job_dir

SWEEP = "https://zenodo.org/record/10271/files/th_8_2.tar.bz2"
N_IMAGES = 300
PDBE = "https://www.ebi.ac.uk/pdbe/entry-files/download"
GAMMA_SEQ = Path(ccp4i2.__file__).parent / "demo_data" / "gamma" / "gamma.pir"


FULL_RUN = []  # the full xia2 run's directory, for the pipeline


def one_letter(model_path):
    """The first polymer's sequence, from the model file."""
    import gemmi
    st = gemmi.read_structure(str(model_path))
    st.setup_entities()
    entity = next(e for e in st.entities if e.entity_type == gemmi.EntityType.Polymer)
    return gemmi.one_letter_code(entity.full_sequence)


def images():
    out = inputs_dir() / "th_8_2"
    if len(list(out.glob("*.cbf"))) < N_IMAGES:
        out.mkdir(exist_ok=True)
        with tarfile.open(fetch(SWEEP, "th_8_2.tar.bz2")) as tar:
            members = [m for m in tar.getmembers() if m.isfile()][:N_IMAGES]
            for m in members:
                m.name = Path(m.name).name
                tar.extract(m, out)
    return out


def sweep(image_dir, first, last):
    return ["--IMAGE_FILE", "imageFile/baseName=th_8_2_0001.cbf",
            f"imageFile/relPath={image_dir}", f"imageStart={first}", f"imageEnd={last}"]


def main(only=()):
    def wanted(task):
        return not only or task in only

    failed = []

    def run(project, task, args, clone=True):
        try:
            i2run(project, task, *args)
        except Exception as err:  # keep going: one task's failure is a finding
            failed.append(f"{project} {task}: {err}")
            return None
        job = last_job_dir(project, task)
        print(task, "->", job)
        if clone:
            clone_last(project, task)
        return job

    if wanted("xia2_dials") or wanted("xia2_multiplex") or wanted("AlternativeImportXIA2"):
        image_dir = images()
        full = run("Thaumatin", "xia2_dials", sweep(image_dir, 1, N_IMAGES))
        FULL_RUN[:] = [full]
        if wanted("xia2_multiplex"):
            half = N_IMAGES // 2
            half1 = run("Thaumatin", "xia2_dials", sweep(image_dir, 1, half), clone=False)
            half2 = run("Thaumatin", "xia2_dials", sweep(image_dir, half + 1, N_IMAGES), clone=False)
            if half1 and half2:
                run("Thaumatin", "xia2_multiplex",
                    ["--XIA2_RUN", f"fullPath={half1}", "--XIA2_RUN", f"fullPath={half2}"])
        if full and wanted("AlternativeImportXIA2"):
            run("Thaumatin", "AlternativeImportXIA2", ["--XIA2_DIRECTORY", f"fullPath={full}"])
    if wanted("xia2_xds"):
        # XDS is not installed here: set up on the same sweep, not run.
        run("Thaumatin", "xia2_xds", ["--delay"] + sweep(images(), 1, N_IMAGES), clone=False)
    for task in ("xia2_ssx_reduce", "import_serial_pipe"):
        if wanted(task):
            run("Thaumatin", task, ["--delay"], clone=False)
    if wanted("dr_mr_modelbuild_pipeline"):
        # From images to a built model: the xia2 run's unmerged reflections,
        # scaled by the pipeline, MR with thaumatin 1rqw, then ModelCraft.
        model = fetch(f"{PDBE}/1rqw.cif", "1rqw.cif")
        sequence = inputs_dir() / "thaumatin.seq"
        sequence.write_text(">thaumatin\n" + one_letter(model) + "\n")
        run("Thaumatin", "ProvideAsuContents", [
            "--ASU_CONTENT", f"sequence={one_letter(model)}", "nCopies=1",
            "name=thaumatin", "description=Thaumatin I", "polymerType=PROTEIN",
            "source/baseName=thaumatin.seq", f"source/relPath={inputs_dir()}"], clone=False)
        unmerged = FULL_RUN[0] / "AUTOMATIC_DEFAULT_NATIVE_SWEEP1_INTEGRATE.mtz"
        run("Thaumatin", "dr_mr_modelbuild_pipeline", [
            "--MERGED_OR_UNMERGED", "UNMERGED", "--UNMERGEDFILES", f"file={unmerged}",
            "--ASUIN", "fileOut=ProvideAsuContents[-1].ASUCONTENTFILE",
            "--XYZINORMRBUMP", "XYZINPUT", "--XYZIN", f"fullPath={model}",
            "--SEARCH_PDB", "False", "--SEARCH_AFDB", "False",
            "--MODELCRAFT_NCYC", "5"])
    # Not installed here (findmyseq; ARP/wARP, licensed separately): set up
    # on the MDM2 placement (job 25), not run, for pages from the interface.
    if wanted("findmyseq"):
        run("MDM2", "findmyseq", ["--delay", "--F_SIGF", "fileOut=[1].HKLOUT[0]",
                                  "--XYZIN", "fileOut=[25].XYZOUT[0]",
                                  "--FPHI", "fileOut=[25].MAPOUT_REFMAC"], clone=False)
    if wanted("arp_warp_classic"):
        run("MDM2", "arp_warp_classic", ["--delay", "--AWA_ARP_MODE", "WARPNTRACEMODEL",
                                         "--AWA_FOBS", "fileOut=[1].HKLOUT[0]",
                                         "--AWA_FREE", "fileOut=[1].FREEROUT",
                                         "--AWA_MODELIN", "fileOut=[25].XYZOUT[0]",
                                         "--AWA_SEQIN", "fileOut=[9].ASUCONTENTFILE"], clone=False)
    if wanted("mrparse"):
        run("Gamma", "mrparse", ["--SEQIN", f"fullPath={GAMMA_SEQ}",
                                 "--DATABASE", "PDB", "--USEAPI", "False"])
    if failed:
        raise SystemExit("Failed:\n" + "\n".join(failed))


if __name__ == "__main__":
    main(sys.argv[1:])
