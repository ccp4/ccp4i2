"""Add to three projects the runs the utility tasks' pages are illustrated
from, each answering the question the task is for:

- ProvideSequence (MDM2): the sequence of the crystallised construct, from
  4hg7.seq.
- ProvideAlignment (Gamma): the gamma demo's alignment, read from a file.
- splitMtz (Gamma): a refinement's full MTZ, as another program would give
  it (refined_gamma.mtz from scenario_import.py), split into CCP4i2's typed
  files in one go: observations, free R set, map coefficients.
- mergeMtz (Gamma): the opposite, for a program outside CCP4i2 that wants
  one file: the observations, the free set the model was refined against,
  and the refinement's two sets of map coefficients, merged.
- pointless_reindexToMatch (BetaBlip): P3(2)21 has two ways to index the
  same data (h,k,l and -h,-k,l). The data are reindexed by -h,-k,l, as a
  second crystal processed independently might have been, and then matched
  back to the ModelCraft model of the first: Pointless must find -h,-k,l.

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_utils.py

Run scenario_refine.py, scenario_mr.py, scenario_mr2.py, scenario_maps.py
and scenario_import.py first. Each task's last job gets an unrun clone.
"""
from pathlib import Path
import xml.etree.ElementTree as ET

from scenario_common import clone_last, i2run, inputs_dir, scratch_home

DEMO = Path(__file__).resolve().parents[3] / "server/ccp4i2/demo_data"


def columns(labels, types, dataset):
    args = []
    for i, (label, kind) in enumerate(zip(labels, types)):
        args += [f"columnList[{i}]/columnLabel={label}",
                 f"columnList[{i}]/columnType={kind}",
                 f"columnList[{i}]/dataset={dataset}"]
    return args


def main():
    # The construct's sequence, as the interface sets it from a sequence
    # file: the task reads only the text (the interface fills it in from a
    # chosen file or model; i2run gives it directly).
    i2run("MDM2", "ProvideSequence",
          "--SEQUENCETEXT", (DEMO / "mdm2" / "4hg7.seq").read_text())

    i2run("Gamma", "ProvideAlignment", "--PASTEORREAD", "ALIGNIN",
          "--ALIGNIN", f"fullPath={DEMO / 'gamma' / 'gamma_self_alignment.fas'}")

    full = inputs_dir() / "refined_gamma.mtz"
    assert full.exists(), "run scenario_import.py first: it makes refined_gamma.mtz"
    i2run("Gamma", "splitMtz", "--HKLIN", f"fullPath={full}",
          "--COLUMNGROUPLIST", "columnGroupType=Obs", "contentFlag=4", "dataset=crystal",
          "selected=True", *columns(["F", "SIGF"], ["F", "Q"], "crystal"),
          "--COLUMNGROUPLIST", "columnGroupType=FreeR", "contentFlag=1", "dataset=HKL_base",
          "selected=True", *columns(["FREER"], ["I"], "HKL_base"),
          "--COLUMNGROUPLIST", "columnGroupType=MapCoeffs", "contentFlag=1",
          "dataset=crystal", "selected=True", *columns(["FWT", "PHWT"], ["F", "P"], "crystal"))

    # Job numbers in the Gamma project: 1 the native data, 2 the job that
    # used the free set the model was refined against, 12 the refinement.
    i2run("Gamma", "mergeMtz",
          "--MINIMTZINLIST", "fileName/fileOut=[1].OBSOUT",
          "--MINIMTZINLIST", "fileName/fileIn=[2].FREERFLAG",
          "--MINIMTZINLIST", "fileName/fileOut=[12].FPHIOUT",
          "--MINIMTZINLIST", "fileName/fileOut=[12].DIFFPHIOUT")

    # Mis-index, then match back.
    data = ["--F_SIGF", "fileOut=import_merged[0].OBSOUT",
            "--FREERFLAG", "fileOut=import_merged[0].FREEOUT"]
    i2run("BetaBlip", "pointless_reindexToMatch", *data,
          "--REFERENCE", "SPECIFY", "--USE_REINDEX", "True",
          "--REINDEX_OPERATOR", "h=-h", "k=-k", "l=l")
    i2run("BetaBlip", "pointless_reindexToMatch",
          "--F_SIGF", "fileOut=pointless_reindexToMatch[-1].F_SIGF_OUT",
          "--FREERFLAG", "fileOut=pointless_reindexToMatch[-1].FREERFLAG_OUT",
          "--REFERENCE", "XYZIN_REF", "--XYZIN_REF", "fileOut=modelcraft[-1].XYZOUT")

    for project, task in (("MDM2", "ProvideSequence"), ("Gamma", "ProvideAlignment"),
                          ("Gamma", "splitMtz"), ("Gamma", "mergeMtz"),
                          ("BetaBlip", "pointless_reindexToMatch")):
        clone_last(project, task)


if __name__ == "__main__":
    main()
