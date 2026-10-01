"""Build the runs the analysis, density-modification and rigid-body pages are
illustrated from, on projects the earlier scenarios made.

- AUSPEX (Gamma): the native data checked for ice rings and other outliers.
- phaser_EP_LLG (GammaXe): an anomalous LLG map phased by crank2's model,
  showing the anomalous scatterers the model does not contain.
- molrep_selfrot (CDK2_CyclinA): the self-rotation function of 1h1s, two
  CDK2/cyclin A complexes related by a non-crystallographic two-fold.
- dm_multidomain (CDK2_CyclinA): density modification from the CDK2-only
  partial model, averaging the two CDK2 copies with the N- and C-lobes as
  separate bodies (kinase lobes move relative to each other).
- phaser_rnp_pipeline_phil (LckKinase): SliceNDice's placement of an Lck
  kinase domain split into its lobes, refined as two rigid bodies.
- shelxeMR (MDM2): SHELXE density modification and tracing from the Phaser
  placement of the pruned MDMX model (job 25).

    env CCP4I2_HOME=/tmp/docs-home ccp4-python ../docs/user/tools/scenario_analysis.py

Needs scenario_import.py (Gamma), scenario_ep.py (GammaXe), scenario_parrot.py
(CDK2_CyclinA), scenario_slicendice.py (LckKinase) and scenario_mr_steps.py
(MDM2) first; shelxeMR needs SHELXDIR. Each task's job gets an unrun clone.
Name tasks on the command line to run only those.
"""
import sys

from scenario_common import clone_last, i2run, last_job_dir

RUNS = [
    ("Gamma", "AUSPEX", ["--F_SIGF", "fileOut=[1].OBSOUT"]),
    ("GammaXe", "phaser_EP_LLG", [
        "--F_SIGF", "fileOut=[1].OBSOUT",
        "--PARTIALMODELORMAP", "MODEL",
        "--XYZIN_PARTIAL", "fileOut=crank2[-1].XYZOUT"]),
    ("CDK2_CyclinA", "molrep_selfrot", ["--F_SIGF", "fileIn=[2].F_SIGF"]),
    ("CDK2_CyclinA", "dm_multidomain", [
        "--F_SIGF", "fileIn=[2].F_SIGF", "--FREERFLAG", "fileIn=[2].FREERFLAG",
        "--XYZIN", "fileOut=[2].XYZOUT", "--PHASE_SOURCE", "model",
        "--ASUIN", "fileOut=[1].ASUCONTENTFILE",
        "--DOMAINS", "segments=0-85", "mode=average",
        "--DOMAINS", "segments=86-296", "mode=average"]),
    ("LckKinase", "phaser_rnp_pipeline_phil", [
        "--F_SIGF", "fileOut=[1].OBSOUT", "--FREERFLAG", "fileOut=[1].FREEOUT",
        "--XYZIN_PARENT", "fileOut=slicendice[-1].XYZOUT",
        "--SELECTIONS", "text=A/", "--SELECTIONS", "text=B/",
        "--COMP_BY", "ASU", "--ASUFILE", "fileOut=[2].ASUCONTENTFILE"]),
    ("MDM2", "shelxeMR", [
        "--F_SIGF", "fileOut=[1].HKLOUT[0]", "--FREERFLAG", "fileOut=[1].FREEROUT",
        "--XYZIN", "fileOut=[25].XYZOUT[0]"]),
]


def main(only=()):
    failed = []
    for project, task, args in RUNS:
        if only and task not in only:
            continue
        try:
            i2run(project, task, *args)
        except Exception as err:  # keep going: one task's failure is a finding
            failed.append(f"{project} {task}: {err}")
            continue
        print(task, "->", last_job_dir(project, task))
        clone_last(project, task)
    if failed:
        raise SystemExit("Failed:\n" + "\n".join(failed))


if __name__ == "__main__":
    main(sys.argv[1:])
