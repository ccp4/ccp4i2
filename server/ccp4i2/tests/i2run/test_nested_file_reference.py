"""A file inside a list item can be named by where it came from.

A Phaser ensemble's structure is a file inside a list item
("pdbItemList/structure"). fileOut= was resolved only for a file argument
itself; given there, the keyword set a dead attribute and left the file unset,
so a scenario had to look up the file's database id instead. The two runs
nest: the harness's test database belongs to the outermost context.
"""
from .utils import demoData, i2run


def test_ensemble_structure_named_by_file_out():
    from ccp4i2.db import models

    with i2run(["ImportCoordinate", "--XYZIN", f"fullPath={demoData('beta_blip', 'beta.pdb')}"]) as imported:
        args = ["phaser_mr_frf_phil",
                "--F_SIGF", f"fullPath={demoData('beta_blip', 'beta_blip_P3221.mtz')}",
                "columnLabels=/*/*/[Fobs,Sigma]",
                "--ENSEMBLES", "label=beta", "use=True", "number=1",
                "pdbItemList/identity_to_target=0.9",
                "pdbItemList/structure/fileOut=ImportCoordinate[-1].XYZOUT",
                "--COMP_BY", "DEFAULT"]
        with i2run(args) as frf:
            assert (frf / "RFILEOUT.phaser_rlist.pkl").is_file()
            source = models.Job.objects.get(number=imported.name.replace("job_", ""))
            job = models.Job.objects.get(number=frf.name.replace("job_", ""))
            used = [u.file for u in models.FileUse.objects.filter(job=job)]
            assert any(f.job_id == source.id and f.job_param_name == "XYZOUT" for f in used), \
                [(f.job.number, f.job_param_name) for f in used]
