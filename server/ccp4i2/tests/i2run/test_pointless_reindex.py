from .utils import demoData, i2run


def test_reindex_to_coords():
    """Test pointless reindexing gamma data to match coordinate reference."""
    args = ["pointless_reindexToMatch"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN_REF", demoData("gamma", "gamma_model.pdb")]
    args += ["--REFERENCE", "XYZIN_REF"]

    with i2run(args) as job:
        mtz_files = list(job.glob("*.mtz"))
        assert len(mtz_files) >= 1, f"No MTZ output from pointless: {list(job.iterdir())}"


def test_reindex_annotations_say_the_operator():
    """Each output says the operator once: the given one, or the one found.

    It appended "Reindexed:[...];SG:...;Resolution:...;Cell:..." to the
    input's annotation, so a second run lengthened the first one's, and a
    match against a reference said "NewSG:..." without the operator found."""
    from ccp4i2.db import models

    def annotations(job):
        record = models.Job.objects.filter(number=job.name.replace("job_", "")).first()
        return {f.job_param_name: f.annotation for f in models.File.objects.filter(job=record)}

    args = ["pointless_reindexToMatch"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--REFERENCE", "SPECIFY", "--USE_REINDEX", "True"]
    args += ["--REINDEX_OPERATOR", "h=k", "k=h", "l=-l"]
    with i2run(args) as job:
        out = annotations(job)["F_SIGF_OUT"]
        assert out.count("reindexed by") == 1 and "reindexed by k,h,-l (" in out, out
        assert "Reindexed:" not in out and "Cell:" not in out, out
