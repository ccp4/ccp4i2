import xml.etree.ElementTree as ET
from .utils import i2run


_GAMMA_FASTA = (
    ">gamma_chymotrypsinogen\n"
    "MKYLLPTAAAGLLLLAAQPAMAMTQKIIMKDGLPSDSKLIQFEWNNPDQ"
    "FQKDRISCGQVSIPNNGDCIYGDIAFKAGAWEALFNAAGATFESAERRQ"
    "ADKEGCYREATTCLPLFTIFCNAFMPQIHGAVEGYHWDKYGGMCDDRYC"
    "GTIMGGMYAGGCQAVKSAAADATFNIKNIFESRKIDGLMSQMMCGAAYPF"
    "YVKDGENRCCQAANFYKSGAVILTHPGSEPNLTYWKF\n"
)


def test_freetext_sequence():
    """Test ProvideSequence with pasted sequence text."""
    args = ["ProvideSequence"]
    args += ["--SEQUENCETEXT", _GAMMA_FASTA]
    with i2run(args) as job:
        asu = job / "CASUCONTENTOUT.asu.xml"
        assert asu.exists(), f"No ASU content file: {list(job.iterdir())}"
        tree = ET.parse(asu)
        seqs = tree.findall(".//sequence")
        assert len(seqs) >= 1, "No sequences found in ASU output"


def test_freetext_sequence_annotations():
    """The outputs are named from the sequence's FASTA title, not left blank."""
    from ccp4i2.db import models
    with i2run(["ProvideSequence", "--SEQUENCETEXT", _GAMMA_FASTA]) as job:
        record = models.Job.objects.filter(number=job.name.replace("job_", "")).first()
        found = {f.job_param_name: f.annotation for f in models.File.objects.filter(job=record)}
        assert "gamma_chymotrypsinogen" in (found.get("CASUCONTENTOUT") or ""), found
