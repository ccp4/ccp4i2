import xml.etree.ElementTree as ET
from .utils import i2run


_GAMMA_SEQ = (
    "MKYLLPTAAAGLLLLAAQPAMAMTQKIIMKDGLPSDSKLIQFEWNNPDQ"
    "FQKDRISCGQVSIPNNGDCIYGDIAFKAGAWEALFNAAGATFESAERRQ"
    "ADKEGCYREATTCLPLFTIFCNAFMPQIHGAVEGYHWDKYGGMCDDRYC"
    "GTIMGGMYAGGCQAVKSAAADATFNIKNIFESRKIDGLMSQMMCGAAYPF"
    "YVKDGENRCCQAANFYKSGAVILTHPGSEPNLTYWKF"
)


def test_freetext_sequence():
    """Test ProvideSequence with pasted sequence text."""
    args = ["ProvideSequence"]
    args += ["--SEQUENCETEXT", _GAMMA_SEQ]
    with i2run(args) as job:
        asu = job / "CASUCONTENTOUT.asu.xml"
        assert asu.exists(), f"No ASU content file: {list(job.iterdir())}"
        tree = ET.parse(asu)
        seqs = tree.findall(".//sequence")
        assert len(seqs) >= 1, "No sequences found in ASU output"
