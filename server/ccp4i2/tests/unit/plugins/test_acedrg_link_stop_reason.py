"""When AceDRG refuses a link, the job says why in AceDRG's own words.

AceDRG writes its reason to <LINK_ID>_errorInfo.txt; the job used to report
only the log lines before it, under a FileNotFoundError from MakeLink going
on to look for the dictionary AceDRG had not written.
"""
from ccp4i2.wrappers.AcedrgLink.script.AcedrgLink import acedrg_stop_reason


def test_the_reason_is_read_and_flattened(tmp_path):
    (tmp_path / "LYS-GLU_errorInfo.txt").write_text(
        "atom C in monomer GLU has a total valence of 3,\n which is not allowed!\n\n")
    assert acedrg_stop_reason(tmp_path, "LYS-GLU") == (
        "atom C in monomer GLU has a total valence of 3, which is not allowed!")


def test_no_file_no_reason(tmp_path):
    assert acedrg_stop_reason(tmp_path, "LYS-GLU") is None


def test_an_empty_file_is_no_reason(tmp_path):
    (tmp_path / "LYS-GLU_errorInfo.txt").write_text("\n")
    assert acedrg_stop_reason(tmp_path, "LYS-GLU") is None
