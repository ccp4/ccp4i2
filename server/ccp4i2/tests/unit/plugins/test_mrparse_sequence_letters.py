"""MrParse is given a sequence it accepts: non-standard letters replaced.

MrParse refuses any letter outside the 20 amino acids ("'X' is not a valid
unambiguous letter for protein"); Haiku's CDK4/cyclin D1 mrparse job failed
on an X (2026-10-13).
"""
from pathlib import Path

from ccp4i2.core.tasks import get_plugin_class


def _plugin(tmp_path, text):
    seq = tmp_path / "target.fasta"
    seq.write_text(text)
    plugin = get_plugin_class("mrparse")(workDirectory=str(tmp_path), name="m")
    plugin.container.inputData.SEQIN.setFullPath(str(seq))
    return plugin, seq


def test_non_standard_letters_are_replaced_keeping_the_length(tmp_path):
    plugin, seq = _plugin(tmp_path, ">cdk4\nMATSRXEPVAEIGBZ*\n")
    used = plugin._sequence_for_mrparse()
    assert used != str(seq)
    body = Path(used).read_text().splitlines()[1]
    assert body == "MATSRAEPVAEIGDE" and len(body) == len("MATSRXEPVAEIGBZ")


def test_a_clean_sequence_is_passed_as_it_is(tmp_path):
    plugin, seq = _plugin(tmp_path, ">gamma\nMKVLAAGIVGLL\n")
    assert plugin._sequence_for_mrparse() == str(seq)
