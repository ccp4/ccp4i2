"""Hydrogens a search model should not carry into refinement."""
from pathlib import Path

import pytest

from ccp4i2.core.tasks import get_plugin_class
from ccp4i2.lib.utils.formats.hydrogens import is_hydrogen_record, strip_hydrogens

DEMO = Path(__file__).resolve().parents[3] / "demo_data" / "beta_blip"


def _record(name, element="", kind="ATOM  "):
    return f"{kind}    1 {name:4s} SER A   1       0.000   0.000   0.000  1.00 20.00          {element:>2s}\n"


@pytest.mark.parametrize("name, element, expected", [
    (" H  ", "", True), (" HA ", "", True), ("1HB ", "", True), ("HD21", "", True),
    (" CA ", "", False), (" N  ", "", False), ("HG  ", "", False),     # mercury, column 13
    (" HG ", "", True),                                               # serine's hydroxyl hydrogen
    (" HG ", "H", True), ("HG  ", "HG", False), (" D  ", "D", True), (" CA ", "C", False),
])
def test_a_hydrogen_record_is_known_by_element_or_name(name, element, expected):
    assert is_hydrogen_record(_record(name, element)) is expected
    assert is_hydrogen_record("REMARK   1 HA2 HA3\n") is False


def test_strip_hydrogens_counts_what_it_removes(tmp_path):
    model = tmp_path / "with_h.pdb"
    model.write_text(_record(" N  ", "N") + _record(" CA ", "C") + _record(" HA2", "H")
                     + _record(" HA3", "H") + "END\n")
    out = tmp_path / "without.pdb"
    assert strip_hydrogens(model, out) == 2
    atoms = [l for l in out.read_text().splitlines() if l.startswith("ATOM")]
    assert len(atoms) == 2 and not any(is_hydrogen_record(l + "\n") for l in atoms)
    assert strip_hydrogens(DEMO / "blip.pdb", tmp_path / "blip.pdb") == 0


def test_chainsaw_drops_hydrogens_and_keeps_the_cell(tmp_path):
    # Opus (2026-10-09): a Chainsaw model kept glycine's HA2/HA3 on residues
    # it mutated, and REFMAC refused it (error 350).
    plugin = get_plugin_class("chainsaw")(workDirectory=str(tmp_path), name="chainsaw")
    out = tmp_path / "XYZOUT.pdb"
    out.write_text(_record(" N  ") + _record(" CA ") + _record(" HA2") + _record("HD21") + "END\n")
    plugin.container.outputData.XYZOUT.setFullPath(str(out))
    plugin.cryst1card = "CRYST1   50.000   60.000   70.000  90.00  90.00  90.00 P 21 21 21\n"
    plugin.processOutputFiles()
    lines = out.read_text().splitlines()
    assert lines[0].startswith("CRYST1   50.000")
    assert [l[12:16] for l in lines if l.startswith("ATOM")] == [" N  ", " CA "]
