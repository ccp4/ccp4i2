"""print() of a character the console encoding lacks must not raise (freerflag, 2026-10-22)."""
import io
import sys

from ccp4i2.lib.utils.stdio import tolerate_unencodable_output


def test_unencodable_character_is_escaped_not_raised(monkeypatch):
    raw = io.BytesIO()
    cp1252 = io.TextIOWrapper(raw, encoding="cp1252", newline="\n")
    monkeypatch.setattr(sys, "stdout", cp1252)
    monkeypatch.setattr(sys, "stderr", object())  # no reconfigure(): skipped

    tolerate_unencodable_output()
    print("Cell: \u03b1=90.00")
    cp1252.flush()

    assert raw.getvalue() == b"Cell: \u03b1=90.00\n"
