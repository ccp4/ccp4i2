"""Shared fixtures for the PanDDA staging tests.

The three-dataset BAZ2B subset in demo_data is source-shaped (dimple names),
so staging is exercised on real files. The ``value_order`` dictionary that
``prepare_dict_for_pandda`` exists for is manufactured here from a ``type``
one, deliberately: no shipped fixture has that spelling, and a test that
only ever saw ``type`` would pass vacuously (design note §15.1).
"""
import os

import pytest

import ccp4i2

MINI = os.path.join(os.path.dirname(ccp4i2.__file__), "demo_data", "pandda_baz2b_mini")
LABELS = ["BAZ2BA-x425", "BAZ2BA-x427", "BAZ2BA-x428"]

_TYPE_TO_PDBX = {"single": "SING", "double": "DOUB", "triple": "TRIP",
                 "aromatic": "AROM", "deloc": "DELO"}


def source_files(label):
    base = os.path.join(MINI, label)
    return (os.path.join(base, f"{label}.dimple.pdb"),
            os.path.join(base, f"{label}.dimple.mtz"),
            os.path.join(base, "ligand.cif"))


def pdbx_spelling(text: str) -> str:
    """A ``type``-style acedrg dictionary re-spelled the PDBx way:
    ``value_order`` + ``pdbx_aromatic_flag``, upper-case four-letter tokens."""
    text = (text.replace("_chem_comp_bond.type", "_chem_comp_bond.value_order")
                .replace("_chem_comp_bond.aromatic", "_chem_comp_bond.pdbx_aromatic_flag"))
    out = []
    for line in text.splitlines():
        fields = line.split()
        if len(fields) == 7 and fields[3] in _TYPE_TO_PDBX:
            fields[3] = _TYPE_TO_PDBX[fields[3]]
            fields[4] = fields[4].upper()
            line = "  ".join(fields)
        out.append(line)
    return "\n".join(out) + "\n"


@pytest.fixture
def mini_specs():
    from ccp4i2.wrappers.pandda_campaign.script.pandda_staging import DatasetSpec
    specs = []
    for label in LABELS:
        pdb, mtz, cif = source_files(label)
        specs.append(DatasetSpec(label=label, xyzin=pdb, hklin=mtz, dictionary=cif,
                                 project_uuid=f"00000000-0000-0000-0000-{LABELS.index(label):012d}"))
    return specs
