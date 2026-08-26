"""
CMiniMtzDataFile.miniMtzType() -- what kind of mini-MTZ is this file?

mergeMtz has called this for years and it did not exist. Nobody saw, because
the task could not get that far: the CLI silently accepted `filename=` where
the field is `fileName`, so both inputs stayed unset, the wrapper skipped them
both, joinMtz was handed an empty list, and the job reported success having
written nothing.

Needs gemmi (available CCP4-free), no CCP4 binaries.
"""
import pytest

gemmi = pytest.importorskip("gemmi")

from ccp4i2.core.CCP4XtalData import (
    CFreeRDataFile,
    CMapCoeffsDataFile,
    CMiniMtzDataFile,
    CObsDataFile,
    CPhsDataFile,
)


def _mtz(tmp_path, name, columns):
    """A minimal MTZ carrying exactly *columns* beside H, K, L."""
    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = gemmi.SpaceGroup("P 1")
    mtz.set_cell_for_all(gemmi.UnitCell(30, 30, 30, 90, 90, 90))
    mtz.add_dataset("data")
    for label in columns:
        mtz.add_column(label, "R")
    import numpy as np
    rows = np.array([[1.0, 2.0, 3.0] + [1.0] * len(columns),
                     [2.0, 3.0, 4.0] + [2.0] * len(columns)], dtype="float32")
    mtz.set_data(rows)
    path = tmp_path / name
    mtz.write_to_file(str(path))

    obj = CMiniMtzDataFile()
    obj.setFullPath(str(path))
    return obj


@pytest.mark.parametrize("columns,expected_class,expected_flag", [
    (["Iplus", "SIGIplus", "Iminus", "SIGIminus"], CObsDataFile, 1),
    (["Fplus", "SIGFplus", "Fminus", "SIGFminus"], CObsDataFile, 2),
    (["I", "SIGI"], CObsDataFile, 3),
    (["F", "SIGF"], CObsDataFile, 4),
    (["HLA", "HLB", "HLC", "HLD"], CPhsDataFile, 1),
    (["PHI", "FOM"], CPhsDataFile, 2),
    (["F", "PHI"], CMapCoeffsDataFile, 1),
    (["FREER"], CFreeRDataFile, 1),
])
def test_each_signature_is_recognised(tmp_path, columns, expected_class, expected_flag):
    obj = _mtz(tmp_path, "mini.mtz", columns)
    cls, flag = obj.miniMtzType()
    assert cls is expected_class
    assert flag == expected_flag


def test_columns_decide_not_the_content_flag(tmp_path):
    """The answer comes from the file, so it is right even when nothing set a
    contentFlag -- which is the case for a file named on the command line."""
    obj = _mtz(tmp_path, "mini.mtz", ["FREER"])
    assert not obj.contentFlag.isSet()
    assert obj.miniMtzType() == (CFreeRDataFile, 1)


def test_unrecognised_columns_are_not_guessed(tmp_path):
    obj = _mtz(tmp_path, "odd.mtz", ["WIBBLE", "SIGWIBBLE"])
    assert obj.miniMtzType() == (None, None)


def test_a_missing_file_is_not_an_exception(tmp_path):
    obj = CMiniMtzDataFile()
    obj.setFullPath(str(tmp_path / "absent.mtz"))
    assert obj.miniMtzType() == (None, None)
