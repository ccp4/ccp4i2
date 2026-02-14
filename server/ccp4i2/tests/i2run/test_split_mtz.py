import gemmi
from .utils import demoData, i2run


def test_gamma_abcd_column_group_list():
    args = ["splitMtz"]
    args += ["--HKLIN", demoData("gamma", "initial_phases.mtz")]
    args += ["--COLUMNGROUPLIST", "columnGroupType=Phs", "contentFlag=1", "dataset=ds1", "selected=True"]
    args += mtzColumnArgs(["HLA", "HLB", "HLC", "HLD"], "A", "ds1")
    with i2run(args) as job:
        checkMtz(job / "ds1_HLA_HLB_HLC_HLD.mtz", ["HLA", "HLB", "HLC", "HLD"])


def test_gamma_abcd_user_column_group():
    args = ["splitMtz"]
    args += ["--HKLIN", demoData("gamma", "initial_phases.mtz")]
    args += ["--USERCOLUMNGROUP", "columnGroupType=Phs", "contentFlag=1", "dataset=ds1"]
    args += mtzColumnArgs(["HLA", "HLB", "HLC", "HLD"], "A", "ds1")
    with i2run(args) as job:
        checkMtz(job / "ds1_HLA_HLB_HLC_HLD.mtz", ["HLA", "HLB", "HLC", "HLD"])


def test_4iid_phi_fom():
    args = ["splitMtz"]
    args += ["--HKLIN", demoData("glyco", "4iid.mtz")]
    args += ["--USERCOLUMNGROUP", "columnGroupType=Phs", "contentFlag=2", "dataset=1"]
    args += mtzColumnArgs(["PHIC", "FOM"], ["P", "W"], "1")
    with i2run(args) as job:
        checkMtz(job / "1_PHIC_FOM.mtz", ["PHI", "FOM"])


def mtzColumnArgs(labels, types, datasets):
    args = []
    types = [types] * len(labels) if isinstance(types, str) else types
    datasets = [datasets] * len(labels) if isinstance(datasets, str) else datasets
    for i, (label, colType, dataset) in enumerate(zip(labels, types, datasets)):
        args += [
            f"columnList[{i}]/columnLabel={label}",
            f"columnList[{i}]/columnType={colType}",
            f"columnList[{i}]/dataset={dataset}",
        ]
    return args


def checkMtz(path, expected):
    mtz = gemmi.read_mtz_file(str(path))
    labels = set(col.label for col in mtz.columns)
    assert labels == {"H", "K", "L"} | set(expected)
