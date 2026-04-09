# Copyright (C) 2026 University of York
import gemmi
from .urls import redo_mtz
from .utils import demoData, download, i2run


# ---------------------------------------------------------------------------
# Existing tests: explicit column group selection
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# New tests: auto-detect mode (no column groups specified)
# ---------------------------------------------------------------------------


def test_auto_detect_redo_1cbs():
    """
    PDB-REDO MTZ files typically contain:
    - Obs (FP/SIGFP or similar)
    - FreeR flags
    - HL phases (HLA..HLD)
    - Map coefficients (FWT/PHWT, DELFWT/PHDELWT)

    When no column groups are specified, splitMtz should auto-detect
    and extract all of them.
    """
    with download(redo_mtz("1cbs")) as mtz:
        # Inspect input
        mtzin = gemmi.read_mtz_file(mtz)
        print(f"\n=== Input MTZ: {len(mtzin.columns)} columns ===")
        for col in mtzin.columns:
            print(f"  {col.label} (type={col.type})")

        # Run with no COLUMNGROUPLIST or USERCOLUMNGROUP — triggers auto-detect
        args = ["splitMtz", "--HKLIN", mtz]
        with i2run(args) as job:
            found_types = check_auto_detect_output(job)
            # PDB-REDO should have at least obs and FreeR
            assert "Obs" in found_types, f"Expected Obs data, found: {found_types}"
            assert "FreeR" in found_types, f"Expected FreeR data, found: {found_types}"


def test_auto_detect_redo_2ceu():
    """Test auto-detect with another PDB-REDO entry."""
    with download(redo_mtz("2ceu")) as mtz:
        args = ["splitMtz", "--HKLIN", mtz]
        with i2run(args) as job:
            found_types = check_auto_detect_output(job)
            assert "Obs" in found_types, f"Expected Obs data, found: {found_types}"


def test_auto_detect_gamma_phases():
    """Test auto-detect on a file with HL phases."""
    args = ["splitMtz", "--HKLIN", demoData("gamma", "initial_phases.mtz")]
    with i2run(args) as job:
        found_types = check_auto_detect_output(job)
        assert "Phs" in found_types, f"Expected phase data, found: {found_types}"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


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


def check_auto_detect_output(job):
    """
    Check auto-detect output: verify each output MTZ is valid.
    Returns set of column group types found (e.g. {"Obs", "FreeR", "Phs", "MapCoeffs"}).
    """
    # Type signatures → group type
    SIG_TO_TYPE = {
        "FQ": "Obs", "JQ": "Obs", "GLGL": "Obs", "KMKM": "Obs",
        "I": "FreeR",
        "AAAA": "Phs", "PW": "Phs",
        "FP": "MapCoeffs", "FQP": "MapCoeffs",
    }

    output_files = list(job.iterdir())
    mtz_files = [f for f in output_files if f.suffix == ".mtz" and f.name != "HKLIN.mtz"]

    print(f"\n=== Output files ({len(output_files)} total) ===")
    for f in sorted(output_files):
        print(f"  {f.name}")

    found_types = set()
    print(f"\n=== Output MTZ files: {len(mtz_files)} ===")
    for mtz_file in mtz_files:
        mtz = gemmi.read_mtz_file(str(mtz_file))
        data_cols = [c for c in mtz.columns if c.label not in ('H', 'K', 'L')]
        col_labels = [c.label for c in data_cols]
        type_sig = "".join(c.type for c in data_cols)
        print(f"  {mtz_file.name}: columns={col_labels} signature={type_sig}")

        # Verify the MTZ is valid
        assert len(mtz.columns) > 3, f"{mtz_file.name} has too few columns"
        assert mtz.nreflections > 0, f"{mtz_file.name} has no reflections"

        group_type = SIG_TO_TYPE.get(type_sig)
        if group_type:
            found_types.add(group_type)

    print(f"\n=== Found types: {found_types} ===")
    assert len(found_types) >= 1, "Should extract at least one column group"
    return found_types
