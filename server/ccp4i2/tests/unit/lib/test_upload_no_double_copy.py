"""An upload that derives nothing must not be stored twice.

Reflection uploads stage into CCP4_DOWNLOADED_FILES first, because splitting an
MTZ down to the columns one parameter wants needs the whole file to read from.
When nothing is split --- an unmerged dataset, or a general CMtzDataFile with
no column selector --- the imported file is byte-identical to the staged one,
and the staged copy is never cleaned up. Every such upload cost twice its size
on disk, forever, and unmerged datasets are the large ones.

Pure Python -- no CCP4 binaries and no database needed. The no-split branch
never reads the reflection data, so the bytes here need not be a real MTZ.
"""
from pathlib import Path

from ccp4i2.core.CCP4XtalData import CUnmergedMtzDataFile
from ccp4i2.lib.utils.files.upload_param import handle_reflections


class _JobIn:
    """Just enough of a Job for the destination path."""

    def __init__(self, project_directory):
        self.project = type("_Project", (), {"directory": str(project_directory)})()


def _staged(project_root, name, payload):
    staged_dir = project_root / "CCP4_DOWNLOADED_FILES"
    staged_dir.mkdir(parents=True, exist_ok=True)
    (project_root / "CCP4_IMPORTED_FILES").mkdir(parents=True, exist_ok=True)
    path = staged_dir / name
    path.write_bytes(payload)
    return path


def test_an_unsplit_upload_leaves_nothing_behind_in_staging(tmp_path):
    payload = b"MTZ\x00" + bytes(range(256)) * 64
    staged = _staged(tmp_path, "hypf_unmerged.mtz", payload)

    result = handle_reflections(
        _JobIn(tmp_path), CUnmergedMtzDataFile(), "hypf_unmerged.mtz", None, staged
    )

    imported = Path(result["path"])
    assert imported.parent == tmp_path / "CCP4_IMPORTED_FILES"
    assert imported.read_bytes() == payload
    assert not staged.exists(), "the staged upload was left behind as a second copy"


def test_the_project_holds_exactly_one_copy_of_the_bytes(tmp_path):
    """The point of the change, stated as the property that matters."""
    payload = b"MTZ\x00unmerged" * 4096
    staged = _staged(tmp_path, "gamma.mtz", payload)

    handle_reflections(
        _JobIn(tmp_path), CUnmergedMtzDataFile(), "gamma.mtz", None, staged
    )

    copies = [
        path for path in tmp_path.rglob("*")
        if path.is_file() and path.read_bytes() == payload
    ]
    assert len(copies) == 1, f"expected one copy, found {[str(c) for c in copies]}"
