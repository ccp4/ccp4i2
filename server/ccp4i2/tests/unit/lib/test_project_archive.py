"""What a project zip has to look like to be importable.

The layout rule used to be implicit: ``import_ccp4_project_zip`` opened
``DATABASE.db.xml`` at the top of the archive and a zip built any other way
died on a bare ``KeyError``. Since ``zip -r proj.zip ./PROJ`` is the obvious
way to zip a project by hand, and produces the other shape, these tests pin
down both layouts as importable and the diagnosis for a zip that is neither.
"""

import zipfile
from pathlib import Path

import pytest

from ccp4i2.db.import_i2xml import (
    ProjectArchiveError,
    _extract_member,
    archive_prefix,
    inspect_ccp4_project_zip,
)

SNAPSHOT = "DATABASE.db.xml"

PROJECT_XML = """<?xml version='1.0' encoding='ASCII'?>
<ccp4:ccp4i2 xmlns:ccp4="http://www.ccp4.ac.uk/ccp4ns">
  <ccp4i2_header>
    <projectName>Demo</projectName>
    <projectId>f48588ee9ca611f1895cacde48001122</projectId>
  </ccp4i2_header>
  <ccp4i2_body>
    <projectTable>
      <project projectid="f48588ee9ca611f1895cacde48001122" projectname="Demo"
               projectdirectory="/elsewhere/CCP4I2_PROJECTS/Demo"/>
    </projectTable>
    <jobTable>
      <job jobid="a" jobnumber="1"/>
      <job jobid="b" jobnumber="2"/>
    </jobTable>
    <fileTable>
      <file fileid="c"/>
    </fileTable>
  </ccp4i2_body>
</ccp4:ccp4i2>
"""


def build_zip(path: Path, prefix: str = "") -> Path:
    """A minimal project archive, with its contents under ``prefix``."""
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr(f"{prefix}{SNAPSHOT}", PROJECT_XML)
        archive.writestr(f"{prefix}CCP4_IMPORTED_FILES/thing.pdb", "ATOM\n")
        archive.writestr(f"{prefix}CCP4_JOBS/job_1/stdout.txt", "hello\n")
    return path


class TestArchivePrefix:
    def test_exported_zip_has_no_prefix(self):
        assert archive_prefix([SNAPSHOT, f"CCP4_JOBS/job_1/x"]) == ""

    def test_hand_rolled_zip_wraps_the_project_in_a_directory(self):
        """``zip -r proj.zip ./PROJ`` — the common mistake, and importable."""
        names = [f"PROJ/{SNAPSHOT}", "PROJ/CCP4_JOBS/job_1/x"]
        assert archive_prefix(names) == "PROJ/"

    def test_leading_dot_slash_is_part_of_the_prefix(self):
        assert archive_prefix([f"./PROJ/{SNAPSHOT}"]) == "./PROJ/"

    def test_shallowest_snapshot_wins(self):
        """A project may hold an older zipped project; the outer one is meant."""
        names = [SNAPSHOT, f"CCP4_IMPORTED_FILES/old/{SNAPSHOT}"]
        assert archive_prefix(names) == ""

    def test_no_snapshot_is_an_explained_refusal(self):
        with pytest.raises(ProjectArchiveError) as excinfo:
            archive_prefix(["notes.txt", "data/thing.mtz"])
        assert SNAPSHOT in str(excinfo.value)

    def test_a_similarly_named_file_is_not_a_snapshot(self):
        with pytest.raises(ProjectArchiveError):
            archive_prefix(["MY_DATABASE.db.xml"])


class TestInspect:
    @pytest.mark.parametrize("prefix", ["", "PROJ/"])
    def test_reads_the_project_through_either_layout(self, tmp_path, prefix):
        archive = build_zip(tmp_path / "p.zip", prefix=prefix)
        summary = inspect_ccp4_project_zip(archive)
        assert summary["prefix"] == prefix
        assert summary["project_name"] == "Demo"
        assert summary["project_uuid"] == "f48588ee9ca611f1895cacde48001122"
        assert summary["recorded_directory"] == "/elsewhere/CCP4I2_PROJECTS/Demo"
        assert summary["jobs"] == 2
        assert summary["files"] == 1

    def test_a_zip_of_something_else_is_rejected_by_name(self, tmp_path):
        archive = tmp_path / "junk.zip"
        with zipfile.ZipFile(archive, "w") as handle:
            handle.writestr("junk/notes.txt", "hello")
        with pytest.raises(ProjectArchiveError):
            inspect_ccp4_project_zip(archive)

    def test_unreadable_snapshot_says_so(self, tmp_path):
        archive = tmp_path / "bad.zip"
        with zipfile.ZipFile(archive, "w") as handle:
            handle.writestr(SNAPSHOT, "<not-xml")
        with pytest.raises(ProjectArchiveError) as excinfo:
            inspect_ccp4_project_zip(archive)
        assert "not readable" in str(excinfo.value)

    def test_snapshot_without_a_project_says_so(self, tmp_path):
        archive = tmp_path / "empty.zip"
        with zipfile.ZipFile(archive, "w") as handle:
            handle.writestr(SNAPSHOT, "<ccp4i2><ccp4i2_body/></ccp4i2>")
        with pytest.raises(ProjectArchiveError) as excinfo:
            inspect_ccp4_project_zip(archive)
        assert "names no project" in str(excinfo.value)

    def test_not_a_zip_at_all(self, tmp_path):
        plain = tmp_path / "notes.txt"
        plain.write_text("hello")
        with pytest.raises(ProjectArchiveError):
            inspect_ccp4_project_zip(plain)


class TestExtractMember:
    """Stripping a prefix by hand is what lets ``../`` out, so guard it."""

    def test_writes_inside_the_project(self, tmp_path):
        archive_path = tmp_path / "p.zip"
        with zipfile.ZipFile(archive_path, "w") as handle:
            handle.writestr("CCP4_JOBS/job_1/stdout.txt", "hello\n")
        root = tmp_path / "project"
        root.mkdir()
        with zipfile.ZipFile(archive_path) as handle:
            _extract_member(
                handle,
                "CCP4_JOBS/job_1/stdout.txt",
                root / "CCP4_JOBS/job_1/stdout.txt",
                root,
            )
        assert (root / "CCP4_JOBS/job_1/stdout.txt").read_text() == "hello\n"

    def test_refuses_to_write_outside_the_project(self, tmp_path):
        archive_path = tmp_path / "evil.zip"
        with zipfile.ZipFile(archive_path, "w") as handle:
            handle.writestr("CCP4_JOBS/../../escaped.txt", "pwned\n")
        root = tmp_path / "project"
        root.mkdir()
        with zipfile.ZipFile(archive_path) as handle:
            with pytest.raises(ProjectArchiveError):
                _extract_member(
                    handle,
                    "CCP4_JOBS/../../escaped.txt",
                    root / "CCP4_JOBS/../../escaped.txt",
                    root,
                )
        assert not (tmp_path.parent / "escaped.txt").exists()
