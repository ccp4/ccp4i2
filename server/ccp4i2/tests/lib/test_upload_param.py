"""Tests for upload_param utilities, especially annotation building."""

from ...lib.utils.files.upload_param import build_file_annotation


class TestBuildFileAnnotation:
    """Tests for build_file_annotation function."""

    def test_no_metadata(self):
        """Basic annotation with no MTZ metadata."""
        result = build_file_annotation("test.mtz")
        assert result == "Imported from test.mtz"

    def test_full_metadata(self):
        """Annotation with all metadata fields populated."""
        metadata = {
            "crystal_name": "xtal1",
            "dataset_name": "native",
            "original_columns": ["F", "SIGF"],
        }
        result = build_file_annotation("test.mtz", metadata)
        assert result == "Imported from test.mtz; Crystal: xtal1; Dataset: native; Columns: F, SIGF"

    def test_partial_metadata_no_crystal(self):
        """Annotation with crystal_name as None."""
        metadata = {
            "crystal_name": None,
            "dataset_name": "peak",
            "original_columns": ["I", "SIGI"],
        }
        result = build_file_annotation("data.mtz", metadata)
        assert result == "Imported from data.mtz; Dataset: peak; Columns: I, SIGI"

    def test_only_columns(self):
        """Annotation with only column information."""
        metadata = {
            "crystal_name": None,
            "dataset_name": None,
            "original_columns": ["FP", "SIGFP"],
        }
        result = build_file_annotation("obs.mtz", metadata)
        assert result == "Imported from obs.mtz; Columns: FP, SIGFP"

    def test_only_crystal(self):
        """Annotation with only crystal name."""
        metadata = {
            "crystal_name": "my_crystal",
            "dataset_name": None,
            "original_columns": None,
        }
        result = build_file_annotation("file.mtz", metadata)
        assert result == "Imported from file.mtz; Crystal: my_crystal"

    def test_empty_metadata_dict(self):
        """Annotation with empty metadata dict (all None/empty)."""
        metadata = {
            "crystal_name": None,
            "dataset_name": None,
            "original_columns": None,
        }
        result = build_file_annotation("file.mtz", metadata)
        # No additional info, just the base annotation
        assert result == "Imported from file.mtz"

    def test_empty_columns_list(self):
        """Annotation with empty columns list."""
        metadata = {
            "crystal_name": "xtal",
            "dataset_name": "data",
            "original_columns": [],
        }
        result = build_file_annotation("file.mtz", metadata)
        # Empty list is falsy, so columns part is omitted
        assert result == "Imported from file.mtz; Crystal: xtal; Dataset: data"

    def test_multiple_columns(self):
        """Annotation with multiple columns."""
        metadata = {
            "crystal_name": "crystal1",
            "dataset_name": "native",
            "original_columns": ["F", "SIGF", "FreeR_flag", "PHI", "FOM"],
        }
        result = build_file_annotation("complex.mtz", metadata)
        assert result == "Imported from complex.mtz; Crystal: crystal1; Dataset: native; Columns: F, SIGF, FreeR_flag, PHI, FOM"
