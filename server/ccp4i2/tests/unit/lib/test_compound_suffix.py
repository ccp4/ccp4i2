"""Compound file extensions must survive an upload intact.

Several CCP4i2 datatypes are identified by a two-part extension -- most
visibly ``.asu.xml`` for ASU contents. ``pathlib.Path.suffix`` sees only the
last part, so the upload path used to slugify the stem "gamma.asu" into
"gamma-asu" and store the file as ``gammaasu.xml``, silently losing the
extension the datatype is named by.
"""
import pathlib

from ccp4i2.lib.utils.files.available_name import (
    COMPOUND_SUFFIXES,
    available_file_name_based_on,
    split_compound_suffix,
)


class TestSplitCompoundSuffix:
    def test_asu_xml_keeps_both_parts(self):
        assert split_compound_suffix("gamma.asu.xml") == ("gamma", ".asu.xml")

    def test_ordinary_extension_is_unchanged(self):
        assert split_compound_suffix("rnase.pdb") == ("rnase", ".pdb")
        assert split_compound_suffix("lib.cif") == ("lib", ".cif")

    def test_every_declared_compound_suffix_round_trips(self):
        for suffix in COMPOUND_SUFFIXES:
            stem, got = split_compound_suffix(f"myfile{suffix}")
            assert (stem, got) == ("myfile", suffix)

    def test_matching_is_case_insensitive_but_preserves_case(self):
        assert split_compound_suffix("X.ASU.XML") == ("X", ".ASU.XML")

    def test_leading_directories_are_ignored(self):
        assert split_compound_suffix("/tmp/d/gamma.asu.xml") == (
            "gamma",
            ".asu.xml",
        )

    def test_name_that_is_only_the_suffix_is_not_left_stemless(self):
        # ".asu.xml" with nothing in front must not produce an empty stem.
        stem, _ = split_compound_suffix(".asu.xml")
        assert stem != ""

    def test_no_extension(self):
        assert split_compound_suffix("noext") == ("noext", "")


class TestAvailableFileName:
    def test_deduplication_preserves_a_compound_suffix(self, tmp_path):
        original = tmp_path / "gamma.asu.xml"
        original.write_text("<xml/>")

        dest = available_file_name_based_on(original)

        assert dest != original, "should have picked a free name"
        assert dest.name == "gamma_1.asu.xml"

    def test_deduplication_of_an_ordinary_suffix(self, tmp_path):
        original = tmp_path / "model.pdb"
        original.write_text("ATOM")

        assert available_file_name_based_on(original).name == "model_1.pdb"

    def test_free_name_is_returned_unchanged(self, tmp_path):
        target = tmp_path / "fresh.asu.xml"
        assert available_file_name_based_on(target).name == "fresh.asu.xml"
