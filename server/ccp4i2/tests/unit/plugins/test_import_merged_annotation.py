"""The imported observations say what they are and where they came from.

An import_merged job's OBSOUT reached the next job's pull-down as a bare
"1: " -- no annotation at all, so a project with several imports offered
several indistinguishable lines (Paul Bond, issue #510). The content type and
the source file were both known by then; neither got as far as the label,
because ``columnthings()`` set the annotation only on the branch where the
interface had *not* already chosen the columns -- the branch its own comment
marks "should not happen". Every real import took the other one.

``annotateObsOut()`` now runs from ``nearlyDone()``, the single point every
import format passes through, just after the final contentFlag is settled.
"""
import pytest

from ccp4i2.core.tasks import get_plugin_class


@pytest.fixture
def plugin(tmp_path):
    cls = get_plugin_class("import_merged")
    return cls(workDirectory=str(tmp_path), name="import_merged")


def test_annotation_names_content_source_crystal_and_dataset(plugin, tmp_path):
    """Paul's case: anomalous intensities from a named crystal and dataset."""
    obs_out = plugin.container.outputData.OBSOUT
    obs_out.contentFlag.set(1)  # CONTENT_FLAG_IPAIR
    plugin.container.inputData.HKLIN.set(
        str(tmp_path / "merged_intensities_native.mtz"))
    plugin.container.inputData.CRYSTALNAME.set("Nat3")
    plugin.container.inputData.DATASETNAME.set("nat3")

    plugin.annotateObsOut()

    assert str(obs_out.annotation) == (
        "Anomalous Is from merged_intensities_native.mtz; "
        "Crystal: Nat3; Dataset: nat3")


def test_annotation_keeps_the_file_extension(plugin, tmp_path):
    """stripedName() would drop the '.mtz'; the annotation keeps it."""
    plugin.container.outputData.OBSOUT.contentFlag.set(3)  # mean Is
    plugin.container.inputData.HKLIN.set(str(tmp_path / "gamma.mtz"))

    plugin.annotateObsOut()

    assert str(plugin.container.outputData.OBSOUT.annotation) == \
        "Mean Is from gamma.mtz"


def test_unset_content_flag_does_not_guess(plugin):
    """int() of an unset contentFlag is 0, and CONTENT_ANNOTATION[-1] would
    then confidently report 'Mean SFs' for data of unknown type."""
    plugin.annotateObsOut()

    annotation = str(plugin.container.outputData.OBSOUT.annotation)
    assert "Mean SFs" not in annotation
    assert annotation == plugin.container.outputData.OBSOUT.qualifiers("guiLabel")


def test_an_existing_annotation_is_left_alone(plugin, tmp_path):
    """A rerun must not overwrite an annotation the user (or an earlier, more
    specific step) has already given the file."""
    obs_out = plugin.container.outputData.OBSOUT
    obs_out.annotation.set("Native data, Feb 2026 beamtime")
    obs_out.contentFlag.set(4)
    plugin.container.inputData.HKLIN.set(str(tmp_path / "gamma.mtz"))

    plugin.annotateObsOut()

    assert str(obs_out.annotation) == "Native data, Feb 2026 beamtime"
