"""What the run check costs per dataset. A 158-dataset campaign on DDU spent
27 s warm in ``runTimeValidity`` (and past the ingress timeout cold), reading
every reflection file in full and three times over. Each file is now opened
once, header-only, for the whole of a check."""
import shutil

import gemmi
import pytest

from ccp4i2.core.tasks import get_plugin_class

from .conftest import LABELS, source_files


@pytest.fixture
def reads(monkeypatch):
    """Every gemmi.read_mtz_file call the plugin makes: (path, with_data)."""
    calls = []
    real = gemmi.read_mtz_file

    def counting(path, *args, **kwargs):
        calls.append((str(path), kwargs.get("with_data", args[0] if args else True)))
        return real(path, *args, **kwargs)

    monkeypatch.setattr(gemmi, "read_mtz_file", counting)
    return calls


@pytest.fixture
def plugin(tmp_path):
    plugin = get_plugin_class("pandda_campaign")(workDirectory=str(tmp_path))
    datasets = plugin.container.inputData.DATASETS
    for label in LABELS:
        pdb, mtz, _cif = source_files(label)
        item = datasets.makeItem()
        item.DTAG.set(label)
        item.XYZIN.setFullPath(pdb)
        item.HKLIN.setFullPath(mtz)
        datasets.append(item)
    return plugin


def test_a_check_opens_each_reflection_file_once_and_header_only(plugin, reads):
    outlook = plugin._resolution_outlook()
    resolutions = plugin._dataset_resolutions()
    hint = plugin._sizing_hint()

    assert len(reads) == len(LABELS), reads
    assert {with_data for _path, with_data in reads} == {False}
    assert outlook[1] == max(resolutions)
    expected = [gemmi.read_mtz_file(source_files(l)[1]).resolution_high() for l in LABELS]
    assert resolutions == pytest.approx(expected)
    assert hint["datasets"] == len(LABELS) and hint["cell_volume_class"]


def test_a_changed_list_is_read_again(plugin, reads, tmp_path):
    plugin._dataset_resolutions()
    assert len(reads) == len(LABELS)

    extra = tmp_path / "extra.mtz"
    shutil.copy(source_files(LABELS[0])[1], extra)
    item = plugin.container.inputData.DATASETS.makeItem()
    item.DTAG.set("extra")
    item.HKLIN.setFullPath(str(extra))
    plugin.container.inputData.DATASETS.append(item)

    assert len(plugin._dataset_resolutions()) == len(LABELS) + 1
    assert len(reads) == 2 * len(LABELS) + 1


def test_an_unreadable_file_is_left_out_not_fatal(plugin, reads, tmp_path):
    bad = tmp_path / "bad.mtz"
    bad.write_bytes(b"not an mtz")
    item = plugin.container.inputData.DATASETS.makeItem()
    item.DTAG.set("bad")
    item.HKLIN.setFullPath(str(bad))
    plugin.container.inputData.DATASETS.append(item)
    assert len(plugin._dataset_resolutions()) == len(LABELS)


def test_the_campaign_summary_reads_headers_only(reads, tmp_path):
    from ccp4i2.lib.pandda_export import dataset_resolution

    class Job:
        directory = tmp_path
    shutil.copy(source_files(LABELS[0])[1], tmp_path / "final.mtz")
    assert dataset_resolution(Job()) == pytest.approx(
        gemmi.read_mtz_file(source_files(LABELS[0])[1]).resolution_high(), abs=0.01)
    assert reads[0][1] is False
