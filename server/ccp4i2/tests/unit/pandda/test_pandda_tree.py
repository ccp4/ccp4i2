"""The receipt's tree reader (design note section 7): what PanDDA declared,
what arrived, and the difference, with the two facts of section 7.3 pinned.
CCP4-free: gemmi for tiny maps, yaml for events.yaml."""
import pytest

pytest.importorskip("gemmi", reason="needs gemmi")
pytest.importorskip("yaml", reason="needs PyYAML")

from ccp4i2.wrappers.pandda_events.script.pandda_tree import (
    DatasetNotFound, find_event_map, read_dataset, read_events_yaml, read_ligand_id)
from .synthetic_tree import event_record, make_tree, write_map


@pytest.fixture
def tree(tmp_path):
    return make_tree(tmp_path / "pandda2_out", {
        "xtal-0004": [event_record(1, bdc=0.80), event_record(2, bdc=0.95, build=False)],
        "xtal-0000": [],
    })


def test_reads_every_declared_thing(tree):
    ds = read_dataset(tree, "xtal-0004")
    assert ds.apo_model is not None and ds.apo_model.is_file()
    assert not ds.apo_model.is_symlink(), "resolved through PanDDA's symlink"
    assert ds.zmap is not None and ds.mean_map is not None
    assert ds.pandda_model is None
    assert [e.idx for e in ds.events] == [1, 2]
    one, two = ds.events
    assert one.bdc == pytest.approx(0.80) and one.centroid == (15.0, 40.0, 30.0)
    assert one.build is not None and one.pose is not None and one.pose.is_file()
    assert one.build.optimal_contour == pytest.approx(0.98)
    assert one.build.rscc == pytest.approx(0.30)
    assert two.build is None and two.pose is None, "no build means no pose expected"
    assert one.event_map is not None and two.event_map is not None
    assert ds.events_table_present
    assert one.site_idx == 1 and two.site_idx == 2
    assert one.hit_probability == pytest.approx(0.9)
    assert ds.shortfalls() == []


def test_event_map_is_matched_on_index_not_bdc_token(tree):
    ddir = tree / "processed_datasets" / "xtal-0004"
    # the token in the name is 1-BDC rounded: 0.95 -> "0.05"; never parse it
    assert find_event_map(ddir, "xtal-0004", 2).name == "xtal-0004-event_2_1-BDC_0.05_map.native.ccp4"
    assert find_event_map(ddir, "xtal-0004", 3) is None


def test_zero_event_dataset_is_complete(tree):
    ds = read_dataset(tree, "xtal-0000")
    assert ds.events == []
    assert ds.n_events == 0 and ds.shortfalls() == []
    assert read_events_yaml(tree / "processed_datasets" / "xtal-0000" / "events.yaml") == {}


def test_missing_dataset_directory_is_the_one_hard_failure(tree):
    with pytest.raises(DatasetNotFound):
        read_dataset(tree, "xtal-9999")


def test_shortfalls_name_what_did_not_arrive(tmp_path):
    tree = make_tree(tmp_path / "out", {
        "xtal-0007": [
            event_record(1),
            (event_record(2), {"event_map": False}),
            (event_record(3), {"pose": False}),
        ],
    }, staged_apo=False)
    ds = read_dataset(tree, "xtal-0007")
    assert ds.apo_model is None, "a dangling -pandda-input.pdb symlink is a missing file"
    assert ds.n_events == 3 and ds.n_event_maps == 2
    assert ds.n_poses_expected == 3 and ds.n_poses == 2
    assert ds.shortfalls() == [
        "apo model (-pandda-input.pdb)",
        "event 2: event map",
        "event 3: candidate pose",
    ]


def test_partial_run_without_events_table(tmp_path):
    tree = make_tree(tmp_path / "out", {"xtal-0001": [event_record(1)]}, events_table=False)
    ds = read_dataset(tree, "xtal-0001")
    assert not ds.events_table_present
    assert ds.events[0].site_idx is None and ds.events[0].hit_probability is None
    assert ds.shortfalls() == []


def test_missing_events_yaml_reads_as_zero_events(tmp_path):
    tree = make_tree(tmp_path / "out", {"xtal-0002": []})
    (tree / "processed_datasets" / "xtal-0002" / "events.yaml").unlink()
    assert read_dataset(tree, "xtal-0002").events == []


def test_pose_falls_back_to_recorded_build_path(tmp_path):
    tree = make_tree(tmp_path / "out", {"xtal-0003": [event_record(1)]})
    ddir = tree / "processed_datasets" / "xtal-0003"
    (ddir / "xtal-0003_event_1_best_autobuild.pdb").unlink()
    ds = read_dataset(tree, "xtal-0003")
    assert ds.events[0].pose == ddir / "autobuild" / "7_1_dict_0.pdb"


def test_ligand_id_is_the_component_pandda_was_given(tmp_path):
    tree = make_tree(tmp_path / "out", {"xtal-0005": [event_record(1)]}, ligand_code="5KX")
    ds = read_dataset(tree, "xtal-0005")
    assert ds.ligand_id == "5KX"
    # the comp_LIG alias staging appends is not the code
    ddir = tree / "processed_datasets" / "xtal-0005"
    cif = ddir / "ligand_files" / "dict.cif"
    cif.write_text(cif.read_text() + "data_comp_LIG\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\nLIG C1 C\n")
    assert read_ligand_id(ddir) == "5KX"
    none = make_tree(tmp_path / "none", {"xtal-0006": []}, ligand_code=None)
    assert read_dataset(none, "xtal-0006").ligand_id is None
