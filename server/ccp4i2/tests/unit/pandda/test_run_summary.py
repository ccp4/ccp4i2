"""The campaign report's analysis: PanDDA's own summary recapitulated from
the tables it writes (events, sites, per-dataset YAML), the way Reinspect
recovered PanDDA 1's HTML. CCP4-free: gemmi for tiny maps, yaml for the
per-dataset files, lxml for the report render."""
import textwrap
from xml.etree import ElementTree as ET

import pytest

pytest.importorskip("gemmi", reason="needs gemmi")
pytest.importorskip("yaml", reason="needs PyYAML")

from ccp4i2.wrappers.pandda_campaign.script import pandda_run_summary as summary
from .synthetic_tree import event_record, make_tree

# Verbatim shape of the bundled PanDDA 2's processed_dataset.yaml: numpy
# scalars serialised as python/object tags with a binary payload, which
# yaml.safe_load refuses outright.
NUMPY_YAML = textwrap.dedent("""\
    Summary:
      Processing Resolution: 2.07
      Comparator Datasets:
      - xtal-0000
      - xtal-0001
      - xtal-0002
      Selected Model: 2
      Selected Model Events:
      - 1
      - 2
    Models:
      1:
        Processed?: true
        Characterization Datasets:
        - xtal-0001
        Model Score: 0.0
        Number of Initial Events: 3
        Number of Size Filtered Events: 1
        Number of Score Filtered Events: 0
        Events: {}
      2:
        Processed?: true
        Characterization Datasets:
        - xtal-0000
        Model Score: 0.0
        Number of Initial Events: 4
        Number of Size Filtered Events: 3
        Number of Score Filtered Events: 2
        Events:
          1:
            Score: !!python/object/apply:numpy.core.multiarray.scalar
            - &id001 !!python/object/apply:numpy.dtype
              args:
              - f4
              - false
              - true
              state: !!python/tuple
              - 3
              - <
              - null
              - null
              - null
              - -1
              - -1
              - 0
            - !!binary |
              WPF8Pw==
            BDC: 0.76
""")


@pytest.fixture
def tree(tmp_path):
    root = make_tree(tmp_path / "pandda2_out", {
        "xtal-0002": [event_record(1, bdc=0.53, score=0.99), event_record(2, bdc=0.91, score=0.62, build=False)],
        "xtal-0000": [],
        "xtal-0001": [],
    })
    (root / "processed_datasets" / "xtal-0002" / "processed_dataset.yaml").write_text(NUMPY_YAML)
    (root / "analyses" / "pandda_analyse_sites.csv").write_text(
        ',site_idx,centroid\n0,1,"(-25.8, -2.5, -1.2)"\n1,2,"(10.0, 20.0, 30.0)"\n')
    return root


def test_numpy_tagged_yaml_loads_and_the_wrapped_values_are_none(tmp_path):
    path = tmp_path / "processed_dataset.yaml"
    path.write_text(NUMPY_YAML)
    data = summary.load_pandda_yaml(path)
    assert data["Summary"]["Selected Model"] == 2
    event = data["Models"][2]["Events"][1]
    assert event["Score"] is None, "a numpy scalar is not reconstructed"
    assert event["BDC"] == pytest.approx(0.76), "plain floats beside it survive"


def test_dataset_summary_reads_the_selected_model_filter_counts(tree):
    record = summary.read_dataset_summary(tree / "processed_datasets" / "xtal-0002")
    assert record["resolution"] == pytest.approx(2.07)
    assert record["n_comparators"] == 3
    assert (record["n_models"], record["selected_model"]) == (2, 2)
    assert (record["n_initial_events"], record["n_size_filtered_events"],
            record["n_score_filtered_events"], record["n_events"]) == (4, 3, 2, 2)
    assert record["analysed"]


def test_dataset_summary_tolerates_a_minimal_or_missing_yaml(tree):
    minimal = summary.read_dataset_summary(tree / "processed_datasets" / "xtal-0000")
    assert minimal["resolution"] == pytest.approx(1.8) and minimal["n_comparators"] is None
    (tree / "processed_datasets" / "xtal-0001" / "processed_dataset.yaml").unlink()
    missing = summary.read_dataset_summary(tree / "processed_datasets" / "xtal-0001")
    assert missing["resolution"] is None and missing["analysed"]


def test_events_join_the_table_with_each_datasets_build_record(tree):
    events = summary.read_events(tree)
    assert [(e["dtag"], e["event_idx"]) for e in events] == [("xtal-0002", 1), ("xtal-0002", 2)]
    one, two = events
    assert one["event_fraction"] == pytest.approx(0.47)
    assert one["site_idx"] == 1 and two["site_idx"] == 2
    assert one["score"] == pytest.approx(0.99), "from events.yaml, not the table"
    assert one["build_score"] == pytest.approx(0.57) and one["rscc"] == pytest.approx(0.30)
    assert two["build_score"] is None, "no Build block, no build columns"
    assert one["interesting"] is False


def test_sites_parse_the_stringified_centroid(tree):
    sites = summary.read_sites(tree)
    assert [s["site_idx"] for s in sites] == [1, 2]
    assert sites[0]["centroid"] == pytest.approx((-25.8, -2.5, -1.2))


def test_no_tables_means_no_events_and_no_sites(tmp_path):
    root = make_tree(tmp_path / "pandda2_out", {"xtal-0000": []}, events_table=False)
    assert summary.read_events(root) == [] and summary.read_sites(root) == []
    run = summary.summarise_run(root)
    assert run["stats"]["n_events"] == 0 and run["stats"]["n_sites"] == 0
    assert run["histograms"]["event_fraction"] == []


def test_bin_values_is_equal_width_and_keeps_the_ends():
    bins = summary.bin_values([0.0, 0.5, 1.0, 1.0], nbins=4)
    assert [b["count"] for b in bins] == [1, 0, 1, 2], "the maximum lands in the last bin, not off the end"
    assert bins[0]["centre"] == pytest.approx(0.125)
    assert summary.bin_values([0.3, 0.3]) == [{"centre": 0.3, "count": 2}]
    assert summary.bin_values([None, None]) == []


def test_run_summary_aggregates_sites_and_datasets(tree):
    run = summary.summarise_run(tree)
    stats = run["stats"]
    assert (stats["n_datasets"], stats["n_analysed"], stats["n_events"], stats["n_sites"]) == (3, 3, 2, 2)
    assert stats["n_datasets_with_events"] == 1
    assert stats["best_score"] == pytest.approx(0.99)
    assert stats["median_resolution"] == pytest.approx(1.8)
    by_site = {s["site_idx"]: s for s in run["sites"]}
    assert by_site[1]["n_events"] == 1 and by_site[1]["n_datasets"] == 1
    assert by_site[1]["centroid"] == pytest.approx((-25.8, -2.5, -1.2))
    assert by_site[2]["best_score"] == pytest.approx(0.62)
    by_dtag = {d["dtag"]: d for d in run["datasets"]}
    assert by_dtag["xtal-0002"]["best_score"] == pytest.approx(0.99)
    assert by_dtag["xtal-0000"]["n_events"] == 0 or by_dtag["xtal-0000"]["n_events"] is None
    assert sum(b["count"] for b in run["histograms"]["event_fraction"]) == 2
    assert sum(b["count"] for b in run["histograms"]["resolution"]) == 3


def test_a_site_only_in_the_events_table_still_counts(tree):
    (tree / "analyses" / "pandda_analyse_sites.csv").unlink()
    run = summary.summarise_run(tree)
    assert [s["site_idx"] for s in run["sites"]] == [1, 2]
    assert all(s["centroid"] is None for s in run["sites"])


def test_xml_carries_every_column_and_blanks_the_absent_ones(tree):
    root = ET.Element("pandda_campaign")
    summary.analysis_to_xml(root, summary.summarise_run(tree))
    node = root.find("analysis")
    assert node.findtext("stats/n_events") == "2"
    events = node.findall("events/event")
    assert len(events) == 2
    assert set(summary.EVENT_ATTRS) <= set(events[0].keys())
    assert events[1].get("build_score") == "", "absent values are empty strings, so the columns line up"
    assert events[0].get("interesting") == "false"
    assert node.find("sites/site").get("centroid") == "-25.80 -2.50 -1.20"
    assert len(node.findall("histograms/histogram")) == 4
    assert node.find("histograms/histogram[@name='resolution']/bin").get("count")


def test_the_report_draws_the_analysis(tree, tmp_path):
    etree = pytest.importorskip("lxml.etree", reason="the report layer needs lxml")
    from ccp4i2.wrappers.pandda_campaign.script.pandda_campaign_report import pandda_campaign_report

    root = ET.Element("pandda_campaign")
    ET.SubElement(root, "state").text = "finished"
    ET.SubElement(root, "n_datasets").text = "3"
    datasets = ET.SubElement(root, "datasets")
    for dtag, label in (("xtal-0000", "Campaign_a"), ("xtal-0001", "Campaign_b"), ("xtal-0002", "Campaign_c")):
        ET.SubElement(datasets, "dataset", xtal=dtag, label=label, dict="yes")
    summary.analysis_to_xml(root, summary.summarise_run(tree))
    report = pandda_campaign_report(xmlnode=etree.fromstring(ET.tostring(root)),
                                    jobInfo={"fileroot": str(tmp_path)}, jobStatus="Finished")
    rendered = ET.tostring(report.as_data_etree(), encoding="unicode")
    assert "Events per site" in rendered
    for title in ("Event fraction (1-BDC)", "Hit-in-site probability", "Processing resolution (A)"):
        assert title in rendered
    assert rendered.count("<barchart") >= 5, "sites (2) and the binned histograms"
    assert "Campaign_c" in rendered, "events are labelled by the dataset's project, not only its xtal"


def test_the_report_without_analysis_still_renders():
    etree = pytest.importorskip("lxml.etree", reason="the report layer needs lxml")
    from ccp4i2.wrappers.pandda_campaign.script.pandda_campaign_report import pandda_campaign_report

    xml = "<pandda_campaign><state>running</state><n_datasets>4</n_datasets></pandda_campaign>"
    report = pandda_campaign_report(xmlnode=etree.fromstring(xml), jobInfo={}, jobStatus="Running")
    assert "PanDDA is running" in ET.tostring(report.as_data_etree(), encoding="unicode")
