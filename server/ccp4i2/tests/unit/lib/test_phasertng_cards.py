"""PhaserTNG Picard's results, read from what it writes, into program.xml.

Fixtures are real: Gamma (one model, one copy) and BetaBlip (beta-lactamase
and BLIP, two components) best.1.dag.cards from Picard runs on
2026-10-08 (ccp4-20260702), and BetaBlip's 'Best Solution:' log line.
"""
from pathlib import Path
from xml.etree import ElementTree as ET

from ccp4i2.lib.utils.formats.phasertng_cards import (
    annotation_components, parse_best_solution, parse_dag_cards, program_xml, same_spacegroup,
)

FIX = Path(__file__).resolve().parent / "fixtures" / "phasertng"


def test_gamma_one_component():
    sols = parse_dag_cards((FIX / "gamma_best.1.dag.cards").read_text())
    assert len(sols) == 1
    best = sols[0]
    assert float(best["zscore"]) == 47.733035
    assert float(best["rfactor"]) == 23.871342 and float(best["llg"]) == 970.011812
    assert best["hermann_mauguin"] == "P 21 21 21"
    assert best["seek_components"] == "1" and best["seek_complete"] == "true"
    assert best["poses"] == [{"id": "2952750000", "tag": "gamma_model", "tfz": 47.733035}]
    assert best["components"] == [{"tag": "gamma_model", "rfz": 32.0, "tfz": 46.0, "packs": "Y"}]


def test_betablip_two_components_in_the_order_placed():
    sols = parse_dag_cards((FIX / "betablip_best.1.dag.cards").read_text())
    best = sols[0]
    assert best["seek_components"] == "2" and best["seek_complete"] == "true"
    assert [(p["tag"], round(p["tfz"], 2)) for p in best["poses"]] == [("beta", 22.08), ("blip", 31.29)]
    assert [(c["tag"], c["rfz"], c["tfz"]) for c in best["components"]] == \
        [("beta", 12.0, 22.0), ("blip", 5.0, 25.0)]
    assert float(best["zscore"]) == 31.292507 and float(best["rfactor"]) == 38.78733


def test_the_best_solution_line_carries_rfree():
    text = (FIX / "betablip_best_solution.txt").read_text()
    assert parse_best_solution(text) == {"rwork": 38.79, "rfree": 46.79, "dmin": 3.0,
                                         "poses": 2, "spacegroup": "P3221"}
    assert parse_best_solution("nothing here") is None


def test_annotation_without_a_find_gives_no_components():
    assert annotation_components("") == []
    assert annotation_components("RF=100%/5z") == []


def test_program_xml_holds_what_a_verdict_reads():
    record = {
        "solutions": parse_dag_cards((FIX / "betablip_best.1.dag.cards").read_text()),
        "best": parse_best_solution((FIX / "betablip_best_solution.txt").read_text()),
        "run": {"resolution": 3.0036, "anomalous": False, "wilson_b": 51.2,
                "matthews_z": 1.0, "matthews_vm": 2.4, "matthews_probability": 0.9,
                "twinning_indicated": False, "tncs_indicated": False, "wall_seconds": 95.3},
    }
    root = program_xml(record)
    assert root.find("Data").get("resolution") == "3.00"
    assert root.find("Composition").get("z") == "1.0"
    sols = root.find("Solutions")
    assert sols.get("count") == "1" and sols.get("components_sought") == "2" and sols.get("complete") == "true"
    best = sols.find("Solution")
    assert best.get("tfz") == "31.29" and best.get("rfactor") == "38.79" and best.get("rfree") == "46.79"
    assert best.get("spacegroup") == "P 32 2 1" and best.get("packs") == "Y"
    assert best.get("scattering_fraction") == "0.784"
    tags = [(c.get("tag"), c.get("tfz"), c.get("search_tfz")) for c in best.find("Components")]
    assert tags == [("beta", "22.08", "22.0"), ("blip", "31.29", "25.0")]
    # the xpaths the judgement file uses
    assert root.find(".//Solutions/Solution[1]").get("tfz") == "31.29"
    assert root.find(".//Solutions/Solution[1]/Components/Component[last()]").get("tfz") == "31.29"
    assert root.find(".//Solutions/Solution[1]/Components/Component[1]").get("tfz") == "22.08"
    ET.tostring(root)  # serialises


def test_two_copies_of_one_model_keep_their_own_search_records():
    cards = (FIX / "gamma_best.1.dag.cards").read_text()
    # a second copy of gamma_model, placed second at a lower search TFZ
    two = cards.replace('phaserdag node annotation " FIND=gamma_model RF=100%/32z',
                        'phaserdag node annotation " FIND=gamma_model RF=100%/32z RFR=616/100%/32z GYRE=272 '
                        'TF=100%/46z/Y FIND=gamma_model RF=100%/9z RFR=1/100%/9z GYRE=1 TF=100%/11z/Y '
                        'POSE=1 PAK=0%/Y RBR=1 XX=1 YY=2 ZZ=3 Z0=0 Z1=1 Z2=2 ZQ=')
    pose = next(l for l in cards.splitlines() if l.startswith("phaserdag node pose"))
    two = two.replace(pose, pose + "\n" + pose.replace("tfz 47.733035", "tfz 11.5"))
    root = program_xml({"solutions": parse_dag_cards(two), "best": None, "run": {}})
    tags = [(c.get("tag"), c.get("tfz"), c.get("search_tfz")) for c in root.find(".//Components")]
    assert tags == [("gamma_model", "47.73", "46.0"), ("gamma_model", "11.50", "11.0")]


def test_the_data_space_group_and_whether_the_solution_changed_it():
    assert same_spacegroup("P 2 21 21", "P22121") is True
    assert same_spacegroup("P 21 21 21", "P22121") is False
    assert same_spacegroup("P 1", None) is None
    sols = parse_dag_cards((FIX / "gamma_best.1.dag.cards").read_text())
    root = program_xml({"solutions": sols, "best": None, "run": {}, "data_spacegroup": "P22121"})
    assert root.find("Data").get("spacegroup") == "P22121"
    assert root.find(".//Solution").get("spacegroup_changed") == "true"
    root = program_xml({"solutions": sols, "best": None, "run": {}, "data_spacegroup": "P212121"})
    assert root.find(".//Solution").get("spacegroup_changed") == "false"
