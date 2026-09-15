"""Hand-selection + confidence logic for molrep_map.

The task recommends a hand from the real-space map-model CC, but must be honest
about how trustworthy that call is -- a weak fit (both hands low CC, e.g. the
CryoMapMR job_1 case at ~0.05/0.05) is reported as not confidently determinable,
not dressed up as a confident pick. Pure logic; no gemmi/CCP4 needed.
"""

from types import SimpleNamespace

from ccp4i2.wrappers.molrep_map.script.molrep_map import molrep_map


def _plugin(cc, placed=("Original", "Flipped"), scores=None):
    p = object.__new__(molrep_map)
    p._cc = dict(cc)
    scores = scores or {}
    p._results = {
        h: SimpleNamespace(placed=(h in placed), score=scores.get(h))
        for h in ("Original", "Flipped")
    }
    return p


def test_confident_when_a_clear_high_cc_winner():
    p = _plugin({"Original": 0.82, "Flipped": 0.10})
    assert p._choose_hand() == "Original"
    assert p._hand_confidence() == "confident"


def test_weak_when_both_hands_fit_poorly():
    # The job_1 situation: both ~0.05 -> not a solvable case.
    p = _plugin({"Original": 0.05, "Flipped": 0.05})
    assert p._hand_confidence() == "weak"


def test_ambiguous_when_good_but_too_close():
    p = _plugin({"Original": 0.42, "Flipped": 0.40})
    assert p._hand_confidence() == "ambiguous"


def test_single_when_only_one_hand_placed():
    p = _plugin({"Original": 0.7}, placed=("Original",))
    assert p._choose_hand() == "Original"
    assert p._hand_confidence() == "single"


def test_falls_back_to_molrep_score_without_cc():
    p = _plugin({"Original": float("nan"), "Flipped": float("nan")},
                scores={"Original": 0.05, "Flipped": 0.09})
    assert p._choose_hand() == "Flipped"   # higher molrep score wins the fallback
