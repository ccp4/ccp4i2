"""The Phaser run record: events recorded, prose read by its fixed sentences.

Both summary fixtures are real r.summary() texts from Phaser 2.8.4 on the
gamma demo data: one copy of the model (a single solution) and two copies
asked of a one-copy crystal (a second component that fails at the search
resolution, is retried at full resolution, and is given up).
"""
from pathlib import Path

from lxml import etree

from ccp4i2.wrappers.phaser_phil.script.phaser_run import (
    PhaserRecorder, annotation_tokens, strategy_attempts, summary_blocks)

FIXTURES = Path(__file__).parent / "fixtures" / "phaser"


def blocks(name):
    return summary_blocks((FIXTURES / name).read_text())


class TestSummaryBlocks:

    def test_one_block_per_module_heading(self):
        one = blocks("summary_gamma_1copy.txt")
        two = blocks("summary_gamma_2copy.txt")
        assert len(one) == 17 and len(two) == 32
        assert one[0][0] == "AUTOMATED MOLECULAR REPLACEMENT"
        assert "Input Search Order" in one[0][1]
        assert "Refinement" in one[-2][1] or "REFINEMENT" in one[-2][0]


class TestStrategyAttempts:

    def test_single_solution_is_one_placed_attempt(self):
        attempts, unparsed = strategy_attempts(blocks("summary_gamma_1copy.txt"))
        assert unparsed == 0
        assert len(attempts) == 1
        # The gamma model is already at the origin: RF*0 TF*0 in the annotation
        assert attempts[0]["outcome"].startswith("placed")
        assert attempts[0]["component"] == "model"

    def test_failed_second_component_is_narrated(self):
        attempts, unparsed = strategy_attempts(blocks("summary_gamma_2copy.txt"))
        assert unparsed == 0
        outcomes = [(a["resolution"], a["outcome"]) for a in attempts]
        assert outcomes[0] == (2.99, "placed")
        assert outcomes[1] == (2.99, "no solution at lower resolution")
        assert outcomes[2] == ("full", "no definite solution")
        assert attempts[0]["llg"] > 3000

    def test_unknown_wording_is_counted_not_guessed(self):
        made_up = [("AUTOMATED MOLECULAR REPLACEMENT", "** Something new and different")]
        attempts, unparsed = strategy_attempts(made_up)
        assert attempts == [] and unparsed == 1


class TestAnnotationTokens:

    def test_documented_tokens(self):
        comps, unknown = annotation_tokens(
            "RFZ=12.3 TFZ=15.1 PAK=0 LLG=350 TFZ==18.0 RFZ=8.0 TFZ=22.5 PAK=2 LLG=900 LLG=1200 TFZ==25.1")
        assert unknown == []
        assert len(comps) == 2
        assert comps[0]["RFZ"] == 12.3 and comps[0]["TFZ"] == 15.1 and comps[0]["TFZeq"] == 18.0
        assert comps[1]["PAK"] == 2 and comps[1]["LLG"] == 1200   # low- then high-res refinement

    def test_flags_and_amalgamation(self):
        comps, unknown = annotation_tokens("RF*0 TF*0 PAK=0 LLG=3047 TFZ==46.3 RF++ TFZ=* ( TFZ=9.1 PAK=1 LLG=400 ) +TNCS *T=2")
        assert unknown == []
        assert "rotation 000" in comps[0]["notes"][0]
        assert comps[0]["LLG"] == 3047
        assert "first molecule in P1" in " ".join(comps[1]["notes"])
        assert comps[1]["amalgamated"] == ["TFZ=9.1", "PAK=1", "LLG=400"]
        assert "tNCS" in " ".join(comps[1]["notes"]) and "template solution 2" in " ".join(comps[1]["notes"])

    def test_the_gamma_annotation(self):
        comps, unknown = annotation_tokens(" RF*0 TF*0 LLG=4650 TFZ==41.2 PAK=0 LLG=4650 TFZ==41.2")
        assert unknown == [] and len(comps) == 1 and comps[0]["TFZeq"] == 41.2

    def test_unknown_tokens_are_reported(self):
        comps, unknown = annotation_tokens("RFZ=5.0 WOBBLE=3 LLG=10")
        assert unknown == ["WOBBLE=3"] and comps[0]["LLG"] == 10


class TestRecorder:

    def make(self):
        root = etree.Element("PhaserMrResults")
        flushes = []
        clock = [0.0]
        rec = PhaserRecorder(root, flush=lambda r: flushes.append(len(r)), min_interval=1.0,
                             clock=lambda: clock[0])
        return root, rec, flushes, clock

    def test_modules_progress_and_verdict_are_recorded(self):
        root, rec, flushes, clock = self.make()
        rec.call_back("phaser_start", "")
        rec.call_back("MODE", "AUTOMATED MOLECULAR REPLACEMENT")
        rec.startProgressBar("Refining solutions", 10)
        for _ in range(3):
            rec.incrementProgressBar()
        rec.warn("No Signal in Translation Function")
        rec.call_back("result", "Single Solution")
        assert [m.get("name") for m in root.findall("Modules/Module")] == ["AUTOMATED MOLECULAR REPLACEMENT"]
        activity = root.find("CurrentActivity")
        assert activity.findtext("label") == "Refining solutions"
        assert activity.findtext("max") == "10" and activity.findtext("value") == "3"
        assert root.findtext("PhaserWarnings/Warning") == "No Signal in Translation Function"
        assert root.findtext("Verdict") == "Single Solution"
        assert root.findtext("CurrentModule") == "AUTOMATED MOLECULAR REPLACEMENT"

    def test_increments_are_throttled_but_events_are_not(self):
        root, rec, flushes, clock = self.make()
        rec.startProgressBar("x", 100)          # forced
        n = len(flushes)
        for _ in range(50):
            rec.incrementProgressBar()          # within the interval: no flush
        assert len(flushes) == n
        clock[0] += 1.5
        rec.incrementProgressBar()
        assert len(flushes) == n + 1
        rec.call_back("MODE", "SOMETHING")      # forced
        assert len(flushes) == n + 2

    def test_attempts_open_on_rotation_and_close_on_placement(self):
        root, rec, flushes, clock = self.make()
        rec.call_back("MODE", "MOLECULAR REPLACEMENT ROTATION FUNCTION")
        rec.call_back("MODE", "MOLECULAR REPLACEMENT TRANSLATION FUNCTION")
        rec.call_back("current best solution", "SOLU SET RF*0 TF*0 LLG=3047")
        rec.call_back("MODE", "MOLECULAR REPLACEMENT ROTATION FUNCTION")
        rec.warn("No Signal in Translation Function")
        rec.call_back("MODE", "MOLECULAR REPLACEMENT ROTATION FUNCTION")
        rec.finish()
        attempts = root.findall("Attempts/Attempt")
        assert [(a.get("component"), a.get("number"), a.get("outcome")) for a in attempts] == [
            ("1", "1", "placed"), ("2", "1", "ended"), ("2", "2", "ended")]
        assert attempts[1].findtext("Warning") == "No Signal in Translation Function"
        assert root.findtext("PhaserCurrentBestSolution").startswith("SOLU SET")


class TestEpCycles:

    def test_the_gamma_sad_block(self):
        cycles, final = __import__("ccp4i2.wrappers.phaser_phil.script.phaser_run", fromlist=["ep_cycles"]).ep_cycles(
            blocks("summary_gamma_ep.txt"))
        assert [c["cycle"] for c in cycles] == list(range(1, len(cycles) + 1))
        assert len(cycles) >= 5
        assert cycles[0]["cycle"] == 1 and cycles[0]["added"] == 1
        assert cycles[0]["llg"] == 968 and cycles[0]["deleted"] == [3, 4]
        assert cycles[0]["converged"] is False
        assert cycles[1]["added"] == 2 and cycles[2]["deleted"] == [6]
        # Phaser prints a figure-of-merit table at the start and the end only
        assert final["llg"] == 1060.85 and final["fom"] == 0.401
