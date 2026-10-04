"""Phaser's expected LLG, read into program.xml as numbers.

It was there only as the text of Phaser's summary: what Phaser expects of
each search model, and how much of the asymmetric unit the models account
for (BetaBlip with beta alone: 0.61; a judge of the job could not see that
a component was missing).
"""
import pytest

pytest.importorskip("lxml")

from lxml import etree  # noqa: E402

from ccp4i2.wrappers.phaser_phil.script import phaser_run  # noqa: E402

# As Phaser 2.8.4 wrote it for BetaBlip (docs scenario, job 4), abridged.
BLOCK = """
   Resolution of Data (Selected):    3.004 (3.00366)
   eLLG Target: 225

-----------------
POLY-ALANINE ELLG
-----------------
   eLLG-target |  0.10  0.20  0.40  0.80  1.60  3.20
           225 |   158   163   181   266 -full--full-

--------------
MONOMERIC ELLG
--------------

   Expected LLG (eLLG)
   -------------------
   eLLG: eLLG of ensemble alone
       eLLG   RMSD frac-scat  Ensemble
      949.5  0.502   0.60978  beta
      377.7  0.461   0.37205  blip

   Resolution for eLLG target
   --------------------------
   eLLG-reso: Resolution to achieve target eLLG (225)
     eLLG-reso  Ensemble
         4.963  beta
         3.693  blip

   Resolution for eLLG target: data collection
   -------------------------------------------
   eLLG-reso: Resolution to achieve target eLLG (225) with perfect data
     eLLG-reso  Ensemble
         4.972  beta
         3.723  blip

   eLLG indicates that placement of a single copy of ensemble "beta" should be easy
   eLLG indicates that placement of a single copy of ensemble "blip" should be easy

   Expected LLG (eLLG): Chains
   ---------------------------
   eLLG: eLLG of chain alone
       eLLG   RMSD frac-scat chain  Ensemble
      949.5  0.502   0.60978  "  "  beta
"""


def test_each_ensemble_is_read():
    found = phaser_run.expected_llg(BLOCK)
    assert found["target"] == 225.0
    assert found["ensembles"] == [
        {"name": "beta", "ellg": 949.5, "rmsd": 0.502, "fraction_scattering": 0.60978,
         "resolution_for_target": 4.963, "call": "easy"},
        {"name": "blip", "ellg": 377.7, "rmsd": 0.461, "fraction_scattering": 0.37205,
         "resolution_for_target": 3.693, "call": "easy"},
    ]


def test_as_xml():
    root = etree.Element("PhaserMrResults")
    phaser_run.expected_llg_xml(BLOCK, root)
    assert root.findtext("ExpectedLLG/Target") == "225.0"
    assert root.findtext("ExpectedLLG/FractionScatteringOfEnsembles") == "0.98183"
    assert [e.findtext("eLLG") for e in root.findall("ExpectedLLG/Ensemble")] == ["949.5", "377.7"]
    assert root.findtext("ExpectedLLG/Ensemble[Name='blip']/ResolutionForTarget") == "3.693"


def test_nothing_without_the_table():
    root = etree.Element("PhaserMrResults")
    assert phaser_run.expected_llg_xml("eLLG Target: 225\n", root) is None
    assert len(root) == 0
