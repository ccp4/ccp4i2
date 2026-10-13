"""A protein entry too short to be the crystallised construct is questioned.

An agent typed "MNNI" for both BetaBlip components (gpt-oss, 2026-10-13);
validation passed it and MR ran on a peptide-sized composition.
"""
from ccp4i2.agent import judgement
from ccp4i2.core.tasks import get_plugin_class


def _codes(plugin):
    return [r["code"] for r in plugin.validity()._reports]


def test_a_four_residue_protein_is_warned_of(tmp_path):
    plugin = get_plugin_class("ProvideAsuContents")(workDirectory=str(tmp_path), name="a")
    items = plugin.container.inputData.ASU_CONTENT
    item = items.makeItem()
    items.append(item)
    item.sequence.set("MNNI")
    item.name.set("beta")
    item.polymerType.set("PROTEIN")
    assert 210 in _codes(plugin)
    item.sequence.set("HPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLL")
    assert 210 not in _codes(plugin)


def test_a_peptide_sized_au_is_told_the_sequences_are_wrong(tmp_path):
    (tmp_path / "program.xml").write_text(
        "<ASUCONTENTMATTHEWS><totalWeight>978.0</totalWeight><matthewsCompositions>"
        "<composition><nMolecules>1</nMolecules><solventPercentage>99.8</solventPercentage>"
        "<matthewsCoeff>500</matthewsCoeff><matthewsProbability>0.0</matthewsProbability>"
        "</composition></matthewsCompositions><cellVolume>651311</cellVolume><polymerMode>P</polymerMode>"
        "</ASUCONTENTMATTHEWS>")
    verdict = judgement.judge("ProvideAsuContents", tmp_path)
    advice = " ".join(step.get("advice", "") for step in verdict.get("next") or [])
    assert verdict.get("outcome") == "too_sparse"
    assert "The sequences are wrong, not the copies" in advice


def test_the_same_chain_entered_twice_is_warned_of(tmp_path):
    # Haiku added gamma by the page's SEQIN and again as an entry's source
    # (2026-10-13): overpacked, and Crank2's Matthews then found no fit
    plugin = get_plugin_class("ProvideAsuContents")(workDirectory=str(tmp_path), name="a")
    items = plugin.container.inputData.ASU_CONTENT
    gamma = "MKVLAAGIVGLLAQPAGAEDHLAAGRLPLMSKAWEAVLKAG" * 3
    for name in ("Chain_gamma", "gamma"):
        item = items.makeItem()
        items.append(item)
        item.sequence.set(gamma)
        item.name.set(name)
        item.polymerType.set("PROTEIN")
    assert 211 in _codes(plugin)
