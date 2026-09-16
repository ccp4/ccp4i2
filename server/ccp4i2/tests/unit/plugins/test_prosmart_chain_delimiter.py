"""prosmart's chain list is comma- *and* space-tolerant.

The React chain selector (CChainSelectElement) stores CHAINLIST_1 comma-
delimited ("A,B"); the wrapper used to split on whitespace only, so "A,B" reached
prosmart as one chain literally named "A,B" -> "could not read chain A,B", which
sank the whole servalcat_pipe refinement (CryoMapMR job 13). The wrapper now
splits on either delimiter; these tests pin that against the built command line.
"""

from ccp4i2.core.tasks import get_plugin_class


def _chains_after_c1(commandLine):
    """The chain tokens prosmart is given: the args after -c1 up to the next flag."""
    assert "-c1" in commandLine, commandLine
    chains = []
    for word in commandLine[commandLine.index("-c1") + 1:]:
        if word.startswith("-"):
            break
        chains.append(word)
    return chains


def _build(tmp_path, chainlist):
    p = get_plugin_class("prosmart")()
    p.workDirectory = tmp_path  # a pathlib.Path; makeFileName does workDir / name
    p.tempFile = str(tmp_path / "RESTRAINTS_TARGET.pdb")
    ref = tmp_path / "ref.pdb"
    ref.write_text("REMARK reference\n")
    p.container.inputData.REFERENCE_MODELS.append(str(ref))
    p.container.inputData.CHAINLIST_1.set(chainlist)
    p.makeCommandAndScript()
    return [str(w) for w in p.commandLine]


def test_comma_delimited_chains_are_split(tmp_path):
    # The React-frontend format that used to break.
    assert _chains_after_c1(_build(tmp_path, "A,B")) == ["A", "B"]


def test_space_delimited_still_works(tmp_path):
    # The legacy Qt format must keep working.
    assert _chains_after_c1(_build(tmp_path, "A B")) == ["A", "B"]


def test_mixed_and_padded(tmp_path):
    assert _chains_after_c1(_build(tmp_path, " A, B ,C ")) == ["A", "B", "C"]


def test_single_chain(tmp_path):
    assert _chains_after_c1(_build(tmp_path, "A")) == ["A"]
