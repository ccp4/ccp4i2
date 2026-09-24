"""Test that the reflection digest carries the content-based `diagnosis` block.

The digest enrichment adds `diagnose_reflection_file()` output under a new
`diagnosis` key *additively* — the legacy `format` (extension-based) and
`merged` (getMerged stub) keys must stay put so existing consumers are
unaffected. CCP4-free (gemmi only), so it runs in the slim CI.
"""
import pytest

from ccp4i2 import I2_TOP
from ccp4i2.core import CCP4XtalData

_GAMMA_MTZ = I2_TOP / "demo_data" / "gamma" / "freeR.mtz"


@pytest.fixture
def gamma_mtz():
    if not _GAMMA_MTZ.exists():
        pytest.skip(f"Demo data not found: {_GAMMA_MTZ}")
    return str(_GAMMA_MTZ)


def test_reflection_digest_has_diagnosis_block(gamma_mtz):
    from ccp4i2.lib.utils.files.digest import (
        digest_cgenericrefldatafile_file_object,
    )

    f = CCP4XtalData.CGenericReflDataFile()
    f.setFullPath(gamma_mtz)
    d = digest_cgenericrefldatafile_file_object(f)

    # Backward-compatible: the legacy keys are still present.
    assert "format" in d
    assert "merged" in d

    # Additive: the new content-based diagnosis block.
    diag = d.get("diagnosis")
    assert diag is not None, "digest should carry a diagnosis block"
    assert diag["format"] == "mtz"
    assert diag["merged"] is True          # real detection, not the getMerged stub
    assert diag["staraniso"] is False      # plain gamma MTZ
    assert diag["spaceGroupNumber"] == 19  # P 21 21 21
    assert diag["cell"] and len(diag["cell"]) == 6
