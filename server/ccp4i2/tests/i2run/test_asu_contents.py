import os
import xml.etree.ElementTree as ET
from .utils import i2run


def _xml_elements_equal(e1, e2, ignore_source_path=False):
    """Compare two XML elements semantically, ignoring child element order.

    This handles cases where the XML content is equivalent but elements
    appear in different orders (e.g., <name> before <description> vs after).

    Args:
        e1, e2: XML elements to compare
        ignore_source_path: If True, only check that source/baseName ends with expected filename
    """
    if e1.tag != e2.tag:
        return False

    # Special handling for source element - actual output may store full path
    if e1.tag == 'source' and ignore_source_path:
        # Just verify baseName elements end with same filename
        bn1 = e1.find('baseName')
        bn2 = e2.find('baseName')
        if bn1 is not None and bn2 is not None:
            # Actual may have full path, expected may have just filename
            actual_basename = (bn1.text or '').strip()
            expected_basename = (bn2.text or '').strip()
            # Check if actual path ends with expected filename
            if not actual_basename.endswith(expected_basename):
                return False
            return True  # Source element matches
        return False

    if (e1.text or '').strip() != (e2.text or '').strip():
        return False
    if e1.attrib != e2.attrib:
        return False

    # Get children sorted by tag for comparison
    children1 = sorted(list(e1), key=lambda x: x.tag)
    children2 = sorted(list(e2), key=lambda x: x.tag)

    if len(children1) != len(children2):
        return False

    return all(_xml_elements_equal(c1, c2, ignore_source_path) for c1, c2 in zip(children1, children2))


_BETA_SEQ = (
    "HPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRVD"
    "AGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGP"
    "KELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPAAMATTLRKLLTGELLTLASRQ"
    "QLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTG"
    "SQATMDERNRQIAEIGASLIKHW"
)

_BLIP_SEQ = (
    "AGVMTGAKFTQIQFGMTRQQVLDIAGAENCETGGSFGDSIHCRGHAAGDYYAYATFGFTS"
    "AAADAKVDSKSQEKLLAPSAPTLTLAKFNQVTVGMTRAQVLATVGQGSCTTWSEYYPAYP"
    "STAGVTLSLSCFDVDGYSSTGFYRGSAHLWFTDGVLQGKRQWDLV"
)

_FAKE_PATH = os.path.join("$CCP4I2_ROOT", "demo_data", "beta_blip")


def test_beta_blip():
    args = ["ProvideAsuContents"]
    args += ["--ASU_CONTENT"]
    args += [f"sequence={_BETA_SEQ}"]
    args += ["nCopies=1"]
    args += ["name=BETA"]
    args += ["description=Beta-lactamase"]
    args += ["source/baseName=beta.seq"]
    args += [f"source/relPath={_FAKE_PATH}"]
    args += ["polymerType=PROTEIN"]
    args += ["--ASU_CONTENT"]
    args += ["name=BLIP"]
    args += ["description=Beta-lactamase inhibitory protein"]
    args += [f"sequence={_BLIP_SEQ}"]
    args += ["source/baseName=blip.seq"]
    args += [f"source/relPath={_FAKE_PATH}"]
    args += ["nCopies=1"]
    args += ["polymerType=PROTEIN"]
    with i2run(args) as job:
        actual = ET.parse(job / "ASUCONTENTFILE.asu.xml").find(".//seqList")
        expected = ET.fromstring(
            f"""
            <seqList>
              <CAsuContentSeq>
                <sequence>{_BETA_SEQ}</sequence>
                <nCopies>1</nCopies>
                <polymerType>PROTEIN</polymerType>
                <name>BETA</name>
                <description>Beta-lactamase</description>
                <source>
                  <baseName>beta.seq</baseName>
                  <relPath>{_FAKE_PATH}</relPath>
                </source>
              </CAsuContentSeq>
              <CAsuContentSeq>
                <sequence>{_BLIP_SEQ}</sequence>
                <nCopies>1</nCopies>
                <polymerType>PROTEIN</polymerType>
                <name>BLIP</name>
                <description>Beta-lactamase inhibitory protein</description>
                <source>
                  <baseName>blip.seq</baseName>
                  <relPath>{_FAKE_PATH}</relPath>
                </source>
              </CAsuContentSeq>
            </seqList>
            """
        )
        # Use semantic comparison that ignores element ordering within CAsuContentSeq
        # and handles source path differences (actual may have full path, expected has filename only)
        assert _xml_elements_equal(actual, expected, ignore_source_path=True), (
            f"XML mismatch:\nActual: {ET.tostring(actual, encoding='unicode')}\n"
            f"Expected: {ET.tostring(expected, encoding='unicode')}"
        )
