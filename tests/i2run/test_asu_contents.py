import os
import xml.etree.ElementTree as ET
from .utils import i2run


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
        canonicalizedActual = ET.canonicalize(ET.tostring(actual), strip_text=True)
        canonicalizedExpected = ET.canonicalize(ET.tostring(expected), strip_text=True)
        assert canonicalizedActual == canonicalizedExpected
