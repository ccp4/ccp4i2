"""The report classes ported off lxml-only APIs actually render.

The lint proves an API is not called; it cannot prove the replacement finds the
same nodes. `phaser_MR_PAK_report` is the one worth checking directly: it called
`self.xmlnode.xpath(...)` fifteen times on what the framework hands it, which is
a standard-library Element from `ET.parse(...).getroot()`. Every one raised
`AttributeError`, so the report could not render at all -- and the task is
registered with that `reportPath`, so this was reachable from the UI.

The XML below is synthetic, since no phaser fixture is in the tree. It carries
the element names the report queries, which is what the port had to preserve.
"""

import xml.etree.ElementTree as ET

import pytest

PHASER_XML = """<PhaserMrResults>
  <PhaserAdvisories><Advisory>an advisory</Advisory></PhaserAdvisories>
  <PhaserWarnings><Warning>a warning</Warning></PhaserWarnings>
  <PhaserCurrentBestSolution>
    <Solution>
      <spaceGroup>P 21 21 21</spaceGroup>
      <overallLLG>123.4</overallLLG>
      <Component><Name>ens1</Name><RFZ>4.1</RFZ><TFZ>9.2</TFZ></Component>
    </Solution>
  </PhaserCurrentBestSolution>
  <PhaserMrSolutions>
    <Solutions><Solution><SPG>P 21 21 21</SPG><TFZ>9.2</TFZ></Solution></Solutions>
  </PhaserMrSolutions>
  <Summary>Phaser Module: MOLECULAR REPLACEMENT PACKING ANALYSIS
 3 accepted of 4 </Summary>
</PhaserMrResults>
"""


@pytest.fixture
def rendered(tmp_path):
    from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR_PAK.script import (
        phaser_MR_PAK_report,
    )

    report = phaser_MR_PAK_report.phaser_MR_PAK_report(
        xmlnode=ET.fromstring(PHASER_XML),
        jobInfo={"fileroot": str(tmp_path)},
        jobStatus="Finished",
    )
    return ET.tostring(report.as_data_etree(), encoding="unicode")


class TestPhaserMrPakReport:
    def test_it_renders_at_all(self, rendered):
        """Before the port this raised AttributeError on the first query."""
        assert rendered

    @pytest.mark.parametrize(
        "content,came_from",
        [
            ("an advisory", "//PhaserAdvisories/Advisory -> .//"),
            ("a warning", "//PhaserWarnings/Warning -> .//"),
            ("P 21 21 21", "xpath0(...) -> find(), then ./spaceGroup"),
            ("123.4", "./overallLLG"),
            ("ens1", "Component/Name, a selector passed as a variable"),
        ],
    )
    def test_each_ported_query_still_finds_its_node(self, rendered, content, came_from):
        assert content in rendered, f"lost by the port of {came_from}"


class TestTheOtherPortedModules:
    """Import-only: these need rdkit, a qtpisa scene, or a refmac verdict."""

    @pytest.mark.parametrize(
        "module",
        [
            "ccp4i2.wrappers.qtpisa.script.qtpisa_report",
            "ccp4i2.pipelines.MakeProjectsAndDoLigandPipeline.script."
            "MakeProjectsAndDoLigandPipeline_report",
            "ccp4i2.pipelines.prosmart_refmac.script.prosmart_refmac_report",
        ],
    )
    def test_it_imports(self, module):
        __import__(module)
