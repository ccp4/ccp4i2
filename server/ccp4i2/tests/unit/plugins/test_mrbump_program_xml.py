"""MrBUMP's final solutions reach program.xml, where a judgement reads them.

Its quick mode wrote them only to results.txt, so the judgement had to read
the rendered report, which existed only once someone had opened it.
"""
import xml.etree.ElementTree as ET

from ccp4i2.core.tasks import get_plugin_class

# The final table of MDM2 job 113 (docs scenario), as MrBUMP wrote it
RESULTS = """
########################################################################################################
#####                            Final MR solution from Phaser...                                  #####
########################################################################################################

                                                       |   Scores                |  Info        |  Phaser                      |  Refmac       
   #                            Model name  Copy  RID  |     eLLG  SeqID  Cover  |   Res  Expt  |     RFZ   TFZ    LLG     SG  |       R  Rfree
   1  loc0_ALL_selected_0_CH_s0.55_r12-100     1    1  |  2928.40   54.5   0.91  |   0.0  UNKN  |     3.4  12.7  106.0  P6522  |    0.51   0.52
"""


def test_the_final_table_is_recorded_in_program_xml(tmp_path):
    plugin = get_plugin_class("mrbump_basic")(workDirectory=str(tmp_path), name="m")
    results = tmp_path / "search_mrbump_1" / "results"
    results.mkdir(parents=True)
    (results / "results.txt").write_text(RESULTS)
    plugin.recordFinalSolutions()
    root = ET.parse(plugin.makeFileName("PROGRAMXML")).getroot()
    (solution,) = root.findall("FinalSolutions/Solution")
    assert {c.tag: c.text for c in solution} == {
        "Model": "loc0_ALL_selected_0_CH_s0.55_r12-100", "Copy": "1", "eLLG": "2928.40",
        "SeqID": "54.5", "Cover": "0.91", "RFZ": "3.4", "TFZ": "12.7", "LLG": "106.0",
        "SpaceGroup": "P6522", "R": "0.51", "Rfree": "0.52"}


def test_mrbumps_own_program_xml_is_added_to_not_replaced(tmp_path):
    plugin = get_plugin_class("mrbump_basic")(workDirectory=str(tmp_path), name="m")
    (tmp_path / "search_mrbump_1" / "results").mkdir(parents=True)
    (tmp_path / "search_mrbump_1" / "results" / "results.txt").write_text(RESULTS)
    with open(plugin.makeFileName("PROGRAMXML"), "w") as f:
        f.write("<MrBUMP><SequenceSearch hits='3'/></MrBUMP>")
    plugin.recordFinalSolutions()
    root = ET.parse(plugin.makeFileName("PROGRAMXML")).getroot()
    assert root.find("SequenceSearch") is not None and root.find("FinalSolutions").get("solutions") == "1"
