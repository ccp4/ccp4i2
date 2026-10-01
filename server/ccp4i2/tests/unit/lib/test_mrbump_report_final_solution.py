"""MrBUMP's results are fixed-width text wider than the report's box, which
cut off TFZ, the space group and R-free. The report now tables the final
solution; this pins the parsing against MrBUMP's real layout."""
from ccp4i2.wrappers.mrbump_basic.script.mrbump_basic_report import final_solutions

RESULTS = """
#####                            Final MR solution from Phaser...                                  #####
########################################################################################################

                                                       |   Scores                |  Info        |  Phaser                      |  Refmac
   #                            Model name  Copy  RID  |     eLLG  SeqID  Cover  |   Res  Expt  |     RFZ   TFZ    LLG     SG  |       R  Rfree
   1  loc0_ALL_selected_0_CH_s0.55_r12-100     1    1  |  2928.40   54.5   0.91  |   0.0  UNKN  |     3.4  12.7  106.0  P6522  |    0.51   0.52
"""


def test_final_solution_row():
    assert final_solutions(RESULTS) == [{
        "model": "loc0_ALL_selected_0_CH_s0.55_r12-100", "rfz": "3.4", "tfz": "12.7",
        "llg": "106.0", "sg": "P6522", "r": "0.51", "rfree": "0.52"}]


def test_no_final_block():
    assert final_solutions("Molecular Replacement results... | 1 | 2 |") == []
    assert final_solutions(None) == []
