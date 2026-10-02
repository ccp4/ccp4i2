"""molrep_selfrot reads molrep.doc. "INFO: pseudo-translation was not detected."
opened the pseudo-translation block, which swallowed every rotation peak (no
SelfRotation in program.xml) and drew a "detected" banner under Molrep's
"not detected"; and the origin Patterson peak, its columns fused
("0.6335E+05281.92"), was dropped. Lines from a real 1h1s run."""
from ccp4i2.core.tasks import get_plugin_class

DOC = """\
 --- Patterson ---
 Number of peaks :      3
      IX  IY  IZ  Xfrac  Yfrac  Zfrac  Xort   Yort    Zort   Dens Dens/sigma
   1   0   0   0  0.000  0.000  0.000   0.00   0.00   0.00  0.6335E+05281.92
   2   3   1   3  0.044  0.005  0.022   3.29   0.69   3.29   2832.     12.60
   3   2   4   0  0.027  0.032  0.000   1.97   4.36   0.00   2361.     10.51
INFO: pseudo-translation was not detected.

 Space group : P 21 21 21
            theta    phi     chi    alpha    beta   gamma      Rf    Rf/sigma
Sol_Rf   1     0.00    0.00    0.00    0.00    0.00    0.00    0.1014E+06 25.16
Sol_Rf   2   128.12   27.86  180.00   27.86 -103.76  152.14     9149.      2.27
"""


def test_not_detected_leaves_rotation_peaks_and_origin(tmp_path):
    from lxml import etree
    (tmp_path / "molrep.doc.txt").write_text(DOC)
    plugin = get_plugin_class("molrep_selfrot")(workDirectory=str(tmp_path), name="molrep_selfrot")
    plugin.xmlnode = etree.Element("MolrepResult")
    plugin.scrapeDocFile()
    root = plugin.xmlnode
    assert root.find(".//PseudoTranslation") is None
    assert len(root.findall(".//SelfRotation/*")) >= 2
    assert [p.findtext("No") for p in root.findall(".//Patterson/Peak")] == ["1", "2", "3"]
