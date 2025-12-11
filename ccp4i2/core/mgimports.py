import ccp4mg  # Modifies sys.path so below imports work

from global_definitions import get_dispobj
from plugins.Sequence import SequenceViewer
from plugins.Sequence.phmmerReport import PhmmerReportNoGui
import global_definitions
import hklfile
import MGApplication
import mmdb2
import mmut
import MolLabel
import PdbView
import pygl_coord
import UtilityThread


def displayTableObjects():
    # If this is imported before a QApplication is created it errors:
    # QWidget: Must construct a QApplication before a QWidget
    import displayTableObjects
    return displayTableObjects


__all__ = [
    "displayTableObjects",
    "get_dispobj",
    "global_definitions",
    "hklfile",
    "MGApplication",
    "mmdb2",
    "mmut",
    "MolLabel",
    "PdbView",
    "PhmmerReportNoGui",
    "pygl_coord",
    "SequenceViewer",
    "UtilityThread",
]
