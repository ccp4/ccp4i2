from ..qtgui.CCP4MessageBox import CMessageBox
from .CCP4Config import GRAPHICAL
from .CCP4ErrorHandling import CErrorReport, Severity


def warningMessage(
    errorReport: CErrorReport,
    windowTitle: str = "",
    message: str = "",
    jobId= None,
    parent=None,
    ifStack: bool = True,
    minSeverity: Severity = Severity.UNDEFINED,
):
    if len(message) > 0:
        message = message.rstrip() + "\n"
    report = errorReport.report(ifStack=ifStack, minSeverity=minSeverity)
    if GRAPHICAL() and parent is not None:
        CMessageBox(parent, windowTitle, message, details=report, jobId=jobId).show()
    else:
        print(report)
