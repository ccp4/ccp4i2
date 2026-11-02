"""The wwPDB validation web service (the ``onedep`` client package)."""
import sys
import time


class ValidationError(RuntimeError):
    pass


def _check(ret, log):
    if ret.get("onedep_error_flag"):
        raise ValidationError("OneDep error: %s" % ret.get("onedep_status_text"))
    if "status" in ret:
        log("OneDep status: %s" % ret["status"])


def run_validation(model_path, sf_path, pdf_out, xml_out, svg_out, log_out, log=print,
                   max_wait_seconds=3600):
    """Submit model + structure factors, wait, and fetch the report (PDF, XML,
    the percentile slider SVG). Returns ``'completed'`` or ``'failed'``.

    Polling backs off quadratically (2, 8, 18, 32 ... seconds) as the service
    asks; a failed run leaves its log at ``log_out``.
    """
    from onedep import __apiUrl__
    from onedep.api.Validate import Validate

    val = Validate(apiUrl=__apiUrl__)
    _check(val.newSession(), log)
    _check(val.inputModelXyzFile(str(model_path)), log)
    _check(val.inputStructureFactorFile(str(sf_path)), log)
    _check(val.run(), log)
    waited = 0
    it = 0
    status = "failed"
    while waited < max_wait_seconds:
        it += 1
        pause = min(it * it * 2, 120)
        time.sleep(pause)
        waited += pause
        ret = val.getStatus()
        if ret.get("onedep_error_flag"):
            ret["status"] = "failed"
        if ret.get("status") in ("completed", "failed"):
            status = ret["status"]
            log("validation %s" % status)
            break
        log("[%4d] waited %d s for the validation service" % (it, waited))
        sys.stdout.flush()
    if status == "completed":
        _check(val.getReport(str(pdf_out)), log)
        _check(val.getReportData(str(xml_out)), log)
        try:
            val.getOutputByType(str(svg_out), contentType="validation-report-slider")
        except Exception as err:
            log("validation slider not downloaded: %s" % err)
    else:
        try:
            val.getReportLog(str(log_out))
        except Exception as err:
            log("validation log not downloaded: %s" % err)
    return status
