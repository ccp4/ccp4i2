import html
import os
import re
import sys
from pathlib import Path

from ccp4i2.report import Report


def read_pairef_results(job_dir):
    """What paired refinement decided, from the files PAIREF leaves in
    <job>/pairef_project and its log: the suggested cutoff, the change in R
    and R-free each added shell made, and PAIREF's own warnings (too few free
    reflections in a shell, intensities missing). Empty where a file is not
    there yet."""
    job_dir = Path(job_dir)
    project = job_dir / "pairef_project"
    results = {"cutoff": None, "shells": [], "warnings": []}
    cutoff = project / "PAIREF_cutoff.txt"
    if cutoff.is_file():
        results["cutoff"] = cutoff.read_text(encoding="utf-8").strip() or None
    values = project / "project_R-values.csv"
    if values.is_file():
        for line in values.read_text(encoding="utf-8").splitlines():
            fields = line.split()
            if len(fields) != 7 or line.lstrip().startswith("#"):
                continue
            shell, *numbers = fields
            try:
                rw0, rw1, drw, rf0, rf1, drf = (float(x) for x in numbers)
            except ValueError:
                continue
            start, end = shell.split("->") if "->" in shell else (shell, "")
            results["shells"].append({
                "from": start.rstrip("A"), "to": end.rstrip("A"),
                "rwork": (rw0, rw1, drw), "rfree": (rf0, rf1, drf)})
    log = job_dir / "log.txt"
    if log.is_file():
        text = log.read_text(encoding="utf-8", errors="replace")
        block = re.search(r"These warning messages appeared during calculation:\n(.*?)"
                          r"(?:\nResults are listed|\Z)", text, re.S)
        if block:
            results["warnings"] = [w.strip()[len("WARNING:"):].strip()
                                   for w in block.group(1).splitlines()
                                   if w.strip().startswith("WARNING:")]
    return results


class pairef_report(Report):
    TASKNAME = 'pairef'
    USEPROGRAMXML = False
    SEPARATEDATA = True
    RUNNING = True
    WATCHED_FILE = os.path.join("pairef_project", "PAIREF_project.html")

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        super().__init__(xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        running = jobStatus in ["Running", "Running remotely"]
        if not running:
            self.addResults()
        self.addPairedRefinement(running)
        pairef_html = Path(self.jobInfo["fileroot"], "pairef_project", "PAIREF_project.html").resolve()
        if pairef_html.exists():
            projectid = self.jobInfo["projectid"]
            jobNumber = self.jobInfo["jobnumber"]
            url = f"/database/projectId/{projectid}/jobNumber/{jobNumber}/fileName/pairef_project/PAIREF_project.html"
            fold = self.addFold(label="PAIREF's own report", initiallyOpen=False)
            fold.append('<span style="font-size:110%">PAIREF\'s graphs and logs for every step: </span>')
            if not running:
                fold.append(f'<a href="{url}">Open Results</a>')
            else:
                if sys.platform == "win32":
                    pairef_html = pairef_html.as_uri()
                fold.append(f'<a href="{pairef_html}">Open Results</a>')
        elif running:
            self.append("The html report is not ready yet")

    def addPairedRefinement(self, running):
        results = read_pairef_results(self.jobInfo["fileroot"])
        if not results["shells"]:
            self.append("Paired refinement is running." if running
                        else "Paired refinement produced no results.")
            return
        fold = self.addFold(label="Paired refinement", initiallyOpen=True)
        if results["cutoff"] and not running:
            fold.append(f"<p><b>Suggested high-resolution cutoff: {results['cutoff']} Å</b></p>")
        fold.append(
            "<p>Each row adds one shell of data and refines again. The change is "
            "measured at the previous resolution, on the same reflections, so the two "
            "R-free values are comparable: a shell is worth including when adding it "
            "lowers R-free there (a negative change). A change of a few ten-thousandths "
            "is within noise, the more so in a shell with few free reflections.</p>")
        table = fold.addTable(transpose=False, id="pairef_shells")
        shells = results["shells"]
        table.addData(title="Shell added (Å)",
                      data=[f"{s['from']} → {s['to']}" for s in shells])
        table.addData(title="R<sub>work</sub> before", data=["%.4f" % s["rwork"][0] for s in shells])
        table.addData(title="after", data=["%.4f" % s["rwork"][1] for s in shells])
        table.addData(title="R<sub>free</sub> before", data=["%.4f" % s["rfree"][0] for s in shells])
        table.addData(title="after", data=["%.4f" % s["rfree"][1] for s in shells])
        table.addData(title="Change in R<sub>free</sub>", data=["%+.4f" % s["rfree"][2] for s in shells])
        if results["warnings"]:
            fold.append("<p>PAIREF warned:</p><ul>%s</ul>" % "".join(
                f"<li>{html.escape(w)}</li>" for w in results["warnings"]))
