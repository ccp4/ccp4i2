import json
from pathlib import Path
from ccp4i2.report.CCP4ReportParser import Report


class modelcraft_report(Report):
    TASKNAME = "modelcraft"
    RUNNING = True
    USEPROGRAMXML = False
    WATCHED_FILE = str(Path("modelcraft/modelcraft.json"))

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(
            self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw
        )
        self.addDiv(style="clear:both;")
        if jobStatus in ["Running", "Running remotely"]:
            self.append("<p><b>The job is currently running.</b></p>")
        json_path = Path(jobInfo["fileroot"], "modelcraft/modelcraft.json")
        if json_path.exists():
            with json_path.open(encoding="utf-8") as text:
                self.json = json.load(text)
            self.add_running_job(parent=self)
            summaryFold = self.addFold(
                label="Summary", initiallyOpen=True, brief="Summary"
            )
            self.add_table(parent=summaryFold)
            self.add_message(jobStatus, parent=summaryFold)
            self.add_graph(parent=summaryFold)
            self.add_picture()
            self.addDiv(style="clear:both;")

    def add_picture(self, parent=None):
        if parent is None:
            parent = self
        parent.addDiv(style="clear:both;")
        pictureFold = parent.addFold(
            label="Picture", initiallyOpen=True, brief="Picture"
        )
        baseScenePath = Path(__file__).resolve().parent / "modelcraft_1.scene.xml"
        pictureFold.addPicture(
            label="Autobuilt structure", sceneFile=str(baseScenePath), id="autobuild_1"
        )

    def add_running_job(self, parent=None):
        if parent is None:
            parent = self
        subjob = self.json.get("running_job")
        if subjob:
            parent.append(f"<p>The sub-job that is currently running is: {subjob}</p>")

    def add_table(self, parent=None):
        if "final" not in self.json:
            return
        if parent is None:
            parent = self
        cycle = self.json["final"]["cycle"]
        residues = self.json["final"]["residues"]
        waters = self.json["final"]["waters"]
        rwork = round(self.json["final"]["r_work"], 2)
        rfree = round(self.json["final"]["r_free"], 2)
        table = parent.addTable(transpose=True)
        table.addData(title="Cycle", data=[cycle])
        table.addData(title="Residues", data=[residues])
        table.addData(title="Waters", data=[waters])
        table.addData(title="R<sub>Work</sub>", data=[rwork])
        table.addData(title="R<sub>Free</sub>", data=[rfree])

    def add_message(self, jobStatus, parent=None):
        if "final" not in self.json:
            return
        if parent is None:
            parent = self
        cycle = self.json["final"]["cycle"]
        rfree = self.json["final"]["r_free"]
        if jobStatus in ["Running", "Running remotely"]:
            message = f"Cycle {cycle} currently has the best model."
        else:
            message = f"The output model was taken from cycle {cycle}."
        message += " On the basis of the refinement statistics,"
        if rfree > 0.50:
            message += " the model is very incomplete or wrong."
        elif rfree > 0.40:
            message += " the model is substantially incomplete"
            message += " and may contain incorrect regions."
        elif rfree > 0.35:
            message += " the model is likely to contain correct regions"
            message += " but requires further work."
        else:
            message += " the model is approaching completion."
        parent.append(f"<p>{message}</p>")

    def add_graph(self, parent=None):
        if len(self.json.get("cycles", [])) == 0:
            return
        if parent is None:
            parent = self
        cycles = []
        residues = []
        rfrees = []
        rworks = []
        for cycle in self.json["cycles"]:
            cycles.append(cycle["cycle"])
            residues.append(cycle["residues"])
            rfrees.append(cycle["r_free"])
            rworks.append(cycle["r_work"])
        title = "Residues and R-factors VS cycle"
        style = "width:600px; height:300px; border:0px;"
        graph = parent.addFlotGraph(title=title, style=style)
        graph.addData(title="Cycle", data=cycles)
        graph.addData(title="Residues", data=residues)
        graph.addData(title="R-free", data=rfrees)
        graph.addData(title="R-work", data=rworks)
        plot = graph.addPlotObject()
        plot.append("title", title)
        plot.append("plottype", "xy")
        plot.append("xlabel", "Cycle")
        plot.append("xintegral", "true")
        plot.append("ylabel", "Residues")
        plot.append("yintegral", "true")
        # plot.append("ylabel", "R-factor", rightaxis="true")
        line = plot.append("plotline", xcol=1, ycol=2)
        line.append("label", "Residues")
        line.append("colour", "blue")
        line = plot.append("plotline", xcol=1, ycol=3, rightaxis="true")
        line.append("label", "R<sub>Free</sub>")
        line.append("colour", "gold")
        line = plot.append("plotline", xcol=1, ycol=4, rightaxis="true")
        line.append("label", "R<sub>Work</sub>")
        line.append("colour", "lightblue")
