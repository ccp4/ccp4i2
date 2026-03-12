import logging
import re

from ccp4i2.report import Report

logger = logging.getLogger(f"ccp4i2:{__name__}")


class phaser_singleMR_report(Report):
    TASKNAME = 'phaser_singleMR'
    RUNNING = False

    def __init__(self, xmlnode=None, jobInfo=None, jobStatus=None, **kw):
        if jobInfo is None:
            jobInfo = {}
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        if jobStatus is None or jobStatus.lower() == 'nooutput':
            return
        self.defaultReport()

    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        self.addResults()
        parent.append("<p>Results for Phaser Single Atom MR Run</p>")
        bda = self.xmlnode.findall('.//SubJob_000/RunDate')[0].text
        parent.append("<p>" + bda + "</p>")
        # KeyText from smartie may be wrapped in <pre> tags; strip them
        # and use addPre which creates a proper CCP4i2ReportPre element.
        bestCC = self.xmlnode.findall('.//SubJob_000/KeyText_003')[0].text
        bestCC = re.sub(r'^\s*<pre>\s*', '', bestCC, flags=re.IGNORECASE)
        bestCC = re.sub(r'\s*</pre>\s*$', '', bestCC, flags=re.IGNORECASE)
        parent.addPre(text=bestCC)
        parent.append("<p>Note : the best available solution was selected for I2</p>")
        self._add_graphs(parent)

    def _add_graphs(self, parent):
        """Build graphs from GraphTable elements embedded in program.xml.

        The wrapper's _convert_log() embeds graph metadata (column titles,
        plot configs, axis scaling) directly in program.xml alongside the
        data, so no separate file is needed.
        """
        graph_tables = self.xmlnode.findall('.//GraphTable')
        if not graph_tables:
            return

        colours = ['gold', 'lightblue', 'red', 'green', 'purple', 'teal', 'blue']
        graph = parent.addFlotGraph(
            title="Phaser Results",
            label="testlabel",
            style="height:300px; width:500px; float:left; border:0px;",
            outputXml=self.outputXml,
            internalId="PhaserGraphs",
        )

        colcount = 0
        for table_info in graph_tables:
            p = graph.addPlotObject()
            p.append('title', table_info.get("title"))
            p.append('plottype', 'xy')
            p.append('xintegral', 'true')

            # Get data from the sibling element referenced by data_ref.
            # ET doesn't support parent traversal, so search from root.
            data_ref = table_info.get("data_ref")
            data_nodes = self.xmlnode.findall(f".//{data_ref}")
            if not data_nodes or not data_nodes[0].text:
                continue
            data_vec = data_nodes[0].text.split()

            dcols = table_info.findall("Column")
            plots = table_info.findall("Plot")

            # Determine which y-columns go on left vs right axis
            yleft = []
            firstset = True
            for plot in plots:
                for _line in plot.findall("PlotLine"):
                    yleft.append(firstset)
                firstset = False

            for i, dcol in enumerate(dcols):
                invec = data_vec[i::len(dcols)]
                colname = dcol.get("title")
                graph.addData(title=colname, data=list(map(float, invec)))
                if i > 0:
                    if i == 1:
                        p.append('ylabel', colname)
                        p.append('xlabel', dcols[0].get("title"))
                    right = 'false'
                    try:
                        right = 'false' if yleft[i - 1] else 'true'
                    except IndexError:
                        pass
                    line = p.append(
                        'plotline',
                        xcol=colcount + 1,
                        ycol=colcount + i + 1,
                        rightaxis=right,
                    )
                    line.append('label', colname)
                    colour = colours[i - 1] if i <= len(colours) else 'black'
                    line.append('colour', colour)

            colcount += len(dcols)
