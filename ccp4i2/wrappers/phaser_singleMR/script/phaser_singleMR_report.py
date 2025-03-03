from __future__ import print_function
from report.CCP4ReportParser import *
from core import CCP4ErrorHandling
import xml.etree.ElementTree as etree

class phaser_singleMR_report(Report):

    TASKNAME = 'phaser_singleMR'
    RUNNING = False

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        if jobStatus is None or jobStatus.lower() == 'nooutput':
            return
        self.defaultReport()

    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        results = self.addResults()
        parent.append("<p>Results for Phaser Single Atom MR Run</p>")
        bda = self.xmlnode.findall('.//SubJob_000/RunDate')[0].text
        bestCC = self.xmlnode.findall('.//SubJob_000/KeyText_003')[0].text
        parent.append("<p>"+bda+"</p>")
        parent.append(bestCC)
        parent.append("<p>Note : the best available solution was selected for I2</p>")
        self.AddGraphsfromXRT(parent)

    def AddGraphsfromXRT(self, parent):
        graph_height = 300
        graph_width = 500
        colour_time = ['gold', 'lightblue', 'red', 'green', 'purple', 'teal', 'blue']
        graph = parent.addFlotGraph(title="Phaser Results", label="testlabel",
                                    style="height:%dpx; width:%dpx; float:left; border:0px;"%(graph_height, graph_width), 
                                    outputXml=self.outputXml, internalId="PhaserGraphs")
        # The xrt file contains graph info (note, this needs to be correctly set in a format specific to the ccp4).
        loadf = etree.parse(os.path.join(self.jobInfo['fileroot'], "program.xrt"))
        if not loadf:
            print("Failed to load xrt file in phaser_singleMR_reports : Report will be missing")
            return
        rootconf = loadf.getroot()
        graphinfo = rootconf.findall(".//xrt:results/xrt:section/xrt:graph/xrt:table", namespaces={"xrt" : "http://www.ccp4.ac.uk/xrt"})
        colcount = 0
        for xrt_inf in graphinfo:
            p = graph.addPlotObject()
            gr_title = xrt_inf.get("title")
            p.append('title', gr_title)
            p.append('plottype', 'xy')
            p.append('xintegral', 'true')
            gr_path = xrt_inf.get("select")  # This is the xml path to the actual data
            gr_path = ".//"+gr_path.lstrip("/Job/")
            dcols = xrt_inf.findall("xrt:data", namespaces={"xrt" : "http://www.ccp4.ac.uk/xrt"})
            xml_data = self.xmlnode.findall(gr_path)[0].text
            data_vec = xml_data.split()
            dsets = xrt_inf.findall("xrt:plot", namespaces={"xrt" : "http://www.ccp4.ac.uk/xrt"})
            yleft = []
            firstset = True
            for dset in dsets:
                # The first group of plotlines goes on the left y-axis, all else on the right.
                lineset = dset.findall("plotline")
                for line in lineset:
                    yleft.append(firstset)
                firstset = False
            for i in range(0, len(dcols)): # Loop over columns in graph (i=0 is the x-axis, the rest y)
                invec = data_vec[i::len(dcols)]
                colname = dcols[i].get("title")
                graph.addData(title=colname, data=list(map(float, invec)))
                if i > 0:
                    if i == 1:
                        p.append('ylabel', colname)
                        p.append('xlabel', dcols[0].get("title"))
                    try:
                        if yleft[i-1]:
                            l = p.append('plotline', xcol=colcount+1, ycol=colcount+i+1, rightaxis='false') # ycol doesn't reset in base classes, hence the counter....
                        else:
                            l = p.append('plotline', xcol=colcount+1, ycol=colcount+i+1, rightaxis='true')
                    except IndexError:
                        l = p.append('plotline', xcol=colcount+1, ycol=colcount+i+1, rightaxis='false')
                        print("PROBLEM: Graphs for singleMR: Issue deciding if y-axis on left or right. Check xml & xrt files.")
                    l.append('label', dcols[i].get("title"))
                    if i < len(colour_time):
                        l.append('colour', colour_time[i-1])
                    else:
                        print("Warning : out of graph colours, increase vector in phaser_singleMR") # Really shouldn't happen, but just in case.
                        l.append('colour', 'black')
            colcount += len(dcols)

