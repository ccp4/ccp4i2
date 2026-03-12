"""
Graph and chart report elements.

Graph, FlotGraph, GraphGroup, FlotGraphGroup, GraphLineChooser, DAGGraph.
"""

import re
import xml.etree.ElementTree as etree
from io import StringIO

from ccp4i2.core.CCP4ErrorHandling import CException
from ccp4i2.report.core import (
    ReportClass, Container,
    XRTNS, CCP4NS,
    applySelect,
)
from ccp4i2.report.elements import BaseTable, Plot


class PictureGroup(Container):
    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super(
            PictureGroup,
            self).__init__(
            xrtnode=xrtnode,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.help = None
        self.launch = None


class GraphGroup(Container):
    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super(
            GraphGroup,
            self).__init__(
            xrtnode=xrtnode,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.help = None

        if kw.get('launcher', None) is not None:
            from ccp4i2.report.actions import Launch
            ele = etree.Element('launch')
            ele.set('label', kw['launcher'])
            ele.set('exe', 'loggraph')
            self.launch = Launch(ele, jobInfo=jobInfo)
        else:
            self.launch = None

        if xrtnode is not None:
            from ccp4i2.report.actions import Help
            help = xrtnode.find(XRTNS + 'help')
            if help is not None:
                self.help = Help(help, mode='graph')
                if xrtnode is not None:
                    xrtnode.remove(help)
        elif 'help' in kw:
            from ccp4i2.report.actions import Help
            self.help = Help(ref=kw['help'])


class FlotGraphGroup(Container):
    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super(
            FlotGraphGroup,
            self).__init__(
            xrtnode=xrtnode,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.help = None
        self.launch = None
        self.launchOnly = False

        if kw.get('launcher', None) is not None:
            ele = etree.Element('launch')
            from ccp4i2.report.actions import Launch
            self.launchOnly = True
            ele.set('label', kw['launcher'])
            ele.set('exe', 'loggraph')
            if kw.get('withLaunch', True):
                self.launch = Launch(ele, jobInfo=jobInfo)

        if xrtnode is not None:
            from ccp4i2.report.actions import Help
            help = xrtnode.find(XRTNS + 'help')
            if help is not None:
                self.help = Help(help, mode='graph')
                if xrtnode is not None:
                    xrtnode.remove(help)
        elif 'help' in kw:
            from ccp4i2.report.actions import Help
            self.help = Help(ref=kw['help'])



class DrawnDiv(Container):
    drawnDivCount = 0

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super(
            DrawnDiv,
            self).__init__(
            xrtnode=xrtnode,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.help = None
        DrawnDiv.drawnDivCount += 1
        self.id = kw.get('id', 'drawnDiv_' + str(DrawnDiv.drawnDivCount))
        self.data_is_urls = kw.get('data_is_urls', False)

        if xrtnode is not None:
            from ccp4i2.report.actions import Help
            help = xrtnode.find(XRTNS + 'help')
            if help is not None:
                self.help = Help(help, mode='gallery')
                if xrtnode is not None:
                    xrtnode.remove(help)
        elif 'help' in kw:
            from ccp4i2.report.actions import Help
            self.help = Help(ref=kw['help'])

        self.height = kw.get('height', '250px')
        self.width = kw.get('width', '250px')
        self.style = kw.get('style', 'margin:0px;padding:0px;display:inline-block;') + \
            'height:' + self.height + ';width:' + self.width + ';'
        self.data_data = kw.get('data', 'None')
        self.data_renderer = kw.get('renderer', 'None')
        self.data_require = kw.get('require', 'None')
        self.data_initially_drawn = kw.get('initiallyDrawn', False)


class ObjectGallery(Container):
    galleryCount = 0

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super(
            ObjectGallery,
            self).__init__(
            xrtnode=xrtnode,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.help = None
        ObjectGallery.galleryCount += 1
        self.id = kw.get('id', 'gallery_' + str(ObjectGallery.galleryCount))

        if xrtnode is not None:
            from ccp4i2.report.actions import Help
            help = xrtnode.find(XRTNS + 'help')
            if help is not None:
                self.help = Help(help, mode='gallery')
                if xrtnode is not None:
                    xrtnode.remove(help)
        elif 'help' in kw:
            from ccp4i2.report.actions import Help
            self.help = Help(ref=kw['help'])

        self.tableWidth = kw.get('tableWidth', '7em')
        self.contentWidth = kw.get('contentWidth', '250px')
        self.height = kw.get('height', '250px')
        self.style = kw.get(
            'style',
            'padding:0px;overflow:auto;display:inline-block;margin:1px;')


class GraphLineChooser(Container):
    graphLineChooserCount = 0

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super(
            GraphLineChooser,
            self).__init__(
            xrtnode=xrtnode,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.help = None
        GraphLineChooser.graphLineChooserCount += 1
        self.id = kw.get('id', 'graphLineChooser_' +
                         str(GraphLineChooser.graphLineChooserCount))

        if xrtnode is not None:
            from ccp4i2.report.actions import Help
            help = xrtnode.find(XRTNS + 'help')
            if help is not None:
                self.help = Help(help, mode='gallery')
                if xrtnode is not None:
                    xrtnode.remove(help)
        elif 'help' in kw:
            from ccp4i2.report.actions import Help
            self.help = Help(ref=kw['help'])

        self.height = kw.get('height', '250px')
        self.tableWidth = kw.get('tableWidth', '200px')
        self.contentWidth = kw.get('contentWidth', '200px')
        if 'px' in self.tableWidth and 'px' in self.contentWidth:
            totalWidth = str(
                int(self.tableWidth[:-2]) + int(self.contentWidth[:-2])) + 'px'
        self.style = kw.get('style', 'margin:0px;padding:0px;display:inline-block;') + \
            'height:' + self.height + '; width:' + totalWidth + ';'


class Graph(ReportClass):

    tableCount = 0
    ERROR_CODES = {
        1: {
            'description': 'Plot definition text unreadable. Maybe invalid XML?'}, 2: {
            'description': 'Plot definition file unreadable. Maybe invalid XML?'}, 3: {
                'description': 'Unable to create RTF file'}}

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super(
            Graph,
            self).__init__(
            xrtnode=xrtnode,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.id = kw.get('id', 'graph_' + str(BaseTable.tableCount))
        BaseTable.tableCount += 1
        self.coldata = []
        self.coltitle = []
        self.plots = []
        self.title = None
        self.tableText = ''
        self.headerText = ''
        self.launch = None
        self.headerSeparator = None
        self.pimpleData = None
        self.outputCsv = kw.get('outputCsv', True)

        # Find title
        if xrtnode is not None:
            if xrtnode.find(XRTNS + 'title') is not None:
                if xrtnode.find(
                        XRTNS + 'title').get('select') is not None and xmlnode is not None:
                    xmlEleList = xmlnode.findall(
                        xrtnode.find(XRTNS + 'title').get('select'))
                    if len(xmlEleList) > 0:
                        if isinstance(xmlEleList[0], str):
                            self.title = xmlEleList[0]
                        else:
                            self.title = xmlEleList[0].text
                    else:
                        self.title = None
                else:
                    self.title = xrtnode.find('title').text
            else:
                self.title = xrtnode.get('title', None)
        elif 'title' in kw:
            self.title = kw['title']

        if kw.get('launcher', None) is not None:
            from ccp4i2.report.actions import Launch
            ele = etree.Element('launch')
            ele.set('label', kw.get('launcher', 'More graphs'))
            ele.set('exe', 'loggraph')
            self.launch = Launch(
                ele,
                jobInfo=jobInfo,
                ccp4_data_id='data_' +
                self.internalId)

        if xrtnode is not None:
            from ccp4i2.report.actions import Help
            help = xrtnode.find(XRTNS + 'help')
        else:
            help = kw.get('help', None)
        if help is not None:
            from ccp4i2.report.actions import Help
            self.help = Help(help, mode='graph')
        else:
            self.help = None

        if xrtnode is not None:
            if xrtnode.get("select") is not None and xmlnode is not None:
                # Make list of selectd xml nodes
                self.xmldata = xmlnode.findall(xrtnode.get("select"))
        elif 'select' in kw and xmlnode is not None:
            self.xmldata = []
            for p in kw['select'].split("|"):
                self.xmldata.extend(xmlnode.findall(p.strip()))
        elif 'selectNodes' in kw and xmlnode is not None:
            self.xmldata = kw['selectNodes']
        else:
            # or put current xml node in a list
            self.xmldata = []
            if xmlnode is not None:
                self.xmldata.append(xmlnode)

         # assemble column data
        if xrtnode is not None:
            for node in xrtnode:
                if node.tag == XRTNS + "data":
                    self.addData(
                        xmldata=self.xmldata,
                        title=node.get('title'),
                        select=node.get('select'),
                        expr=node.get('expr'))
            # Handle table as a block in the xml file
                elif node.tag == XRTNS + "table":
                    self.addTable(
                        xmldata=self.xmldata,
                        select=node.get('select'),
                        headers=node.find(
                            XRTNS + 'headers'))

        # get plot definitions
        if xrtnode is not None:
            for node in xrtnode:
                if node.tag == XRTNS + "plot":
                    self.addPlot(
                        xrtnode=node,
                        xmlnode=xmlnode,
                        select=node.get('select'))

    def addPimpleData(self, xmlnode=None, select=None, usePlotly=False):
        if xmlnode is None:
            xmlnode = self.xmlnode
        if select is not None:
            xmlnode = xmlnode.findall(select)
        self.pimpleData = etree.Element('pimple_data')
        if usePlotly:
            self.pimpleData.set('usePlotly', 'True')
        from copy import deepcopy
        for key in ['headers', 'data', 'plot']:
            eleList = xmlnode.findall(key)
            for ele in eleList:
                self.pimpleData.append(deepcopy(ele))
        for attr in list(xmlnode.attrib.keys()):
            self.pimpleData.set(attr, xmlnode.get(attr))

    def addData(
            self,
            xmldata=None,
            title=None,
            select=None,
            expr=None,
            data=[]):
        colvals = []
        if len(data) > 0:
            colvals.extend(data)
        elif select:
            if xmldata is None:
                xmldata = self.xmldata
            for x in xmldata:
                selectedList = x.findall(select)
                for selitem in selectedList:
                    if selitem.text is None:
                        val = '-'
                    else:
                        val = selitem.text.strip()
                        if expr is not None:  # allow an expression to be applied to the data
                            try:
                                val = eval(expr, {"x": float(val)})
                            except BaseException:
                                val = '-'
                    colvals.append(val)
        self.coldata.append(colvals)
        self.coltitle.append(title)

    def addTable(self, xmlnode=None, **kw):
        if xmlnode is None:
            if len(self.xmldata) > 0:
                xmlnode = self.xmldata[0]
            else:
                return
        if kw.get('select', None) is not None and xmlnode is not None:
            tableEleList = xmlnode.findall(kw['select'])
            if len(tableEleList) > 0:
                self.tableText = tableEleList[0].text

        if kw.get('headers', None) is not None and xmlnode is not None:
            headersEleList = xmlnode.findall(kw['headers'])
            if len(headersEleList) > 0:
                self.headerText = headersEleList[0].text
                self.headerSeparator = headersEleList[0].get('separator', None)

    def addPlot(self, xrtnode=None, xmlnode=None, **kw):
        if xmlnode is None:
            if len(self.xmldata) > 0:
                xmlnode = self.xmldata[0]
        if xrtnode is None:
            if kw.get('plotFile', None):
                from ccp4i2.core import CCP4Utils
                try:
                    text = CCP4Utils.readFile(kw['plotFile'])

                    text = re.sub('<xrt:', '<', text)
                    text = re.sub('</xrt:', '</', text)
                    xrtnode = etree.fromstring(text)
                except BaseException:
                    raise CException(self.__class__, 1, str(kw['plotFile']))
            if kw.get('plot', None) is not None:
                if hasattr(kw['plot'], "decode"):
                    text = re.sub('<xrt:', '<', kw['plot'].decode())
                else:
                    text = re.sub('<xrt:', '<', kw['plot'])
                text = re.sub('</xrt:', '</', text)
                xrtnode = etree.fromstring(text)

        if xrtnode is not None and kw.get('select', None) is None:
            # Just copy plot directives from xrt - substitute infor any
            # 'select' attribute
            xrtnode.tag = 'plot'
            xrtnode = applySelect(xrtnode, xmlnode)
            self.plots.append(xrtnode)
        elif xmlnode is not None and kw.get('select', None) is not None:
            # Copy plot directives from xml
            plotEleList = xmlnode.findall(kw['select'])
            for plotEle in plotEleList:
                plotEle.tag = 'plot'
                self.plots.append(plotEle)

    def addPlotObject(self):
        plot = Plot()
        self.plots.append(plot)
        return plot

    def makeTableText(self):
        # pad data arrays to a uniform length
        if len(self.coldata) > 0:
            maxlen = max([len(data) for data in self.coldata])
            for i in range(len(self.coldata)):
                while len(self.coldata[i]) < maxlen:
                    self.coldata[i].append(None)
            # Convert to single text block
            text = ''
            for i in range(len(self.coldata[0])):
                for col in self.coldata:
                    if col[i] is None:
                        text += '- '
                    else:
                        text += str(col[i]) + ' '
                text += "\n"
            self.tableText = text
            # Headers to single text block
            for head in self.coltitle:
                self.headerText += head + ' '
        self.coldata = []
        self.coltitle = []

    def data_as_etree(self):
        eleTree = etree.parse(
            StringIO(
                '<ccp4:ccp4_data xmlns:ccp4="' +
                CCP4NS +
                '"></ccp4:ccp4_data>'))
        ccp4_data = eleTree.getroot()
        if self.title is not None:
            ccp4_data.set('title', self.title)
        ccp4_data.set('id', 'data_' + self.internalId)
        ccp4_data.set('style', 'display:none;')

        if self.pimpleData is not None:
            import copy

            for item in self.pimpleData:
                ccp4_data.append(copy.deepcopy(item))
            if self.pimpleData.get('title') is not None:
                ccp4_data.set('title', self.pimpleData.get('title'))
            if self.pimpleData.get('usePlotly') is not None:
                ccp4_data.set('usePlotly', self.pimpleData.get('usePlotly'))
        else:
            header = etree.Element('headers')
            header.text = self.headerText
            if self.headerSeparator is not None:
                header.set('separator', self.headerSeparator)
            ccp4_data.append(header)

            data = etree.Element('data')
            data.text = self.tableText
            ccp4_data.append(data)

            for plot in self.plots:
                if isinstance(plot, Plot):
                    ccp4_data.append(plot.as_etree())
                else:
                    ccp4_data.append(plot)

        return ccp4_data



class FlotGraph(Graph):

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super(
            FlotGraph,
            self).__init__(
            xrtnode=xrtnode,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.initiallyDrawn = kw.get('initiallyDrawn', True)

        self.launch = None
        self.launchOnly = False
        self.launcherLabel = None
        self.flot_id = kw.get('internalId', None)
        ele = etree.Element('launch')
        if kw.get('launcher', None) is not None:
            from ccp4i2.report.actions import Launch
            self.launchOnly = True
            self.launcherLabel = kw.get('launcher', 'More graphs')
            ele.set('label', self.launcherLabel)
            ele.set('exe', 'loggraph')
            if kw.get('withLaunch', True):
                self.launch = Launch(
                    ele,
                    jobInfo=jobInfo,
                    ccp4_data_id='data_' +
                    self.internalId)

    def as_data_etree(self):
        self.makeTableText()
        root = super().as_data_etree()
        # Add launcher attribute if this graph should be launched in a separate window
        if self.launchOnly and self.launcherLabel:
            root.set('launcher', self.launcherLabel)
        root.append(self.data_as_etree())
        return root


class DAGGraph(Container):
    """Report element for rendering a directed acyclic graph.

    The elements string contains vis-network JSON ({nodes, edges})
    built from PhaserTNG dag.html solution pathway files.
    The frontend CCP4i2ReportDAG component renders it using
    vis-network with a hierarchical layout matching pyvis.
    """

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super().__init__(xrtnode=xrtnode, xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        self.title = kw.get('title', '')
        self.elements = kw.get('elements', '[]')
        self.layout = kw.get('layout', 'dagre')

    def as_data_etree(self):
        root = super().as_data_etree()
        root.set('title', self.title)
        root.set('layout', self.layout)
        root.text = self.elements
        return root
