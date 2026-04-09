"""
Graph and chart report elements.

Graph, FlotGraph, GraphGroup, FlotGraphGroup, GraphLineChooser, DAGGraph.

These elements render data as interactive charts in the report frontend.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as etree
from io import StringIO
from typing import Any

from ccp4i2.core.CCP4ErrorHandling import CException
from ccp4i2.report.core import ReportClass, Container, CCP4NS
from ccp4i2.report.elements import BaseTable, Plot


class PictureGroup(Container):
    """Container for a group of Picture elements."""

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        super(PictureGroup, self).__init__(
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.help = None
        self.launch = None


class GraphGroup(Container):
    """Container for a group of Graph elements with optional loggraph launcher."""

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        super(GraphGroup, self).__init__(
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

        if 'help' in kw:
            from ccp4i2.report.actions import Help
            self.help = Help(ref=kw['help'])


class FlotGraphGroup(Container):
    """Container for a group of FlotGraph elements.

    Supports a ``launcher`` kwarg to open all graphs in loggraph,
    and ``withLaunch=False`` to suppress the launch button.
    """

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        super(FlotGraphGroup, self).__init__(
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.help = None
        self.launch = None
        self.launchOnly: bool = False

        if kw.get('launcher', None) is not None:
            ele = etree.Element('launch')
            from ccp4i2.report.actions import Launch
            self.launchOnly = True
            ele.set('label', kw['launcher'])
            ele.set('exe', 'loggraph')
            if kw.get('withLaunch', True):
                self.launch = Launch(ele, jobInfo=jobInfo)

        if 'help' in kw:
            from ccp4i2.report.actions import Help
            self.help = Help(ref=kw['help'])



class DrawnDiv(Container):
    """A container whose content is rendered by a JavaScript renderer.

    Used for custom visualizations (e.g. Ramachandran plots) that require
    client-side drawing with a specified ``renderer`` and ``data``.
    """

    drawnDivCount: int = 0

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        super(DrawnDiv, self).__init__(
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.help = None
        DrawnDiv.drawnDivCount += 1
        self.id: str = kw.get('id', 'drawnDiv_' + str(DrawnDiv.drawnDivCount))
        self.data_is_urls: bool = kw.get('data_is_urls', False)

        if 'help' in kw:
            from ccp4i2.report.actions import Help
            self.help = Help(ref=kw['help'])

        self.height: str = kw.get('height', '250px')
        self.width: str = kw.get('width', '250px')
        self.style: str = kw.get('style', 'margin:0px;padding:0px;display:inline-block;') + \
            'height:' + self.height + ';width:' + self.width + ';'
        self.data_data: str = kw.get('data', 'None')
        self.data_renderer: str = kw.get('renderer', 'None')
        self.data_require: str = kw.get('require', 'None')
        self.data_initially_drawn: bool = kw.get('initiallyDrawn', False)


class ObjectGallery(Container):
    """Gallery container showing selectable objects (e.g. molecules) with a sidebar table."""

    galleryCount: int = 0

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        super(ObjectGallery, self).__init__(
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.help = None
        ObjectGallery.galleryCount += 1
        self.id: str = kw.get('id', 'gallery_' + str(ObjectGallery.galleryCount))

        if 'help' in kw:
            from ccp4i2.report.actions import Help
            self.help = Help(ref=kw['help'])

        self.tableWidth: str = kw.get('tableWidth', '7em')
        self.contentWidth: str = kw.get('contentWidth', '250px')
        self.height: str = kw.get('height', '250px')
        self.style: str = kw.get(
            'style',
            'padding:0px;overflow:auto;display:inline-block;margin:1px;')


class GraphLineChooser(Container):
    """Side-by-side table + graph where selecting a table row highlights a graph line."""

    graphLineChooserCount: int = 0

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        super(GraphLineChooser, self).__init__(
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.help = None
        GraphLineChooser.graphLineChooserCount += 1
        self.id: str = kw.get('id', 'graphLineChooser_' +
                         str(GraphLineChooser.graphLineChooserCount))

        if 'help' in kw:
            from ccp4i2.report.actions import Help
            self.help = Help(ref=kw['help'])

        self.height: str = kw.get('height', '250px')
        self.tableWidth: str = kw.get('tableWidth', '200px')
        self.contentWidth: str = kw.get('contentWidth', '200px')
        if 'px' in self.tableWidth and 'px' in self.contentWidth:
            totalWidth = str(
                int(self.tableWidth[:-2]) + int(self.contentWidth[:-2])) + 'px'
        self.style: str = kw.get('style', 'margin:0px;padding:0px;display:inline-block;') + \
            'height:' + self.height + '; width:' + totalWidth + ';'


class Graph(ReportClass):
    """Data graph with column data, plot definitions, and optional pimple data.

    Data columns are added via ``addData()``, plots via ``addPlot()`` or
    ``addPlotObject()``. The ``data_as_etree()`` method builds the CCP4
    data element consumed by the frontend charting component.
    """

    tableCount: int = 0
    ERROR_CODES: dict = {
        1: {
            'description': 'Plot definition text unreadable. Maybe invalid XML?'}, 2: {
            'description': 'Plot definition file unreadable. Maybe invalid XML?'}, 3: {
                'description': 'Unable to create RTF file'}}

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        super(Graph, self).__init__(
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.id: str = kw.get('id', 'graph_' + str(BaseTable.tableCount))
        BaseTable.tableCount += 1
        self.coldata: list[list] = []
        self.coltitle: list[str | None] = []
        self.plots: list[etree.Element | Plot] = []
        self.title: str | None = kw.get('title', None)
        self.tableText: str = ''
        self.headerText: str = ''
        self.launch = None
        self.headerSeparator: str | None = None
        self.pimpleData: etree.Element | None = None
        self.outputCsv: bool = kw.get('outputCsv', True)

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

        help = kw.get('help', None)
        if help is not None:
            from ccp4i2.report.actions import Help
            self.help = Help(help, mode='graph')
        else:
            self.help = None

        if 'select' in kw and xmlnode is not None:
            self.xmldata: list[etree.Element] = []
            for p in kw['select'].split("|"):
                self.xmldata.extend(xmlnode.findall(p.strip()))
        elif 'selectNodes' in kw and xmlnode is not None:
            self.xmldata = kw['selectNodes']
        else:
            # or put current xml node in a list
            self.xmldata = []
            if xmlnode is not None:
                self.xmldata.append(xmlnode)

    def addPimpleData(
        self,
        xmlnode: etree.Element | None = None,
        select: str | None = None,
        usePlotly: bool = False,
    ) -> None:
        """Add pimple (matplotlib) data from XML."""
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
        xmldata: list[etree.Element] | None = None,
        title: str | None = None,
        select: str | None = None,
        expr: str | None = None,
        data: list | None = None,
    ) -> None:
        """Add a column of data to the graph."""
        if data is None:
            data = []
        colvals: list = []
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

    def addTable(self, xmlnode: etree.Element | None = None, **kw: Any) -> None:
        """Add table data from XML (headers + data block)."""
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

    def addPlot(
        self,
        xmlnode: etree.Element | None = None,
        **kw: Any,
    ) -> None:
        """Add a plot definition from file, text, or XML select."""
        if xmlnode is None:
            if len(self.xmldata) > 0:
                xmlnode = self.xmldata[0]

        plot_node = None
        if kw.get('plotFile', None):
            from ccp4i2.core import CCP4Utils
            try:
                text = CCP4Utils.readFile(kw['plotFile'])
                text = re.sub('<xrt:', '<', text)
                text = re.sub('</xrt:', '</', text)
                plot_node = etree.fromstring(text)
            except BaseException:
                raise CException(self.__class__, 1, str(kw['plotFile']))
        if kw.get('plot', None) is not None:
            if hasattr(kw['plot'], "decode"):
                text = re.sub('<xrt:', '<', kw['plot'].decode())
            else:
                text = re.sub('<xrt:', '<', kw['plot'])
            text = re.sub('</xrt:', '</', text)
            plot_node = etree.fromstring(text)

        if plot_node is not None and kw.get('select', None) is None:
            plot_node.tag = 'plot'
            self.plots.append(plot_node)
        elif xmlnode is not None and kw.get('select', None) is not None:
            # Copy plot directives from xml
            plotEleList = xmlnode.findall(kw['select'])
            for plotEle in plotEleList:
                plotEle.tag = 'plot'
                self.plots.append(plotEle)

    def addPlotObject(self) -> Plot:
        """Create and add an empty Plot element, returning it for configuration."""
        plot = Plot()
        self.plots.append(plot)
        return plot

    def makeTableText(self) -> None:
        """Convert column data arrays into a single text block for CCP4 data format."""
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

    def data_as_etree(self) -> etree.Element:
        """Build the CCP4 data element tree for the frontend charting component."""
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
    """Graph variant with Flot/Plotly rendering and optional external launcher."""

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        super(FlotGraph, self).__init__(
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.initiallyDrawn: bool = kw.get('initiallyDrawn', True)

        self.launch = None
        self.launchOnly: bool = False
        self.launcherLabel: str | None = None
        self.flot_id: str | None = kw.get('internalId', None)
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

    def as_data_etree(self) -> etree.Element:
        self.makeTableText()
        root = super().as_data_etree()
        # Add launcher attribute if this graph should be launched in a separate window
        if self.launchOnly and self.launcherLabel:
            root.set('launcher', self.launcherLabel)
        root.append(self.data_as_etree())
        return root


class DAGGraph(Container):
    """Directed acyclic graph visualization (e.g. PhaserTNG solution pathways).

    The ``elements`` string contains vis-network JSON ({nodes, edges}).
    The frontend ``CCP4i2ReportDAG`` component renders it using
    vis-network with a hierarchical layout.
    """

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        super().__init__(xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        self.title: str = kw.get('title', '')
        self.elements: str = kw.get('elements', '[]')
        self.layout: str = kw.get('layout', 'dagre')

    def as_data_etree(self) -> etree.Element:
        root = super().as_data_etree()
        root.set('title', self.title)
        root.set('layout', self.layout)
        root.text = self.elements
        return root
