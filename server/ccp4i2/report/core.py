"""
Core report classes and utilities.

This module contains the foundation classes that all report elements
build upon: ReportClass (base), Container (children management),
Report (top-level), and supporting utilities.
"""

import logging
import os
import xml.etree.ElementTree as etree

from ccp4i2.core.CCP4ErrorHandling import SEVERITY_WARNING, CErrorReport, CException
from ccp4i2 import I2_TOP

# Import error handling (lazy to avoid circular imports)
_diagnostics_module = None


def _get_diagnostics():
    """Lazy import of diagnostics module to avoid circular imports."""
    global _diagnostics_module
    if _diagnostics_module is None:
        from ccp4i2.report import errors as _diagnostics_module
    return _diagnostics_module


logger = logging.getLogger(f"ccp4i2:{__name__}")

XRTNS = "{http://www.ccp4.ac.uk/xrt}"
CCP4NS = "http://www.ccp4.ac.uk/ccp4ns"
XHTMLNS = "{http://www.w3.org/1999/xhtml}"
CURRENT_CSS_VERSION = '0.1.0'



def htmlBase():
    """Return base path for report static files (images, CSS, etc.).

    Static files are served from Next.js public directory.
    Path: /report_files/{version}/
    """
    return '/report_files/' + CURRENT_CSS_VERSION


def toBoolean(text):
    try:
        i = int(text)
        return bool(i)
    except BaseException:
        return bool(text)


def getChildObject(child, xmlnode, jobInfo={}, report=None):
    # Lazy imports to avoid circular dependencies
    from ccp4i2.report.elements import (
        Text, Copy, Generic, Fold, Results, Status,
        BaseTable, Table,
    )
    from ccp4i2.report.graphs import Graph, GraphGroup
    from ccp4i2.report.graphs import ObjectGallery
    from ccp4i2.report.pictures import Picture
    from ccp4i2.report.io_data import InputData, OutputData, ImportedFiles
    from ccp4i2.report.metadata import Title, JobDetails, JobLogFiles
    from ccp4i2.report.actions import Launch

    if child.tag == XRTNS + "table":
        obj = Table(child, xmlnode)
    elif child.tag == XRTNS + "graphgroup":
        obj = GraphGroup(child, xmlnode, jobInfo)
    elif child.tag == XRTNS + "graph":
        obj = Graph(child, xmlnode, jobInfo)
    elif child.tag == XRTNS + "objectgallery":
        obj = ObjectGallery(child, xmlnode, jobInfo)
    elif child.tag == XRTNS + "results":
        obj = Results(child, xmlnode, jobInfo)
    elif child.tag == XRTNS + "fold":
        obj = Fold(child, xmlnode)
    elif child.tag == XRTNS + "text":
        obj = Text(child, xmlnode)
    elif child.tag == XRTNS + "status":
        obj = Status(child, xmlnode)
    elif child.tag == XRTNS + "if":
        obj = IfContainer(child, xmlnode, jobInfo)
    elif child.tag == XRTNS + "ifnot":
        obj = IfContainer(child, xmlnode, jobInfo, False)
    elif child.tag == XRTNS + "loop":
        obj = Loop(child, xmlnode, jobInfo)
    elif child.tag == XRTNS + "copy":
        obj = Copy(child, xmlnode)
    elif child.tag == XRTNS + "inputdata":
        obj = InputData(child, xmlnode, jobInfo)
    elif child.tag == XRTNS + "outputdata":
        obj = OutputData(child, xmlnode, jobInfo)
    elif child.tag == XRTNS + "importedfiles":
        obj = ImportedFiles(child, xmlnode, jobInfo)
    elif child.tag == XRTNS + "title":
        obj = Title(child, xmlnode, jobInfo)
    elif child.tag == XRTNS + "jobdetails":
        obj = JobDetails(child, xmlnode, jobInfo)
    elif child.tag == XRTNS + "logfiles":
        obj = JobLogFiles(child, xmlnode, jobInfo)
    elif child.tag == XRTNS + "picture":
        obj = Picture(child, xmlnode, jobInfo)
    elif child.tag == XRTNS + "launch":
        obj = Launch(child, xmlnode, jobInfo)
    elif child.tag == XRTNS + "comment":
        obj = None
    else:
        obj = Generic(child, xmlnode)
    return obj


def findChildren(obj, cls, id=None):
    retList = []
    if not hasattr(obj, 'children'):
        return retList
    for child in obj.children:
        if isinstance(child, cls):
            retList.append(child)
        elif isinstance(child, Container):
            retList.extend(findChildren(child, cls, id))
    return retList


def applySelect(xrtnode, xmlnode, jobInfo={}):
    # Test all nodes in xrtnode for a 'select' attribute
    # Find the value in xmlnode and copy into xrtnode node.text and remove the select attribute
    # Currently this used to process the plotting input
    if xrtnode.get('select') is not None:
        if xmlnode is not None:
            values = xmlnode.findall(xrtnode.get('select'))
            if len(values) > 0:
                xrtnode.text = str(values[0].text)
            else:
                xrtnode.text = 'NO DATA AVAILABLE'
        del xrtnode.attrib['select']
    elif xrtnode.get('database') is not None:
        pathList = xrtnode.get('database').split('/')
        while pathList.count('') > 0:
            pathList.remove('')
        data = jobInfo.get(pathList[0], {})
        for path in pathList[1:]:
            data = data.get(path, {})
        if not isinstance(data, dict):
            if isinstance(data, list):
                # I think the length should always be greater than 0, but a
                # Keith refmac job in 09/02/2022 created a situation where it
                # was not.
                if len(data) > 0:
                    xrtnode.text = str(data[0])
            else:
                xrtnode.text = str(data)
        del xrtnode.attrib['database']
    for child in xrtnode:
        applySelect(child, xmlnode, jobInfo)
    return xrtnode


class ReportClass(object):
    """
    Base class for all report elements.

    Provides common attributes and methods for report elements including
    ID generation, parent traversal, and XML serialization.

    Attributes:
        id: Unique identifier for this element
        label: Human-readable label
        title: Title/tooltip text
        class_: CSS class(es) for styling
        style: Inline CSS styles
        parent: Parent element in the hierarchy
        internalId: Auto-generated internal ID for XML output

    Modern Usage:
        Elements can access the diagnostics system for error reporting:

            self.add_diagnostic('warning', 'TABLE_EMPTY', 'No data rows')
    """

    def __init__(self, *arg, **kw):
        super(ReportClass, self).__init__()
        self.id = kw.get('id', None)
        self.label = kw.get('label', None)
        self.title = kw.get('title', None)
        self.class_ = kw.get('class_', None)
        self.style = kw.get('style', None)
        self.parent = kw.get('parent', None)
        self.outputXml = kw.get('outputXml', False)
        self.internalId = kw.get('internalId', None)

        # Modern diagnostics collector (lazy initialized)
        self._diagnostics = None

        if self.internalId is None:
            report = self.getReport()
            if report is not None:
                className = type(self).__name__
                instanceCount = len(findChildren(report, type(self)))
                self.internalId = className + '_' + str(instanceCount)
            else:
                className = type(self).__name__
                self.internalId = className + '_' + str(Report.elementCount)
                Report.elementCount += 1

    @property
    def diagnostics(self):
        """Get or create the diagnostics collector for this element."""
        if self._diagnostics is None:
            errors = _get_diagnostics()
            self._diagnostics = errors.DiagnosticCollector()
        return self._diagnostics

    def add_diagnostic(self, level: str, code: str, message: str, **kwargs):
        """
        Add a diagnostic to this element.

        Args:
            level: 'debug', 'info', 'warning', 'error', or 'critical'
            code: Machine-readable error code (e.g., 'TABLE_EMPTY')
            message: Human-readable description
            **kwargs: Additional context (location, details, exception)
        """
        errors = _get_diagnostics()
        level_enum = errors.DiagnosticLevel(level.lower())
        location = kwargs.pop('location', type(self).__name__)
        return self.diagnostics.add(level_enum, code, message, location=location, **kwargs)

    def warning(self, code: str, message: str, **kwargs):
        """Add a warning diagnostic."""
        return self.add_diagnostic('warning', code, message, **kwargs)

    def error(self, code: str, message: str, **kwargs):
        """Add an error diagnostic."""
        return self.add_diagnostic('error', code, message, **kwargs)

    def getReport(self):
        """
        Get the root Report object.

        Traverses up the parent hierarchy to find the Report.
        """
        if hasattr(self, 'parent') and self.parent is not None:
            return self.parent.getReport()
        return None

    # Modern alias
    def get_report(self):
        """Modern alias for getReport()."""
        return self.getReport()

    def data_id(self):
        """Get the data ID for this element (used in XML output)."""
        return 'data_' + self.internalId

    def data_url(self):
        """Get the data URL for this element."""
        return './tables_as_xml_files/' + self.data_id() + '.xml'

    def as_data_etree(self):
        """
        Generate XML element for frontend consumption.

        Returns:
            ET.Element with tag 'CCP4i2Report{ClassName}'
        """
        root = etree.Element(
            'CCP4i2Report{}'.format(type(self).__name__),
            key=self.internalId,
        )
        root.set('class', self.class_ if self.class_ is not None else '')
        root.set('style', self.style if self.style is not None else '')
        root.text = self.text if (
            hasattr(self, 'text') and self.text is not None) else ''
        root.tail = self.tail if (
            hasattr(self, 'tail') and self.tail is not None) else ''
        return root

# Container class - base class for reports and folds


class Container(ReportClass):
    tag = 'container'
    ERROR_CODES = {
        1: {
            'severity': SEVERITY_WARNING, 'description': 'Failed to interpret xrt node tag'}, 2: {
            'description': 'Can not create IfContainer instance without xrtnode argument (it should not be used in Python mode)'}, 3: {
                'description': 'Can not create Loop instance without xrtnode argument (it should not be used in Python mode)'}}

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super(
            Container,
            self).__init__(
            xrtnode=xrtnode,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        if getattr(self, 'errReport', None) is None:
            self.errReport = CException()
        self.children = []
        self.jobInfo = {}
        self.jobInfo.update(jobInfo)
        self.xrtnode = xrtnode
        self.xmlnode = xmlnode
        self.text = ''
        self.tail = ''
        self._scriptErrorReport = None
        if xrtnode is not None:
            self.interpretXrt(xrtnode=xrtnode)

    def interpretXrt(self, xrtnode=None):
        if xrtnode is not None:
            if xrtnode.tag in [XRTNS + 'remove', 'remove']:
                self.tag = XRTNS + 'remove'
            else:
                self.tag = XHTMLNS + 'root'
                # Initialising 'ifContainer' may skip initialising children if
                # the tag is changed to dummy
                for child in xrtnode:
                    obj = getChildObject(child, self.xmlnode, self.jobInfo)
                    if obj is None:
                        self.errReport.append(
                            self.__class__, 1, str(child.tag))
                    else:
                        self.children.append(obj)
            self.text = xrtnode.text
            self.tail = xrtnode.tail
        else:
            self.text = kw.get('text', None)
            self.tail = kw.get('tail', None)
        if self.text is None:
            self.text = ''
        if self.tail is None:
            self.tail = ''

    def append(self, child):
        if isinstance(child, str):
            from ccp4i2.report.elements import Generic
            self.children.append(Generic(xmlnode=self.xmlnode, text=child))
        else:
            self.children.append(child)

    def insert(self, indx, child):
        if isinstance(child, str):
            from ccp4i2.report.elements import Generic
            self.children.insert(
                indx, Generic(
                    xmlnode=self.xmlnode, text=child))
        else:
            self.children.insert(indx, child)

    def __len__(self):
        return len(self.children)

    def addObjectOfClass(
            self,
            classOfObject,
            xrtnode=None,
            xmlnode=None,
            jobInfo=None,
            **kw):
        if xrtnode is None:
            xrtnode = self.xrtnode
        if xmlnode is None:
            xmlnode = self.xmlnode
        if jobInfo is None:
            jobInfo = self.jobInfo
        obj = classOfObject(
            xrtnode=xrtnode,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            parent=self,
            **kw)
        self.children.append(obj)
        return obj

    def addText(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        """This adds plain text to a report. Tags are converted to references. append is also available."""
        from ccp4i2.report.elements import Text
        return self.addObjectOfClass(Text, xrtnode, xmlnode, jobInfo, **kw)

    def addPre(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        from ccp4i2.report.elements import Pre
        return self.addObjectOfClass(Pre, xrtnode, xmlnode, jobInfo, **kw)

    def addFetchPre(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        from ccp4i2.report.elements import FetchPre
        return self.addObjectOfClass(FetchPre, xrtnode, xmlnode, jobInfo, **kw)

    def addGraph(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        from ccp4i2.report.graphs import Graph
        return self.addObjectOfClass(Graph, xrtnode, xmlnode, jobInfo, **kw)

    def addProgress(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        from ccp4i2.report.elements import Progress
        return self.addObjectOfClass(Progress, xrtnode, xmlnode, jobInfo, **kw)

    def addFlotGraph(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        from ccp4i2.report.graphs import FlotGraph
        return self.addObjectOfClass(
            FlotGraph, xrtnode, xmlnode, jobInfo, **kw)

    def addGraphGroup(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        from ccp4i2.report.graphs import GraphGroup
        return self.addObjectOfClass(
            GraphGroup, xrtnode, xmlnode, jobInfo, **kw)

    def addFlotGraphGroup(
            self,
            xrtnode=None,
            xmlnode=None,
            jobInfo=None,
            **kw):
        from ccp4i2.report.graphs import FlotGraphGroup
        return self.addObjectOfClass(
            FlotGraphGroup, xrtnode, xmlnode, jobInfo, **kw)

    def addDrawnDiv(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        from ccp4i2.report.graphs import DrawnDiv
        return self.addObjectOfClass(DrawnDiv, xrtnode, xmlnode, jobInfo, **kw)

    def addObjectGallery(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        from ccp4i2.report.graphs import ObjectGallery
        return self.addObjectOfClass(
            ObjectGallery, xrtnode, xmlnode, jobInfo, **kw)

    def addGraphLineChooser(
            self,
            xrtnode=None,
            xmlnode=None,
            jobInfo=None,
            **kw):
        from ccp4i2.report.graphs import GraphLineChooser
        return self.addObjectOfClass(
            GraphLineChooser, xrtnode, xmlnode, jobInfo, **kw)

    def addPictureGroup(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        from ccp4i2.report.graphs import PictureGroup
        return self.addObjectOfClass(
            PictureGroup, xrtnode, xmlnode, jobInfo, **kw)

    def addDiv(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        from ccp4i2.report.elements import Div
        return self.addObjectOfClass(Div, xrtnode, xmlnode, jobInfo, **kw)

    def addTable(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        from ccp4i2.report.elements import Table
        return self.addObjectOfClass(Table, xrtnode, xmlnode, jobInfo, **kw)

    def addPicture(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        from ccp4i2.report.pictures import Picture
        return self.addObjectOfClass(Picture, xrtnode, xmlnode, jobInfo, **kw)

    def addTitle(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        from ccp4i2.report.metadata import Title
        return self.addObjectOfClass(Title, xrtnode, xmlnode, jobInfo, **kw)

    def addLaunch(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        from ccp4i2.report.actions import Launch
        return self.addObjectOfClass(Launch, xrtnode, xmlnode, jobInfo, **kw)

    def addDownload(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        from ccp4i2.report.actions import Download
        return self.addObjectOfClass(Download, xrtnode, xmlnode, jobInfo, **kw)

    def addCopyToClipboard(
            self,
            xrtnode=None,
            xmlnode=None,
            jobInfo=None,
            text="",
            label="",
            **kw):
        from ccp4i2.report.actions import CopyToClipboard
        return self.addObjectOfClass(
            CopyToClipboard, text=text, label=label, **kw)

    def addCopyUrlToClipboard(
            self,
            xrtnode=None,
            xmlnode=None,
            jobInfo=None,
            text="",
            label="",
            **kw):
        from ccp4i2.report.actions import CopyUrlToClipboard
        return self.addObjectOfClass(
            CopyUrlToClipboard, text=text, label=label, **kw)

    def addDAGGraph(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        from ccp4i2.report.graphs import DAGGraph
        return self.addObjectOfClass(DAGGraph, xrtnode, xmlnode, jobInfo, **kw)

    def addResults(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        from ccp4i2.report.elements import Results
        return self.addObjectOfClass(Results, xrtnode, xmlnode, jobInfo, **kw)

    def addFold(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        from ccp4i2.report.elements import Fold
        return self.addObjectOfClass(Fold, xrtnode, xmlnode, jobInfo, **kw)

    def addCopy(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        from ccp4i2.report.elements import Copy
        return self.addObjectOfClass(Copy, xrtnode, xmlnode, jobInfo, **kw)

    def addHelp(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        from ccp4i2.report.actions import Help
        return self.addObjectOfClass(Help, xrtnode, xmlnode, jobInfo, **kw)

    # Modern grid layout methods

    def addGrid(self, spacing=2, direction='row', **kw):
        """
        Add a grid container for responsive layouts.

        Creates a MUI Grid container that arranges its children in a
        responsive grid layout.

        Args:
            spacing: Gap between grid items (0-10)
            direction: 'row', 'column', 'row-reverse', 'column-reverse'
            **kw: Additional arguments passed to GridContainer

        Returns:
            GridContainer that can have GridItems added

        Example:
            grid = self.addGrid(spacing=2)
            left = grid.addItem(xs=12, md=8)
            left.append(self.addTable(...))
            right = grid.addItem(xs=12, md=4)
            right.append(self.addGraph(...))
        """
        from ccp4i2.report.grid import GridContainer, GridDirection
        if isinstance(direction, str):
            direction = GridDirection(direction)
        grid = GridContainer(spacing=spacing, direction=direction, parent=self, **kw)
        self.children.append(grid)
        return grid

    def addGridItem(self, xs=12, sm=None, md=None, lg=None, xl=None, **kw):
        """
        Add a grid item directly to this container.

        Usually you would add items to a GridContainer instead, but this
        method allows adding a single grid item for simple layouts.

        Args:
            xs: Columns at extra-small (0-600px) - default 12 (full width)
            sm: Columns at small (600-900px)
            md: Columns at medium (900-1200px)
            lg: Columns at large (1200-1536px)
            xl: Columns at extra-large (1536px+)
            **kw: Additional arguments passed to GridItem

        Returns:
            GridItem that can have content added
        """
        from ccp4i2.report.grid import GridItem, GridSpan
        span = GridSpan(xs=xs, sm=sm, md=md, lg=lg, xl=xl)
        item = GridItem(span=span, parent=self, **kw)
        self.children.append(item)
        return item

    def addTwoColumnLayout(self, left_span=8, right_span=4, spacing=2, **kw):
        """
        Add a two-column layout (convenience method).

        Creates a grid with two columns. Returns a tuple of (left, right)
        containers that you can add content to.

        Args:
            left_span: Column width for left column (1-12, default 8)
            right_span: Column width for right column (1-12, default 4)
            spacing: Gap between columns

        Returns:
            Tuple of (left_item, right_item) GridItems

        Example:
            left, right = self.addTwoColumnLayout()
            left.append(self.addTable(...))
            right.append(self.addGraph(...))
        """
        from ccp4i2.report.grid import GridContainer, GridItem, GridSpan
        grid = GridContainer(spacing=spacing, parent=self, **kw)
        self.children.append(grid)

        left = GridItem(span=GridSpan(xs=12, md=left_span), parent=grid)
        grid.append(left)

        right = GridItem(span=GridSpan(xs=12, md=right_span), parent=grid)
        grid.append(right)

        return left, right

    def getPictures(self):
        from ccp4i2.report.pictures import Picture
        rv = []
        for child in self.children:
            if isinstance(child, Picture):
                rv.append(child.picDefFile.__str__())
            elif isinstance(child, (Container, Loop)):
                rv.extend(child.getPictures())
        return rv

    def as_data_etree(self):
        root = super().as_data_etree()
        for child in self.children:
            if hasattr(child, 'as_data_etree'):
                root.append(child.as_data_etree())
        return root

    def errorReport(self):
        err = CException()
        err.extend(self.errReport)
        for child in self.children:
            if getattr(child, 'errReport', None) is not None:
                err.extend(child.errorReport())
        return err

    def loadXmlFile(self, fileName):
        from ccp4i2.core import CCP4Utils
        text = CCP4Utils.readFile(fileName)
        ele = etree.fromstring(text)
        return ele



# If class - container for conditional content
class IfContainer(Container):
    def __init__(
            self,
            xrtnode=None,
            xmlnode=None,
            jobInfo={},
            target=True,
            **kw):
        if xrtnode is None:
            raise CException(self.__class__, 2)
        if xmlnode is not None and len(
                xmlnode.findall(xrtnode.get("select"))) != 1:
            # Failed to find the 'if' parameter in the output so do not put
            # anything in report
            Container.__init__(
                self,
                etree.Element(
                    XRTNS +
                    "remove"),
                xmlnode,
                jobInfo,
                **kw)
            return
        try:
            value = toBoolean(xmlnode.findall(xrtnode.get("select"))[0].text)
        except BaseException:
            Container.__init__(
                self,
                etree.Element(
                    XRTNS +
                    "remove"),
                xmlnode,
                jobInfo,
                **kw)
            return
        if value == target:
            Container.__init__(self, xrtnode, xmlnode, jobInfo)
        else:
            Container.__init__(
                self,
                etree.Element(
                    XRTNS +
                    "remove"),
                xmlnode,
                jobInfo,
                **kw)


class Loop(ReportClass):

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super(
            Loop,
            self).__init__(
            xrtnode=xrtnode,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        if xrtnode is None:
            raise CException(self.__class__, 3)
        import copy
        self.children = []
        self.text = xrtnode.text
        self.tail = xrtnode.tail
        selectStr = xrtnode.get('select')
        if selectStr is not None and xmlnode is not None:
            selectedNodes = xmlnode.findall(selectStr)
            for node in selectedNodes:
                for child in xrtnode:
                    childCopy = copy.deepcopy(child)
                    obj = getChildObject(childCopy, node, jobInfo)
                    if obj is not None:
                        self.children.append(obj)

    def getPictures(self):
        from ccp4i2.report.pictures import Picture
        rv = []
        for child in self.children:
            if isinstance(child, Picture):
                rv.append(child.picDefFile.__str__())
            elif isinstance(child, (Container, Loop)):
                rv.extend(child.getPictures())
        return rv

# Report class - container for ccp4 html report


class Report(Container):
    TASKNAME = ''
    RUNNING = False
    WATCHED_FILE = None
    FAILED = False
    USEPROGRAMXML = True
    SEPARATEDATA = False
    elementCount = 1
    ERROR_CODES = {101: {'description': 'Import xrt file does not exist'},
                   102: {'description': 'Failed to read import xrt file'},
                   103: {'description': 'Failed to find insert xrt element in program file'},
                   104: {'description': 'Failed to find insert element in imported xrt file'},
                   105: {'description': 'Failed loading job info for jobId'},
                   106: {'description': 'Failed loading XML file'},
                   107: {'description': 'Failed creating csv file'},
                   108: {'description': 'Failed creating xml data file'}, }

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        Container.__init__(
            self,
            xrtnode=xrtnode,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.tag = 'report'
        self.pictureQueue = None
        self.standardise = kw.get('standardise', False)
        self.jobNumber = kw.get('jobNumber', None)
        self.projectId = kw.get('projectId', None)
        self.jobStatus = kw.get('jobStatus', None)
        self.cssFile = kw.get('cssFile', None)
        self.jsFile = kw.get('jsFile', None)
        self.additionalCssFiles = kw.get('additionalCssFiles', None)
        self.additionalJsFiles = kw.get('additionalJsFiles', None)
        self.additionalScript = kw.get('additionalScript', None)
        self.requireDataMain = kw.get('requireDataMain', 'mosflmApp')
        self.htmlBase = kw.get('htmlBase', None)
        self.cssVersion = kw.get('cssVersion', None)
        if 'jobId' in kw:
            try:
                from ccp4i2.report import CCP4ReportGenerator
                jobInfo = CCP4ReportGenerator.getReportJobInfo(kw['jobId'])
                self.jobInfo.update(jobInfo)
            except BaseException:
                self.errReport.append(self.__class__, 105, kw['jobId'])
                return
        if xmlnode is None and 'xmlFile' in kw:
            try:
                text = open(kw['xmlFile']).read()
                self.xmlnode = etree.fromstring(text, etree.XMLParser(encoding='utf-8'))
            except Exception as e:
                self.errReport.append(
                    self.__class__, 106, 'Reading file: ' + str(kw['xmlFile']) + '\n' + str(e))
                return
        if xrtnode is not None:
            xrtnode = self.removeComments(xrtnode)
            if xmlnode is not None:
                xrtnode = self.makeInserts(xrtnode, xmlnode)
            if self.standardise and str(xrtnode.tag) == 'report':
                xrtnode = self.standardiseXrtReport(xrtnode)
            Container.interpretXrt(self, xrtnode=xrtnode)

    def getJobFolder(self):
        return self.jobInfo.get('fileroot')

    def removeComments(self, xrtnode):
        for ele in xrtnode.iterfind('.//' + XRTNS + 'comment'):
            ele.getparent().remove(ele)
        return xrtnode

    def makeInserts(self, xrtnode, xmlnode):
        for insertEle in xrtnode.iterfind('.//' + XRTNS + 'insertXrt'):
            ele = None
            if insertEle.get('filename') is not None:
                from ccp4i2.core import CCP4Utils
                fileName = insertEle.get('filename')
                if fileName[0:8] == '$CCP4I2/':
                    fileName = os.path.join(
                        CCP4Utils.getCCP4I2Dir(), fileName[8:])
                if not os.path.exists(fileName):
                    self.errReport.append(self.__class__, 101, fileName)
                else:
                    try:
                        ele = etree.fromstring(open(fileName).read(), etree.XMLParser(encoding='utf-8'))
                    except BaseException:
                        self.errReport.append(self.__class__, 102, fileName)
                    else:
                        # Look for sub-element of the inserted file
                        if insertEle.get('select') is not None:
                            path = insertEle.get('select')
                            ele0 = ele.find(path)
                            if ele0 is None:
                                self.errReport.append(
                                    self.__class__, 104, 'file: ' + str(fileName) + ' findall: ' + str(path))
                            else:
                                ele = ele0
            # This is looking in the program output for xrt - dubious.
            elif insertEle.get('select') is not None:
                path = insertEle.get('select')
                ele = xmlnode.find(path)
                if ele is None:
                    self.errReport.append(self.__class__, 103, path)
            if ele is not None:
                for child in ele.getchildren():
                    insertEle.addnext(child)
                insertEle.getparent().remove(insertEle)
        return xrtnode

    def standardisePythonReport(self):
        from ccp4i2.report.metadata import Title, JobDetails, JobLogFiles, ReferenceGroup, Reference
        from ccp4i2.report.io_data import InputData, OutputData
        from ccp4i2.report.elements import Fold
        self.children.insert(0, Title(jobInfo=self.jobInfo))
        for cls in (InputData, OutputData, JobDetails, JobLogFiles):
            add = True
            for child in self.children:
                if isinstance(child, cls):
                    add = False
                    break
            if add:
                self.children.append(cls(jobInfo=self.jobInfo))
        # If there is not list of references try adding the default for the
        # task
        fold = Fold(
            label='Bibliographic references',
            brief='Biblio',
            jobInfo=self.jobInfo)
        referenceGroup = ReferenceGroup()
        referenceGroup.loadFromMedLine(self.TASKNAME)
        fold.children.append(referenceGroup)
        if isinstance(self.children[-1], JobDetails):
            self.children.insert(-1, fold)
        else:
            self.children.append(fold)

    def standardiseXrtReport(self, report):
        if report.find(XRTNS + 'inputdata') is None:
            inEle = etree.Element(XRTNS + 'inputdata')
            report.append(inEle)
        if report.find(XRTNS + 'outputdata') is None:
            outEle = etree.Element(XRTNS + 'outputdata')
            report.append(outEle)
        if report.find(XRTNS + 'jobdetails') is None:
            jobEle = etree.Element(XRTNS + 'jobdetails')
            report.append(jobEle)
        return report

    def as_data_etree(self):
        if self.standardise:
            self.standardisePythonReport()
        root = super().as_data_etree()
        return root

    def getReport(self):
        return self

    def showWarnings(self):
        """
        Deal with any Warnings blocks created by CPluginScript.createWarningsXML
        """
        xmlnode = self.xmlnode

        warningsNodes = xmlnode.findall('Warnings')
        if len(warningsNodes) > 0:

            haveWarnings = False
            for wsN in warningsNodes:
                logFileNodes = wsN.findall("logFile")
                for lfN in logFileNodes:
                    fileWarningsNode = lfN.findall("warning")
                    if len(fileWarningsNode) > 0:
                        haveWarnings = True
                        break
                if haveWarnings:
                    break

            if haveWarnings:
                warningsFold = self.addFold(
                    label='Warnings', initiallyOpen=False, brief='Warnings')
                topDiv = warningsFold.addDiv(
                    style='border:0px solid blue; overflow:auto;')

                for wsN in warningsNodes:
                    logFileNodes = wsN.findall("logFile")
                    for lfN in logFileNodes:
                        fileWarningsNode = lfN.findall("warning")
                        if len(fileWarningsNode) > 0:
                            logFileDiv = topDiv.addDiv(
                                style='border:0px solid black; position: relative; float: left;')
                            logFileDiv.append(
                                "<h4>" + lfN.findall("fileName")[0].text + "</h4>")
                            for wN in fileWarningsNode:
                                textNode = wN.findall("text")[0]
                                typeNode = wN.findall("type")[0]
                                if typeNode.text == "ERROR":
                                    logFileDiv.addDiv(
                                        style="color:red;")
                                    logFileDiv.append(
                                        "<pre>" + textNode.text + "</pre>")
                                elif typeNode.text == "WARNING":
                                    logFileDiv.addDiv()
                                    logFileDiv.append(
                                        "<pre>" + textNode.text + "</pre>")
                                elif typeNode.text == "ATTENTION":
                                    logFileDiv.addDiv()
                                    logFileDiv.append(
                                        "<pre>" + textNode.text + "</pre>")
