"""
Core report classes and utilities.

This module contains the foundation classes that all report elements
build upon: ReportClass (base), Container (children management),
Report (top-level), and supporting utilities.
"""

from __future__ import annotations

import logging
import xml.etree.ElementTree as etree
from typing import Any, TYPE_CHECKING

from ccp4i2.core.CCP4ErrorHandling import SEVERITY_WARNING, CErrorReport, CException
from ccp4i2 import I2_TOP

if TYPE_CHECKING:
    from ccp4i2.report.grid import GridContainer, GridItem

# Import error handling (lazy to avoid circular imports)
_diagnostics_module = None


def _get_diagnostics() -> Any:
    """Lazy import of diagnostics module to avoid circular imports."""
    global _diagnostics_module
    if _diagnostics_module is None:
        from ccp4i2.report import errors as _diagnostics_module
    return _diagnostics_module


logger = logging.getLogger(f"ccp4i2:{__name__}")

CCP4NS: str = "http://www.ccp4.ac.uk/ccp4ns"
XHTMLNS: str = "{http://www.w3.org/1999/xhtml}"
CURRENT_CSS_VERSION: str = '0.1.0'



def htmlBase() -> str:
    """Return base path for report static files (images, CSS, etc.).

    Static files are served from Next.js public directory.
    Path: /report_files/{version}/
    """
    return '/report_files/' + CURRENT_CSS_VERSION


def toBoolean(text: str) -> bool:
    """Convert a string to a boolean, trying integer conversion first."""
    try:
        i = int(text)
        return bool(i)
    except BaseException:
        return bool(text)


def findChildren(obj: Any, cls: type, id: str | None = None) -> list:
    """Recursively find all children of a given class."""
    retList: list = []
    if not hasattr(obj, 'children'):
        return retList
    for child in obj.children:
        if isinstance(child, cls):
            retList.append(child)
        elif isinstance(child, Container):
            retList.extend(findChildren(child, cls, id))
    return retList


class ReportClass(object):
    """Base class for all report elements.

    Provides common attributes and methods for report elements including
    ID generation, parent traversal, and XML serialization via
    ``as_data_etree()``.

    Attributes:
        id: Unique identifier for this element
        label: Human-readable label
        title: Title/tooltip text
        class_: CSS class(es) for styling
        style: Inline CSS styles
        parent: Parent element in the hierarchy
        internalId: Auto-generated internal ID for XML output
    """

    def __init__(self, *arg: Any, **kw: Any) -> None:
        super(ReportClass, self).__init__()
        self.id: str | None = kw.get('id', None)
        self.label: str | None = kw.get('label', None)
        self.title: str | None = kw.get('title', None)
        self.class_: str | None = kw.get('class_', None)
        self.style: str | None = kw.get('style', None)
        self.parent: ReportClass | None = kw.get('parent', None)
        self.outputXml: bool = kw.get('outputXml', False)
        self.internalId: str | None = kw.get('internalId', None)

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
    def diagnostics(self) -> Any:
        """Get or create the diagnostics collector for this element."""
        if self._diagnostics is None:
            errors = _get_diagnostics()
            self._diagnostics = errors.DiagnosticCollector()
        return self._diagnostics

    def add_diagnostic(self, level: str, code: str, message: str, **kwargs: Any) -> Any:
        """Add a diagnostic to this element.

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

    def warning(self, code: str, message: str, **kwargs: Any) -> Any:
        """Add a warning diagnostic."""
        return self.add_diagnostic('warning', code, message, **kwargs)

    def error(self, code: str, message: str, **kwargs: Any) -> Any:
        """Add an error diagnostic."""
        return self.add_diagnostic('error', code, message, **kwargs)

    def getReport(self) -> Report | None:
        """Traverse up the parent hierarchy to find the root Report."""
        if hasattr(self, 'parent') and self.parent is not None:
            return self.parent.getReport()
        return None

    # Modern alias
    def get_report(self) -> Report | None:
        """Modern alias for getReport()."""
        return self.getReport()

    def data_id(self) -> str:
        """Get the data ID for this element (used in XML output)."""
        return 'data_' + self.internalId

    def data_url(self) -> str:
        """Get the data URL for this element."""
        return './tables_as_xml_files/' + self.data_id() + '.xml'

    def as_data_etree(self) -> etree.Element:
        """Generate XML element for frontend consumption.

        Returns an ``etree.Element`` with tag ``CCP4i2Report{ClassName}``.
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


class Container(ReportClass):
    """Base class for report elements that can contain children.

    Provides child management (append, insert), the ``addX()`` convenience
    methods for building reports programmatically, and recursive
    ``as_data_etree()`` serialization.
    """

    tag: str = 'container'
    ERROR_CODES: dict = {}

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        super(Container, self).__init__(
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        if getattr(self, 'errReport', None) is None:
            self.errReport = CException()
        self.children: list[ReportClass] = []
        self.jobInfo: dict[str, Any] = {}
        self.jobInfo.update(jobInfo)
        self.xmlnode = xmlnode
        self.text: str = kw.get('text', '') or ''
        self.tail: str = kw.get('tail', '') or ''
        self._scriptErrorReport = None

    def append(self, child: ReportClass | str) -> None:
        """Append a child element. Strings are wrapped in Generic."""
        if isinstance(child, str):
            from ccp4i2.report.elements import Generic
            self.children.append(Generic(xmlnode=self.xmlnode, text=child))
        else:
            self.children.append(child)

    def insert(self, indx: int, child: ReportClass | str) -> None:
        """Insert a child element at position. Strings are wrapped in Generic."""
        if isinstance(child, str):
            from ccp4i2.report.elements import Generic
            self.children.insert(
                indx, Generic(
                    xmlnode=self.xmlnode, text=child))
        else:
            self.children.insert(indx, child)

    def __len__(self) -> int:
        return len(self.children)

    def addObjectOfClass(
        self,
        classOfObject: type,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> Any:
        """Create an instance of classOfObject and append it as a child."""
        if xmlnode is None:
            xmlnode = self.xmlnode
        if jobInfo is None:
            jobInfo = self.jobInfo
        obj = classOfObject(
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            parent=self,
            **kw)
        self.children.append(obj)
        return obj

    def addText(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add plain text to the report."""
        from ccp4i2.report.elements import Text
        return self.addObjectOfClass(Text, xmlnode, jobInfo, **kw)

    def addPre(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add a preformatted text block."""
        from ccp4i2.report.elements import Pre
        return self.addObjectOfClass(Pre, xmlnode, jobInfo, **kw)

    def addAlignment(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add a sequence alignment viewer (ClustalW format)."""
        from ccp4i2.report.elements import Alignment
        return self.addObjectOfClass(Alignment, xmlnode, jobInfo, **kw)

    def addFetchPre(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add a preformatted text block (fetched from file)."""
        from ccp4i2.report.elements import FetchPre
        return self.addObjectOfClass(FetchPre, xmlnode, jobInfo, **kw)

    def addGraph(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add a data graph."""
        from ccp4i2.report.graphs import Graph
        return self.addObjectOfClass(Graph, xmlnode, jobInfo, **kw)

    def addProgress(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add a progress bar."""
        from ccp4i2.report.elements import Progress
        return self.addObjectOfClass(Progress, xmlnode, jobInfo, **kw)

    def addFlotGraph(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add a Flot/Plotly graph."""
        from ccp4i2.report.graphs import FlotGraph
        return self.addObjectOfClass(FlotGraph, xmlnode, jobInfo, **kw)

    def addGraphGroup(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add a group of graphs."""
        from ccp4i2.report.graphs import GraphGroup
        return self.addObjectOfClass(GraphGroup, xmlnode, jobInfo, **kw)

    def addFlotGraphGroup(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add a group of Flot graphs."""
        from ccp4i2.report.graphs import FlotGraphGroup
        return self.addObjectOfClass(FlotGraphGroup, xmlnode, jobInfo, **kw)

    def addDrawnDiv(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add a JavaScript-rendered div."""
        from ccp4i2.report.graphs import DrawnDiv
        return self.addObjectOfClass(DrawnDiv, xmlnode, jobInfo, **kw)

    def addObjectGallery(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add an object gallery."""
        from ccp4i2.report.graphs import ObjectGallery
        return self.addObjectOfClass(ObjectGallery, xmlnode, jobInfo, **kw)

    def addGraphLineChooser(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add a graph line chooser (table + graph)."""
        from ccp4i2.report.graphs import GraphLineChooser
        return self.addObjectOfClass(GraphLineChooser, xmlnode, jobInfo, **kw)

    def addPictureGroup(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add a group of pictures."""
        from ccp4i2.report.graphs import PictureGroup
        return self.addObjectOfClass(PictureGroup, xmlnode, jobInfo, **kw)

    def addDiv(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add a div container."""
        from ccp4i2.report.elements import Div
        return self.addObjectOfClass(Div, xmlnode, jobInfo, **kw)

    def addTable(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add a data table."""
        from ccp4i2.report.elements import Table
        return self.addObjectOfClass(Table, xmlnode, jobInfo, **kw)

    def addPicture(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add a picture/scene visualization."""
        from ccp4i2.report.pictures import Picture
        return self.addObjectOfClass(Picture, xmlnode, jobInfo, **kw)

    def addTitle(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add a job title element."""
        from ccp4i2.report.metadata import Title
        return self.addObjectOfClass(Title, xmlnode, jobInfo, **kw)

    def addLaunch(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add a viewer launch button."""
        from ccp4i2.report.actions import Launch
        return self.addObjectOfClass(Launch, xmlnode, jobInfo, **kw)

    def addDownload(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add a download button."""
        from ccp4i2.report.actions import Download
        return self.addObjectOfClass(Download, xmlnode, jobInfo, **kw)

    def addFileLink(self, label: str, relativePath: str, fileType: str = 'text', **kw: Any) -> Any:
        """Add a clickable link to a file in the job directory.

        Args:
            label: Display text (e.g. "Show Pointless logfile")
            relativePath: Path relative to job dir (e.g. "job_1/log.txt")
            fileType: "text" for Monaco preview, "html" for browser tab
        """
        from ccp4i2.report.actions import FileLink
        obj = FileLink(
            jobInfo=self.jobInfo,
            label=label,
            relativePath=relativePath,
            fileType=fileType,
            projectId=kw.get('projectId', self.jobInfo.get('project_pk')),
            **{k: v for k, v in kw.items() if k != 'projectId'},
        )
        self.children.append(obj)
        return obj

    def addCopyToClipboard(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        text: str = "",
        label: str = "",
        **kw: Any,
    ) -> Any:
        """Add a copy-to-clipboard button."""
        from ccp4i2.report.actions import CopyToClipboard
        return self.addObjectOfClass(
            CopyToClipboard, text=text, label=label, **kw)

    def addCopyUrlToClipboard(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        text: str = "",
        label: str = "",
        **kw: Any,
    ) -> Any:
        """Add a copy-URL-to-clipboard button."""
        from ccp4i2.report.actions import CopyUrlToClipboard
        return self.addObjectOfClass(
            CopyUrlToClipboard, text=text, label=label, **kw)

    def addDAGGraph(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add a directed acyclic graph visualization."""
        from ccp4i2.report.graphs import DAGGraph
        return self.addObjectOfClass(DAGGraph, xmlnode, jobInfo, **kw)

    def addResults(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add a results container."""
        from ccp4i2.report.elements import Results
        return self.addObjectOfClass(Results, xmlnode, jobInfo, **kw)

    def addFold(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add a collapsible fold section."""
        from ccp4i2.report.elements import Fold
        return self.addObjectOfClass(Fold, xmlnode, jobInfo, **kw)

    def addCopy(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add a copy element (copies XML content into the report)."""
        from ccp4i2.report.elements import Copy
        return self.addObjectOfClass(Copy, xmlnode, jobInfo, **kw)

    def addHelp(self, xmlnode: etree.Element | None = None, jobInfo: dict[str, Any] | None = None, **kw: Any) -> Any:
        """Add a help button."""
        from ccp4i2.report.actions import Help
        return self.addObjectOfClass(Help, xmlnode, jobInfo, **kw)

    # Modern grid layout methods

    def addGrid(self, spacing: int = 2, direction: str = 'row', **kw: Any) -> GridContainer:
        """Add a grid container for responsive layouts.

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

    def addGridItem(self, xs: int = 12, sm: int | None = None, md: int | None = None, lg: int | None = None, xl: int | None = None, **kw: Any) -> GridItem:
        """Add a grid item directly to this container.

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

    def addTwoColumnLayout(self, left_span: int = 8, right_span: int = 4, spacing: int = 2, **kw: Any) -> tuple[GridItem, GridItem]:
        """Add a two-column layout (convenience method).

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

    def getPictures(self) -> list[str]:
        """Recursively collect all picture file paths from children."""
        from ccp4i2.report.pictures import Picture
        rv: list[str] = []
        for child in self.children:
            if isinstance(child, Picture):
                rv.append(child.picDefFile.__str__())
            elif isinstance(child, Container):
                rv.extend(child.getPictures())
        return rv

    def as_data_etree(self) -> etree.Element:
        root = super().as_data_etree()
        for child in self.children:
            if hasattr(child, 'as_data_etree'):
                root.append(child.as_data_etree())
        return root

    def errorReport(self) -> CException:
        """Collect error reports from this container and all children."""
        err = CException()
        err.extend(self.errReport)
        for child in self.children:
            if getattr(child, 'errReport', None) is not None:
                err.extend(child.errorReport())
        return err

    def loadXmlFile(self, fileName: str) -> etree.Element:
        """Load and parse an XML file."""
        from ccp4i2.core import CCP4Utils
        text = CCP4Utils.readFile(fileName)
        ele = etree.fromstring(text)
        return ele


class Report(Container):
    """Top-level report container.

    Handles report standardization (adding Title, InputData, OutputData,
    JobDetails, references).
    """

    TASKNAME: str = ''
    RUNNING: bool = False
    WATCHED_FILE: str | None = None
    FAILED: bool = False
    USEPROGRAMXML: bool = True
    SEPARATEDATA: bool = False
    elementCount: int = 1
    ERROR_CODES: dict = {105: {'description': 'Failed loading job info for jobId'},
                   106: {'description': 'Failed loading XML file'},
                   107: {'description': 'Failed creating csv file'},
                   108: {'description': 'Failed creating xml data file'}, }

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        Container.__init__(
            self,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.tag: str = 'report'
        self.pictureQueue = None
        self.standardise: bool = kw.get('standardise', False)
        self.jobNumber: int | None = kw.get('jobNumber', None)
        self.projectId: str | None = kw.get('projectId', None)
        self.jobStatus: str | None = kw.get('jobStatus', None)
        self.cssFile: str | None = kw.get('cssFile', None)
        self.jsFile: str | None = kw.get('jsFile', None)
        self.additionalCssFiles = kw.get('additionalCssFiles', None)
        self.additionalJsFiles = kw.get('additionalJsFiles', None)
        self.additionalScript: str | None = kw.get('additionalScript', None)
        self.requireDataMain: str = kw.get('requireDataMain', 'mosflmApp')
        self.htmlBase: str | None = kw.get('htmlBase', None)
        self.cssVersion: str | None = kw.get('cssVersion', None)
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

    def getJobFolder(self) -> str | None:
        """Return the job's file root directory."""
        return self.jobInfo.get('fileroot')

    def standardisePythonReport(self) -> None:
        """Add standard sections (Title, I/O data, references, JobDetails) to a Python-built report."""
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

    def as_data_etree(self) -> etree.Element:
        if self.standardise:
            self.standardisePythonReport()
        root = super().as_data_etree()
        return root

    def getReport(self) -> Report:
        return self

    def showWarnings(self) -> None:
        """Process Warnings blocks created by CPluginScript.createWarningsXML."""
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
