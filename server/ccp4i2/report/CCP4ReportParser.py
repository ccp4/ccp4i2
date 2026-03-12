"""
CCP4i2 Report Parser and Element Classes.

This module provides the core report element classes for generating
XML reports from CCP4 job output. Reports are serialized to XML and
consumed by the React frontend for display.

Architecture:
    Report (Container)
        ├── Title
        ├── Results (Container)
        │   ├── Fold (Container)
        │   │   ├── Table
        │   │   ├── Graph/FlotGraph
        │   │   └── Text/Pre
        │   └── ...
        ├── InputData
        ├── OutputData
        └── ReferenceGroup

Key Classes:
    Report - Main report container, subclassed by plugin reports
    Container - Base for elements that contain children
    Fold - Collapsible section
    Table - Data table
    Graph/FlotGraph - Chart/graph visualization
    Text/Pre - Text content

Backward Compatibility:
    This module maintains full backward compatibility with existing
    plugin report classes in wrappers/, pipelines/, etc. Legacy code
    that subclasses Report continues to work unchanged.

Modern Usage:
    New reports can use the grid layout system for responsive layouts:

        from ccp4i2.report import GridContainer, GridItem, GridSpan

        grid = GridContainer(spacing=2)
        left = grid.item(span=GridSpan(xs=12, md=8))
        left.append(self.addTable(...))
        self.append(grid)
"""

import logging
import os
import re
import sys
import xml.etree.ElementTree as etree
from io import StringIO

from lxml import html as lxml_html

from ccp4i2.core.CCP4ErrorHandling import SEVERITY_OK, SEVERITY_WARNING, CErrorReport, CException
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
# From http://www.w3.org/QA/2002/04/valid-dtd-list.html
DOCTYPE = '''<?xml version="1.0"?>
    <!DOCTYPE html>
    <html xmlns:xsi="http://www.w3.org/1999/xhtml"></html>
    '''
# DOCTYPE = '''<?xml version="1.0" encoding="utf-8"?>\n<!DOCTYPE html>'''

XHTMLNS = "{http://www.w3.org/1999/xhtml}"
CURRENT_CSS_VERSION = '0.1.0'



def htmlBase():
    """Return base path for report static files (images, CSS, etc.).

    Static files are served from Next.js public directory.
    Path: /report_files/{version}/
    """
    return '/report_files/' + CURRENT_CSS_VERSION


# htmlDoc function removed - HTML generation no longer used
# The data path uses as_data_etree() which returns XML for frontend rendering


def testPathExists(self, relPath):
    return True


# getLastestMinorVersion removed - CSS version management no longer needed


def toBoolean(text):
    try:
        i = int(text)
        return bool(i)
    except BaseException:
        return bool(text)


def getChildObject(child, xmlnode, jobInfo={}, report=None):
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
    elif child.tag == XRTNS + "drilldown":
        obj = DrillDown(jobInfo)
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
        # print 'findChildren',child,
        # isinstance(child,cls),isinstance(child,Container)
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


def saveToFile(tree, fileName):
    text = etree.tostring(tree, xml_declaration=True)
    f = open(fileName, 'w')
    f.write(text)
    f.close()


def findallEval(srcexpr, xmlnode):
    """Evaluate a python expression containing findall elements deliminated by {}. The findall expressions are evaluated and return strings which are substituted into the expression. The expression is evaluated and the result returned."""
    tgtexpr = ""
    for i in range(srcexpr.count('{')):
        i0 = srcexpr.find('{')
        if i0 < 0:
            break
        i1 = srcexpr.find('}', i0)
        if i1 < 0:
            break
        tgtexpr += srcexpr[:i0]
        xp = srcexpr[i0 + 1:i1]
        nodes = xmlnode.findall(xp)
        tgtexpr += '"""'
        for node in nodes:
            tgtexpr += node.text
        tgtexpr += '"""'
        srcexpr = srcexpr[i1 + 1:]
    tgtexpr += srcexpr
    return eval(tgtexpr)


def PARSER():
    I2XmlParser.insts = etree.XMLParser(encoding='utf-8')
    return I2XmlParser.insts


class I2XmlParser:
    insts = None


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
                'description': 'Can not create Loop instance without xrtnode argument (it should not be used in Python mode)'}, 4: {
                    'description': 'Unable to create RTF file - unable to import PyQt4'}}

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
            self.children.append(Generic(xmlnode=self.xmlnode, text=child))
        else:
            self.children.append(child)

    def insert(self, indx, child):
        if isinstance(child, str):
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
        return self.addObjectOfClass(Text, xrtnode, xmlnode, jobInfo, **kw)

    def addPre(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        return self.addObjectOfClass(Pre, xrtnode, xmlnode, jobInfo, **kw)

    def addFetchPre(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        return self.addObjectOfClass(FetchPre, xrtnode, xmlnode, jobInfo, **kw)

    def addGraph(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        return self.addObjectOfClass(Graph, xrtnode, xmlnode, jobInfo, **kw)

    def addProgress(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        return self.addObjectOfClass(Progress, xrtnode, xmlnode, jobInfo, **kw)

    def addFlotGraph(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        return self.addObjectOfClass(
            FlotGraph, xrtnode, xmlnode, jobInfo, **kw)

    def addGraphGroup(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        return self.addObjectOfClass(
            GraphGroup, xrtnode, xmlnode, jobInfo, **kw)

    def addFlotGraphGroup(
            self,
            xrtnode=None,
            xmlnode=None,
            jobInfo=None,
            **kw):
        return self.addObjectOfClass(
            FlotGraphGroup, xrtnode, xmlnode, jobInfo, **kw)

    def addDrawnDiv(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        return self.addObjectOfClass(DrawnDiv, xrtnode, xmlnode, jobInfo, **kw)

    def addObjectGallery(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        return self.addObjectOfClass(
            ObjectGallery, xrtnode, xmlnode, jobInfo, **kw)

    def addGraphLineChooser(
            self,
            xrtnode=None,
            xmlnode=None,
            jobInfo=None,
            **kw):
        return self.addObjectOfClass(
            GraphLineChooser, xrtnode, xmlnode, jobInfo, **kw)

    def addPictureGroup(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        return self.addObjectOfClass(
            PictureGroup, xrtnode, xmlnode, jobInfo, **kw)

    def addDiv(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        return self.addObjectOfClass(Div, xrtnode, xmlnode, jobInfo, **kw)

    def addTable(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        return self.addObjectOfClass(Table, xrtnode, xmlnode, jobInfo, **kw)

    def addPicture(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        return self.addObjectOfClass(Picture, xrtnode, xmlnode, jobInfo, **kw)

    def addTitle(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        return self.addObjectOfClass(Title, xrtnode, xmlnode, jobInfo, **kw)

    def addLaunch(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        return self.addObjectOfClass(Launch, xrtnode, xmlnode, jobInfo, **kw)

    def addDownload(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        return self.addObjectOfClass(Download, xrtnode, xmlnode, jobInfo, **kw)

    def addCopyToClipboard(
            self,
            xrtnode=None,
            xmlnode=None,
            jobInfo=None,
            text="",
            label="",
            **kw):
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
        return self.addObjectOfClass(
            CopyUrlToClipboard, text=text, label=label, **kw)

    def addDAGGraph(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        return self.addObjectOfClass(DAGGraph, xrtnode, xmlnode, jobInfo, **kw)

    def addResults(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        return self.addObjectOfClass(Results, xrtnode, xmlnode, jobInfo, **kw)

    def addFold(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        return self.addObjectOfClass(Fold, xrtnode, xmlnode, jobInfo, **kw)

    def addCopy(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
        return self.addObjectOfClass(Copy, xrtnode, xmlnode, jobInfo, **kw)

    def addHelp(self, xrtnode=None, xmlnode=None, jobInfo=None, **kw):
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
        rv = []
        for child in self.children:
            if isinstance(child, Picture):
                rv.append(child.picDefFile.__str__())
            elif isinstance(child, (Container, Loop)):
                rv.extend(child.getPictures())
        # print 'Container.getPictures',rv
        return rv

    def as_data_etree(self):
        root = super().as_data_etree()
        for child in self.children:
            if hasattr(child, 'as_data_etree'):
                root.append(child.as_data_etree())
        return root

    def graph_data_as_rtf(self, fileName=None):
        """
        RTF/ODF export functionality - DEPRECATED.

        This method previously used Qt (PySide2) for RTF document generation.
        RTF export is no longer supported in the Qt-free version.
        Use data_as_csv() for data export instead.
        """
        import logging
        logger = logging.getLogger(f"ccp4i2:{__name__}")
        logger.warning(
            "graph_data_as_rtf() is deprecated - RTF export requires Qt which is no longer available")
        return None

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

    def scriptErrorReport(self):
        from ccp4i2.core import CCP4ErrorHandling, CCP4File
        if self._scriptErrorReport is None:
            self._scriptErrorReport = CCP4ErrorHandling.CErrorReport()
            root = self.jobInfo.get('fileroot', None)
            print('scriptErrorReport root', root)
            if root is not None:
                filePath = os.path.join(root, 'diagnostic.xml')
                if os.path.exists(filePath):
                    i2xml = CCP4File.CI2XmlDataFile(filePath)
                    self._scriptErrorReport.setEtree(i2xml.getBodyEtree())
        return self._scriptErrorReport


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
        # print 'ReportParser IfContainer',xmlnode.findall( xrtnode.get(
        # "select" ) ),type(xmlnode.findall( xrtnode.get( "select" ) ))
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
        # print 'ReportParser IfContainer from program
        # output:',value,'target:',target
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
        # print 'Loop.init',selectStr
        if selectStr is not None and xmlnode is not None:
            selectedNodes = xmlnode.findall(selectStr)
            # print 'Loop.init selectedNodes',selectedNodes
            for node in selectedNodes:
                for child in xrtnode:
                    childCopy = copy.deepcopy(child)
                    obj = getChildObject(childCopy, node, jobInfo)
                    if obj is not None:
                        self.children.append(obj)
                    # self.children[-1].text = xrtnode.text
                    # self.children[-1].tail = xrtnode.tail

    def getPictures(self):
        rv = []
        for child in self.children:
            if isinstance(child, Picture):
                rv.append(child.picDefFile.__str__())
            elif isinstance(child, (Container, Loop)):
                rv.extend(child.getPictures())
        # print 'Loop.getPictures',rv
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
    ERROR_CODES = {0: {'severity': SEVERITY_OK, 'description': 'OK'},
                   101: {'description': 'Import xrt file does not exist'},
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
                self.xmlnode = etree.fromstring(text, PARSER())
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
        # self._makeMgPicture = None

    def getJobFolder(self):
        return self.jobInfo.get('fileroot')

    def removeComments(self, xrtnode):
        for ele in xrtnode.iterfind('.//' + XRTNS + 'comment'):
            # print 'Removing comment',ele.text
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
                # print 'Report.makeInserts fileName',fileName
                if not os.path.exists(fileName):
                    self.errReport.append(self.__class__, 101, fileName)
                else:
                    try:
                        ele = etree.fromstring(open(fileName).read(), PARSER())
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
            # print etree.tostring(xrtnode,pretty_print=True)
        return xrtnode

    def containsPictures(self):
        # print 'Report.containsPictures',self._makeMgPicture
        # if self._makeMgPicture is None or
        # len(self._makeMgPicture.pictureQueue)==0:
        if self.pictureQueue is None:
            self.pictureQueue = self.getPictures()
        if len(self.pictureQueue) == 0:
            return False
        else:
            return True

    def makeMgPicture(self):
        if self._makeMgPicture is None:
            import CCP4MakeMgPicture
            self._makeMgPicture = CCP4MakeMgPicture.CMakeMgPicture(
                jobInfo=self.jobInfo)
        return self._makeMgPicture

    def standardisePythonReport(self):
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
        # if report.find(XRTNS+'title') is None:
        #  titleEle = etree.Element(XRTNS+'title')
        #  report.insert(0,titleEle)
        if report.find(XRTNS + 'inputdata') is None:
            inEle = etree.Element(XRTNS + 'inputdata')
            report.append(inEle)
        if report.find(XRTNS + 'outputdata') is None:
            outEle = etree.Element(XRTNS + 'outputdata')
            report.append(outEle)
        # if report.find(XRTNS+'drilldown') is None:
        #  jobEle = etree.Element(XRTNS+'drilldown')
        #  report.append(jobEle)
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
                                    clearDiv = logFileDiv.addDiv(
                                        style="color:red;")
                                    logFileDiv.append(
                                        "<pre>" + textNode.text + "</pre>")
                                elif typeNode.text == "WARNING":
                                    clearDiv = logFileDiv.addDiv()
                                    logFileDiv.append(
                                        "<pre>" + textNode.text + "</pre>")
                                elif typeNode.text == "ATTENTION":
                                    clearDiv = logFileDiv.addDiv()
                                    logFileDiv.append(
                                        "<pre>" + textNode.text + "</pre>")


class Results(Container):
    """Results container - uses Container.as_data_etree() for data path."""
    pass


class DrillDown(ReportClass):

    def __init__(self, jobInfo, **kw):
        super(DrillDown, self).__init__(**kw)
        self.jobInfo = jobInfo
        # print 'DrillDown jobInfo',jobInfo

    def getSubJobs(self, jobDir):
        import glob

        from ccp4i2.core import CCP4File
        if jobDir is None:
            return []
        globDirs = glob.glob(os.path.join(jobDir, 'job_*'))

        def sortSubJobsKey(inp):
            return int(inp.split('_')[-1])
        sortedDirs = sorted(globDirs, key=sortSubJobsKey)
        retSubJobs = []
        for path in sortedDirs:
            subJobs = self.getSubJobs(path)
            parsFile = glob.glob(os.path.join(path, '*.params.xml'))
            if len(parsFile) > 0:
                # Get taskName from file
                f = CCP4File.CI2XmlDataFile(fullPath=parsFile[0])
                pluginName = str(f.header.get('pluginName'))
                jobId = int(f.header.get('jobId'))
            else:
                pluginName = None
                jobId = None
            '''
      reportFile =  glob.glob(os.path.join(path,'*.report.html'))
      if len(reportFile) == 0:
        reportFile = None
      else:
        reportFile = reportFile[0]
      '''
            reportFile = parsFile[0][0:-10] + 'report.html'

            retSubJobs.append([os.path.split(path)[-1].split('_')
                              [-1], jobId, pluginName, reportFile, subJobs])

        return retSubJobs


def foldTitleLine(label, initiallyOpen, brief=None):
    root = etree.Element('root')
    anchor = etree.Element('a')
    anchor.set('name', label)
    root.append(anchor)
    span = etree.Element('span')
    span.set('class', 'folder')
    if brief is not None:
        span.set('title', brief)
    span.set('onclick', "toggleview(this)")
    if initiallyOpen:
        span.text = u"\u25BC" + " " + label
    else:
        span.text = u"\u25B6" + " " + label
    root.append(span)
    div = etree.Element('div')
    if initiallyOpen:
        div.set('class', 'hidesection displayed')
    else:
        div.set('class', 'hidesection hidden')
    root.append(div)
    return root


def foldLinkLine(label, href, id):
    root = etree.Element('root')
    span = etree.Element('span')
    span.set('class', 'folder')
    span.set('onclick', "togglesubjob(this)")
    span.text = 'Show ' + label
    root.append(span)
    div = etree.Element('div')
    div.set('id', id)
    div.set('class', 'subjob')
    root.append(div)
    return root

# Fold class - container for a hidden report within a report


class Fold(Container):
    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super(
            Fold,
            self).__init__(
            xrtnode=xrtnode,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        if xrtnode is not None:
            self.label = xrtnode.get('label', 'Show details')
            self.brief = xrtnode.get('brief', None)
        else:
            self.label = kw.get('label', 'Fold')
            self.brief = kw.get('brief', None)
        self.initiallyOpen = kw.get('initiallyOpen', False)

    def as_data_etree(self):
        if self.brief is None or len(self.brief) == 0:
            self.brief = self.label
        root = super().as_data_etree()
        root.set('label', self.label)
        root.set('initiallyOpen', str(self.initiallyOpen))
        root.set('brief', self.brief)
        return root


class Text(ReportClass):
    tag = 'span'

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super(
            Text,
            self).__init__(
            xrtnode=xrtnode,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.text = ""
        self.tail = ""
        # print 'Text.init',xrtnode.get( "if" ),xrtnode.get( "select"
        # ),xrtnode.text,xrtnode.tail
        if xrtnode is not None and len(xrtnode) > 0:
            if xrtnode.get("if") is not None and (
                xmlnode is None or not bool(
                    xmlnode.findall(
                        xrtnode.get("if")))):
                return
            if xrtnode.text:
                self.text += xrtnode.text
            # print 'Text.init self.text',self.text
            if xrtnode.tail:
                self.tail += xrtnode.tail
            if xrtnode.get("select") and xmlnode is not None:
                nodes = xmlnode.findall(xrtnode.get("select"))
                for node in nodes:
                    self.text += node.text
        elif kw.get('text', None) is not None:
            self.text += kw['text']
            # print 'Text.__init__',text
        elif kw.get('select', None) is not None and xmlnode is not None:
            nodes = xmlnode.findall(kw['select'])
            for node in nodes:
                self.text += node.text

    def as_data_etree(self):
        root = super().as_data_etree()
        return root

    def data_as_xml(self, fileName=None):
        # Tag as dummy so parent container will strip the tags

        if self.text is not None:
            root = etree.Element("div")
        else:
            root = etree.Element(self.tag)

        if self.style is not None:
            root.set('style', self.style)
        if self.class_ is not None:
            root.set('class', self.class_)
        try:
            if self.text is not None:
                # This is a hack to deal with libxml2 crash with "some" strings
                # on macOS. SJM 01/09/2023
                ts = self.text.split("\n")
                for t in ts:
                    el = etree.Element(self.tag)
                    if self.style is None:
                        el.set(
                            'style', "line-height: normal; margin: 0px 0px 0px 0px;")
                    el.text = t
                    root.append(el)
        except BaseException:
            print('Failed creating XML for:', self.text)
        try:
            if self.tail is not None:
                root.tail = self.tail
        except BaseException:
            print('Failed creating XML for:', self.tail)
        # print 'Text.as_etree',root.text,root.tail

        if fileName is not None:
            with open(fileName, 'w') as outputFile:
                outputFile.write(etree.tostring(root).decode())

        return root

    def getText(self):
        if len(self.text) == 0 or len(
                self.tail) == 0 or self.text[-1] == ' ' or self.tail[0] == ' ':
            return self.text + self.tail
        else:
            return self.text + ' ' + self.tail


class Pre(Text):
    tag = 'pre'


class FetchPre(Text):
    tag = 'pre'

    def data_as_xml(self, fileName=None):
        dataXml = Text.data_as_xml(self, fileName)
        if (hasattr(dataXml, "findall")):
            pres = dataXml.findall("pre")
            if len(pres) > 0:
                thePre = pres[0]
                theClass = thePre.get("class")
                if not theClass:
                    thePre.set("class", "urlfetchtag")
                else:
                    thePre.set("class", theClass + " urlfetchtag")
        return dataXml

# Status class


class Status(Container):
    """Status container - uses Container.as_data_etree() for data path."""
    pass


class Copy(ReportClass):
    def __init__(self, xrtnode=None, xmlnode=None, **kw):
        super(Copy, self).__init__(xrtnode=xrtnode, xmlnode=xmlnode, **kw)
        self.text = ""
        self.root = etree.Element('root')
        if xmlnode is not None:
            if xrtnode is not None:
                self.root.extend(xmlnode.findall(xrtnode.get("select")))
            elif kw.get('select', False):
                self.root.extend(xmlnode.findall(kw.get('select', False)))

        # print 'Copy __init__ root',self.root.tag


class Generic(ReportClass):

    ERROR_CODES = {1: {'description': 'Can not interpret text'}, 2: {
        'description': 'Error substituting values into generic item'}}

    def __init__(
            self,
            xrtnode=None,
            xmlnode=None,
            jobInfo={},
            text=None,
            defaultTag='p',
            **kw):
        super().__init__()
        self.id = kw.get('id', None)
        self.class_ = kw.get('class_', None)
        self.xmltree = None

        if xrtnode is None and text is not None:
            try:
                xrtnode = etree.fromstring(text.encode('utf-8'))
            except BaseException:
                try:
                    if isinstance(text, bytes):
                        xrtnode = etree.fromstring(
                            '<' + defaultTag + '>' + text.encode('utf-8') + '</' + defaultTag + '>')
                    else:
                        xrtnode = etree.fromstring(
                            '<' + defaultTag + '>' + text + '</' + defaultTag + '>')
                except BaseException:
                    raise
                    raise CException(self.__class__, 1, str(text))

        if xrtnode is not None:
            try:
                self.xmltree = applySelect(xrtnode, xmlnode, jobInfo)
            except BaseException:
                if text is not None:
                    self.errReport.append(self.__class__, 2, str(text))

    def as_data_etree(self):
        root = super().as_data_etree()
        root.append(self.xmltree)
        return root


class BaseTable(ReportClass):
    tableCount = 0

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super(BaseTable, self).__init__(xrtnode=xrtnode, xmlnode=xmlnode, **kw)
        BaseTable.tableCount += 1
        self.id = kw.get('id', 'table_' + str(BaseTable.tableCount))
        self.coldata = []
        self.coltitle = []
        self.colsubtitle = []
        self.colTips = []
        self.xmldata = []
        self.help = None
        self.download = None
        self.outputCsv = kw.get('outputCsv', True)
        downloadable = kw.get('downloadable', False)

        if downloadable:
            if self.title is not None:
                self.download = Download(
                    jobInfo=jobInfo, dataName=self.id + '_' + self.title)
            else:
                self.download = Download(jobInfo=jobInfo, dataName=self.id)
        if xrtnode is not None:
            self.excludeIfDataMissing = xrtnode.get(
                'excludeIfDataMissing', False)
            if xmlnode is not None:
                self.xmldata = xmlnode.findall(xrtnode.get("select"))
            self.transpose = (xrtnode.get("transpose") is not None)
            if xrtnode.get("help") is not None:
                self.help = Help(ref=xrtnode.get("help"))
        else:
            self.excludeIfDataMissing = kw.get('excludeIfDataMissing', False)
            if xmlnode is not None:
                if kw.get('select', None) is not None:
                    self.xmldata = []
                    tableEleList = []
                    for p in kw['select'].split("|"):
                        self.xmldata.extend(xmlnode.findall(p.strip()))
                elif kw.get('selectNodes', None) is not None:
                    self.xmldata = kw['selectNodes']
                else:
                    self.xmldata = [xmlnode]
            self.transpose = kw.get('transpose', False)
            if kw.get('help', None) is not None:
                self.help = Help(ref=kw['help'])
        # assemble column data

        if xrtnode is not None:
            for col in xrtnode:
                if col.tag == XRTNS + "data":
                    colttl1 = col.get("title")
                    colttl2 = col.get("subtitle")
                    colexpr = col.get("expr")
                    colsel = col.get("select")
                    self.addData(
                        xmldata=self.xmldata,
                        title=colttl1,
                        subtitle=colttl2,
                        expr=colexpr,
                        select=colsel)

    def addData(
            self,
            xmldata=None,
            title=None,
            subtitle=None,
            expr=None,
            function=None,
            select=None,
            data=[],
            tip=None):
        if tip is not None:
            self.colTips.append(tip)
        else:
            self.colTips += [None]
        colvals = []
        if xmldata is not None:
            if not isinstance(xmldata, list):
                xmldata = [xmldata]
        else:
            xmldata = self.xmldata
        if select is not None:         # select the values from the xml
            for x in xmldata:
                for selitem in x.findall(select):
                    if selitem.text is None:
                        val = '-'
                    else:
                        val = selitem.text.strip()
                        if expr is not None:  # allow an expression to be applied to the data
                            try:
                                val = eval(expr, {"x": float(val)})
                            except BaseException:
                                val = '-'
                        if function is not None:
                            val = function(val)
                    colvals.append(val)
        else:                      # allow a list of values to be given
            colvals.extend(data)
        if self.excludeIfDataMissing and len(colvals) == 0:
            pass
        else:
            self.coldata.append(colvals)
            maxlen = max([len(col) for col in self.coldata])
            for i in range(len(self.coldata)):
                while len(self.coldata[i]) < maxlen:
                    self.coldata[i].append(None)
            self.coltitle.append(title)
            self.colsubtitle.append(subtitle)

    def data_as_csv(self, fileName=None):
        import csv
        if len(self.coldata) == 0:
            return CErrorReport()
        try:
            f = open(fileName, 'w', newline='')
        except BaseException:
            return CErrorReport(Report, 107, str(fileName))
        try:
            writer = csv.writer(f)
            if self.transpose:
                iRow = 0
                for row in self.coldata:
                    if len(self.coltitle) > 0:
                        if iRow < len(self.coltitle):
                            outputRow = [self.coltitle[iRow]] + row
                        else:
                            outputRow = [" "] + row
                    else:
                        outputRow = row
                    writer.writerow(outputRow)
                    iRow += 1
            else:
                if len(self.coltitle) > 0:
                    writer.writerow(self.coltitle)
                nc = len(self.coldata)
                nr = len(self.coldata[0])
                for n in range(nr):
                    row = []
                    for col in self.coldata:
                        row.append(col[n])
                    writer.writerow(row)
        except BaseException:
            return CErrorReport(Report, 107, str(fileName))
        try:
            f.close()
        except BaseException:
            return CErrorReport(Report, 107, str(fileName))
        return CErrorReport()

    def data_as_xml(self, fileName=None):
        import csv
        if len(self.coldata) == 0:
            return CErrorReport()
        try:
            f = open(fileName, 'wb')
        except BaseException:
            return CErrorReport(Report, 107, str(fileName))
        try:
            dataAsEtree = self.data_as_etree()
            # As it comes out here, the data tag is in ccp4 namespace (i.e.
            # ccp4_data:data)
            rootNode = etree.Element('CCP4ApplicationOutput')
            for childTable in dataAsEtree:
                childTable.tag = 'CCP4Table'
                rootNode.append(childTable)
            dataAsText = etree.tostring(rootNode)
            f.write(dataAsText)
        except BaseException:
            print('Failed to select data')
            return CErrorReport(Report, 107, str(fileName))
        try:
            f.close()
        except BaseException:
            return CErrorReport(Report, 107, str(fileName))
        return CErrorReport()

# JavaScript-Enhances Table class


def _set_cell_content(cell_element, content):
    """Set table cell content, preserving inline HTML if present.

    If content contains valid XML/HTML markup (e.g. <i>, <b>, <br/>),
    it is embedded as child elements. Otherwise it is set as plain text,
    letting ElementTree handle entity encoding automatically.
    """
    if not content:
        cell_element.text = content or ""
        return
    try:
        parsed = etree.fromstring('<_w>' + content + '</_w>')
        cell_element.text = parsed.text
        for child in parsed:
            cell_element.append(child)
    except etree.ParseError:
        cell_element.text = content


class Table(BaseTable):

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super(
            Table,
            self).__init__(
            xrtnode=xrtnode,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        # self.internalId = kw.get("internalId", "tabletable_{0}_{1}".format(jobInfo.get('jobid','xxx'),BaseTable.tableCount))

    def as_data_etree(self):
        root = super().as_data_etree()
        root.set('transpose', 'True' if self.transpose else 'False')
        for child in self.data_as_etree().findall('.//table'):
            root.append(child)
        return root

    def data_as_etree(self, fileName=None):

        hassubhead = sum([x is not None for x in self.colsubtitle]) > 0
        eleTree = etree.parse(StringIO('<script></script>'))
        ccp4_data = eleTree.getroot()
        ccp4_data.set('id', 'data_' + self.internalId)
        ccp4_data.set('type', 'application/xml')

        table = etree.Element('table')
        if not self.transpose:
            head = etree.Element('thead')
            table.append(head)
            headtr = etree.Element('tr')
            for i in range(len(self.coltitle)):
                if self.coltitle[i] is not None:
                    span = 1
                    while i + span < len(self.coltitle):
                        if self.coltitle[i + span] is not None:
                            break
                        span += 1
                    th = etree.Element('th')
                    if i >= 0 and i < len(
                            self.colTips) and self.colTips[i] is not None:
                        th.set('title', self.colTips[i])
                    _set_cell_content(th, self.coltitle[i] if self.coltitle[i] else "")
                    if span > 1:
                        th.set('colspan', str(span))
                    headtr.append(th)
            head.append(headtr)
            tbody = etree.Element('tbody')
            if len(self.coldata) > 0:
                for i in range(len(self.coldata[0])):
                    tr = etree.Element('tr')
                    for col in self.coldata:
                        td = etree.Element('td')
                        if i >= 0 and i < len(
                                self.colTips) and self.colTips[i] is not None:
                            td.set('title', self.colTips[i])
                        _set_cell_content(td, str(col[i]))
                        tr.append(td)
                    tbody.append(tr)
            table.append(tbody)
            if hasattr(self, "class_") and self.class_ is not None:
                table.set('class', self.class_)
        else:
            tbody = etree.Element('tbody')
            for i in range(len(self.coldata)):
                tr = etree.Element('tr')
                tbody.append(tr)
                th = etree.Element('th')
                if self.coltitle[i] is not None:
                    _set_cell_content(th, str(self.coltitle[i]))
                if i >= 0 and i < len(
                        self.colTips) and self.colTips[i] is not None:
                    th.set('title', self.colTips[i])
                tr.append(th)
                if hassubhead:
                    th1 = etree.Element('th')
                    if self.colsubtitle[i] is not None:
                        _set_cell_content(th1, self.colsubtitle[i])
                    tr.append(th1)
                for j in range(len(self.coldata[i])):
                    td = etree.Element('td')
                    if i >= 0 and i < len(
                            self.colTips) and self.colTips[i] is not None:
                        td.set('title', self.colTips[i])
                    _set_cell_content(td, str(self.coldata[i][j]))
                    tr.append(td)
            table.append(tbody)
            if hasattr(self, "class_") and self.class_ is not None:
                table.set('class', self.class_ + ' transpose')
            else:
                table.set('class', 'transpose')

        ccp4_data.append(table)

        # print etree.tostring(ccp4_data,encoding='utf-8')
        return ccp4_data


class Div(Container):
    tag = 'div'

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super(
            Div,
            self).__init__(
            xrtnode=xrtnode,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)


class Progress(ReportClass):
    tag = 'progress'

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super(
            Progress,
            self).__init__(
            xrtnode=xrtnode,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.value = 0
        self.label = kw.get('label', '')
        self.outputXml = kw.get('outputXml', False)
        self.initiallyDrawn = kw.get('initiallyDrawn', True)
        if 'value' in kw:
            self.value = kw['value']
        self.max = 100
        if 'max' in kw:
            self.max = kw['max']

    def as_data_etree(self):
        root = super().as_data_etree()
        root.set('value', str(self.value))
        root.set('max', str(self.max))
        root.set('label', str(self.label))
        return root

    def data_as_xml(self, fileName=None):
        dataAsEtree = etree.Element('span')
        dataAsEtree.text = self.label
        dataAsEtree.append(
            etree.Element(
                'progress', value=str(
                    self.value), max=str(
                    self.max), style=self.style))
        rootNode = etree.Element('CCP4ApplicationOutput')
        rootNode.append(dataAsEtree)
        dataAsText = etree.tostring(rootNode)
        with open(fileName, 'wb') as f:
            f.write(dataAsText)
        return CErrorReport()


class PictureGroup(Container):
    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super(
            PictureGroup,
            self).__int__(
            xrtnode=xrtnode,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.help = None
        self.launch = None
        Container.__init__(
            self,
            xrtnode=xrtnode,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)


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
        # print 'GraphGroup.__init__',self.children

        if kw.get('launcher', None) is not None:
            ele = etree.Element('launch')
            ele.set('label', kw['launcher'])
            ele.set('exe', 'loggraph')
            self.launch = Launch(ele, jobInfo=jobInfo)
        else:
            self.launch = None

        if xrtnode is not None:
            help = xrtnode.find(XRTNS + 'help')
            if help is not None:
                self.help = Help(help, mode='graph')
                if xrtnode is not None:
                    xrtnode.remove(help)
        elif 'help' in kw:
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
        # print 'GraphGroup.__init__',self.children

        self.launch = None
        ele = etree.Element('launch')
        if kw.get('launcher', None) is not None:
            self.launchOnly = True
            ele.set('label', kw['launcher'])
            ele.set('exe', 'loggraph')
            if kw.get('withLaunch', True):
                self.launch = Launch(ele, jobInfo=jobInfo)

        if xrtnode is not None:
            help = xrtnode.find(XRTNS + 'help')
            if help is not None:
                self.help = Help(help, mode='graph')
                if xrtnode is not None:
                    xrtnode.remove(help)
        elif 'help' in kw:
            self.help = Help(ref=kw['help'])

    def data_as_xmls(self, fileRoot=None):
        # Beware children could include Loop or other non-Graph elements
        graphObjList = []
        for child in self.children:
            if isinstance(child, FlotGraph):
                graphObjList.append(child)


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
            help = xrtnode.find(XRTNS + 'help')
            if help is not None:
                self.help = Help(help, mode='gallery')
                if xrtnode is not None:
                    xrtnode.remove(help)
        elif 'help' in kw:
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
            help = xrtnode.find(XRTNS + 'help')
            if help is not None:
                self.help = Help(help, mode='gallery')
                if xrtnode is not None:
                    xrtnode.remove(help)
        elif 'help' in kw:
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
            help = xrtnode.find(XRTNS + 'help')
            if help is not None:
                self.help = Help(help, mode='gallery')
                if xrtnode is not None:
                    xrtnode.remove(help)
        elif 'help' in kw:
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
                'description': 'Unable to create RTF file - unable to import PyQt4 '}}

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super(
            Graph,
            self).__init__(
            xrtnode=xrtnode,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.id = kw.get('id', 'graph_' + str(BaseTable.tableCount))
        # self.internalId = kw.get("internalId","table_{0}_{1}".format(jobInfo.get('jobid','xxx'),BaseTable.tableCount))
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
                # print
                # 'Graph.__init__',xrtnode.find(XRTNS+'title'),xrtnode.find(XRTNS+'title').get('select')
                if xrtnode.find(
                        XRTNS + 'title').get('select') is not None and xmlnode is not None:
                    xmlEleList = xmlnode.findall(
                        xrtnode.find(XRTNS + 'title').get('select'))
                    # print 'Graph.__init__ xmlEleList',xmlEleList
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
            ele = etree.Element('launch')
            ele.set('label', kw.get('launcher', 'More graphs'))
            ele.set('exe', 'loggraph')
            self.launch = Launch(
                ele,
                jobInfo=jobInfo,
                ccp4_data_id='data_' +
                self.internalId)

        if xrtnode is not None:
            help = xrtnode.find(XRTNS + 'help')
        else:
            help = kw.get('help', None)
        if help is not None:
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
            xmlnode = xmlnode.findall0(select)
        self.pimpleData = etree.Element('pimple_data')
        if usePlotly:
            self.pimpleData.set('usePlotly', 'True')
        from copy import deepcopy
        for key in ['headers', 'data', 'plot']:
            eleList = xmlnode.findall(key)
            for ele in eleList:
                # print('addPimpleData',key,ele.get('id'))
                self.pimpleData.append(deepcopy(ele))
        for attr in list(xmlnode.attrib.keys()):
            # print('addPimpleData',attr,xmlnode.get(attr))
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
            # print 'Graph.addTable headersEleList',headersEleList
            if len(headersEleList) > 0:
                self.headerText = headersEleList[0].text
                self.headerSeparator = headersEleList[0].get('separator', None)

    def addPlot(self, xrtnode=None, xmlnode=None, **kw):
        if xmlnode is None:
            if len(self.xmldata) > 0:
                xmlnode = self.xmldata[0]
        if xrtnode is None:
            # try:
            if kw.get('plotFile', None):
                from ccp4i2.core import CCP4Utils
                try:
                    text = CCP4Utils.readFile(kw['plotFile'])

                    text = re.sub('<xrt:', '<', text)
                    text = re.sub('</xrt:', '</', text)
                    # print 'addPlot plot',text
                    xrtnode = etree.fromstring(text)
                except BaseException:
                    raise CException(self.__class__, 1, str(kw['plotFile']))
            if kw.get('plot', None) is not None:
                print(kw['plot'])
                

                if hasattr(kw['plot'], "decode"):
                    text = re.sub('<xrt:', '<', kw['plot'].decode())
                else:
                    text = re.sub('<xrt:', '<', kw['plot'])
                text = re.sub('</xrt:', '</', text)
                xrtnode = etree.fromstring(text)
            # except:
            #  raise CException(self.__class__,1,str(text))

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
        # ccp4_data = etree.Element('ccp4_data')
        ccp4_data.set('id', 'data_' + self.internalId)
        ccp4_data.set('style', 'display:none;')

        if self.pimpleData is not None:
            import copy

            # print
            # 'Graph.data_as_etree',self.pimpleData,len(self.pimpleData),self.pimpleData.get('title')
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

            # print 'data_as_etree',self.plots
            for plot in self.plots:
                if isinstance(plot, Plot):
                    ccp4_data.append(plot.as_etree())
                else:
                    ccp4_data.append(plot)

        return ccp4_data

    def getListOfRows(self):
        rowdatalist = []
        if len(self.coldata) > 0:
            rowdatalist = self.coldata
        elif len(self.tableText) > 0:
            for row in self.tableText.split('\n'):
                rowdata = []
                for item in row.split():
                    rowdata.append(item)
                if len(rowdata) > 0:
                    rowdatalist.append(rowdata)
        elif self.pimpleData is not None:
            dataEleList = self.pimpleData.findall('./data')
            # print 'Graph.getListOfRows dataEleList',dataEleList
            if len(dataEleList) > 0:
                for row in dataEleList[0].text.split('\n'):
                    rowdata = []
                    for item in row.split():
                        rowdata.append(item)
                    if len(rowdata) > 0:
                        rowdatalist.append(rowdata)
            # Do we have some column headers - should be space separated words
            headerEleList = self.pimpleData.findall('./headers')
            if len(headerEleList) > 0:
                separator = headerEleList[0].get('separator', ' ')
                headerList = headerEleList[0].text.split(separator)
                if len(rowdatalist) > 0 and len(
                        rowdatalist[0]) == len(headerList):
                    rowdatalist.insert(0, headerList)
        return rowdatalist

    def data_as_rtf(self, document=None):
        """
        RTF table export functionality - DEPRECATED.

        This method previously used Qt (PySide2) for RTF table generation.
        RTF export is no longer supported in the Qt-free version.
        Use data_as_csv() for data export instead.
        """
        import logging
        logger = logging.getLogger(f"ccp4i2:{__name__}")
        logger.warning(
            "data_as_rtf() is deprecated - RTF export requires Qt which is no longer available")
        return None

    def data_as_csv(self, fileName=None):
        rowList = self.getListOfRows()
        if len(rowList) == 0:
            return CErrorReport(Report, 107, fileName)
        import csv
        try:
            f = open(fileName, 'wb')
        except BaseException:
            return CErrorReport(Report, 107, fileName)
        try:
            writer = csv.writer(f)
            if len(self.coltitle) > 0:
                writer.writerow(self.coltitle)
            for row in rowList:
                writer.writerow(row)
        except Exception as e:
            return CErrorReport(Report, 107, fileName)
        try:
            f.close()
        except BaseException:
            return CErrorReport(Report, 107, fileName)

        return CErrorReport()

    def data_as_xml(self, fileName=None):
        rowList = self.getListOfRows()
        if len(rowList) == 0:
            return CErrorReport(Report, 108, fileName)
        try:
            f = open(fileName, 'wb')
        except BaseException:
            return CErrorReport(Report, 108, fileName)
        try:
            dataAsEtree = self.data_as_etree()
            # As it comes out here, the data tag is in ccp4 namespace (i.e.
            # ccp4_data:data)
            dataAsEtree.tag = 'CCP4Table'
            rootNode = etree.Element('CCP4ApplicationOutput')
            rootNode.append(dataAsEtree)
            dataAsText = etree.tostring(rootNode)
            f.write(dataAsText)
        except Exception as e:
            return CErrorReport(Report, 108, fileName)
        try:
            f.close()
        except BaseException:
            return CErrorReport(Report, 108, fileName)

        return CErrorReport()

# Graph class


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

        # if kw.get('launcher',None) is not None:
        self.launch = None
        self.launchOnly = False
        self.launcherLabel = None
        self.flot_id = kw.get('internalId', None)
        ele = etree.Element('launch')
        if kw.get('launcher', None) is not None:
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


class GenericElement(ReportClass):
    def __init__(self, tag=None, text=None, **kw):
        super().__init__()
        self.tag = tag
        self.id = kw.get('id', None)
        self.class_ = kw.get('class_', None)
        self.text = text
        self.attributes = kw
        self.children = []
        # print 'GenericElement.init', self.tag , self.text

    def append(self, name, text=None, **kw):
        if isinstance(name, GenericElement):
            self.children.append(name)
        else:
            self.children.append(GenericElement(name, text, **kw))
        return self.children[-1]

    def set(self, key, value):
        self.attributes[key] = value

    def as_data_etree(self):
        root = super().as_data_etree()
        root.append(self.as_etree())
        return root

    def as_etree(self):
        # print 'GenericElement.as_etree',self.tag,self.text,self.attributes
        ele = etree.Element(self.tag)
        if self.id is not None:
            ele.set('id', self.id)
        if self.class_ is not None:
            ele.set('class', self.class_)
        if self.text is not None:
            ele.text = self.text
        for key, value in list(self.attributes.items()):
            ele.set(key, str(value).strip())
        for obj in self.children:
            ele.append(obj.as_etree())
        return ele


def parse_from_unicode(unicode_str):
    from lxml import etree as lxml_etree
    utf8_parser = lxml_etree.XMLParser(encoding='utf-8')
    s = unicode_str.encode('utf-8')
    return lxml_etree.fromstring(s, parser=utf8_parser)


class Plot(GenericElement):

    ERROR_CODES = {
        101: {
            'severity': SEVERITY_WARNING,
            'description': 'Failed loading Plot schema'},
        102: {
            'description': 'Failed validation'},
    }

    def __init__(self, text=None, **kw):
        GenericElement.__init__(self, tag='plot', text=text, **kw)

    def validate(self):
        import os

        from lxml import etree as lxml_etree

        from ccp4i2.core import CCP4Utils
        try:
            schemafile = os.path.join(
                CCP4Utils.getCCP4I2Dir(),
                'pimple',
                'CCP4ApplicationOutput.xsd')
            schematree = parse_from_unicode(open(schemafile).read())
            schema = lxml_etree.XMLSchema(schematree)
        except BaseException:
            raise
            return CErrorReport(self.__class__, 101, schemafile)

        tree = self.as_etree()

        table = etree.Element('CCP4Table')
        table.append(tree)
        output = etree.Element('CCP4ApplicationOutput')
        output.append(table)
        valid = schema.validate(lxml_etree.fromstring(etree.tostring(output)))
        if not valid:
            log = str(schema.error_log)
            # print 'validate log',log
            err = CErrorReport()
            err.append(self.__class__, 102, log)
            return err
        else:
            return CErrorReport()

    def as_etree(self):
        tree = GenericElement.as_etree(self)
        return tree


class IODataList(ReportClass):

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super().__init__()
        self.jobInfo = {}
        self.jobInfo.update(jobInfo)

    def _file_ids_for_role(self, role):
        """Return list of file IDs for the given role."""
        if role not in self.jobInfo:
            return []
        return [str(fi['fileId']) for fi in self.jobInfo[role]]

    def _build_data_etree(self, title, roles):
        """Build minimal as_data_etree with just file IDs and title.

        The frontend (CCP4i2ReportInputOutputData) only needs:
        - <h5> elements for the accordion title
        - <div id="input_file_{UUID}"> elements to extract file UUIDs
        All file metadata is fetched client-side from the project's file list.
        """
        root = super().as_data_etree()
        inner = etree.Element('root')
        head = etree.Element('h5')
        head.text = title
        inner.append(head)
        for role in roles:
            for file_id in self._file_ids_for_role(role):
                div = etree.Element('div')
                div.set('id', 'input_file_' + file_id)
                inner.append(div)
        root.append(inner)
        return root


class InputData(IODataList):

    def as_data_etree(self):
        return self._build_data_etree('Input Data', ['inputfiles', 'importedfiles'])


class OutputData(IODataList):

    def as_data_etree(self):
        return self._build_data_etree('Output Data', ['outputfiles'])


class ImportedFiles(IODataList):

    def as_data_etree(self):
        return self._build_data_etree('Imported Files', ['importedfiles'])


class Title(ReportClass):
    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super().__init__()
        import time

        # print 'Title.__init__ jobInfo',jobInfo
        self.title0 = 'Job ' + \
            str(jobInfo['jobnumber']) + ': ' + jobInfo['tasktitle']
        if jobInfo.get(
                'jobtitle',
                None) is not None and len(
                jobInfo['jobtitle']) > 0:
            self.title1 = jobInfo['jobtitle']
        else:
            self.title1 = None

        self.title2 = time.strftime(
            '%H:%M %d-%b-%Y',
            time.localtime(
                jobInfo['creationtime']))

    def as_data_etree(self):
        root = super().as_data_etree()
        title1 = getattr(self, 'title1', '')
        if title1 is None:
            title1 = ''
        root.set('title1', title1)
        title2 = getattr(self, 'title2', '')
        if title2 is None:
            title2 = ''
        root.set('title2', title2)
        return root


class JobDetails(ReportClass):
    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super().__init__()
        self.id = kw.get('id', None)
        self.class_ = kw.get('class_', None)
        self.jobInfo = {}
        self.jobInfo.update(jobInfo)

    def as_data_etree(self):
        import time
        root = super().as_data_etree()
        root.set(
            'creationtime',
            time.strftime(
                '%H:%M %d-%b-%Y',
                time.localtime(
                    self.jobInfo['creationtime'])))
        root.set(
            'finishtime',
            time.strftime(
                '%H:%M %d-%b-%Y',
                time.localtime(
                    self.jobInfo['finishtime'])))
        root.set('status', self.jobInfo.get('status', 'Unknown'))
        return root

    def makeRow(self, key, value):
        tr = etree.Element('tr')
        th = etree.Element('th')
        th.text = str(key)
        tr.append(th)
        td = etree.Element('td')
        td.text = str(value)
        tr.append(td)
        return tr


class JobLogFiles(ReportClass):
    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super().__init__()
        self.id = kw.get('id', None)
        self.class_ = kw.get('class_', None)
        self.jobInfo = {}
        self.jobInfo.update(jobInfo)

    def as_data_etree(self):
        import time
        root = super().as_data_etree()
        root.set(
            'creationtime',
            time.strftime(
                '%H:%M %d-%b-%Y',
                time.localtime(
                    self.jobInfo['creationtime'])))
        root.set(
            'finishtime',
            time.strftime(
                '%H:%M %d-%b-%Y',
                time.localtime(
                    self.jobInfo['finishtime'])))
        root.set('status', self.jobInfo.get('status', 'Unknown'))
        return root

    def makeRow(self, key, value):
        tr = etree.Element('tr')
        th = etree.Element('th')
        th.text = str(key)
        tr.append(th)
        td = etree.Element('td')
        td.text = str(value)
        tr.append(td)
        return tr


class Help:
    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        self.id = kw.get('id', None)
        if xrtnode is not None:
            self.ref = xrtnode.get('ref', None)
        else:
            self.ref = kw.get('ref', None)
        if self.ref is not None and self.ref[0] == '$':

            from ccp4i2.core import CCP4Utils
            if sys.platform == "win32":
                # This had better be sane.
                tweak = CCP4Utils.getCCP4I2Dir().replace('\\', '/')
                tweakref = self.ref.replace('\\', '/')              # Ditto
                self.ref = re.sub(r'\$CCP4I2', tweak, tweakref)
                self.ref = os.path.normpath(self.ref)
            else:
                self.ref = re.sub(
                    r'\$CCP4I2', CCP4Utils.getCCP4I2Dir(), self.ref)
        if xrtnode is not None:
            self.label = xrtnode.get(
                'label',
                'About this ' +
                kw.get(
                    'mode',
                    ''))
        else:
            self.label = kw.get('label', 'About this ' + kw.get('mode', ''))


class Launch:

    counter = 0

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):

        Launch.counter += 1
        self.id = kw.get('id', None)
        # self.internalId = 'launcher_'+str(jobInfo.get('jobid','xxx'))+'_'+str(Launch.counter)
        self.jobId = jobInfo.get('jobid', None)
        self.exe = None
        self.label = None
        # This is a list - could be more than one graph
        self.ccp4_data_id = []
        self.sceneFile = None

        if xrtnode is not None:
            self.exe = xrtnode.get('exe', None)
            self.label = xrtnode.get('label', None)
            if xrtnode.get('ccp4_data_id', None) is not None:
                self.ccp4_data_id.append(xrtnode.get('ccp4_data_id'))
            self.sceneFile = xrtnode.get('sceneFile', None)
        self.exe = kw.get('exe', self.exe)
        self.label = kw.get('label', self.label)
        if kw.get('ccp4_data_id', None) is not None:
            self.ccp4_data_id.append(kw.get('ccp4_data_id'))
        # Use relative path in case project moved - Launcher widget will know
        # jobId and be able to find file
        self.sceneFile = kw.get('sceneFile', self.sceneFile)
        if self.sceneFile is not None:
            self.sceneFile = './' + os.path.split(self.sceneFile)[-1]

    def appendDataId(self, ccp4_data_id=None):
        if self.ccp4_data_id.count(ccp4_data_id) == 0:
            self.ccp4_data_id.append(ccp4_data_id)


class CopyToClipboard:
    def __init__(self, text="", label="Copy to clipboard", **kw):
        self.text = text
        self.label = label


class CopyUrlToClipboard:
    def __init__(self, text="", label="Copy to clipboard", **kw):
        self.text = text
        self.label = label
        self.projectId = kw.get("projectId")
        self.jobnumber = kw.get("jobnumber")


class Download:

    counter = 0

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        Launch.counter += 1
        self.id = kw.get('id', None)
        # self.internalId = 'download_'+str(jobInfo.get('jobid','xxx'))+'_'+str(Download.counter)
        self.jobId = jobInfo.get('jobid', None)
        self.dataName = None
        self.label = None

        if xrtnode is not None:
            self.dataName = xrtnode.get('dataName', None)
            self.label = xrtnode.get('label', None)
        self.dataName = kw.get('dataName', self.dataName)
        self.label = kw.get('label', self.label)


class LaunchTask:

    counter = 0

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        Launch.counter += 1
        self.id = kw.get('id', None)
        # self.internalId = 'launcher_'+str(jobInfo.get('jobid','xxx'))+'_'+str(Launch.counter)
        self.jobId = jobInfo.get('jobid', None)
        if xrtnode is not None:
            self.taskName = xrtnode.get('taskName', None)
            self.label = xrtnode.get('label', None)
            self.ccp4_data_id = xrtnode.get('ccp4_data_id', ccp4_data_id)

        self.taskName = kw.get('taskName', self.taskName)
        self.label = kw.get('label', self.label)
        self.ccp4_data_id = kw.get('ccp4_data_id', self.ccp4_data_id)


class Picture:
    """
    Picture element for CCP4mg/Moorhen scene visualization.

    When a sceneFile is provided, it is used directly without copying.
    When building a scene from scratch (via xrtnode or scene param),
    a new scene file is created with a deterministic name based on label.
    """

    ERROR_CODES = {
        101: {
            'description': 'Error reading picture definition'}, 102: {
            'description': 'Error parsing xml from scene file'}, 103: {
                'description': 'Error parsing xml from scene description'}, 104: {
                    'description': 'No scene description provided'}}

    # Class-level counter for generating unique scene file names within a session
    _scene_counter = {}

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        import copy
        self.id = kw.get('id', None)
        self.class_ = kw.get('class_', None)
        self.launchList = []

        # Track whether we're using an externally-provided scene file
        external_scene_file = None

        bodyEle = etree.Element('ccp4i2_body')
        sceneEle = etree.Element('scene')
        bodyEle.append(sceneEle)

        sceneRoot = None
        if xrtnode is not None:
            sceneRoot = xrtnode
        elif kw.get('scene', None) is not None:
            try:
                sceneRoot = etree.fromstring(kw['scene'], PARSER())
            except BaseException:
                raise CException(self.__class__, 103, kw['scene'])
        elif kw.get('sceneFile', None) is not None:
            import os

            from ccp4i2.core import CCP4Utils
            fileName = kw.get('sceneFile')
            if fileName[0:8] == '$CCP4I2/':
                fileName = os.path.join(CCP4Utils.getCCP4I2Dir(), fileName[8:])
            # print 'Report.Pictures fileName',fileName
            if not os.path.exists(fileName):
                raise CException(self.__class__, 101, fileName)
            else:
                try:
                    sceneRoot = etree.fromstring(
                        open(fileName).read(), PARSER())
                    # Remember that we have an external scene file - no need to copy
                    external_scene_file = fileName
                except BaseException:
                    raise CException(self.__class__, 102, fileName)

        if sceneRoot is None:
            raise CException(self.__class__, 104, fileName)

        for child in sceneRoot:
            sceneEle.append(copy.deepcopy(child))
            bodyEle = applySelect(bodyEle, xmlnode, jobInfo)

        # print etree.tostring(bodyEle,pretty_print=True)

        import os

        from ccp4i2.core import CCP4File

        # If an external scene file was provided, use it directly without copying
        if external_scene_file is not None:
            self.picDefFile = external_scene_file
        else:
            # Generate a deterministic scene file name based on job and label
            # Use label if available, otherwise use a counter per job
            job_id = jobInfo.get('jobid', 'unknown')
            label = kw.get('label', None)
            if xrtnode is not None:
                label = xrtnode.get('label', label)

            if label:
                # Create filename from label (sanitize for filesystem)
                safe_label = "".join(c if c.isalnum() or c in '_-' else '_' for c in label)
                scene_name = f'scene_{safe_label}.scene.xml'
            else:
                # Use counter for this job
                if job_id not in Picture._scene_counter:
                    Picture._scene_counter[job_id] = 0
                Picture._scene_counter[job_id] += 1
                scene_name = f'scene_{Picture._scene_counter[job_id]}.scene.xml'

            if 'fileroot' in jobInfo and jobInfo['fileroot'] is not None:
                scene_path = jobInfo['fileroot'] + scene_name
            else:
                scene_path = os.path.join(os.getcwd(), scene_name)

            self.picDefFile = CCP4File.CI2XmlDataFile(fullPath=scene_path)
            self.picDefFile.header.setCurrent()
            self.picDefFile.header.function.set('MGSCENE')
            self.picDefFile.header.projectId.set(jobInfo.get('projectid', None))
            self.picDefFile.header.projectName.set(
                jobInfo.get('projectname', None))
            self.picDefFile.header.jobId.set(jobInfo.get('jobid', None))
            self.picDefFile.header.jobNumber.set(jobInfo.get('jobnumber', None))
            self.picDefFile.saveFile(bodyEtree=bodyEle, useLXML=False)

        # report.pictureQueue.append(self.picDefFile.__str__())
        if xrtnode is not None:
            self.label = xrtnode.get('label', None)
        else:
            self.label = kw.get('label', None)

        launchNode = etree.Element('launch')
        launchNode.set('exe', 'CCP4mg')
        launchNode.set('label', 'View in CCP4mg')
        if self.picDefFile is not None:
            launchNode.set('sceneFile', str(self.picDefFile))
        self.launchList.append(Launch(launchNode, jobInfo=jobInfo))

        launchNode = etree.Element('launch')
        launchNode.set('exe', 'Coot')
        launchNode.set('label', 'View in Coot')
        # launchNode.set('cootScript',self.cootScript)
        self.launchList.append(Launch(xrtnode=launchNode, jobInfo=jobInfo))


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


class GenericReport(Report):
    def __init__(self, xmlnode=None, jobInfo={}, **kw):
        Report. __init__(self, xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        title = jobInfo.get('tasktitle', '')
        self.addText(text=title)


class Reference(ReportClass):
    ERROR_CODES = {}

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super().__init__()
        self.id = kw.get('id', None)
        if xrtnode is not None:
            data = xrtnode
        else:
            data = kw
        self.href = data.get('href', None)
        self.authorList = data.get('authorList', [])
        if data.get('author', None) is not None:
            self.authorList.append(data.get('author', None))
        self.source = data.get('source', None)
        self.articleTitle = data.get('articleTitle', None)
        self.articleLink = data.get('articleLink', None)

    def as_data_etree(self):
        root = super().as_data_etree()
        if self.articleTitle is not None:
            root.set('articleTitle', self.articleTitle)
        if self.articleLink is not None:
            root.set('articleLink', self.articleLink)
        if self.source is not None:
            root.set('source', self.source)
        if len(self.authorList) > 0:
            root.set('authorList', str(self.authorList))
        return root


class ReferenceGroup(Container):
    ERROR_CODES = {
        100: {
            'description': 'Failed attempting to load MedLine file - file not found'}}

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        Container.__init__(
            self,
            xrtnode=xrtnode,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.label = 'References'
        self.tag = 'div'
        self._class = 'bibreference_group'
        self.taskName = kw.get('taskName', None)

    def loadFromMedLine(self, taskName):
        path = I2_TOP / "references" / f"{taskName}.medline.txt"
        if not path.exists():
            self.errReport.append(
                self.__class__,
                100,
                f'Taskname: {taskName} Filename: {path}')
            return
        self.taskName = taskName

        from ccp4i2.core import CCP4Utils
        try:
            text = CCP4Utils.readFile(fileName=path)
        except CException as e:
            self.errReport.extend(e)
            return
        textList = text.split('\nPMID- ')
        for text in textList:
            ref = Reference()
            m = re.search(r'TI  -(.*)', text)
            if m is not None:
                ref.articleTitle = m.groups()[0].strip()
            m = re.search(r'SO  -(.*)', text)
            if m is not None:
                ref.source = m.groups()[0].strip()
            m = re.findall(r'AU  -(.*)', text)
            for item in m:
                ref.authorList.append(item.strip())
            m = re.search(r'URL -(.*)', text)
            if m is not None:
                ref.articleLink = m.groups()[0].strip()
            if ref.source is not None:
                self.append(ref)
