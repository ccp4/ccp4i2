# Copyright (C) 2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
"""
Report element classes.

Text, tables, folds, and other basic report building blocks.
These are the most commonly used elements in task report scripts.
"""

from __future__ import annotations

import xml.etree.ElementTree as etree
from io import StringIO
from typing import Any, Callable

from ccp4i2.core.CCP4ErrorHandling import SEVERITY_WARNING, CErrorReport, CException
from ccp4i2.report.core import ReportClass, Container


class Results(Container):
    """Results container — uses Container.as_data_etree() for data path."""
    pass


class Fold(Container):
    """Collapsible section in the report.

    Wraps child elements in an expandable/collapsible fold with a label.
    The ``brief`` attribute provides a short label for the fold tab bar.
    """

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        super(Fold, self).__init__(
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.label: str = kw.get('label', 'Fold')
        self.brief: str | None = kw.get('brief', None)
        self.initiallyOpen: bool = kw.get('initiallyOpen', False)

    def as_data_etree(self) -> etree.Element:
        if self.brief is None or len(self.brief) == 0:
            self.brief = self.label
        root = super().as_data_etree()
        root.set('label', self.label)
        root.set('initiallyOpen', str(self.initiallyOpen))
        root.set('brief', self.brief)
        return root


class Text(ReportClass):
    """Inline text element, optionally populated from XML via ``select``."""

    tag: str = 'span'

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        super(Text, self).__init__(
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.text: str = ''
        self.tail: str = ''
        if kw.get('text', None) is not None:
            self.text = kw['text']
        elif kw.get('select', None) is not None and xmlnode is not None:
            nodes = xmlnode.findall(kw['select'])
            for node in nodes:
                self.text += node.text



class Pre(Text):
    """Preformatted text block."""

    tag: str = 'pre'


class Alignment(Text):
    """Sequence alignment viewer element.

    Displays ClustalW-format alignment text in an interactive viewer
    with colour-coded residues, conservation bars, and hover tooltips.
    """

    tag: str = 'alignment'


class FetchPre(Text):
    """Preformatted text block (fetched from file)."""

    tag: str = 'pre'


class Status(Container):
    """Status container — uses Container.as_data_etree() for data path."""
    pass


class Copy(ReportClass):
    """Copies XML content from the program output into the report."""

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        **kw: Any,
    ) -> None:
        super(Copy, self).__init__(xmlnode=xmlnode, **kw)
        self.text: str = ""
        self.root: etree.Element = etree.Element('root')
        if xmlnode is not None and kw.get('select', False):
            self.root.extend(xmlnode.findall(kw.get('select', False)))


class Generic(ReportClass):
    """Generic HTML element, parsed from a text string."""

    ERROR_CODES: dict = {1: {'description': 'Can not interpret text'}, 2: {
        'description': 'Error substituting values into generic item'}}

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        text: str | None = None,
        defaultTag: str = 'p',
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        super().__init__()
        self.id: str | None = kw.get('id', None)
        self.class_: str | None = kw.get('class_', None)
        self.xmltree: etree.Element | None = None

        parsed = None
        if text is not None:
            try:
                parsed = etree.fromstring(text.encode('utf-8'))
            except BaseException:
                try:
                    if isinstance(text, bytes):
                        parsed = etree.fromstring(
                            '<' + defaultTag + '>' + text.encode('utf-8') + '</' + defaultTag + '>')
                    else:
                        parsed = etree.fromstring(
                            '<' + defaultTag + '>' + text + '</' + defaultTag + '>')
                except BaseException:
                    raise

        if parsed is not None:
            self.xmltree = parsed

    def as_data_etree(self) -> etree.Element:
        root = super().as_data_etree()
        root.append(self.xmltree)
        return root


class BaseTable(ReportClass):
    """Base class for data tables with column-oriented storage.

    Columns are added via ``addData()`` with XML xpath selectors or
    explicit data lists. Supports transposed layout and column tooltips.
    """

    tableCount: int = 0

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        super(BaseTable, self).__init__(xmlnode=xmlnode, **kw)
        BaseTable.tableCount += 1
        self.id: str = kw.get('id', 'table_' + str(BaseTable.tableCount))
        self.coldata: list[list] = []
        self.coltitle: list[str | None] = []
        self.colsubtitle: list[str | None] = []
        self.colTips: list[str | None] = []
        self.xmldata: list[etree.Element] = []
        self.help = None
        self.download = None
        self.outputCsv: bool = kw.get('outputCsv', True)
        downloadable: bool = kw.get('downloadable', False)

        if downloadable:
            from ccp4i2.report.actions import Download
            if self.title is not None:
                self.download = Download(
                    jobInfo=jobInfo, dataName=self.id + '_' + self.title)
            else:
                self.download = Download(jobInfo=jobInfo, dataName=self.id)
        self.excludeIfDataMissing = kw.get('excludeIfDataMissing', False)
        if xmlnode is not None:
            if kw.get('select', None) is not None:
                self.xmldata = []
                for p in kw['select'].split("|"):
                    self.xmldata.extend(xmlnode.findall(p.strip()))
            elif kw.get('selectNodes', None) is not None:
                self.xmldata = kw['selectNodes']
            else:
                self.xmldata = [xmlnode]
        self.transpose = kw.get('transpose', False)
        if kw.get('help', None) is not None:
            from ccp4i2.report.actions import Help
            self.help = Help(ref=kw['help'])

    def addData(
        self,
        xmldata: list[etree.Element] | None = None,
        title: str | None = None,
        subtitle: str | None = None,
        expr: str | None = None,
        function: Callable | None = None,
        select: str | None = None,
        data: list | None = None,
        tip: str | None = None,
    ) -> None:
        """Add a column of data to the table."""
        if data is None:
            data = []
        if tip is not None:
            self.colTips.append(tip)
        else:
            self.colTips += [None]
        colvals: list = []
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



# JavaScript-Enhanced Table class


def _set_cell_content(cell_element: etree.Element, content: str) -> None:
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
    """HTML table rendered from column data.

    Supports both normal and transposed layouts. The ``data_as_etree()``
    method builds the ``<table>`` element with ``<thead>``/``<tbody>``.
    """

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        super(Table, self).__init__(
            xmlnode=xmlnode, jobInfo=jobInfo, **kw)

    def as_data_etree(self) -> etree.Element:
        root = super().as_data_etree()
        root.set('transpose', 'True' if self.transpose else 'False')
        for child in self.data_as_etree().findall('.//table'):
            root.append(child)
        return root

    def data_as_etree(self, fileName: str | None = None) -> etree.Element:
        """Build the HTML table element tree from column data."""
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

        return ccp4_data


class Div(Container):
    """Generic ``<div>`` container."""

    tag: str = 'div'

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        super(Div, self).__init__(
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)


class Progress(ReportClass):
    """Progress bar element (0–max, default 0–100)."""

    tag: str = 'progress'

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        super(Progress, self).__init__(
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.value: int | float = 0
        self.label: str = kw.get('label', '')
        self.outputXml: bool = kw.get('outputXml', False)
        self.initiallyDrawn: bool = kw.get('initiallyDrawn', True)
        if 'value' in kw:
            self.value = kw['value']
        self.max: int | float = 100
        if 'max' in kw:
            self.max = kw['max']

    def as_data_etree(self) -> etree.Element:
        root = super().as_data_etree()
        root.set('value', str(self.value))
        root.set('max', str(self.max))
        root.set('label', str(self.label))
        return root



class GenericElement(ReportClass):
    """Programmatically-built HTML element with children.

    Unlike ``Generic`` (which parses a text string), ``GenericElement``
    builds an element tree via ``append()`` and ``set()`` calls.
    """

    def __init__(self, tag: str | None = None, text: str | None = None, **kw: Any) -> None:
        super().__init__()
        self.tag: str | None = tag
        self.id: str | None = kw.get('id', None)
        self.class_: str | None = kw.get('class_', None)
        self.text: str | None = text
        self.attributes: dict[str, Any] = kw
        self.children: list[GenericElement] = []

    def append(self, name: str | GenericElement, text: str | None = None, **kw: Any) -> GenericElement:
        """Append a child element (by tag name or GenericElement instance)."""
        if isinstance(name, GenericElement):
            self.children.append(name)
        else:
            self.children.append(GenericElement(name, text, **kw))
        return self.children[-1]

    def set(self, key: str, value: Any) -> None:
        """Set an attribute on this element."""
        self.attributes[key] = value

    def as_data_etree(self) -> etree.Element:
        root = super().as_data_etree()
        root.append(self.as_etree())
        return root

    def as_etree(self) -> etree.Element:
        """Convert to a plain ``etree.Element`` tree."""
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


def parse_from_unicode(unicode_str: str) -> Any:
    """Parse a unicode string as XML using lxml with UTF-8 encoding."""
    from lxml import etree as lxml_etree
    utf8_parser = lxml_etree.XMLParser(encoding='utf-8')
    s = unicode_str.encode('utf-8')
    return lxml_etree.fromstring(s, parser=utf8_parser)


class Plot(GenericElement):
    """Plot definition element, validated against the CCP4 schema."""

    ERROR_CODES: dict = {
        101: {
            'severity': SEVERITY_WARNING,
            'description': 'Failed loading Plot schema'},
        102: {
            'description': 'Failed validation'},
    }

    def __init__(self, text: str | None = None, **kw: Any) -> None:
        GenericElement.__init__(self, tag='plot', text=text, **kw)

    def validate(self) -> CErrorReport:
        """Validate this plot against the CCP4ApplicationOutput XSD schema."""
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
            err = CErrorReport()
            err.append(self.__class__, 102, log)
            return err
        else:
            return CErrorReport()

    def as_etree(self) -> etree.Element:
        tree = GenericElement.as_etree(self)
        return tree
