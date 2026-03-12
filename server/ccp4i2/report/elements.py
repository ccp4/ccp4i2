"""
Report element classes.

Text, tables, folds, and other basic report building blocks.
"""

import xml.etree.ElementTree as etree
from io import StringIO

from ccp4i2.core.CCP4ErrorHandling import SEVERITY_WARNING, CErrorReport, CException
from ccp4i2.report.core import (
    ReportClass, Container,
    XRTNS,
    applySelect,
)


class Results(Container):
    """Results container - uses Container.as_data_etree() for data path."""
    pass


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
        self.text = ''
        self.tail = ''
        if xrtnode is not None:
            if xrtnode.text is not None:
                self.text = xrtnode.text
            if xrtnode.tail is not None:
                self.tail = xrtnode.tail
            if xrtnode.get('select') is not None and xmlnode is not None:
                self.text = ''
                nodes = xmlnode.findall(xrtnode.get('select'))
                for node in nodes:
                    if node.text is not None:
                        self.text += node.text
        elif kw.get('text', None) is not None:
            self.text = kw['text']
        elif kw.get('select', None) is not None and xmlnode is not None:
            nodes = xmlnode.findall(kw['select'])
            for node in nodes:
                self.text += node.text

    def as_data_etree(self):
        root = super().as_data_etree()
        return root



class Pre(Text):
    tag = 'pre'


class FetchPre(Text):
    tag = 'pre'


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
            from ccp4i2.report.actions import Download
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
                from ccp4i2.report.actions import Help
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
                from ccp4i2.report.actions import Help
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



# JavaScript-Enhanced Table class


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
        super(Table, self).__init__(
            xrtnode=xrtnode, xmlnode=xmlnode, jobInfo=jobInfo, **kw)

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



class GenericElement(ReportClass):
    def __init__(self, tag=None, text=None, **kw):
        super().__init__()
        self.tag = tag
        self.id = kw.get('id', None)
        self.class_ = kw.get('class_', None)
        self.text = text
        self.attributes = kw
        self.children = []

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
            err = CErrorReport()
            err.append(self.__class__, 102, log)
            return err
        else:
            return CErrorReport()

    def as_etree(self):
        tree = GenericElement.as_etree(self)
        return tree
