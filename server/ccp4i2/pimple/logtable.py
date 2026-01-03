"""
CCP4 log table parsing utilities (Qt-free).

This module provides functions to parse CCP4 log file tables into XML/etree format.
Extracted from MGQTmatplotlib.py to avoid Qt dependencies in Django mode.
"""

from lxml import etree

# Resolution scaling functions for graph axes
scaling_functions = {
    "oneoversqrt": lambda x: 1.0 / (x ** 0.5) if x > 0 else 0,
}


class CCP4Table:
    """Parse and represent a CCP4 $TABLE log entry."""

    def __str__(self):
        string = ''
        string = string + self.title + '\n'
        string = string + str(self.graphs_list) + '\n'
        string = string + str(self.headers) + '\n'
        string = string + str(self.data_lines) + '\n'
        return string

    def __init__(self, table):
        tname_start = table.find('$TABLE') + len('$TABLE')
        tname_start = table[tname_start:].find(':') + tname_start + 1

        tname_end = table[tname_start:].find(':')
        if table[tname_start:].find('\n') < tname_end:
            tname_end = table[tname_start:].find('\n')
        self.title = table[tname_start:tname_end + tname_start].strip()

        rest = table[tname_end + tname_start:]

        self.graphs_list = []

        graphs_start = rest.find('$GRAPHS') + len('$GRAPHS')
        graphs_end = rest[graphs_start:].find('$$')
        graphs = rest[graphs_start:graphs_end + graphs_start].strip()
        colCount = 0
        currentHeader = ''

        for g in graphs:
            if g == ':':
                colCount += 1
            if colCount == 4:
                s, name, n, cols = currentHeader.strip().split(':')
                theRange = None
                if n.startswith('N'):
                    theRange = ((None, None), (0, None))
                if not n.startswith('XD') and not n.startswith('N') and not n.startswith('A') and n.find('x') > 0 and n.find('|') > 0:
                    try:
                        xRange, yRange = n.split('x')
                        xRange = xRange.split('|')
                        yRange = yRange.split('|')
                        theRange = (float(xRange[0]), float(xRange[1])), (float(yRange[0]), float(yRange[1]))
                    except:
                        pass
                self.graphs_list.append((name.strip(), cols, theRange))
                colCount = 0
                currentHeader = ''
            else:
                currentHeader += g

        rest2 = rest[graphs_end + graphs_start:]

        headers_start = rest2.find('$$') + 2
        headers_end = rest2[headers_start:].find('$$')
        self.headers = rest2[headers_start:headers_end + 2].replace('\n', ' ').split()

        rest3 = rest2[headers_end + 2:]
        rest3 = rest3[rest3.find('$$') + 2:]
        data_start = rest3.find('$$') + 2
        rest3 = rest3[data_start:].lstrip()
        data_end = rest3.find('$$')
        data_start = 0
        data_lines_str = rest3[data_start:data_end].strip().split('\n')
        self.data_lines = []
        for line_str in data_lines_str:
            vals = []
            vals_str = line_str.split()
            for val_str in vals_str:
                try:
                    vals.append(float(val_str))
                except:
                    vals.append('NA')
            self.data_lines.append(vals)
        self.x_scaling = False
        self.y_scaling = False
        self.custom_x_label = False
        # Check if first column needs special scaling
        if len(self.headers) > 0 and (
            self.headers[0] == 'M(4SSQ/LL)' or
            self.headers[0] == '1/d^2' or
            self.headers[1] == '1/d^2' or
            self.headers[0] == '1/resol^2' or
            self.headers[1] == '1/resol^2'
        ):
            self.custom_x_label = u'Resolution / \xc5'
            self.x_scaling = scaling_functions["oneoversqrt"]

    def toEtree(self):
        """Convert table to lxml etree Element."""
        t = etree.Element("CCP4Table")
        t.attrib["title"] = self.title
        theData = etree.Element("data")
        theDataText = ''
        for tdl in self.data_lines:
            theDataText += ''.join([("%s " % n) for n in tdl[:]]).strip() + '\n'
        theData.text = theDataText
        t.append(theData)
        theHeaders = etree.Element("headers", separator=' ')
        theHeaders.text = ''.join([("%s " % n) for n in self.headers[:]]).strip() + '\n'
        t.append(theHeaders)
        for g in self.graphs_list:
            p = etree.Element('plot')
            theTitle = etree.Element("title")
            theTitle.text = g[0]
            p.append(theTitle)
            t.append(p)
            for s in g[1].split(',')[1:]:
                plotline = etree.Element('plotline', xcol=g[1].split(',')[0], ycol=s)
                p.append(plotline)
        return t


def CCP4LogToEtree(b):
    """
    Parse CCP4 log text containing $TABLE entries into an lxml etree.

    Args:
        b: Log file content as string

    Returns:
        lxml.etree.Element: CCP4ApplicationOutput element containing parsed tables
    """
    splits = b[b.find("$TABLE"):].split("$TABLE")
    bigtree = etree.Element("CCP4ApplicationOutput")
    for ns in splits[1:]:
        ns = "$TABLE" + ns
        table = CCP4Table(ns)
        tree = table.toEtree()
        bigtree.append(tree)
    return bigtree


def CCP4LogToXML(b):
    """
    Parse CCP4 log text containing $TABLE entries into XML string.

    Args:
        b: Log file content as string

    Returns:
        str: XML string representation
    """
    splits = b[b.find("$TABLE"):].split("$TABLE")
    header = """<?xml version="1.0" encoding="UTF-8" ?>\n"""

    NSMAP = {'xsi': "http://www.w3.org/2001/XMLSchema-instance"}
    NS = NSMAP['xsi']
    location_attribute = '{%s}noNamespaceSchemaLocation' % NS
    bigtree = etree.Element(
        "CCP4ApplicationOutput",
        nsmap=NSMAP,
        attrib={location_attribute: 'http://www.ysbl.york.ac.uk/~mcnicholas/schema/CCP4ApplicationOutput.xsd'}
    )
    for ns in splits[1:]:
        ns = "$TABLE" + ns
        table = CCP4Table(ns)
        tree = table.toEtree()
        bigtree.append(tree)

    status_xml = header + etree.tostring(bigtree, encoding='unicode', pretty_print=True)
    return status_xml


def CCP4LogFileNameToXML(f):
    """
    Parse CCP4 log file into XML string.

    Args:
        f: Path to log file

    Returns:
        str: XML string representation
    """
    if f.endswith('.log') or f.endswith('.txt'):
        with open(f) as fo:
            b = fo.read()
        return CCP4LogToXML(b)
    else:
        return CCP4LogToXML("")


def CCP4LogFileNameToEtree(f):
    """
    Parse CCP4 log file into lxml etree.

    Args:
        f: Path to log file

    Returns:
        lxml.etree.Element: CCP4ApplicationOutput element containing parsed tables
    """
    if f.endswith('.log') or f.endswith('.txt'):
        with open(f) as fo:
            b = fo.read()
        return CCP4LogToEtree(b)
    else:
        return CCP4LogToEtree("")
