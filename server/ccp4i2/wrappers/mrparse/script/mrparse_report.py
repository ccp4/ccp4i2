
import os
import xml.etree.ElementTree as etree
from pathlib import Path

from ccp4i2.report import Report
from ccp4i2.report.embedded_assets import localise_report_assets, vendored_asset


def parse_from_unicode(unicode_str):
    utf8_parser = etree.XMLParser(encoding='utf-8')
    s = unicode_str.encode('utf-8')
    return etree.fromstring(s, parser=utf8_parser)

class mrparse_report(Report):

    TASKNAME = 'mrparse'
    USEPROGRAMXML = False
    SEPARATEDATA = True

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        if jobStatus is None or jobStatus.lower() == 'nooutput':
            return
        self.outputXml = self.jobStatus is not None and self.jobStatus.lower().count('running')
        if self.jobStatus is not None and not self.jobStatus.lower().count('running'):
            self.outputXml = False
        self.defaultReport()
        return

    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        parent.append("<p>Finished running MrParse</p>")
        basepath = self.jobInfo['fileroot']
        mrparse_rep = os.path.join(basepath, "mrparse_0", 'mrparse.html')

        self.complexTemplates(parent)
        ResultsI2Folder = parent.addFold(label='MrParse Reports', initiallyOpen=True)
        if not os.path.exists(mrparse_rep):
            ResultsI2Folder.append('<p>MrParse report not found</p>')
            return

        # MrParse writes its stylesheets and scripts as absolute paths into its
        # own installation. Copy them in beside the report and rewrite the
        # references relative to it, so the report renders when served over
        # HTTP and survives project export, import and relocation.
        #
        # This replaces a rewrite that ran on Windows only (elsewhere the Qt
        # app opened the file straight from disk, where absolute paths worked)
        # and pointed at /database/projectid/N/jobnumber/M/file/... — a route
        # that no longer exists, carrying a project identity that changes when
        # a project is imported somewhere else.
        localised = localise_report_assets(
            Path(mrparse_rep),
            marker='mrparse/html',
            subdirectory='mrparse_html',
            # MrParse loads d3 from d3js.org. Serve our own copy: the report
            # is an archival record that should still draw its feature viewer
            # offline, and on a machine with no route to the internet.
            extra_assets={
                'https://d3js.org/d3.v5.min.js': vendored_asset('d3.v5.min.js'),
            },
        )
        report_name = localised.name if localised else 'mrparse.html'

        ResultsI2Folder.append('<span style="font-size:110%">Click on the '
                               'following link to display the browser report '
                               'for the MrParse job</span>')
        ResultsI2Folder.addFileLink(
            label='Open MrParse Results',
            relativePath=f'mrparse_0/{report_name}',
            fileType='html',
        )

    def complexTemplates(self, parent):
        """The entries whose hits match more than one of the sequences
        searched for: each written as one model holding those chains, to be
        placed as a single rigid body (docs/multi-component-mr.md)."""
        try:
            found = self.xmlnode.findall('.//Complexes/Complex')
        except Exception:
            found = []
        if not found:
            return
        fold = parent.addFold(label='Complex templates', initiallyOpen=True)
        fold.append('<p>Hits from one PDB entry match more than one of the sequences '
                    'searched for. Each entry below is also written as one model holding '
                    'those chains, in the entry\'s own frame, to be placed as a single '
                    'rigid body: the number of copies is then the number of copies of the '
                    'complex. Use it when the subunits are expected to sit as they do in '
                    'the entry; otherwise search for the single-chain models as separate '
                    'components in Expert molecular replacement.</p>')
        table = fold.addTable()
        entries, files, components = [], [], []
        for element in found:
            entries.append(element.get('entry', ''))
            files.append(element.get('file', ''))
            components.append('; '.join(
                f"chain {c.get('chain')}: {c.get('target')} "
                f"({100 * float(c.get('identity') or 0):.0f}% identity, hit {c.get('hit')})"
                for c in element.findall('Component')))
        table.addData(title='Entry', data=entries)
        table.addData(title='Model', data=files)
        table.addData(title='Chains and the sequences they match', data=components)
