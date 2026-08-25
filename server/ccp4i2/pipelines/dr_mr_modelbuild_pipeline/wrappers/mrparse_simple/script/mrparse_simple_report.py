import os
from pathlib import Path

from ccp4i2.report import Report
from ccp4i2.report.embedded_assets import localise_report_assets, vendored_asset


class mrparse_simple_report(Report):

    TASKNAME = 'mrparse_simple'
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

        ResultsI2Folder = parent.addFold(label='MrParse Reports', initiallyOpen=True)
        if os.path.exists(mrparse_rep):
            # MrParse writes its stylesheets and scripts as absolute paths into
            # its own installation, which 404 when the report is served over
            # HTTP and do not survive export or relocation. Copy them in beside
            # the report and rewrite the references to be relative.
            localised = localise_report_assets(
                Path(mrparse_rep),
                marker='mrparse/html',
                subdirectory='mrparse_html',
                # MrParse loads d3 from d3js.org. Serve our own copy: the
                # report is an archival record that should still draw its
                # feature viewer offline, and on a machine with no route to
                # the internet.
                extra_assets={
                    'https://d3js.org/d3.v5.min.js': vendored_asset('d3.v5.min.js'),
                },
            )
            report_name = localised.name if localised else 'mrparse.html'
            ResultsI2Folder.addFileLink(
                label='Open MrParse Results',
                relativePath=f'mrparse_0/{report_name}',
                fileType='html',
            )
        else:
            ResultsI2Folder.append('<p>MrParse report not found</p>')
        return
