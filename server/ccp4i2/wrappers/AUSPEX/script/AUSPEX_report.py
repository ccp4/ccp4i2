import os

from ccp4i2.report import Report


class AUSPEX_report(Report):

    TASKNAME = 'AUSPEX'
    USEPROGRAMXML = False
    SEPARATEDATA = True

    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,jobStatus=jobStatus,**kw)
        if jobStatus is None or jobStatus.lower() == 'nooutput':
            return
        self.outputXml = self.jobStatus is not None and self.jobStatus.lower().count('running')
        if self.jobStatus is not None and not self.jobStatus.lower().count('running'):
            self.outputXml = False
        self.defaultReport()

    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        results = self.addResults()
        fileroot = self.jobInfo.get('fileroot', '')

        # Map of plot filenames to their display labels
        plot_labels = [
            ("I_plot.png", "Plot of intensities against resolution"),
            ("SigI_plot.png", "Plot of sigma (intensity) against resolution"),
            ("ISigI_plot.png", "Plot of intensity / sigma against resolution"),
            ("F_plot.png", "Plot of amplitudes against resolution"),
            ("SigF_plot.png", "Plot of sigma (intensity) against resolution"),
            ("FSigF_plot.png", "Plot of amplitude / sigma against resolution"),
            ("score.png", "Plot of icefinder score and intensities vs resolution"),
            ("intensities.png", "Plots for intensities"),
            ("amplitudes.png", "Plots for amplitudes"),
        ]

        im_out = self.jobInfo.get('filenames', {}).get('IM_OUT', [])

        found_any = False
        for plot_name, label in plot_labels:
            for img in im_out:
                if not img or os.path.basename(img) != plot_name:
                    continue
                # Convert absolute path to relative (to job directory)
                rel_path = os.path.relpath(img, fileroot) if fileroot else img
                found_any = True
                results.addFileLink(
                    label=label,
                    relativePath=rel_path,
                    fileType='image',
                )

        if not found_any:
            results.append("<p>No AUSPEX plot files found</p>")

