from ccp4i2.report import Report


class coot_rebuild_report(Report):
    TASKNAME = 'coot_rebuild'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        self.addText(text='Happily finished')
        pdbsWrittenPath = './/coot_rebuild/number_output_files'
        pdbsWrittenStringS = xmlnode.findall(pdbsWrittenPath)
        if len(pdbsWrittenStringS) > 0:
            pdbsWrittenString = xmlnode.findall(pdbsWrittenPath)[0].text
            self.append('<br/>')
            self.addText(text='Number of PDBs written: ' + pdbsWrittenString)
        dictsWrittenPath = './/coot_rebuild/number_output_dicts'
        dictsWrittenStringS = xmlnode.findall(dictsWrittenPath)
        if len(dictsWrittenStringS) > 0:
            dictsWrittenString = xmlnode.findall(dictsWrittenPath)[0].text
            self.append('<br/>')
            self.addText(text='Number of DICTs written: ' + dictsWrittenString)
        self.drawWarnings()

    def drawWarnings(self, parent=None):
        if parent is None: parent = self
        warnings = self.xmlnode.findall('.//Warnings/Warning')
        warningsFolder = parent.addFold(label='Warnings', initiallyOpen=True)
        if len(warnings)>0:
            for warning in warnings:
                warningsFolder.addText(text=warning.text,style='color:red;')
                warningsFolder.append('<br/>')
        else:
            warningsFolder.addText(text='No warnings from coot_rebuild')
        return
