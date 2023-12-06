from core import CCP4Modules
from PySide2 import QtCore
from qtgui import CCP4TaskWidget

class CTaskDials_image(CCP4TaskWidget.CTaskWidget):

    TASKNAME = 'dials_image'
    TASKVERSION = 0.0
    TASKMODULE ='data_processing'
    TASKTITLE = 'DIALS Image Viewer'
    SHORTTASKTITLE = "DIALS Image Viewer"
    DESCRIPTION = 'DIALS Image Viewer'
    WHATNEXT = []

    def __init__(self,parent):
        CCP4TaskWidget.CTaskWidget.__init__(self,parent)

    def drawContents(self):
        self.setProgramHelpFile('shelxeMR')
        folder = self.openFolder(folderFunction='inputData', title='Input Data and Run Parameters')
        self.openSubFrame(title='Select input data', tip='A JSON file from DIALS is required', frame=[True])
        self.createLine(['widget', 'JSON_IN'])
        self.container.inputData.JSON_IN.dataChanged.connect(self.populatePickleFile)
        self.createLine(['widget', 'PICKLE_IN'])
        self.closeSubFrame()
        self.closeFolder()

    @QtCore.Slot()
    def populatePickleFile(self):
        pm = CCP4Modules.PROJECTSMANAGER()
        if self.container.inputData.JSON_IN.dbFileId.isSet():
            fileId = self.container.inputData.JSON_IN.dbFileId
            jobparamname = pm.db().getFileInfo(fileId=fileId)['jobparamname']
            jobId =  pm.db().getFileInfo(fileId=fileId)['jobid']
            if jobparamname.split('.')[0] == 'DIALSJOUT':
                fileId =  pm.db().getJobFilesInfo(jobId=jobId,jobParamName='DIALSPOUT.%s' % jobparamname.split('.')[1])[0]['fileId']
                self.getWidget('PICKLE_IN').model.setDbFileId(fileId)
            else:
                self.container.inputData.PICKLE_IN.unSet()
        else:
            self.container.inputData.PICKLE_IN.unSet()

