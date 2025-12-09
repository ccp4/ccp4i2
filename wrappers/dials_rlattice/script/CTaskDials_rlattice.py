
from ccp4i2.core import CCP4Modules
from qtgui import CCP4TaskWidget
from ccp4i2.baselayer import QtCore

class CTaskDials_rlattice(CCP4TaskWidget.CTaskWidget):

    TASKNAME = 'dials_rlattice'
    TASKVERSION = 0.0
    TASKMODULE ='data_processing'
    TASKTITLE = 'DIALS Reciprocal Lattice Viewer'
    SHORTTASKTITLE = "DIALS rLattice Viewer"
    DESCRIPTION = 'DIALS Reciprocal Lattice Viewer'
    WHATNEXT = []

    def __init__(self,parent):
        CCP4TaskWidget.CTaskWidget.__init__(self,parent)

    def drawContents(self):
        self.setProgramHelpFile('shelxeMR')
        folder = self.openFolder(folderFunction='inputData', title='Input Data and Run Parameters')
        self.openSubFrame(title='Select input data', tip='A JSON file and Pickle file from DIALS are required', frame=[True])
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
                self.getWidget('PICKLE_IN').updateViewFromModel()
            else:
                self.container.inputData.PICKLE_IN.unSet()
        else:
            self.container.inputData.PICKLE_IN.unSet()
