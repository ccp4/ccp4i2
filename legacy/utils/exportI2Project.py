from __future__ import print_function
import sys
import os
import time
import functools

from baselayer import QtCore

if __name__ == "__main__":
    sys.path.append(os.path.join(os.path.dirname(__file__),".."))

class CompressClass(QtCore.QObject):

    doneSignal = QtCore.Signal()

    def __init__(self,parent=None,projectId=None,after=None,jobList=None,excludeI2files=False,fileName=None,projectName=None):
        QtCore.QObject.__init__(self,parent)
        from core import CCP4Modules
        from qtcore import CCP4Export
        print('compressProject',after,jobList,excludeI2files,fileName)

        projectInfo =  CCP4Modules.PROJECTSMANAGER().db().getProjectInfo(projectName=projectName)
        projectId = projectInfo['projectid']

        # If there is a limited set of jobs then find the input jobs that are not output by jobs on that list
        inputFilesList,inputFileIdList,fromJobList,errReport =  CCP4Modules.PROJECTSMANAGER().getJobInputFiles(projectDir=projectInfo['projectdirectory'],jobIdList=jobList,useDb=True,excludeI2files=excludeI2files)
        fromJobIdList = []
        fromJobNumberList = []
        for item in fromJobList:
          fromJobIdList.append(item[0])
          fromJobNumberList.append(item[1])
          
        dbxml = os.path.join( projectInfo['projectdirectory'],'CCP4_TMP','DATABASE'+str(int(time.time()))+'.db.xml')
        print('Creating XML database:'+dbxml)
        jobNumberList,errReport = CCP4Modules.PROJECTSMANAGER().db().exportProjectXml(projectId,fileName=dbxml,recordExport=True,status='exportable',after=after,jobList=jobList,inputFileList=inputFileIdList,inputFileFromJobList=fromJobIdList)

        if jobList is not None:
            directoriesList = []
        else:
            directoriesList = ['CCP4_IMPORTED_FILES','CCP4_PROJECT_FILES']

        self.setObjectName("ExportThreadObjectParent")
        self.exportThread = CCP4Export.ExportProjectThread(self,projectDir=projectInfo['projectdirectory'],dbxml=dbxml,target=fileName,jobList=jobNumberList,inputFilesList=inputFilesList,directoriesList=directoriesList,)
        @QtCore.Slot()
        def updateSavingJobData():
            print("updateSavingJobData")
        @QtCore.Slot()
        def progressSavingJobData():
            print("progressSavingJobData")
        @QtCore.Slot(str)
        def doneSavingJobData(name,fname):
            print("Exported project",name,"to",fname)
            self.doneSignal.emit()
        self.exportThread.savingJobData.connect(updateSavingJobData)
        self.exportThread.startSavingJobData.connect(progressSavingJobData)
        self.exportThread.finished.connect(functools.partial(doneSavingJobData,projectInfo['projectname'],fileName))

    def run(self):
        self.exportThread.start()

def startDb(parent=None, fileName=None, mode='sqlite', userName=None, userPassword=None,**kw):
    from core import CCP4Utils
    from core import CCP4ErrorHandling
    from dbapi import CCP4DbApi
    from qtgui import CCP4DbManagerGui
    db = CCP4DbApi.CDbApi(parent=parent, fileName=fileName, mode=mode, createDb=True, userName=userName, userPassword=userPassword, loadDiagnostic=kw.get('loadDiagnostic',True))
    return db

if __name__ == "__main__":

    from core import CCP4Config
    from utils.QApp import CGuiApplication # KJS : does this need to go in here (is fn used ?)
    app = CGuiApplication(sys.argv)
    CCP4Config.CONFIG().set('graphical', False)
    CCP4Config.CONFIG().set('qt', True)
    from core import CCP4Modules
    from core import CCP4Utils

    pm = CCP4Modules.PROJECTSMANAGER()
    db = startDb(pm, mode='sqlite')
    pm.setDatabase(db)

    #projects = CCP4Modules.PROJECTSMANAGER().db().getProjectDirectoryList()
    projects = sys.argv[1].split(",")

    @QtCore.Slot(dict)
    def exportNextProject(**args):
        try:
            name = projects.pop()
            compClass = CompressClass(projectName=name,fileName=name+".ccp4_project.zip")
            compClass.doneSignal.connect(exportNextProject)
            compClass.run()
        except:
            app.quit()

    exportNextProject()

    sys.exit(app.exec_())
