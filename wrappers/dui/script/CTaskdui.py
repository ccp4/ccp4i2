from pathlib import Path
from core.CCP4Modules import PROJECTSMANAGER, TASKMANAGER
from dbapi.CCP4DbApi import JOB_STATUS_FINISHED
from qtgui.CCP4TaskWidget import CTaskWidget

class CTaskdui(CTaskWidget):
    TASKMODULE = 'data_processing'
    TASKNAME = 'dui'
    TASKVERSION = 0.01
    SHORTTASKTITLE = 'DIALS User Interface'
    TASKTITLE = 'Integrate images with DIALS'
    DESCRIPTION = 'Launch DUI and capture output'
    WHATNEXT = ['aimless_pipe']

    def drawContents(self):
        self.setProgramHelpFile('dui')
        self.openFolder(folderFunction='inputData', title='Input Data and Run Parameters', followFrom=False)
        self.createLine(['subtitle', 'Start the DIALS User Interface (DUI) and capture data on output'])

        menuText = ["Do not continue from a previous DUI session"]
        enumerators = [""]
        jobList = PROJECTSMANAGER().db().getProjectJobListInfo(
            projectId=self.projectId(),
            jobStatus=[JOB_STATUS_FINISHED],
            mode=['jobid', 'jobnumber', 'taskname', 'jobtitle'],
            topLevelOnly=True,
            order='DESC',
        )
        for job in jobList:
            if job["taskname"] == "dui":
                jobDir = PROJECTSMANAGER().db().jobDirectory(jobId=job["jobid"])
                path = Path(jobDir, "run_dui2_nodes", "run_data")
                if path.exists():
                    title = job.get("jobtitle") or TASKMANAGER().getTitle(job["taskname"])
                    text = f'{job["jobnumber"]} {title}'
                    menuText.append(text)
                    enumerators.append(str(path))

        qualifiers = {'enumerators': enumerators, 'menuText': menuText}
        self.container.inputData.DUI2_RUN_DATA.setQualifiers(qualifiers)
        self.container.inputData.DUI2_RUN_DATA.set("")
        self.createLine(['label', ''])
        self.createLine(['label', 'If you wish to continue from a previous DUI session from this project'])
        self.createLine(['label', 'then please select the job below'])
        self.createLine(['label', ''])
        self.createLine(['widget', 'DUI2_RUN_DATA'])
        self.createLine(['label', ''])
        self.createLine(['label', 'Click "Run" to start'])
        self.createLine(['label', 'When you have finished integration in DIALS, exit DUI'])
        self.createLine(['label', 'Then use the "Data reduction" follow-on task'])
