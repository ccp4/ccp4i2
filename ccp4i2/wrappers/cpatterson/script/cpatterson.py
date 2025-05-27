"""
Copyright (C) 2014 The University of York, 2024 STFC
"""

from ....core.CCP4PluginScript import CPluginScript


class cpatterson(CPluginScript):

    TASKMODULE = 'wrappers'
    TASKTITLE = 'Prepare map coefficients'
    TASKNAME = 'cpatterson'
    TASKCOMMAND = 'cpatterson'
    TASKVERSION= 0.0
    MAINTAINER = 'stuart.mcnicholas@york.ac.uk'

    def processInputFiles ( self ):
        from ....core import CCP4XtalData
        self.hklin,error = self.makeHklin([['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]])

    def makeCommandAndScript(self):

        out = self.container.outputData

        self.appendCommandLine([ '-stdin' ])

        self.appendCommandScript([ '-mtzin ' + self.hklin ])
        self.appendCommandScript([ '-mapout %s'%(str(out.MAPOUT.fullPath)) ])
        self.appendCommandScript ([ '-colin-fo F,SIGF' ])

        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        self.container.outputData.MAPOUT.annotation = 'Computed using ' + str(self.container.inputData.F_SIGF.annotation)
        return CPluginScript.SUCCEEDED
