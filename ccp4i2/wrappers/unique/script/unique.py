"""
     unique.py: CCP4 GUI Project
     Author: Martyn Winn
     Copyright (C) 2011 STFC
"""

from core.CCP4PluginScript import CPluginScript

class unique(CPluginScript):

    TASKTITLE = 'Create dummy dataset' # A short title for gui menu
    TASKNAME = 'unique'   # Task name - should be same as class name
    TASKVERSION= 0.1               # Version of this plugin

    # used by the base class startProcess()
    TASKCOMMAND = 'unique'   # The command to run the executable
    # used by the base class makeCommandAndScript()
    COMLINETEMPLATE = '''1 HKLOUT $HKLOUT'''
    COMTEMPLATE = '''1 SYMM $SPACEGROUPCELL.spaceGroup.quote
1 CELL $SPACEGROUPCELL.cell.a $SPACEGROUPCELL.cell.b $SPACEGROUPCELL.cell.c $SPACEGROUPCELL.cell.alpha $SPACEGROUPCELL.cell.beta $SPACEGROUPCELL.cell.gamma  
1 RESO $RESOLUTION
$TITLE TITLE $TITLE
1 END'''

    def postProcess(self,processId=-1,data={}):
      rv,exitStatus,exitCode = self.postProcessCheck(processId)
      if rv>0:
        self.reportStatus(rv)
        return

      # sanity check that unique has produced something
      # -1 means we have not managed to get value out of log file
      self.container.outputData.NO_REFLECTIONS = -1
      logText = self.logFileText()
      pyListLogLines = logText.split("\n")
      for j, pyStrLine in enumerate(pyListLogLines):
         if "reflections within resolution limits" in pyStrLine:
            num_ref = int(pyStrLine.split()[0])
            if num_ref < 1:
               rv = 1
            else:
               self.container.outputData.NO_REFLECTIONS  = num_ref

      # Copy params, including quantities derived here in post-processing, to
      # the PROGRAMXML file which .xrt expects.
      self.container.saveDataToXml(self.makeFileName('PROGRAMXML'))

      self.reportStatus(rv)

