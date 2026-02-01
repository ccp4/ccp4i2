from ccp4i2.wrappers.x2mtz.script import x2mtz


class cif2mtz(x2mtz.x2mtz):
    TASKMODULE = 'test'
    TASKTITLE = 'Import mmCIF reflection file'
    TASKNAME = 'cif2mtz'
    TASKCOMMAND = 'cif2mtz'

    ERROR_CODES = { 301 : { 'description' : 'No output file found after cif2mtz conversion' },
                    302 : { 'description' : 'Output from cif2mtz does not contain recognised reflection or FreeR set data' }
                    }

    def makeCommandAndScript(self):
      inputData = self.container.inputData

      self.appendCommandLine(['HKLIN',inputData.HKLIN.fullPath])
      self.appendCommandLine(['HKLOUT',self.container.outputData.HKLOUT.fullPath])
      print('cif2mtz.makeCommandAndScript',inputData.SPACEGROUPCELL.spaceGroup,inputData.SPACEGROUPCELL.spaceGroup.isSet())
      if inputData.SPACEGROUPCELL.spaceGroup.isSet():
          self.appendCommandScript("SYMM '%s'"%str(inputData.SPACEGROUPCELL.spaceGroup))
      if inputData.SPACEGROUPCELL.cell.isSet():
          # angles default to 90 etc..
          inputData.SPACEGROUPCELL.cell.set(inputData.SPACEGROUPCELL.cell.fix(inputData.SPACEGROUPCELL.cell.get()))
          self.appendCommandScript("CELL %f %f %f %f %f %f"%
              (inputData.SPACEGROUPCELL.cell.a.__float__(),inputData.SPACEGROUPCELL.cell.b.__float__(),inputData.SPACEGROUPCELL.cell.c.__float__(),
              inputData.SPACEGROUPCELL.cell.alpha.__float__(),inputData.SPACEGROUPCELL.cell.beta.__float__(),inputData.SPACEGROUPCELL.cell.gamma))
      if inputData.CRYSTALNAME.isSet():
          self.appendCommandScript("XNAME '%s'"%str(inputData.CRYSTALNAME))
      if inputData.DATASETNAME.isSet():
          self.appendCommandScript("DNAME '%s'"%str(inputData.DATASETNAME))
      if inputData.MMCIF_SELECTED_BLOCK.isSet():
          self.appendCommandScript("BLOCK '%s'"%str(inputData.MMCIF_SELECTED_BLOCK))
      self.appendCommandScript('END')
      return 0
