from ccp4i2.wrappers.x2mtz.script import x2mtz


class scalepack2mtz(x2mtz.x2mtz):

    TASKMODULE = 'test'
    TASKTITLE = 'Convert scalepack merged reflection file to MTZ'
    TASKNAME = 'scalepack2mtz'  
    TASKVERSION= 0.1    
    TASKCOMMAND = 'scalepack2mtz'  

    def makeCommandAndScript(self):

      inp = self.container.inputData
      # input parameters
      self.appendCommandLine(["HKLIN", inp.HKLIN.fullPath])

      if self.container.inputData.SPACEGROUPCELL.cell.isSet():
         inp.SPACEGROUPCELL.cell.set(inp.SPACEGROUPCELL.cell.fix(inp.SPACEGROUPCELL.cell.get()))
         self.appendCommandScript("CELL %f %f %f %f %f %f" %
           (inp.SPACEGROUPCELL.cell.a.__float__(), inp.SPACEGROUPCELL.cell.b.__float__(), inp.SPACEGROUPCELL.cell.c.__float__(),
            inp.SPACEGROUPCELL.cell.alpha.__float__(), inp.SPACEGROUPCELL.cell.beta.__float__(), inp.SPACEGROUPCELL.cell.gamma.__float__()))
      if inp.SPACEGROUPCELL.spaceGroup.isSet():
         self.appendCommandScript("SYMMETRY %d" % inp.SPACEGROUPCELL.spaceGroup.number())     
      if inp.WAVELENGTH.isSet():
          self.appendCommandScript("WAVE %f" % inp.WAVELENGTH.__float__() )                                                
      if inp.CRYSTALNAME.isSet():self.appendCommandScript("XNAME '%s'"%str(inp.CRYSTALNAME))
      if inp.DATASETNAME.isSet():self.appendCommandScript("DNAME '%s'"%str(inp.DATASETNAME))

      # control parameters
      if self.container.controlParameters.ANOMALOUS:
          self.appendCommandScript("anomalous")

      s = self.resolutionRangeCommand()
      if s != "":
          self.appendCommandScript(s)
          
      # output parameters
      if self.container.outputData.HKLOUT.isSet():
         self.appendCommandLine(["HKLOUT", self.container.outputData.HKLOUT.fullPath])
      self.appendCommandScript('END')

      return 0

    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def resolutionRangeCommand(self):
        r1 = self.container.inputData.RESOLUTION_RANGE.start
        r2 = self.container.inputData.RESOLUTION_RANGE.end
        print("S2M resolutionRangeCommand",
              self.container.inputData.RESOLUTION_RANGE,
              r1, r2)
        
        high = 0.0
        low  = 99999.0
        if not r1.isSet() and not r2.isSet():
            # nothing set
            s = ""
        else:
            if r1.isSet():
                low = float(r1)
            if r2.isSet():
                high = float(r2)
            if low < high:
                low, high = high, low
            s = "RESOLUTION %f  %f" % (low, high)

        return s
