import os

import gemmi
from lxml import etree

from ccp4i2.core import CCP4ErrorHandling, CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript


class freerflag(CPluginScript):
    TASKMODULE = 'test'      # Where this plugin will appear on the gui
    TASKTITLE = 'Add a freeR flag' # A short title for gui menu
    TASKNAME = 'freerflag'   # Task name - should be same as class name
    TASKVERSION= 0.0               # Version of this plugin
    MAINTAINER = 'liz.potterton@york.ac.uk'
    TASKCOMMAND = 'freerflag'   # The command to run the executable
    COMLINETEMPLATE = None 
    COMTEMPLATE = None  

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def processInputFiles(self):
        print("freerflag wrapper, processInputFiles")

        # We have two possible resolution cutoffs
        #  1) global cutoff if RESMAX is set
        #  2) if COMPLETE and CUTRESOLUTION and FreeR is higher resolution
        #     than the data, then FreeR is cut to the data resolution
        self.freerCutoff = 0.0   # will be set >0 if Freer cutoff is done
        self.globalCutoff = 0.0
       
        self.highResD = float(self.container.inputData.F_SIGF.fileContent.resolutionRange.high)
        self.highResF = 0.0
        if self.container.inputData.FREERFLAG.isSet():
            self.highResF = float(self.container.inputData.FREERFLAG.fileContent.resolutionRange.high)

        if self.container.controlParameters.GEN_MODE == 'COMPLETE':
            if self.container.controlParameters.CUTRESOLUTION:
                #  cut the resolution of the FreeR set if it is higher than the data
                self.cutResolution()
        
            # ignoreErrorCodes to say makeHklin can ignore cmtzjoin gives exitCode 101 for incomplete freerflag
            self.hklin,error = self.makeHklin(['F_SIGF','FREERFLAG'],ignoreErrorCodes=[36])
            print('freerflag.processInputFiles',self.hklin,error)
            if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
                return CPluginScript.FAILED
            else:
                # Optional global cutoff
                if self.container.controlParameters.RESMAX:
                    resmax = self.container.controlParameters.RESMAX.get()
                    self.globalResolutionCutoff(resmax)

                return CPluginScript.SUCCEEDED

        else:
            self.hklin = self.container.inputData.F_SIGF.__str__()

        return CPluginScript.SUCCEEDED

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def cutResolution(self):

        if self.highResD <= self.highResF:
            # Data go to higher resolution, no need to do anything
            return

        # File names
        FSIGF_file = str(self.container.inputData.F_SIGF)
        FREER_file = str(self.container.inputData.FREERFLAG)

        # Cut resolution of FreeR file to match data file
        outfile = os.path.join(self.workDirectory, 'FREEOUT.mtz')        
        highRes= self.highResD-0.001
        print("* Cutting Freer to ", highRes)
        mtz = gemmi.read_mtz_file(FREER_file)
        mtz.set_data(mtz.array[mtz.make_d_array() >= highRes])
        mtz.write_to_file(outfile)
        self.container.inputData.FREERFLAG.setFullPath(outfile)
        self.freerCutoff = highRes
    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def globalResolutionCutoff(self, resmax):
        # after mtzjoin, if needed. Joined file is in self.hklin
        mtz = gemmi.read_mtz_file(self.hklin)
        highRes = mtz.resolution_high()
        if highRes < resmax-0.001:
            # yes cut the data
            print(">>**>> cutting all data", highRes, resmax)
            mtz.set_data(mtz.array[mtz.make_d_array() >= resmax])
            mtz.write_to_file(self.hklin)
            self.globalCutoff = resmax

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def makeCommandAndScript(self):
      self.hklout = os.path.join(self.workDirectory,"hklout.mtz")

      self.appendCommandLine(['HKLIN', self.hklin])
      self.appendCommandLine(['HKLOUT', self.hklout])

      if self.container.controlParameters.GEN_MODE == 'COMPLETE':
          if self.container.inputData.FREERFLAG.isSet():
            self.appendCommandScript("COMPLETE FREE=FREER")
          else:
            self.appendErrorReport(101)
      else:
          # FREERFLAG keyword only applies if not completing existing freeR set
          if self.container.controlParameters.FRAC.isSet():
              self.appendCommandScript("FREERFRAC %s"%(str(self.container.controlParameters.FRAC)))

      if self.container.controlParameters.UNIQUEIFY:
          self.appendCommandScript("UNIQUE")
          if self.container.controlParameters.RESMAX:
              self.appendCommandScript("RESOL %s"%(str(self.container.controlParameters.RESMAX)))

      self.appendCommandScript('END')

      return 0

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def processOutputFiles(self):
      print("freerflag - processOutputFiles")
      self.xmlRoot = etree.Element('FREERFLAGINFO')
      annotation = ''

      if self.container.controlParameters.GEN_MODE == 'COMPLETE':
         annotation = 'Extended freeR;'
         self.addElement(self.xmlRoot, 'Mode', 'Complete')
         if self.container.controlParameters.CUTRESOLUTION:
             self.addElement(self.xmlRoot, 'CutFreerResolution', 'True')
         else:
             self.addElement(self.xmlRoot, 'CutFreerResolution', 'False')

         self.addElement(self.xmlRoot, 'ObservedDataResolution',
                         '{:6.2f}'.format(self.highResD))
         self.addElement(self.xmlRoot, 'FreeR_Resolution',
                         '{:6.2f}'.format(self.highResF))

         if self.freerCutoff > 0.0:
             self.addElement(self.xmlRoot, 'FreerCutResolution',
                             '{:6.2f}'.format(self.freerCutoff))
      else:
          annotation = 'New freeR'
          self.addElement(self.xmlRoot, 'Mode', 'New')
          fraction = 0.05  # default
          if self.container.controlParameters.FRAC.isSet():
              fraction = self.container.controlParameters.FRAC
          self.addElement(self.xmlRoot, 'Fraction', str(fraction))

      if self.container.controlParameters.UNIQUEIFY:
         self.addElement(self.xmlRoot, 'Unique', 'True')

      if self.container.controlParameters.RESMAX:
          resmax = self.globalCutoff
          if resmax > 0.0:
              self.addElement(self.xmlRoot, 'GlobalResolutionLimit',
                              '{:6.2f}'.format(resmax))

      with open ( self.makeFileName('PROGRAMXML'),'w' ) as xmlFile:
         xmlString = etree.tostring (self.xmlRoot, pretty_print=True )
         CCP4Utils.writeXML(xmlFile,xmlString)

      if self.container.controlParameters.GEN_MODE == 'COMPLETE':
          error = self.splitHklout(['FREEROUT'],['FREER'])
      else:          
          error = self.splitHklout(['FREEROUT'],['FreeR_flag'])

      if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
        return CPluginScript.FAILED
     
      # Annotation cf aimless_pipe
      if self.container.outputData.FREEROUT.fileContent.spaceGroup.isSet():
          sgname = self.container.outputData.FREEROUT.fileContent.spaceGroup.__str__()
      else:
          sgname = 'Unk'

      highresFRformatted = "%7.2f" % float(self.container.outputData.FREEROUT.fileContent.resolutionRange.high)
      title =' Spg:'+str(sgname).strip()+';Resln:'+highresFRformatted.strip() + "A;"
      try:
          title = title + "Cell:"+self.container.outputData.FREEROUT.fileContent.cell.guiLabel()
      except Exception as e:
          print('Error writing cell parameters',e)

      annotation += title
      print("Annotation", annotation)
      self.container.outputData.FREEROUT.annotation.set(annotation)

      return CPluginScript.SUCCEEDED

    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def getXML(self):
        try:
            return self.xmlRoot
        except:
            return None
    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def addElement(self, containerXML, elementname, elementtext):
        #print 'addElement', elementname, type(elementtext), elementtext 
        e2 = etree.Element(elementname)
        e2.text = elementtext
        containerXML.append(e2)
