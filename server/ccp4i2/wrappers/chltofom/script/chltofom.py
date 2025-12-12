from __future__ import print_function


from ccp4i2.core import CCP4PluginScript
from ccp4i2.core import CCP4XtalData

  
class chltofom(CCP4PluginScript.CPluginScript):


    TASKNAME = 'chltofom'                                  # Task name - should be same as class name
    TASKCOMMAND = 'chltofom'                             # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    MAINTAINER = 'liz.potterton@york.ac.uk'
    PURGESEARCHLIST = [  [ 'hklout.mtz' , 0 ],    
                         [ 'HKLOUT.mtz' , 0 ]
                       ]
    
    def makeCommandAndScript(self):
      self.appendCommandLine( ['-mtzin',self.container.inputData.HKLIN.__str__() ] )
      if self.container.inputData.HKLIN.contentFlag == CCP4XtalData.CPhsDataFile.CONTENT_FLAG_PHIFOM:
        self.container.controlParameters.DIRECTION.set('FOMTOHL') 
        self.appendCommandLine ( ['-colin-phifom', '/*/*/[PHI,FOM]' ])
        self.appendCommandLine ( ['-colout', '/*/*/[HLA,HLB,HLC,HLD]' ])
      else:
        self.appendCommandLine ( ['-colin-hl', '/*/*/[HLA,HLB,HLC,HLD]' ])
        self.appendCommandLine ( ['-colout','/*/*/[PHI,FOM]' ])
      if self.container.controlParameters.OUTPUTMINIMTZ:
          from ccp4i2.core import CCP4Utils
          self.tmpHklout = CCP4Utils.makeTmpFile(extension='mtz')
          self.appendCommandLine( [ '-mtzout',self.tmpHklout] )
      else:
        self.appendCommandLine( [ '-mtzout',self.container.outputData.HKLOUT.__str__()] )
      return CCP4PluginScript.CPluginScript.SUCCEEDED

    def processOutputFiles(self):
      #print 'chltofom.processOutputFiles',self.container.controlParameters.OUTPUTMINIMTZ,self.__dict__.get('tmpHklout','NONE')
      if self.container.controlParameters.OUTPUTMINIMTZ:
        import os
        logFile = os.path.splitext(self.tmpHklout)[0]+'.log'

        if self.container.inputData.HKLIN.annotation.isSet():
          anno = str(self.container.inputData.HKLIN.annotation)
        else:
          anno = self.container.outputData.HKLOUT.qualifiers('guiLabel')
       
        if self.container.controlParameters.DIRECTION.isSet() and self.container.controlParameters.DIRECTION.__str__() == 'FOMTOHL':
          status = self.splitMtz(self.tmpHklout,[[self.container.outputData.HKLOUT.__str__(),'HLA,HLB,HLC,HLD']],logFile)
          self.container.outputData.HKLOUT.annotation = anno + ' as HLcoeffs'
        else:
          status = self.splitMtz(self.tmpHklout,[[self.container.outputData.HKLOUT.__str__(),'PHI,FOM']],logFile)
          self.container.outputData.HKLOUT.annotation = anno  + ' as phi/FOM'
          self.container.outputData.HKLOUT.contentFlag = CCP4XtalData.CPhsDataFile.CONTENT_FLAG_PHIFOM
        if status != CCP4PluginScript.CPluginScript.SUCCEEDED: return status
        if not os.path.exists(self.container.outputData.HKLOUT.__str__()):
          try:
            shutil.copyfile(self.tmpHklfile,self.container.outputData.HKLOUT.__str__())
          except:
            pass
          return CCP4PluginScript.CPluginScript.FAILED
      return CCP4PluginScript.CPluginScript.SUCCEEDED
