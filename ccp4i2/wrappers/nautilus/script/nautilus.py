import os

from ....core import CCP4ErrorHandling
from ....core import CCP4XtalData
from ....core.CCP4PluginScript import CPluginScript


class nautilus(CPluginScript):
    TASKMODULE="developer_tools"
    TASKNAME = 'nautilus'   # Task name - should be same as class name and match pluginTitle in the .def.xml file
    TASKVERSION= 0.1               # Version of this plugin
    MAINTAINER = 'kevin.cowtan@york.ac.uk'
    ERROR_CODES = { 201 : {'description' : 'Failed to analyse output files' } } 
    PURGESEARCHLIST = [ [ 'hklin.mtz' , 0 ], ['log_mtzjoin.txt', 0] ]
    TASKCOMMAND="cnautilus"

    def __init__(self, *args, **kws):
        super(nautilus, self).__init__(*args, **kws)

    def processInputFiles(self):
        #Preprocess reflections to generate an "HKLIN" file
        if self.container.inputData.FWT_PHWT_IN.isSet():
          if self.container.inputData.FREERFLAG.isSet():
            self.hklin,columns,error = self.makeHklin0([['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN], 'ABCD', 'FREERFLAG', 'FWT_PHWT_IN' ])
            print('FREERFLAG is set, so joining all data objects')
          else :
            print('FREERFLAG is not set, so joining the rest of the data objects')
            self.hklin,columns,error = self.makeHklin0([['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN], 'ABCD', 'FWT_PHWT_IN' ])

          if error.maxSeverity()>CCP4ErrorHandling.Severity.WARNING:
            print('ERROR creating input HKLIN with FWT_PHWT_IN')
            print(error.report())
            return CPluginScript.FAILED
        else:
          if self.container.inputData.FREERFLAG.isSet():
            self.hklin,columns,error = self.makeHklin0([['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN],'ABCD','FREERFLAG'])
            print('FREERFLAG is set, so joining all data objects')
          else :
            print('FREERFLAG is not set, so joining the rest of the data objects')
            self.hklin,columns,error = self.makeHklin0([['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN],'ABCD' ])
          
          if error.maxSeverity()>CCP4ErrorHandling.Severity.WARNING:
            print('ERROR creating input HKLIN')
            print(error.report())
            return CPluginScript.FAILED

        # convert the sequence
        self.seqin = os.path.join(self.getWorkDirectory(),'seqin.fasta')
        self.container.inputData.ASUIN.writeFasta(self.seqin, polymerTypes=["RNA", "DNA"])

        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self,**kw):
        self.appendCommandLine(['-stdin'])

        # INPUT DATA
        self.appendCommandScript("mtzin "+self.hklin)
        self.appendCommandScript("colin-fo F_SIGF_F,F_SIGF_SIGF")
        if self.container.inputData.ABCD.contentFlag == CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL:
          self.appendCommandScript("colin-hl ABCD_HLA,ABCD_HLB,ABCD_HLC,ABCD_HLD")
        else:
          self.appendCommandScript("colin-phifom ABCD_PHI,ABCD_FOM")
        if self.container.inputData.FREERFLAG.isSet( ) :
          self.appendCommandScript("colin-free FREERFLAG_FREER")
        if self.container.inputData.FWT_PHWT_IN.isSet( ) :
          self.appendCommandScript("colin-fc FWT_PHWT_IN_F,FWT_PHWT_IN_PHI")

        self.appendCommandScript("seqin %s"%(self.seqin))

        if self.container.inputData.XYZIN.isSet():
          self.appendCommandScript("pdbin %s"%(str(self.container.inputData.XYZIN.fullPath)))

        # OUTPUT DATA
        if self.container.outputData.XYZOUT.isSet():
          self.appendCommandScript("pdbout %s"%(str(self.container.outputData.XYZOUT.fullPath)))

        # CONTROL PARAMETERS
        if self.container.controlParameters.CYCLES.isSet():
          self.appendCommandScript("cycles %s"%(str(self.container.controlParameters.CYCLES)))
        if self.container.controlParameters.ANISOTROPY_CORRECTION:
          self.appendCommandScript("anisotropy-correction")
        if self.container.controlParameters.RESOLUTION.isSet():
          self.appendCommandScript("resolution %s"%(str(self.container.controlParameters.RESOLUTION)))
        else:
          self.appendCommandScript("resolution 2.0")  # I've added the default value here

        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
      self.container.outputData.XYZOUT.annotation = 'Model built with Nautilus'
      self.container.outputData.XYZOUT.subType=1
      nf,nr,lf,ncr = 0,0,0,0
      try:
        lastc, lastr = "", -999999
        for l in open(str(self.container.outputData.XYZOUT.fullPath)):
          if l[0:4] == "ATOM":
            if l[12:16] == " C1'":
              thisc, thisr = l[21:22], int(l[22:26])
              if thisc != lastc or thisr > lastr+1:
                nf += 1
                ncr = 0
              nr  += 1
              ncr += 1
              lf = max(lf,ncr)
              lastc,lastr = thisc,thisr
        f = open( self.makeFileName('PROGRAMXML'), "w" )
        f.write("<NautilusResult>\n <Final>\n")
        f.write("  <FragmentsBuilt>%d</FragmentsBuilt>\n"%nf)
        f.write("  <NucletidesBuilt>%d</NucletidesBuilt>\n"%nr)
        f.write("  <NucletidesLongestFragment>%d</NucletidesLongestFragment>\n"%lf)
        f.write(" </Final>\n</NautilusResult>\n")
        f.close()
      except Exception as e:
        print(str(e))
      return CPluginScript.SUCCEEDED
