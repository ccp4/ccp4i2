from __future__ import print_function

"""
     ctruncate.py: CCP4 GUI Project
     Copyright (C) 2012 STFC
"""

from core.CCP4PluginScript import CPluginScript
from core.CCP4ErrorHandling import *
class ctruncate(CPluginScript):

    TASKMODULE = 'expt_data_utility'      # Where this plugin will appear on the gui
    TASKTITLE = 'Intensities to amplitudes' # A short title for gui menu
    DESCRIPTION = 'Convert reflection intensities to structure factors (ctruncate)'
    TASKNAME = 'ctruncate'   # Task name - should be same as class name
    TASKVERSION= 0.0               # Version of this plugin

    # used by the base class startProcess()
    TASKCOMMAND = 'ctruncate'   # The command to run the executable
    # used by the base class makeCommandAndScript()
    COMLINETEMPLATE = None 
    COMTEMPLATE = None
    MAINTAINER = 'charles.ballard@stfc.ac.uk'

    ERROR_CODES = { 201 : { 'severity' : SEVERITY_WARNING , 'description' : 'Error creating XML output' } }


    def makeCommandAndScript(self):

      par = self.container.controlParameters
      inp = self.container.inputData
      self.ifAnom = False

      if inp.OBSIN.isSet() :
        #Input from gui - need to poulate HKLIN,FSIGF etc as expected from scripted input
        # OBSIN is a mini-MTZ so columns are know and dependent on contentFlag
        inp.HKLIN.set(inp.OBSIN)
        inp.OBSIN.setContentFlag(reset=True)
        par.AMPLITUDES = int(inp.OBSIN.contentFlag) in (2,4)
        if  int(inp.OBSIN.contentFlag) == 1:
          inp.ISIGIanom.set({ 'Ip' : 'Iplus', 'SIGIp': 'SIGIplus', 'Im' :'Iminus', 'SIGIm' :'SIGIminus'}  )
        elif  int(inp.OBSIN.contentFlag) == 2:
          inp.FSIGFanom.set ( { 'Fp' :'Fplus', 'SIGFp' : 'SIGFplus', 'Fm' :'Fminus', 'SIGFm' :'SIGFminus'}  )
        elif  int(inp.OBSIN.contentFlag) == 3:
          inp.ISIGI.set( { 'I' : 'I', 'SIGI' : 'SIGI' } )
        elif  int(inp.OBSIN.contentFlag) == 4:
          inp.FSIGF.set( { 'F': 'F' , 'SIGF' : 'SIGF' } )
        
      ### files
      self.appendCommandLine(['-hklin',inp.HKLIN.fullPath])
      if inp.SEQIN.isSet():
         self.appendCommandLine(['-seqin',inp.SEQIN.fullPath])
      #print 'ctruncate.makeCommandAndScript OUTPUTMINIMTZ',self.container.controlParameters.OUTPUTMINIMTZ,type(self.container.controlParameters.OUTPUTMINIMTZ)
      #if self.container.controlParameters.OUTPUTMINIMTZ:
        #import os
        #self.tmpHklout = os.path.join(self.workDirectory,'ctruncate_output.mtz')
        #self.appendCommandLine(['-hklout',self.tmpHklout])
       # else:
      self.appendCommandLine(['-hklout',self.container.outputData.HKLOUT.fullPath])

      ### column assignments
      #print 'ctruncate.makeCommandAndScript FSIGF.isSet',inp.FSIGF.isSet()
      self.colinIMEAN = False  # true if input file has an IMEAN column (eg from Aimless)
      if par.AMPLITUDES:
         if inp.FSIGFanom.isSet():
           colarg = '/*/*/['+str(inp.FSIGFanom.Fp)+','+str(inp.FSIGFanom.SIGFp)+','+str(inp.FSIGFanom.Fm)+','+str(inp.FSIGFanom.SIGFm)+']'
           self.appendCommandLine(['-colano',colarg])
           self.ifAnom = True
         if inp.FSIGF.isSet():
           colarg = '/*/*/['+str(inp.FSIGF.F)+','+str(inp.FSIGF.SIGF)+']'
           self.appendCommandLine(['-colin',colarg])
      else:
         if inp.ISIGIanom.isSet():
           colarg = '/*/*/['+str(inp.ISIGIanom.Ip)+','+str(inp.ISIGIanom.SIGIp)+','+str(inp.ISIGIanom.Im)+','+str(inp.ISIGIanom.SIGIm)+']'
           self.appendCommandLine(['-colano',colarg])
           self.ifAnom = True
         if inp.ISIGI.isSet():
           colarg = '/*/*/['+str(inp.ISIGI.I)+','+str(inp.ISIGI.SIGI)+']'
           self.appendCommandLine(['-colin',colarg])
           if str(inp.ISIGI.I) == 'IMEAN':
               self.colinIMEAN = True

      self.appendCommandLine(['-xmlout',self.makeFileName('PROGRAMXML')])

      if par.NRES.isSet():
         self.appendCommandLine(['-nres',par.NRES])
      if par.NO_ANISO:
         self.appendCommandLine(['-no-aniso'])
      if par.AMPLITUDES:
         self.appendCommandLine(['-amplitudes'])
      if par.OUTPUT_INTENSITIES:
         if not self.colinIMEAN:
             # only if input file does not already have an IMEAN column
             self.appendCommandLine(['-Imean'])

      return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
      import os,shutil
      #print 'ctruncate.processOutputFiles',self.container.controlParameters.OUTPUTMINIMTZ,self.container.controlParameters.OUTPUTMINIMTZCONTENTFLAG
      #print 'ctruncate.processOutputFiles HKLOUT',self.container.outputData.HKLOUT.__str__(),os.path.exists(self.container.outputData.HKLOUT.__str__())
              
      if self.container.controlParameters.OUTPUTMINIMTZ:
        # Output one miniMtz of content type determined by OUTPUTMINIMTZCONTENTFLAG - beware column names not i2 standard
        # outputContent should be a CCP4XtalData.CONTENT_FLAG_? value
        if self.container.controlParameters.OUTPUTMINIMTZCONTENTFLAG > 0:
          outputContent = int(self.container.controlParameters.OUTPUTMINIMTZCONTENTFLAG)
        else:
          outputContent = 1
          if not self.ifAnom: outputContent += 2
          if not self.container.controlParameters.OUTPUT_INTENSITIES: outputContent += 1
        #print 'ctruncate.processOutputFiles',outputContent
          
        logFile = os.path.join(self.workDirectory,'cmtzsplit.log')
        # ***** Check ctruncate column names
        from core import CCP4XtalData
        if outputContent == 4:
          # MN Kludge here..*FIXME*.looks to me like the column labels output by ctruncate for Fmean, SIGFmean
          # depend on the type of data it started with (ISIGI, vs ISIGIanom)
          if self.container.inputData.ISIGI.isSet(): colin = 'F,SIGF'
          if self.container.inputData.ISIGIanom.isSet(): colin = 'FMEAN,SIGFMEAN'
          colout = 'F,SIGF'
        elif outputContent == 3:
          if self.container.inputData.ISIGI.isSet(): colin = 'IMEAN,SIGIMEAN'
          if self.container.inputData.ISIGIanom.isSet(): colin = '_MEAN.I_sigI.I,_MEAN.I_sigI.sigI'
          colout = 'I,SIGI'
        elif outputContent == 2:
          colin = 'F(+),SIGF(+),F(-),SIGF(-)'
          colout = 'Fplus,SIGFplus,Fminus,SIGFminus'
        elif outputContent == 1:
          colin = 'I(+),SIGI(+),I(-),SIGI(-)'
          colout = 'Iplus,SIGIplus,Iminus,SIGIminus'
        else:
          colin = CCP4XtalData.CObsDataFile.CONTENT_SIGNATURE_LIST[outputContent-1]
        outputLst = [[self.container.outputData.OBSOUT.__str__(),colin,colout]]
        #print('\n*processOutputFiles OUTPUTMINIMTZ',outputLst,outputContent)
        if outputContent == 1:
          self.container.outputData.OBSOUT1.set(os.path.splitext(str(self.container.outputData.OBSOUT))[0]+'_asIMEAN.mtz')
          if self.colinIMEAN:
              outputLst.append ([ str(self.container.outputData.OBSOUT1) ,
                                  'IMEAN,SIGIMEAN' , 'I,SIGI' ])
          else:
              outputLst.append ([ str(self.container.outputData.OBSOUT1) ,
                                  '_MEAN.I_sigI.I,_MEAN.I_sigI.sigI' , 'I,SIGI' ])
          #print("\n*outputContent == 1, outputLst", outputLst)
        elif outputContent == 2:
          # The usual output of Intensities to SFs tasks
          self.container.outputData.OBSOUT1.set(os.path.splitext(str(self.container.outputData.OBSOUT))[0]+'_asFMEAN.mtz')
          self.container.outputData.OBSOUT1.contentFlag.set(4)
          outputLst.append ([ str(self.container.outputData.OBSOUT1) , 'FMEAN,SIGFMEAN' , 'F,SIGF' ])
        #print('\n**splitMtz output',outputLst)
        status = self.splitMtz(str(self.container.outputData.HKLOUT),outputLst,logFile)
        self.container.outputData.OBSOUT.contentFlag = outputContent
        dName = self.container.outputData.OBSOUT.datasetName()
        if dName == '': dName = 'Reflections'
        self.container.outputData.OBSOUT.annotation.set(dName +' as '+self.container.outputData.OBSOUT.CONTENT_ANNOTATION[outputContent-1])
        self.container.outputData.OBSOUT.contentFlag.set(outputContent)
        if outputContent == 1:
          self.container.outputData.OBSOUT1.contentFlag.set(3)
        elif outputContent == 2:
          self.container.outputData.OBSOUT1.annotation.set(dName + ' as '+self.container.outputData.OBSOUT.CONTENT_ANNOTATION[3])
          self.container.outputData.OBSOUT1.contentFlag.set(4)
        #print('\n***ctruncate.processOutputFiles after splitMtz status',status,'contentFlag',self.container.outputData.OBSOUT.contentFlag)

        '''
        if status != CPluginScript.SUCCEEDED: return status
        if os.path.exists(self.container.outputData.HKLOUT.__str__()):
          bakup,ext = os.path.splitext(self.container.outputData.HKLOUT.__str__())
          bakup = bakup + '_bak'+ext
          shutil.move(self.container.outputData.HKLOUT.__str__(),bakup)
          print 'bakup',bakup,os.path.exists(bakup)
        try:
          shutil.move(self.tmpHklfile,self.container.outputData.HKLOUT.__str__())
        except:
          print 'Failed ctruncATE.processOutputFiles'
          return CPluginScript.FAILED
        '''  

      return CPluginScript.SUCCEEDED

    
        

#======================================================
# PLUGIN TESTS
# See Python documentation on unittest module

import unittest

class testctruncate(unittest.TestCase):

   def setUp(self):
    from core import CCP4Modules
    self.app = CCP4Modules.QTAPPLICATION()
    # make all background jobs wait for completion
    # this is essential for unittest to work
    CCP4Modules.PROCESSMANAGER().setWaitForFinished(10000)

   def tearDown(self):
    from core import CCP4Modules
    CCP4Modules.PROCESSMANAGER().setWaitForFinished(-1)

   def test_1(self):
     from core import CCP4Modules
     from core import CCP4Utils
     import os

     workDirectory = os.path.join(CCP4Utils.getTestTmpDir(),'test1')
     if not os.path.exists(workDirectory): os.mkdir(workDirectory)

     self.wrapper = ctruncate(parent=CCP4Modules.QTAPPLICATION(),name='test1',workDirectory=workDirectory)
     self.wrapper.container.loadDataFromXml(os.path.join(CCP4Utils.getCCP4I2Dir(),'wrappers','ctruncate','test_data','test1.data.xml'))

     self.wrapper.setWaitForFinished(1000000)
     pid = self.wrapper.process()
     self.wrapper.setWaitForFinished(-1)
     if len(self.wrapper.errorReport)>0: print(self.wrapper.errorReport.report())

def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testctruncate)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
