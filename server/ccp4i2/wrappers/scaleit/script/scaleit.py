import os

from lxml import etree

from ccp4i2.core import CCP4XtalData
from ccp4i2.core.CCP4ErrorHandling import CErrorReport
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.wrappers.scaleit.script.scaleit_logscraper import addElement, scaleitLogScraper
from ccp4i2.wrappers.scaleit.script.scaleit_utils import DatalistCheck


class scaleit(CPluginScript):
    TASKMODULE = 'test'
    TASKTITLE = 'Compare two or more datafiles'
    TASKNAME = 'scaleit'
    MAINTAINER = 'pre@mrc-lmb.cam.ac.uk'
    TASKCOMMAND = 'scaleit'

    ERROR_CODES = {  200 : { 'description' : 'Too few datasets defined' },
                     201 : { 'description' : 'Inconsistent cells between datasets'}
                     }

    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def processInputFiles(self):
      print("scaleit wrapper, processInputFiles")

      if len(self.container.inputData.MERGEDFILES) < 2:
          self.appendErrorReport(200)
          self.makeFailXML('Less than two datasets specified')
          self.reportStatus(CPluginScript.FAILED)
          return CPluginScript.FAILED

      # Check compatibility
      datalistcheck =  DatalistCheck(self.container.inputData.MERGEDFILES)
      valid = datalistcheck.checkAll()
      if not valid:
          print("valid fail")
          failmessage = datalistcheck.formatFail()
          self.appendErrorReport(201)
          self.makeFailXML(failmessage)
          self.reportStatus(CPluginScript.FAILED)
          return CPluginScript.FAILED

      mtzNameBase = 'MF'
      self.hklin, self.coloutlist, error = \
                  self.myMakeHklInput(self.container.inputData.MERGEDFILES,
                                      mtzNameBase,
                               CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN)
      print("\n**** self.coloutlist", self.coloutlist)

      return CPluginScript.SUCCEEDED

    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def myMakeHklInput(self, miniMtzList, mtzNameBase, targetContent=None,
                       hklin='hklin', ignoreErrorCodes=[],
                       extendOutputColnames=True, useInputColnames=False):
        # Modfied from CCP4PluginScript to take files in a CList PRE Oct 2022
        #  miniMtzList is CMiniMtzDataFileList of CObsDataFile
        #  mtzNameBase is the prefix to make unique column names
        
        # This function takes a list of mini-mtz files and runs mtzjoin to
        # merge them into one file (the two flags are used to 
        # adjust the input to mtzjoin, the first adds the mtz file type to
        # the output column names, the second uses fixed input 
        # column names). False, False is equiv. to the old makeHlin ftn
        # and True, True is the same as makeHklin0. 
        # It's advised to use the defaults for new interfaces.
        # The input miniMtzsIn is list of either CDataFile object names 
        # in inputData or list of [CDataFile object names in
        # inputData, target content type]
            
        error = CErrorReport()
        infiles = []
        allColout = ','
        outfile = os.path.join(self.workDirectory, hklin + '.mtz')
        count = 1
        self.fileSignatures = []   # column types for each input file
        for miniMtz in miniMtzList:
            mtzName = '{}-{}'.format(mtzNameBase, count)
            count += 1
            obj = miniMtz
            if obj is None:
                self.appendErrorReport(31, mtzName)
                error.append(self.__class__, 31, mtzName)
            elif not obj.isSet():
                self.appendErrorReport(30, mtzName)
                error.append(self.__class__, 30, mtzName)
            else:
                signature = ''
                for col in obj.fileContent.listOfColumns:
                    signature += str(col.columnType)
                self.fileSignatures.append(signature)
                self.appendMakeHklinInput(obj) # comments for log file
                self._buildInputVector(obj, mtzName, targetContent, error,
                                       extendOutputColnames, useInputColnames,
                                       infiles, allColout)
        #print('/\ outfile, infiles', outfile, infiles)
        #print("/\/\ ", self.fileSignatures)
        # Combine files
        status, ret = self.joinMtz(outfile, infiles)
        if status != CPluginScript.SUCCEEDED and ret not in ignoreErrorCodes:
            error.append(CPluginScript, ret, hklin)
            self.appendErrorReport(ret, hklin, cls=CPluginScript)

        allcoloutlist = []
        for infile in infiles:
            allcoloutlist.append(infile[1])
        return outfile, allcoloutlist, error

    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def makeCommandAndScript(self):

      print("makeCommandAndScript")
      # Input file
      self.appendCommandLine(['HKLIN', self.hklin])
      # Output file
      self.hklout = os.path.join(self.workDirectory,"hklout.mtz")
      self.appendCommandLine(['HKLOUT', self.hklout])

      # LABIN: first dataset is "native"
      s = ''
      #  self.coloutlist has labels for each dataset
      for i, collabel in enumerate(self.coloutlist):
          labels = collabel.split(',')
          if i == 0:
              s += ' FP='+labels[0]
              s += ' SIGFP='+labels[1]
          else:
              sn = str(i)
              s += ' FPH'+sn+'='+labels[0]
              s += ' SIGFPH'+sn+'='+labels[1]
      
      self.appendCommandScript('LABIN '+s)

      if self.container.controlParameters.RESOLUTION_MAX.isSet():
          maxres = str(self.container.controlParameters.RESOLUTION_MAX)
          print("** Maxres ", maxres)
          self.appendCommandScript('RESOLUTION ' + maxres)
      
    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def processOutputFiles(self):
        print("processOutputFiles")

        logfile = self.makeFileName('LOG')
        logscrape = scaleitLogScraper()
        logscrape.process(logfile)
        logxml = logscrape.makeXML()  # returns XML block SCALEITLOG

        # Complete XML file, including information about input files
        nativeDname, derivativeDnames = logscrape.fileInfo()
        xmlroot = etree.Element('SCALEIT')

        inputfiletypes = etree.Element('InputFileTypes')
        label = 'Native'
                  
        for signature in self.fileSignatures:
            filetype = self.interpretFileSignature(signature)
            addElement(inputfiletypes, label, filetype)
            label = 'Derivative'

        xmlroot.append(inputfiletypes)
        xmlroot.append(logxml)

        self.writeXML(xmlroot)
        return CPluginScript.SUCCEEDED

    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def interpretFileSignature(self, signature):
        print('interpretFileSignature', signature)
        s = ''
        if signature == 'FQ':
            s = 'Fmean'
        elif signature == 'JQ':
            s = 'Imean'
        elif signature == 'GLGL': 
            s = 'F+/F-'
        elif signature == 'KMKM':
            s = 'I+/I-'
        else:
            s = 'Unrecognised column types: '+signature
        return s

    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def makeFailXML(self, message):
        print("makeFailXML", message)
        xmlfail = etree.Element('SCALEIT')
        addElement(xmlfail, 'Fail', message)
        self.writeXML(xmlfail)
    # - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
    def writeXML(self, xml):
        print("*writeXML")
        et = etree.ElementTree(xml)
        # and write out the XML
        et.write(self.makeFileName('PROGRAMXML'), pretty_print=True)
