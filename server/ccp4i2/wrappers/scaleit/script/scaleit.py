import os

from lxml import etree

from ccp4i2.core import CCP4XtalData
from ccp4i2.core.CCP4ErrorHandling import CErrorReport
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.wrappers.scaleit.script.scaleit_logscraper import addElement, scaleitLogScraper
from ccp4i2.wrappers.scaleit.script.scaleit_utils import DatalistCheck


class scaleit(CPluginScript):
    TASKNAME = 'scaleit'
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
                       hklin='hklin'):
        # Takes files from a CMiniMtzDataFileList (CList of CObsDataFile)
        # and merges them into a single HKLIN file using makeHklinGemmi.
        # Returns (outfile, allcoloutlist, error) where allcoloutlist is
        # a list of comma-separated column name strings, one per input file.

        error = CErrorReport()
        file_objects = []
        self.fileSignatures = []

        for count, miniMtz in enumerate(miniMtzList, start=1):
            mtzName = f'{mtzNameBase}-{count}'
            obj = miniMtz

            if obj is None:
                self.appendErrorReport(31, mtzName)
                error.append(self.__class__, 31, mtzName)
                continue
            elif not obj.isSet():
                self.appendErrorReport(30, mtzName)
                error.append(self.__class__, 30, mtzName)
                continue

            # Build file signature (for processOutputFiles)
            signature = ''
            for col in obj.fileContent.listOfColumns:
                signature += str(col.columnType)
            self.fileSignatures.append(signature)

            # Pass file object directly via file_obj key
            spec = {
                'file_obj': obj,
                'name': mtzName,
                'display_name': mtzName,
            }
            if targetContent is not None:
                spec['target_contentFlag'] = targetContent
            file_objects.append(spec)

        # Merge files using makeHklinGemmi
        outfile = os.path.join(str(self.workDirectory), hklin + '.mtz')
        try:
            output_path = self.makeHklinGemmi(
                file_objects=file_objects,
                output_name=hklin,
                merge_strategy='rename'
            )
            outfile = str(output_path)
        except Exception as e:
            error.append(self.__class__, 200, str(e))

        # Build allcoloutlist: list of comma-separated column names per file
        allcoloutlist = []
        for count, miniMtz in enumerate(miniMtzList, start=1):
            obj = miniMtz
            if obj is None or not obj.isSet():
                continue
            mtzName = f'{mtzNameBase}-{count}'
            # Determine which content flag applies
            if targetContent is not None:
                cf = targetContent
            else:
                cf = int(obj.contentFlag)
            columns = obj.CONTENT_SIGNATURE_LIST[cf - 1]
            prefixed = [f'{mtzName}_{col}' for col in columns]
            allcoloutlist.append(','.join(prefixed))

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
