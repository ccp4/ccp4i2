from __future__ import print_function

import os
import tempfile

from core.CCP4PluginScript import CPluginScript
from core import CCP4Utils
import pathlib
from smartie import smartie

class areaimol(CPluginScript):

    TASKTITLE='Solvent accessible surface - Areaimol'
    TASKNAME = 'areaimol'
    TASKMODULE= 'model_data_utility'
    TASKCOMMAND = 'areaimol'
    TASKVERSION= 0.0
    COMLINETEMPLATE = None
    COMTEMPLATE = None

    def makeCommandAndScript(self):

      inp = self.container.inputData
      out = self.container.outputData

      if inp.XYZIN.fullPath.isSet():
          xyzin_target_file = str( inp.XYZIN.fullPath )
          self.appendCommandLine( [ "XYZIN",xyzin_target_file ] )

      if str(self.container.controlParameters.DIFFMODE) == "COMPARE" and inp.XYZIN2.fullPath.isSet():
          xyzin2_target_file = str( inp.XYZIN2.fullPath )
          self.appendCommandLine( [ "XYZIN2",xyzin2_target_file ] )

      if out.XYZOUT.fullPath.isSet():
          xyzout_target_file = str( out.XYZOUT.fullPath )
          self.appendCommandLine( [ "XYZOUT",xyzout_target_file ] )

      s = ""
      if self.container.controlParameters.EXTRA_AREAIMOL_KEYWORDS.isSet():
           for kwLine in str(self.container.controlParameters.EXTRA_AREAIMOL_KEYWORDS).split('\n'):
              kw = kwLine.lstrip().rstrip()
              if len(kw)>0:
                 if str(kw)[0] != '#':
                    if kw == "END":
                        break
                    s += kw + "\n"

      print("Keyword script:",kw)

      self.appendCommandScript( s+"\nDIFFMODE " + str(self.container.controlParameters.DIFFMODE) + "\n" )
      if str(self.container.controlParameters.DIFFMODE) == "COMPARE" and inp.XYZIN2.fullPath.isSet():
          self.appendCommandScript( s+"\nMATCHUP NOCOORDS\n" )
      if str(self.container.controlParameters.DIFFMODE) == "IMOL" and str(self.container.controlParameters.SYMMETRY) != "":
          self.appendCommandScript( s+"\nMODE NOHOH\n")
          self.appendCommandScript( s+"\nSYMMETRY 1\n")
          self.appendCommandScript( s+"\nSYMMETRY " + str(self.container.controlParameters.SYMMETRY) + "\n" )
          self.appendCommandScript( s+"\nTRANS 2\n")

      if str(self.container.controlParameters.DIFFMODE) == "OFF" and inp.XYZIN2.fullPath.isSet():
          self.appendCommandScript( s+"\nOUTPUT "+str(self.container.controlParameters.OUTPUT_MODE)+"\n" )
      else:
          self.appendCommandScript( s+"\nOUTPUT "+str(self.container.controlParameters.OUTPUT_MODE_COMPARE)+"\n" )
      self.appendCommandScript( s+"\nEND\n" )

      return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        from lxml import etree
        import sys
        import base64
        with open(self.makeFileName("PROGRAMXML"),"w") as programXMLFile:
            xmlStructure = etree.Element("areaimol")
            logText = etree.SubElement(xmlStructure,"LogText")
            with open(self.makeFileName("LOG"),"rb") as logFile:
                logText.text = base64.b64encode(logFile.read())
            try:
                smartie_text = ""
                smartie_logfile = smartie.parselog(self.makeFileName("LOG"))
                smartieText = etree.SubElement(xmlStructure,"SummaryText")
                for i in range(2,smartie_logfile.nsummaries()-1):
                    summary = smartie_logfile.summary(i).retrieve().split("\n")
                    smartie_text += ("\n".join(summary[1:-2])) + "\n"
                smartieText.text = base64.b64encode(bytes(smartie_text,"utf-8"))
            except:
                print("Failed to extract summaries with smartie")
                sys.stderr.write("Failed to extract summaries with smartie")
                exc_type, exc_value,exc_tb = sys.exc_info()[:3]
                sys.stderr.write(str(exc_type)+'\n')
                sys.stderr.write(str(exc_value)+'\n')
                print(str(exc_type)+'\n')
                print(str(exc_value)+'\n')

            CCP4Utils.writeXML(programXMLFile,etree.tostring(xmlStructure))
        self.container.outputData.XYZOUT.annotation = "Coordinates with atom surface area in 'B-factor' column"

        return CPluginScript.SUCCEEDED

