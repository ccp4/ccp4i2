from __future__ import print_function

"""
    mrbump_model_prep.py: CCP4 GUI Project
     Copyright (C) 2020 STFC

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the
     license to address the requirements of UK law.

     You should have received a copy of the modified GNU Lesser General
     Public License along with this library.  If not, copies may be
     downloaded from http://www.ccp4.ac.uk/ccp4license.php

     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
"""

import sys, os, shutil, copy
import json
from collections import OrderedDict

from core.CCP4PluginScript import CPluginScript
from core import CCP4Utils
from core import CCP4ErrorHandling
#from lxml import etree
from xml.etree import ElementTree as ET
class mrbump_model_prep(CPluginScript):

    TASKNAME = 'mrbump_model_prep'
    TASKCOMMAND = 'mrbump'
    TASKTITLE = 'MrBUMP model preparation'
    MAINTAINER = 'stuart.mcnicholas@york.ac.uk'
    WHATNEXT = []
    ERROR_CODES = { 301 : { 'description' : 'Failed somehow in MrBUMP model prep' }
                    }

    def writeProgramXML(self,logDir):

        MRBUMPFIELDS = ['chainSource', 'coverage', 'eLLG', 'evalue', 'experiment', 'mgName', 'modelName', 'modelPDBfile', 'rank', 'resolution', 'score', 'seqID', 'source', 'sourceChainID', 'tarEnd', 'tarGroupEnd', 'tarGroupStart', 'tarStart', 'type']

        from mrbump.output import modelout
        mjson=modelout.Json()
        modelsJsonFile=os.path.join(logDir, "models.json")
        if os.path.isfile(modelsJsonFile):
             model_dict=mjson.readJson(modelsJsonFile)
             xmlroot = ET.Element('mrbump_model_prep')

             for k,v in model_dict.items():
                 model = ET.SubElement(xmlroot,"model")
                 for field in MRBUMPFIELDS:
                     ele = ET.SubElement(model,field)
                     ele.text = str(getattr(v,field))

             with open(modelsJsonFile) as f:
                 j = json.load(f, object_pairs_hook=OrderedDict)

             for k,v in j.items():
                 theKey = k
                 break
                 
             bestFile = model_dict[theKey].modelPDBfile
             bestModel = ET.SubElement(xmlroot,"bestModel")
             bestModel.text = str(bestFile)

             with open(str(self.makeFileName('PROGRAMXML')), 'w') as ostream:
                 ET.indent(xmlroot)
                 CCP4Utils.writeXML(ostream,ET.tostring(xmlroot))

    def findOutputFileFromLog(self,logDir):

        from mrbump.output import modelout
        mjson=modelout.Json()
        modelsJsonFile=os.path.join(logDir, "models.json")
        if os.path.isfile(modelsJsonFile):
             model_dict=mjson.readJson(modelsJsonFile)
             with open(modelsJsonFile) as f:
                 j = json.load(f, object_pairs_hook=OrderedDict)

             for k,v in j.items():
                 theKey = k
                 break
                 
             newFile = model_dict[theKey].modelPDBfile
             if sys.version_info < (3,0):
                 return str(newFile)
             else:
                 return newFile
        else:
                 sys.write("Error: Can't find MrBUMP models json file:\n %s\n" % modelsJsonFile)

    def processInputFiles(self):
        from core import CCP4XtalData
        error = None
        self.hklin = None
        dataObjects = []
        if not self.container.inputData.F_SIGF.isSet():
            return CPluginScript.SUCCEEDED
        #Append Observation with representation dependent on whether we are detwining on Is or not
        dataObjects += [['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]]
        #Include FreeRflag if called for
        if self.container.inputData.FREERFLAG.isSet():
            dataObjects += ['FREERFLAG']
        self.hklin,error = self.makeHklin(dataObjects)
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        else:
            return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self):

      inp = self.container.inputData
      out = self.container.outputData

      from core import CCP4Utils
      import os

      keyin = "GESMAX 1\n" 
      keyin += "PICKLE False\n" 
      keyin += "MDLC True\n" 
      keyin += "MDLM False\n" 
      keyin += "MDLS False\n" 

      if inp.SEARCH_AFDB and inp.AFDBLEVEL != 'None':
          keyin += "AFLEVEL "+str(inp.AFDBLEVEL)+"\n" 
      else:
          keyin += "RLEVEL all\n" 

      if inp.SEARCH_PDB and inp.REDUNDANCYLEVEL:
          if str(inp.REDUNDANCYLEVEL) == '110':
             keyin += "RLEVEL ALL\n" 
          elif str(inp.REDUNDANCYLEVEL) == '100':
             keyin += "RLEVEL 100\n" 
          elif str(inp.REDUNDANCYLEVEL) == '95':
             keyin += "RLEVEL 95\n" 
          elif str(inp.REDUNDANCYLEVEL) == '90':
             keyin += "RLEVEL 90\n" 
          elif str(inp.REDUNDANCYLEVEL) == '70':
             keyin += "RLEVEL 70\n" 
          elif str(inp.REDUNDANCYLEVEL) == '50':
             keyin += "RLEVEL 50\n" 
      elif inp.SEARCH_PDB:
         keyin += "RLEVEL 95\n" 

      if inp.MRMAX:
          keyin += "MRNUM %s\n" % str(inp.MRMAX)
      else:
          keyin += "MRNUM 10\n" 

      keyin += "EXIT True\n" 

      keyfile=os.path.join(self.getWorkDirectory(), "keywords.txt")
      kf=open(keyfile, "w")
      kf.write(keyin)
      kf.close()

      if self.container.inputData.F_SIGF.isSet():
          self.appendCommandLine( [ 'HKLIN', self.hklin ] )
      
      seqFile = os.path.join(self.workDirectory,'SEQIN.fasta')
      inp.ASUIN.writeFasta(fileName=seqFile)
      self.appendCommandLine( [ 'SEQIN', seqFile ] )
      self.appendCommandLine( [ 'KEYIN', str( keyfile ) ] )

      return CPluginScript.SUCCEEDED

    def processOutputFiles(self):

        self.writeProgramXML(os.path.join(self.getWorkDirectory(),"search_mrbump_1","logs"))
        xyzout = self.findOutputFileFromLog(os.path.join(self.getWorkDirectory(),"search_mrbump_1","logs"))
        if os.path.exists(xyzout):
            self.container.outputData.XYZOUT.set(xyzout)
            self.container.outputData.XYZOUT.annotation.set('Model from MrBump model search')

        return CPluginScript.SUCCEEDED
