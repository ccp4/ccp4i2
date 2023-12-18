from __future__ import print_function
"""
    MakeLink_report.py: CCP4 GUI Project
    
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

from report.CCP4ReportParser import Report
import sys
#from lxml import etree
from xml.etree import ElementTree as ET
from core import CCP4Utils

class MakeLink_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'MakeLink'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        clearingDiv = self.addDiv(style="clear:both;")
        self.addDefaultReport(self)
        clearingDiv = self.addDiv(style="clear:both;")

    def addDefaultReport(self, parent=None):
        if parent is None: parent=self
        self.picture(self)
        for AcedrgLinkNode in self.xmlnode.findall("."):
            try:
                cycleNode = AcedrgLinkNode.findall("Cycle")[0]
            except:
                print("Missing cycle")
            try:
                logTextNode = AcedrgLinkNode.findall("LogText")[0]
            except:
                print("Missing logText")
            try:
                newFold = parent.addFold(label="Log text for iteration "+cycleNode.text, initiallyOpen=True)
            except:
                print("Unable to make fold")
            try:
                newFold.addPre(text = logTextNode.text)
            except:
                print("Unable to add Pre")

    def addLigandToGallery(self,pictureGallery,baseScenePath,scenePath,pdbPath,dictPath,tlc,annotation):
       with open(baseScenePath,'r') as baseScene:
           baseSceneText = baseScene.read()
           specializedText = baseSceneText.replace('SUBSTITUTEME',pdbPath)
           rootNode = ET.fromstring(specializedText)
           molDataNode = rootNode.findall('./scene/data/MolData')[0]
           customResNode = ET.fromstring('''<customResCIFFiles> <cifmonomer> <name>'''+tlc+'''</name> <filename>'''+dictPath+'''</filename> </cifmonomer> </customResCIFFiles>''')
           molDataNode.append(customResNode)
           with open(scenePath,'w') as specializedScene:
               CCP4Utils.writeXML(specializedScene,ET.tostring(rootNode))
           pic = pictureGallery.addPicture(label=annotation,sceneFile=scenePath)
       return

    def picture(self,parent=None) :
      ccp4i2_root = CCP4Utils.getCCP4I2Dir()
      import os
      
      baseScenePath = os.path.join(ccp4i2_root,'pipelines','MakeLink','script','MakeLink.scene.xml')

      scenePath = os.path.join(self.jobInfo['fileroot'],'my.scene1.xml')
      pdbPath = self.jobInfo['filenames'].get('UNL_PDB', None)
      dictPath = self.jobInfo['filenames'].get('UNL_CIF', None)
      tlc = "UNL"  # This must match the dictionary
      annotation = self.jobInfo['filenames'].get('ANNOTATION', None)
#FIXME - XML PICTURE
      #self.addLigandToGallery(parent,baseScenePath,scenePath,pdbPath,dictPath,tlc,annotation)

      #pictureGallery = parent.addObjectGallery(style='float:left;border:1px solid black;',height='400px', tableWidth='170px', contentWidth='450px')
      #self.addLigandToGallery(baseScenePath,pdbPath,dictPath,tlc,scenePath,annotation,pictureGallery)

      return
