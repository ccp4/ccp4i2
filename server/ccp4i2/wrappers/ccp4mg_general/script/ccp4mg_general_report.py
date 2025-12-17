import xml.etree.ElementTree as etree

from ccp4i2.core import CCP4Utils
from ccp4i2.report.CCP4ReportParser import *


class ccp4mg_general_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'ccp4mg_general'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        self.addResults()
        self.addText(text='The CCP4MG session is finished.')
        pdbsWrittenPath = './/ccp4mg_general/number_output_files'
        pdbsWrittenStringS = xmlnode.findall(pdbsWrittenPath)
        if len(pdbsWrittenStringS) > 0:
          pdbsWrittenString = xmlnode.findall(pdbsWrittenPath)[0].text
          self.append('<br/>')
          self.addText(text=pdbsWrittenString+' PDB files were written to the CCP4i2 database during the session.')
        else:
          self.append('<br/>')
          self.addText(text='No PDB files were written to the CCP4i2 database during the session.')
        self.drawPictures()
        #self.drawWarnings()

    def drawWarnings(self, parent=None):
        if parent is None: parent = self
        warnings = self.xmlnode.findall('.//Warnings/Warning')
        warningsFolder = parent.addFold(label='Warnings', initiallyOpen=True)
        if len(warnings)>0:
            for warning in warnings:
                warningsFolder.addText(text=warning.text,style='color:red;')
                warningsFolder.append('<br/>')
        else:
            warningsFolder.addText(text='No warnings from ccp4mg_general')
        return

    def drawPictures(self, xmlnode=None, jobInfo=None, parent=None, objectNameMap={}, initiallyOpen=True):
        if parent is None: parent = self
        if xmlnode is None: xmlnode = self.xmlnode
        if jobInfo is None: jobInfo = self.jobInfo
        ccp4i2_root = CCP4Utils.getCCP4I2Dir()
        baseScenePath = os.path.join(ccp4i2_root,'wrappers','ccp4mg_general','script','ccp4mg_general.scene.xml')

        pictureFold = parent.addFold(label='Pictures', initiallyOpen=initiallyOpen,brief='Picture')
        pictureGallery = pictureFold.addObjectGallery(style='float:left;',height='450px', tableWidth='260px', contentWidth='450px')
        jobDirectory = jobInfo['fileroot']

        if self.jobInfo and "filenames" in self.jobInfo and "XYZOUT" in self.jobInfo["filenames"]:
          i = 0
          for fname in self.jobInfo['filenames']["XYZOUT"]:
             baseSceneXML = CCP4Utils.openFileToEtree(baseScenePath)
             et = etree.ElementTree(baseSceneXML)
             filename_element = et.findall(".//data/MolData/filename")[0]
             del filename_element.attrib["database"]
             filename_element.text = fname
             sceneFilePath = os.path.join(jobDirectory,'ccp4mg_general_scene'+str(i)+'.scene.xml')
             et.write(sceneFilePath)
             pic = pictureGallery.addPicture(sceneFile=sceneFilePath,label='Picture of structure '+str(i+1))
             i = i + 1
