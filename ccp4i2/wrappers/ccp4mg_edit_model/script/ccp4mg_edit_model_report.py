import os
import xml.etree.ElementTree as ET

from ....core import CCP4Utils
from ....report.CCP4ReportParser import Report


class ccp4mg_edit_model_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'ccp4mg_edit_model'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        self.addResults()
        self.addText(text='The CCP4MG session is finished.')
        pdbsWrittenPath = './/number_output_files'
        pdbsWrittenStringS = xmlnode.findall(pdbsWrittenPath)
        if len(pdbsWrittenStringS) > 0:
          pdbsWrittenString = xmlnode.findall(pdbsWrittenPath)[0].text
          self.append('<br/>')
          self.addText(text=pdbsWrittenString+' PDB files were written to the CCP4i2 database during the session.')
        else:
          self.append('<br/>')
          self.addText(text='No PDB files were written to the CCP4i2 database during the session.')

        self.append('<div style="clear:both;"/>')
        stdoutFileNames = './/ccp4mg_edit_model/stdout.txt'
        stdoutStrings = xmlnode.findall(stdoutFileNames)
        if len(stdoutStrings)>0:
            self.addText(text="MrBUMP started with the following PDB files/chains from the template model search")
            self.append('<div style="clear:both;"/>')
            table = self.addTable(title='MrBUMP template search summary',label='MrBUMP template search summary')
            jobDirectory = jobInfo['fileroot']
            logName = os.path.basename(stdoutStrings[0].text)
            logText = ""
            model_chain   = []
            score         = []
            loc_seqid     = []
            overall_seqid = []
            source        = []
            with open(os.path.join(jobDirectory,"mrbump_"+logName)) as f:
                logTextLines = f.readlines()

                haveResults1 = False
                haveResults2 = False

                for l in logTextLines:
                    if "Template Model Search Results" in l and "####" in l:
                        haveResults1 = True
                    if "CHAIN_ID" in l and "SCORE" in l and "localSEQID" in l and "overallSEQID" in l:
                        haveResults2 = True
                    if haveResults1 and haveResults2:
                        if l.strip() == "":
                            break
                        if not l.strip().startswith("CHAIN_ID"):
                            split = l.strip().split()
                            model_chain.append(split[0])
                            score.append(split[1])
                            loc_seqid.append(split[2])
                            overall_seqid.append(split[3])
                            source.append(" ".join(split[4:]))
                            logText += l
            table.addData(title="PDB/Chain ID",data=model_chain)
            table.addData(title="Score",data=score)
            table.addData(title="Local Seq. Identity",data=loc_seqid)
            table.addData(title="Overall Seq. Identity",data=overall_seqid)
            table.addData(title="Source",data=source)
    
        self.append('<div style="clear:both;"/>')
            
#FIXME - XML PICTURE
        """    
        if len(self.jobInfo['filenames']["XYZOUT"]) > 0:
            self.drawPictures()
            self.append('<div style="clear:both;"/>')
        """    
        self.addMrBUMPLogs()
        self.append('<div style="clear:both;"/>')
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
            warningsFolder.addText(text='No warnings from ccp4mg_edit_model')
        return

    def addMrBUMPLogs(self, xmlnode=None, jobInfo=None, parent=None, objectNameMap={}, initiallyOpen=True):
        if parent is None: parent = self
        if xmlnode is None: xmlnode = self.xmlnode
        if jobInfo is None: jobInfo = self.jobInfo
        ccp4i2_root = CCP4Utils.getCCP4I2Dir()
        mrBumpDirPath = './/ccp4mg_edit_model/MGMRBUMPDIR'
        mrBumpDirStringS = xmlnode.findall(mrBumpDirPath)
        jobDirectory = jobInfo['fileroot']
        if len(mrBumpDirStringS)>0:
            mrBumpDir = xmlnode.findall(mrBumpDirPath)[0].text
            logFold = parent.addFold(label='MrBUMP Log Files', initiallyOpen=initiallyOpen,brief='Log files')
            logFold.addText(text="MrBUMP log files copied from"+mrBumpDir)

            logFileNames = ".//ccp4mg_edit_model/*[substring(name(), string-length(name()) - 3) = '.log']"
            txtFileNames = ".//substring(name(), string-length(name()) - 3) = '.txt']"
            logFilesStrings = xmlnode.findall(logFileNames)
            logFilesStrings.extend(xmlnode.findall(txtFileNames))
            for logNameFull in logFilesStrings:
                logName = os.path.basename(logNameFull.text)
                if os.path.exists(os.path.join(jobDirectory,"mrbump_"+logName)):
                    logFoldSub = logFold.addFold(label=logName, initiallyOpen=False,brief='Log files')
                    logText = ""
                    with open(os.path.join(jobDirectory,"mrbump_"+logName)) as f:
                        logText = f.read()
                    logFoldSub.addPre(text=logText)
                    logFoldSub.append('<div style="clear:both;"/>')
            logFold.append('<div style="clear:both;"/>')

    def drawPictures(self, xmlnode=None, jobInfo=None, parent=None, objectNameMap={}, initiallyOpen=True):
        if parent is None: parent = self
        if xmlnode is None: xmlnode = self.xmlnode
        if jobInfo is None: jobInfo = self.jobInfo
        ccp4i2_root = CCP4Utils.getCCP4I2Dir()
        baseScenePath = os.path.join(ccp4i2_root,'wrappers','ccp4mg_edit_model','script','ccp4mg_edit_model.scene.xml')

        pictureFold = parent.addFold(label='Pictures', initiallyOpen=initiallyOpen,brief='Picture')
        pictureGallery = pictureFold.addObjectGallery(style='float:left;',height='450px', tableWidth='260px', contentWidth='450px')
        jobDirectory = jobInfo['fileroot']

        if self.jobInfo and "filenames" in self.jobInfo and "XYZOUT" in self.jobInfo["filenames"]:
          i = 0
          for fname in self.jobInfo['filenames']["XYZOUT"]:

             annot = ""
             if self.jobInfo and "outputfiles" in self.jobInfo:
                 for f in self.jobInfo['outputfiles']:
                     if os.path.abspath(os.path.join(jobDirectory,f['filename'])) == os.path.abspath(fname):
                         annot = f['annotation']
                         break

             baseSceneXML = ET.parse(baseScenePath).getroot()
             et = ET.ElementTree(baseSceneXML)
             filename_element = et.findall(".//data/MolData/filename")[0]
             MolData_element = et.findall(".//data/MolData")[0]
             if len(annot)>0:
                 name_element = ET.Element("name")
                 name_element.text = annot
                 MolData_element.append(name_element)
             del filename_element.attrib["database"]
             filename_element.text = fname
             sceneFilePath = os.path.join(jobDirectory,'ccp4mg_edit_model_scene'+str(i)+'.scene.xml')
             et.write(sceneFilePath)
             pic = pictureGallery.addPicture(sceneFile=sceneFilePath,label='Picture of structure '+str(i+1))
             i = i + 1
