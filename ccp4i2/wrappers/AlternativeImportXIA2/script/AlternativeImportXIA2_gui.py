"""
Martin Noble
"""

import glob
import os

from lxml import etree
from PySide2 import QtCore

from ....core.CCP4ErrorHandling import CException
from ....qtgui import CCP4TaskWidget


class CTaskAlternativeImportXIA2(CCP4TaskWidget.CTaskWidget):
    
    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'AlternativeImportXIA2'
    TASKVERSION = 0.0
    TASKMODULE='data_entry'
    TASKTITLE='Import Xia2 results'
    DESCRIPTION='Harvest merged and unmerged files from selected XIA2 data reduction protocols'
    WHATNEXT = ['aimless_pipe']
    AUTOPOPULATEINPUT = False
    
    def __init__(self,parent):
        CCP4TaskWidget.CTaskWidget.__init__(self,parent)
    
    def drawContents(self):
        
        self.setProgramHelpFile('AlternativeImportXIA2')
        
        folder = self.openFolder(folderFunction='inputData',title='Xia2 runs',followFrom=False)
        
        self.createLine( [ 'advice','Xia2 directory to import'] )
        
        runSummaries = self.container.controlParameters.runSummaries
        
        self.createLine( [ 'widget', '-browseDb', True, 'XIA2_DIRECTORY', 'tip', 'Browse to the directory "xia2"' ] )
        self.container.inputData.XIA2_DIRECTORY.dataChanged.connect(self.handleSelectXia2Directory)
    
        self.createLine( [ 'advice','When a valid top-level XIA2 directory is selected above, the list of successful '] )
        self.createLine( [ 'advice','data reduction protocols that XIA2 performed will be summarised below.'] )
        self.createLine( [ 'advice','Select and delete ("-") the protocols you do not wish to import'] )
        self.createLine( [ 'advice',''] )
        self.createLine( [ 'advice','Protocols to keep from this XIA2 directory'] )
        self.createLine(['widget','runSummaries'])
        self.getWidget('runSummaries').setListVisible(True)
        try:
            self.getWidget('runSummaries').setEditorVisible(False)
        except:
            print('No editor stack')

    @QtCore.Slot()
    def handleSelectXia2Directory(self): 
        pattern = os.path.join(self.container.inputData.XIA2_DIRECTORY.__str__(),"")+"*"
        self.container.controlParameters.directoryPath = self.container.inputData.XIA2_DIRECTORY.__str__()
        candidateJobs = [os.path.split(path)[1] for path in glob.glob(pattern)]
        
        runSummaries = self.container.controlParameters.runSummaries
        runSummaries.setQualifiers({'listMinLength' : 0})
        while len(runSummaries) > 0:
            try:
                runSummaries.remove(runSummaries[-1])
            except CException as e:
                #Here should catch the exception that limits list lenght to >0
                print(e)
        runSummaries.setQualifiers({'listMinLength' : 1})
    
        for candidateJob in candidateJobs:
            xmlFilePath = os.path.join(self.container.inputData.XIA2_DIRECTORY.__str__(),candidateJob,'ispyb.xml')
            if os.path.isfile(xmlFilePath):
                with open(xmlFilePath) as xmlFile:
                    text = xmlFile.read()
                    xmlOfFile = etree.fromstring(text)
                    try:
                        spaceGroupText = (' '+xmlOfFile.xpath('//spaceGroup')[0].text)
                        resolutionText = (' Res: '+xmlOfFile.xpath('//resolutionLimitHigh')[0].text)
                        rMeasText = (' Rmeas: '+xmlOfFile.xpath('//rMeasAllIPlusIMinus')[0].text)
                        runSummaries.append(runSummaries.makeItem())
                        runSummaries[-1] = candidateJob + ':'
                        runSummaries[-1] += spaceGroupText
                        runSummaries[-1] += resolutionText
                        runSummaries[-1] += rMeasText
                    except:
                        print('Unable to find fields')
        self.validate()
        self.getWidget('runSummaries').updateViewFromModel()
