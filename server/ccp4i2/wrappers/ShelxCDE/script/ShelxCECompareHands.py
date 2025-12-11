from __future__ import print_function
from ccp4i2.baselayer import QtCore
from ccp4i2.core.CCP4PluginScript import CPluginScript
from lxml import etree
import os
from ccp4i2.wrappers.ShelxCDE.script import ShelxCE
import functools
from ccp4i2.core import CCP4Utils

class ShelxCECompareHands(ShelxCE.ShelxCE):
    TASKNAME = 'ShelxCECompareHands'                     # Task name - should be same as class name
    TASKVERSION= 0.0                                     # Version of this plugin
    ASYNCHRONOUS = False
    TIMEOUT_PERIOD = 9999999.9
    PERFORMANCECLASS = 'CExpPhasPerformance'
    ERROR_CODES = {  301 : { 'description' : 'Failed in first hand run' }, 302 : { 'description' : 'No MAPOUT created' }, 303:{'description':'Message from mars'}}
    
    def process(self):
        
        invalidFiles = self.checkInputData()
        if len(invalidFiles)>0:
            self.reportStatus(CPluginScript.FAILED)
        self.checkOutputData()
        
        self.xmlroot = etree.Element('ShelxCECompareHands')
        
        self.firstHandPlugin = self.makePluginObject('ShelxCE')
        self.firstHandPlugin.container.inputData.copyData(self.container.inputData)
        self.firstHandPlugin.container.controlParameters.copyData(self.container.controlParameters)
        self.firstHandPlugin.container.keywords.copyData(self.container.keywords)
        self.firstHandPlugin.container.keywords.i.set('False')
        self.firstHandPlugin.doAsync=True
        self.firstHandPlugin.finished.connect(functools.partial(self.pluginFinished, self.firstHandPlugin))
        self.watchFile(self.firstHandPlugin.makeFileName('PROGRAMXML'), functools.partial(self.subjobXMLChanged,'FirstHand'))
        self.firstHandPlugin.process()
        self.runningJobs=[self.firstHandPlugin]

        self.secondHandPlugin = self.makePluginObject('ShelxCE')
        self.secondHandPlugin.container.inputData.copyData(self.container.inputData)
        self.secondHandPlugin.container.keywords.copyData(self.container.keywords)
        self.secondHandPlugin.container.keywords.i.set('True')
        self.secondHandPlugin.container.controlParameters.copyData(self.container.controlParameters)
        self.secondHandPlugin.doAsync=True
        self.secondHandPlugin.finished.connect(functools.partial(self.pluginFinished, self.secondHandPlugin))
        self.watchFile(self.secondHandPlugin.makeFileName('PROGRAMXML'), functools.partial(self.subjobXMLChanged,'SecondHand'))
        self.secondHandPlugin.process()
        self.runningJobs.append(self.secondHandPlugin)

    @QtCore.Slot('CPluginScript',dict)
    def pluginFinished(self, plugin, statusDict={}):
        if plugin == self.firstHandPlugin: tag = 'FirstHand'
        elif plugin == self.secondHandPlugin: tag = 'SecondHand'
        self.checkFinishStatus(statusDict, 301, plugin.container.outputData.MAPOUT, 302)
        self.subjobXMLChanged(tag, plugin.makeFileName('PROGRAMXML'))
        self.runningJobs.remove(plugin)
        if len(self.runningJobs) == 0:
            self.finishUp()
    
    def finishUp(self):

        #Pick a winner !
        firstHandCC = float(self.firstHandPlugin.container.outputData.PERFORMANCE.CC)
        secondHandCC = float(self.secondHandPlugin.container.outputData.PERFORMANCE.CC)
        choiceNode = etree.SubElement(self.xmlroot,'HandChosen')
        if secondHandCC > firstHandCC:
            pluginOutput = self.secondHandPlugin.container.outputData
            choiceNode.text='Inverted'
        else:
            pluginOutput = self.firstHandPlugin.container.outputData
            choiceNode.text='Original'
        self.flushXML()
        
        pipelineOutput = self.container.outputData
        pipelineOutput.PERFORMANCE.CC.set(pluginOutput.PERFORMANCE.CC)
        pipelineOutput.PERFORMANCE.FOM.set(pluginOutput.PERFORMANCE.FOM)
        self.harvestFile(pluginOutput.MAPOUT, pipelineOutput.MAPOUT)
        self.harvestFile(pluginOutput.ANOMMAPOUT, pipelineOutput.ANOMMAPOUT)
        if pluginOutput.XYZOUT.exists():
            self.harvestFile(pluginOutput.XYZOUT, pipelineOutput.XYZOUT)
        if pluginOutput.HAOUT.exists():
            self.harvestFile(pluginOutput.HAOUT, pipelineOutput.HAOUT)
        self.harvestFile(pluginOutput.PHSOUT, pipelineOutput.PHSOUT)
        
        self.reportStatus(CPluginScript.SUCCEEDED)

    def subjobXMLChanged(self, nodeName, changedFile):
            with open(changedFile,'r') as changedXML:
                newText = changedXML.read()
            if True:
                from lxml import etree
                changedRoot = etree.fromstring(newText)
                #Remove old copies of XML
                oldNodes = self.xmlroot.xpath(nodeName)
                for oldNode in oldNodes:
                    oldNode.getparent().remove(oldNode)
                newNode = etree.SubElement(self.xmlroot,nodeName)
                newNode.append(changedRoot)
                import datetime
                timeNow = datetime.datetime.now()
                if not hasattr(self,"lastFlushTime"): self.lastFlushTime = timeNow
                deltaTime = timeNow - self.lastFlushTime
                if deltaTime.seconds > 0 or deltaTime.days > 0:
                    self.lastFlushTime = timeNow
                    self.flushXML()
            else:
                print('Failed to copyXML across')

    def flushXML(self):
        import tempfile
        import sys
        try:
            from lxml import etree
            xmlPath = self.makeFileName('PROGRAMXML')
            tfile = tempfile.NamedTemporaryFile()
            xmlTmpPath = tfile.name
            tfile.close()
            with open(xmlTmpPath,'w') as tmpFile:
                CCP4Utils.writeXML(tmpFile,etree.tostring(self.xmlroot, pretty_print=True))
            import shutil
            shutil.move(xmlTmpPath, xmlPath)
            
        except BaseException as e:
            print("########################################")
            print(e); sys.stdout.flush()
            print("########################################")

    def harvestFile(self, pluginOutputItem, pipelineOutputItem):
        import shutil
        try:
            shutil.copyfile(str(pluginOutputItem.fullPath), str(pipelineOutputItem.fullPath))
            pipelineOutputItem.annotation = pluginOutputItem.annotation
            pipelineOutputItem.contentFlag = pluginOutputItem.contentFlag
            pipelineOutputItem.subType = pluginOutputItem.subType
        except:
            self.appendErrorReport(202,str(pluginOutputItem.fullPath)+' '+str(pipelineOutputItem.fullPath))
            self.reportStatus(CPluginScript.FAILED)

    def checkFinishStatus( self, statusDict,failedErrCode,outputFile = None,noFileErrCode= None):
        import os
        if len(statusDict)>0 and statusDict['finishStatus'] == CPluginScript.FAILED:
            self.appendErrorReport(failedErrCode)
            self.reportStatus(statusDict['finishStatus'])
        try:
            assert outputFile.exists(),'Entity provided is not CDataFile or does not exist'
        except:
            self.appendErrorReport(noFileErrCode,'Expected file: '+str(outputFile))
            self.reportStatus(CPluginScript.FAILED)

