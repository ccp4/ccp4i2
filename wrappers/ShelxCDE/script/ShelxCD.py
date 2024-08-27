from core.CCP4PluginScript import CPluginScript
from xml.etree import ElementTree as ET
import os
from wrappers.ShelxCDE.script import ShelxCDEBase

class ShelxCD(ShelxCDEBase.ShelxCDEBase):
    TASKMODULE = 'test'                               # Where this plugin will appear on the gui
    TASKNAME = 'ShelxCD'                                  # Task name - should be same as class name
    TASKCOMMAND = 'shelxd'                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 9999999.9
    #WHATNEXT = ['ShelxCE','phaser_EP_AUTO']
    PERFORMANCECLASS = 'CExpPhasPerformance'

    def processInputFiles(self):
    
        result = self.prepareShelxc()
        if result != CPluginScript.SUCCEEDED: return result
    
        result = self.runShelxc()
        if result != CPluginScript.SUCCEEDED: return result

        self.xmlroot = ET.Element('ShelxCD')
        
        logFile =  os.path.normpath(os.path.join(self.getWorkDirectory(), 'shelxc.log'))
        self.scrapeShelxcLog(self.xmlroot, logFile)
        self.flushXML()
        
        self.watchFile(self.makeFileName('LOG'), self.handleShelxdLogChanged)
        return CPluginScript.SUCCEEDED
    
    def makeCommandAndScript(self):
        self.appendCommandLine('result_fa')
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        result = self.scrapeShelxdLog(self.xmlroot)
        if result != CPluginScript.SUCCEEDED: return result
        self.flushXML()
        
        resultFaPath = os.path.join(self.getWorkDirectory(),'result_fa.hkl')
        if os.path.exists(resultFaPath):
            self.container.outputData.FA.setFullPath(resultFaPath)
            self.container.outputData.FA.annotation = 'Derived SAD FAs'
        
        resultPdbPath = os.path.join(self.getWorkDirectory(),'result_fa.pdb')
        if os.path.exists(resultPdbPath):
            reElementedPdbPath = os.path.join(self.getWorkDirectory(),'reElemented.pdb')
            outputElement = self.container.controlParameters.SFAC.__str__().strip()
            with open(resultPdbPath,'r') as resultFile:
                with open(reElementedPdbPath,'w') as reElementedFile:
                    lines = resultFile.readlines()
                    for line in lines:
                        outputLine = line
                        if line.startswith('HETATM'):
                            if len(outputElement) == 1: outputLine = line[0:13] + outputElement.upper() + line[14:]
                            else: outputLine = line[0:12] + outputElement.upper() + line[14:]
                        reElementedFile.write(outputLine)
                    self.container.outputData.XYZOUT.setFullPath(reElementedPdbPath)
                    self.container.outputData.XYZOUT.annotation = 'SHELXD HA sites'
                    self.container.outputData.XYZOUT.subType = 4
        self.txtOutputFiles()

        unsortedList = []
        tries = self.xmlroot.findall('.//Try')
        for aTry in tries:
            if len(aTry.findall('CFOM')) > 0 and len(aTry.findall('CCAll'))>0:
                unsortedList.append({'CFOM':float(aTry.findall('CFOM')[-1].text), 'CC':float(aTry.findall('CCAll')[-1].text)})
        if len(unsortedList) > 0:
            sortedList = sorted(unsortedList, key=lambda aTry: aTry['CFOM'])
            self.container.outputData.PERFORMANCE.CFOM.set(sortedList[-1]['CFOM'])
            self.container.outputData.PERFORMANCE.CC.set(sortedList[-1]['CC'])
        
        return CPluginScript.SUCCEEDED


