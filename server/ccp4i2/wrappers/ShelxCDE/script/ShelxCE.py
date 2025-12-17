import os

from lxml import etree

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.wrappers.ShelxCDE.script import ShelxCDEBase


class ShelxCE(ShelxCDEBase.ShelxCDEBase):
    TASKMODULE = 'test'                               # Where this plugin will appear on the gui
    TASKTITLE = 'ShelxE'     # A short title for gui menu
    DESCRIPTION = 'Phasing, density modification and autobuilding'
    TASKNAME = 'ShelxCE'                                  # Task name - should be same as class name
    TASKCOMMAND = 'shelxe'                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 9999999.9
    WHATNEXT = ['coot_rebuild','parrot',['modelcraft','$CCP4I2/wrappers/modelcraft/script/experimental.params.xml']]
    PERFORMANCECLASS = 'CExpPhasPerformance'
    ERROR_CODES = {  300 : { 'description' : 'No overall FOM found' },}
    
    def processInputFiles(self):
        
        result = self.prepareShelxc()
        if result != CPluginScript.SUCCEEDED: return result
        
        result = self.runShelxc()
        if result != CPluginScript.SUCCEEDED: return result
        
        self.xmlroot = etree.Element('ShelxCE')
        
        logFile =  os.path.normpath(os.path.join(self.getWorkDirectory(), 'shelxc.log'))
        self.scrapeShelxcLog(self.xmlroot, logFile)
        self.flushXML()
        
        insFilePath  = os.path.normpath(os.path.join(self.getWorkDirectory(),'result_fa.ins'))
        result = self.pdbToRes(self.container.inputData.HAIN, insFilePath)
        if result != CPluginScript.SUCCEEDED: return result
        
        self.watchFile(self.makeFileName('LOG'), self.handleShelxeLogChanged)
        return CPluginScript.SUCCEEDED
    
    def makeCommandAndScript(self):
        self.appendCommandLine('result')
        self.appendCommandLine('result_fa')
        for parameterName in ['i','h','z','sX','aN','n','q','mN']:
            parameter = getattr(self.container.keywords,parameterName)
            if parameter.isSet():
                if len(parameterName) == 1:
                    if parameter: self.appendCommandLine('-'+parameterName)
                else:
                    self.appendCommandLine('-'+parameterName[0:1]+parameter.__str__())
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        print('#shelxcd processOutputFiles')
        processId = self.getProcessId()
        from ccp4i2.core import CCP4Modules
        exitStatus = CCP4Modules.PROCESSMANAGER().getJobData(processId,'exitStatus')
        exitCode = CCP4Modules.PROCESSMANAGER().getJobData(processId,'exitCode')
        if exitStatus != CPluginScript.SUCCEEDED:
            print('ShelxCE ended with non success status')
            self.appendErrorReport(210)
            return CPluginScript.FAILED
        if exitCode != 0:
            print('ShelxCE exited with non zero code')
            self.appendErrorReport(211)
            return CPluginScript.FAILED
        self.scrapeShelxeLog(self.xmlroot)
        self.flushXML()
        if len(self.xmlroot.xpath('//BETAEXPIRED')) != 0:
            print('ShelxCE exited with beta expired')
            self.appendErrorReport(212)
            return CPluginScript.FAILED
        if len(self.xmlroot.xpath('//OverallFOMs/FOM')) > 0:
            self.container.outputData.PERFORMANCE.FOM.set(float(self.xmlroot.xpath('//OverallFOMs/FOM')[-1].text))
        if len(self.xmlroot.xpath('//OverallFOMs/PseudoFreeCC')) > 0:
            self.container.outputData.PERFORMANCE.CC.set(float(self.xmlroot.xpath('//OverallFOMs/PseudoFreeCC')[-1].text))
        
        phsFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'result.phs'))
        if self.container.keywords.i.isSet() and self.container.keywords.i:
            phsFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'result_i.phs'))
       
        if os.path.isfile(phsFilePath):
            mapFilePath = self.container.outputData.MAPOUT.fullPath.__str__()
            phsOutFilePath = self.container.outputData.PHSOUT.fullPath.__str__()
            result = self.harvestPhsFile(phsFilePath, mapFilePath, phsOutFilePath)
            
            if result != CPluginScript.SUCCEEDED: return result
            
            if self.container.keywords.i.isSet() and self.container.keywords.i:
                self.container.outputData.PHSOUT.annotation.set("ShelxE phases - reversed hand")
                self.container.outputData.MAPOUT.annotation.set("ShelxE map - reversed hand")
            else:
                self.container.outputData.PHSOUT.annotation.set("ShelxE phases")
                self.container.outputData.MAPOUT.annotation.set("ShelxE map")

        phaFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'result.pha'))
        if self.container.keywords.i.isSet() and self.container.keywords.i:
            phaFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'result_i.pha'))
        
        if os.path.isfile(phaFilePath):
            mapFilePath = self.container.outputData.ANOMMAPOUT.fullPath.__str__()
            result = self.harvestPhsFile(phaFilePath, mapFilePath)
            
            if result != CPluginScript.SUCCEEDED: return result
            
            if self.container.keywords.i.isSet() and self.container.keywords.i:
                self.container.outputData.ANOMMAPOUT.annotation.set("ShelxE anom. map - reversed hand")
            else:
                self.container.outputData.ANOMMAPOUT.annotation.set("ShelxE anom. map")

        pdbFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'result.pdb'))
        if self.container.keywords.i.isSet() and self.container.keywords.i:
            pdbFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),'result_i.pdb'))
        
        if os.path.isfile(pdbFilePath):
            self.container.outputData.XYZOUT.setFullPath(pdbFilePath)
            if self.container.keywords.i.isSet() and self.container.keywords.i:
                self.container.outputData.XYZOUT.annotation.set('Autobuild results - reversed hand')
            else: self.container.outputData.XYZOUT.annotation.set('Autobuild results')
        
        
        if self.container.keywords.i.isSet() and self.container.keywords.i:
            self.invertCoordinates(self.container.inputData.HAIN.fullPath.__str__(),
                                   self.container.outputData.HAOUT.fullPath.__str__())
            self.container.outputData.HAOUT.annotation.set('HA coords - reversed hand')

        self.txtOutputFiles()
        
        return CPluginScript.SUCCEEDED

