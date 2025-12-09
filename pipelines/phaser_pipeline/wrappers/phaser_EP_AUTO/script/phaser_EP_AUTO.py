from __future__ import print_function

from ccp4i2.core.CCP4PluginScript import CPluginScript
import sys, os
from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core import CCP4Modules
from pipelines.phaser_pipeline.wrappers.phaser_MR.script import phaser_MR
from lxml import etree
from ccp4i2.core import CCP4Utils

class EPAUTOCallbackObject(phaser_MR.CallbackObject):
    def __init__(self, xmlroot=None, xmlResponders = [],workDirectory=None):
        super(EPAUTOCallbackObject, self).__init__(xmlroot, xmlResponders)
        #self.summaryFile = open(os.path.join(workDirectory,'summary.txt'),'w')
        self._summary_buffer = ""
    
    def call_back(self, label, text):
        if label == 'current best solution':
            try:
                for oldNode in self.xmlroot.xpath('//PhaserCurrentBestSolution'):
                    oldNode.getparent().remove(oldNode)
                bestSolNode =etree.SubElement(self.xmlroot,'PhaserCurrentBestSolution')
                phaser_MR.xmlFromSol(text, bestSolNode)
                self.notifyResponders()
            except:
                print('\n\n Exception in analysing current best solution')
        elif label == 'summary':
            print('***Summary callback:[', text.replace('\n','\n***********'),'\n******]\n')
            if text.startswith("**********") and not self._summary_buffer.strip().endswith("***"):
                self.flushSummary()
            self._summary_buffer += text
        else:
            print('***',label,' callback:', text)


    def flushSummary(self):
        summaryNode = etree.SubElement(self.xmlroot,'Summary')
        summaryNode.text = self._summary_buffer
        self.notifyResponders()
        self._summary_buffer = ""

class phaser_EP_AUTO(phaser_MR.phaser_MR):

    TASKNAME = 'phaser_EP_AUTO'                                  # Task name - should be same as class name
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template
    ASYNCHRONOUS = False
    RUNEXTERNALPROCESS=False

    ERROR_CODES = { 201 : { 'description' : 'Failed to find file' }, 202 : { 'description' : 'Failed to interpret searches from Ensemble list' },}
    requiredDefaultList = ['PART_VARI', 'PART_DEVI']

    '''
    def __init__(self,parent=None,name=None,workDirectory=''):
        phaser_MR.phaser_MR. __init__(self,parent=parent,name=name,workDirectory=workDirectory)
    '''
    def __init__(self, *args, **kw):
        super(phaser_EP_AUTO, self).__init__(*args, **kw)
        
        #Create a callback Object that will respond to callbacks from Phaser, principally by putting information
        #intp the outputXML of this plugin

        self.xmlroot = etree.Element('PhaserEpResults')
        self.callbackObject = EPAUTOCallbackObject(self.xmlroot, [self.flushXML])
    
    def startProcess(self, command, **kw):
        
        import phaser
        outputObject = phaser.Output()
        outputObject.setPhenixCallback(self.callbackObject)

        xtalid = 'Scatterers'
        #Note, in phaser distributed with ccp4-7, waveid is nolonger "space" toleratnt
        waveid = 'Wave1'

        inputObject = phaser.InputEP_DAT()
        inputObject.setHKLI(str(self.hklin))
        inputObject.addCRYS_ANOM_LABI(xtalid,waveid,"Fplus","SIGFplus","Fminus","SIGFminus")
        
        inp = self.container.inputData
        if inp.RESOLUTION_LOW.isSet():
            if inp.RESOLUTION_HIGH.isSet():
                inputObject.setRESO(float(inp.RESOLUTION_LOW), float(inp.RESOLUTION_HIGH))
            else:
                inputObject.setHIRES(float(inp.RESOLUTION_LOW))
        elif inp.RESOLUTION_HIGH.isSet():
            inputObject.setHIRES(float(inp.RESOLUTION_HIGH))

        inputObject.setMUTE(False)
        resultObject = phaser.runEP_DAT(inputObject, outputObject)
        self.inputHall = resultObject.getSpaceGroupHall()
        
        if not resultObject.Success():
            self.appendErrorReport(105, resultObject.ErrorName() + '-' + resultObject.ErrorMessage())
            return CPluginScript.FAILED
        
        with open (self.makeFileName('LOG'),'w') as logfile:
            logfile.write(resultObject.logfile())

        hkl = resultObject.getMiller()
        Fpos = resultObject.getFpos(xtalid,waveid)
        Spos = resultObject.getSIGFpos(xtalid,waveid)
        Ppos = resultObject.getPpos(xtalid,waveid)
        Fneg = resultObject.getFneg(xtalid,waveid)
        Sneg = resultObject.getSIGFneg(xtalid,waveid)
        Pneg = resultObject.getPneg(xtalid,waveid)

        i = phaser.InputEP_AUTO()
        i.setSPAC_HALL(resultObject.getSpaceGroupHall())
        cell = resultObject.getUnitCell()
        i.setCELL(cell[0],cell[1],cell[2],cell[3],cell[4],cell[5])
        i.setCRYS_MILLER(hkl)
        i.addCRYS_ANOM_DATA(xtalid,waveid,Fpos,Spos,Ppos,Fneg,Sneg,Pneg)

        if self.setKeywords(i) == CPluginScript.FAILED:
            self.appendErrorReport(102, 'Failed to set phaser keywords in phaser_EP_AUTO')
            return CPluginScript.FAILED

        if self.container.inputData.PARTIALMODELORMAP.__str__() == 'MODEL':
            if self.container.inputData.XYZIN_PARTIAL.isSet() and os.path.isfile(self.container.inputData.XYZIN_PARTIAL.fullPath.__str__()):
                i.setPART_PDB(self.container.inputData.XYZIN_PARTIAL.fullPath.__str__())
        elif self.container.inputData.PARTIALMODELORMAP.__str__() == 'MAP':
            if self.container.inputData.MAPCOEFF_PARTIAL.isSet() and os.path.isfile(self.container.inputData.MAPCOEFF_PARTIAL.fullPath.__str__()):
                i.setPART_HKLI(self.container.inputData.MAPCOEFF_PARTIAL.fullPath.__str__())
                i.setPART_LABI_FC('F')
                i.setPART_LABI_PHIC('PHI')
        elif self.container.inputData.PARTIALMODELORMAP.__str__() == 'NONE' and self.container.inputData.XYZIN_HA.isSet() and os.path.isfile(self.container.inputData.XYZIN_HA.fullPath.__str__()):
            i.setATOM_PDB(xtalid, self.container.inputData.XYZIN_HA.fullPath.__str__())
        
        i.setLLGC_COMP(True)
        i.setLLGC_NCYC(int(self.container.inputData.LLGC_CYCLES))
        
        # Search for "Generic" anomalous scatterer
        if self.container.inputData.PURE_ANOMALOUS.isSet() and self.container.inputData.PURE_ANOMALOUS:
            i.addLLGC_ANOM(True)
        else:
            i.addLLGC_ANOM(False)
        
        for element in self.container.inputData.ELEMENTS:
            i.addLLGC_ELEM(str(element))
        
        i.setWAVE(float(self.container.inputData.WAVELENGTH))

        if self.parseContent(i) == CPluginScript.FAILED:
            self.appendErrorReport(104, 'Failed to parse ASU content in phaser_EP_AUTO')
            return CPluginScript.FAILED

        i.setMUTE(False)
        i.setVERB(True)
        i.setOUTP_LEVE('SUMMARY')
        self.resultObject = phaser.runEP_AUTO(i, outputObject)
        if not self.resultObject.Success():
            self.appendErrorReport(105, self.resultObject.ErrorName() + '-' + self.resultObject.ErrorMessage())
            return CPluginScript.FAILED
        
        with open (self.makeFileName('LOG'),'a') as logfile:
            logfile.write(self.resultObject.logfile())
            jobInfo = CCP4Modules.PROJECTSMANAGER().db().getJobInfo(jobId=self.jobId)
            if "jobtitle" in jobInfo and jobInfo["jobtitle"]:
                logfile.write(str(jobInfo["jobtitle"])+"\n")
            while "parentjobid" in jobInfo and jobInfo["parentjobid"]:
                jobInfo = CCP4Modules.PROJECTSMANAGER().db().getJobInfo(jobId=jobInfo["parentjobid"])
                if "jobtitle" in jobInfo and jobInfo["jobtitle"]:
                    logfile.write(str(jobInfo["jobtitle"])+"\n")

        return CPluginScript.SUCCEEDED

    def processInputFiles(self):
        from ccp4i2.core import CCP4XtalData
        self.hklin,error = self.makeHklin([['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FPAIR]])
        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            for report in error._reports:
                if report['code'] == 32:
                    report['details'] = 'F+ and F- cannot be derived from data. Check file import.'
            self.appendErrorReport(201, 'Failed to prepare input MTZ file (makeHklin failed). ' + str(error))
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED
    
    def subElementWithNameAndText(self, parentNode, name, text):
        newNode = etree.SubElement(parentNode, name)
        newNode.text = text
        return newNode

    def flushXML(self, xml):
        from lxml import etree
        import os
        tmpFilename = self.makeFileName('PROGRAMXML')+'_tmp'
        with open(tmpFilename,'w') as tmpFile:
            CCP4Utils.writeXML(tmpFile,etree.tostring(xml, pretty_print=True))
        self.renameFile(tmpFilename, self.makeFileName('PROGRAMXML'))


    def processOutputFiles(self):
        import os,shutil
        resultObject = self.resultObject

        from ccp4i2.core import CCP4XtalData

        for hand in ['','.hand']:
            possibleCoords = os.path.join(self.getWorkDirectory(),'PHASER.1'+hand+'.pdb')
            if os.path.isfile(possibleCoords):
                self.container.outputData.XYZOUT.append(self.container.outputData.XYZOUT.makeItem())
                self.container.outputData.XYZOUT[-1].setFullPath(possibleCoords)
                self.container.outputData.XYZOUT[-1].annotation.set('Sites')
        if len(self.container.outputData.XYZOUT)>0:
            self.container.outputData.XYZOUT[0].annotation.set(str(self.container.outputData.XYZOUT[0].annotation) + ' - original hand')
            self.container.outputData.XYZOUT[0].subType = 4
        if len(self.container.outputData.XYZOUT)>1:
            self.container.outputData.XYZOUT[1].annotation.set(str(self.container.outputData.XYZOUT[1].annotation) + ' - reversed hand')
            self.container.outputData.XYZOUT[1].subType = 4

        for hand in ['','.hand']:
            possibleHKLOUT = os.path.join(self.getWorkDirectory(),'PHASER.1'+hand+'.mtz')
            if os.path.isfile(possibleHKLOUT):
                self.container.outputData.HKLOUT.append(self.container.outputData.HKLOUT.makeItem())
                self.container.outputData.HKLOUT[-1].setFullPath(possibleHKLOUT)

        if len(self.container.outputData.HKLOUT) > 0:
            self.splitHkloutList(miniMtzsOut=['ABCDOUT','MAPOUT'],programColumnNames=['HLA,HLB,HLC,HLD','FWT,PHWT'],outputBaseName=['ABCDOUT','MAPOUT'],outputContentFlags=[CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL,1],infileList=self.container.outputData.HKLOUT)
        for hlCoeffs in self.container.outputData.ABCDOUT:
            hlCoeffs.annotation.set('Phase estimates')
        if len(self.container.outputData.ABCDOUT)>0:
            self.container.outputData.ABCDOUT[0].annotation.set(str(self.container.outputData.ABCDOUT[0].annotation) + ' - original hand')
        if len(self.container.outputData.ABCDOUT)>1:
            self.container.outputData.ABCDOUT[1].annotation.set(str(self.container.outputData.ABCDOUT[1].annotation) + ' - reversed hand')
        
        for map in self.container.outputData.MAPOUT:
            map.annotation.set('Phased map')
            map.subType = 1
        if len(self.container.outputData.MAPOUT)>0:
            self.container.outputData.MAPOUT[0].annotation.set(str(self.container.outputData.MAPOUT[0]) + ' - original hand')
        if len(self.container.outputData.MAPOUT)>1:
            self.container.outputData.MAPOUT[1].annotation.set(str(self.container.outputData.MAPOUT[1]) + ' - reversed hand')

        while len(self.container.outputData.HKLOUT) > 0:
            self.container.outputData.HKLOUT.remove(self.container.outputData.HKLOUT[-1])

        for hand in ['','.hand']:
            possibleHKLOUT = os.path.join(self.getWorkDirectory(),'PHASER.1'+hand+'.llgmaps.mtz')
            if os.path.isfile(possibleHKLOUT):
                self.container.outputData.HKLOUT.append(self.container.outputData.HKLOUT.makeItem())
                self.container.outputData.HKLOUT[-1].setFullPath(possibleHKLOUT)
        if len(self.container.outputData.HKLOUT) > 0:
            self.splitHkloutList(miniMtzsOut=['LLGMAPOUT'],programColumnNames=['FLLG_AX,PHLLG_AX'],outputBaseName=['LLGMAPOUT'],infileList=self.container.outputData.HKLOUT)

        for llgMap in self.container.outputData.LLGMAPOUT:
            llgMap.annotation.set('Anomalous LLG map')
            llgMap.contentFlag = 1
            llgMap.subType = 2
        if len(self.container.outputData.LLGMAPOUT)>0:
            self.container.outputData.LLGMAPOUT[0].annotation.set(str(self.container.outputData.LLGMAPOUT[0]) + ' - original hand')
        if len(self.container.outputData.LLGMAPOUT)>1:
            self.container.outputData.LLGMAPOUT[1].annotation.set(str(self.container.outputData.LLGMAPOUT[1]) + ' - reversed hand')
        
        return CPluginScript.SUCCEEDED

