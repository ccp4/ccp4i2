import os
import pickle
import sys

from lxml import etree

from ccp4i2.core import CCP4ErrorHandling, CCP4Modules, CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR.script import phaser_MR


class MRAUTOCallbackObject(phaser_MR.CallbackObject):
    def __init__(self, xmlroot=None, xmlResponders = [],workDirectory=None):
        super(MRAUTOCallbackObject, self).__init__(xmlroot, xmlResponders)
        #self.summaryFile = open(os.path.join(workDirectory,'summary.txt'),'w')
        self._summary_buffer = ""
        self.minimumDelayInSeconds = 5
    
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
            #self.summaryFile.write ("{"+text+"}")
            #self.summaryFile.flush()
            if text.startswith("**********") and not self._summary_buffer.strip().endswith("***"):
                self.flushSummary()
            self._summary_buffer += text
        else:
            #self.summaryFile.write ("["+text+"]")
            #self.summaryFile.flush()
            pass


    def flushSummary(self):
        summaryNode = etree.SubElement(self.xmlroot,'Summary')
        summaryNode.text = self._summary_buffer
        self.notifyResponders()
        self._summary_buffer = ""

    # Here I override notifyResponders so as to notify reponders maximally once per 5 seconds.
    def notifyResponders(self):
        import datetime
        if not hasattr(self,'lastNotification'):
            self.lastNotification = datetime.datetime.now()
        datetimeNow = datetime.datetime.now()
        timeSinceNotification = datetimeNow - self.lastNotification
        #self.summaryFile.write(str( timeSinceNotification.microseconds)+"\n")
        if timeSinceNotification.days > 0 or timeSinceNotification.seconds > self.minimumDelayInSeconds:
            #self.summaryFile.write('*******'+str( timeSinceNotification.microseconds)+"\n")
            super(MRAUTOCallbackObject,self).notifyResponders()
            self.lastNotification = datetimeNow

class phaser_MR_AUTO(phaser_MR.phaser_MR):

    TASKNAME = 'phaser_MR_AUTO'
    WHATNEXT = ['prosmart_refmac','modelcraft','coot_rebuild']

    ERROR_CODES = { 201 : { 'description' : 'Failed to find file' }, 202 : { 'description' : 'Failed to interpret searches from Ensemble list' },}

    def __init__(self, *args, **kw):
        super(phaser_MR_AUTO, self).__init__(*args, **kw)
        #Create a callback Object that will respond to callbacks from Phaser, principally by putting information
        #intp the outputXML of this plugin
        self.xmlroot = etree.Element('PhaserMrResults')
        self.callbackObject = MRAUTOCallbackObject(self.xmlroot, [self.flushXML], self.workDirectory)

    def runMR_DAT(self, outputObject):
        import phaser
        inputObject = phaser.InputMR_DAT()
        self.inputObject = inputObject
        inputObject.setHKLI(str(self.hklin))
        #print '\n\n\n****Dir of data input object'
        #print dir(inputObject)
        if self.container.inputData.F_OR_I.isSet() and self.container.inputData.F_OR_I.__str__() == 'I':
            inputObject.setLABI_I_SIGI('I','SIGI')
        else:
            inputObject.setLABI_F_SIGF('F','SIGF')
        inp = self.container.inputData
        if inp.RESOLUTION_LOW.isSet():
            if inp.RESOLUTION_HIGH.isSet():
                inputObject.setRESO(float(inp.RESOLUTION_LOW), float(inp.RESOLUTION_HIGH))
            else:
                inputObject.setHIRES(float(inp.RESOLUTION_LOW))
        elif inp.RESOLUTION_HIGH.isSet():
            inputObject.setHIRES(float(inp.RESOLUTION_HIGH))
        inputObject.setMUTE(False)
        resultObject = phaser.runMR_DAT(inputObject, outputObject)
        # Note: Don't write logfile() here - it truncates the LOG while C++ stdout
        # is still redirected to it. The logfile() write is done after
        # finishCaptureCPlusPlusStdout() in startProcess().

        if not resultObject.Success():
            self.appendErrorReport(105, resultObject.ErrorName() + '-' + resultObject.ErrorMessage())
            return CPluginScript.FAILED
        return resultObject
    
    def startProcess(self):
        import phaser
        outputObject = phaser.Output()
        outputObject.setPhenixCallback(self.callbackObject)
        self.prepareCaptureCPlusPlusStdoutToLog()
        resultObject = self.runMR_DAT(outputObject)
        self.finishCaptureCPlusPlusStdout()
        # Write phaser's internal log buffer AFTER restoring stdout (not before!)
        # The first call uses 'w' to start fresh, subsequent calls use 'a' to append
        if resultObject != CPluginScript.FAILED:
            with open(self.makeFileName('LOG'), 'w') as logfile:
                logfile.write(resultObject.logfile())

        if resultObject == CPluginScript.FAILED: return CPluginScript.FAILED
        self.inputHall = resultObject.getSpaceGroupHall()
        self.inputSpaceGroup = resultObject.getSpaceGroupName()
        inputObject = phaser.InputMR_AUTO()
        inputObject.setKILL_FILE(os.path.join(self.getWorkDirectory(),'INTERRUPT'))
        inputObject.setSPAC_HALL(resultObject.getSpaceGroupHall())
        inputObject.setCELL6(resultObject.getUnitCell())
        #print '\n\n\nresutObject',dir(resultObject)
        if self.container.inputData.F_OR_I.isSet() and self.container.inputData.F_OR_I.__str__() == 'I':
            inputObject.setREFL_I_SIGI(resultObject.getMiller(),resultObject.getIobs(),resultObject.getSigIobs())
        else:
            inputObject.setREFL_F_SIGF(resultObject.getMiller(),resultObject.getFobs(),resultObject.getSigFobs())
        if self.setKeywords(inputObject) == CPluginScript.FAILED:
            return CPluginScript.FAILED
        if self.parseContent(inputObject) == CPluginScript.FAILED:
            return CPluginScript.FAILED
        if self.parseEnsembles(inputObject) == CPluginScript.FAILED:
            return CPluginScript.FAILED
        if self.addSearches(inputObject) == CPluginScript.FAILED:
            return CPluginScript.FAILED
        if self.parseSolutions(inputObject) == CPluginScript.FAILED:
            return CPluginScript.FAILED
        if self.container.inputData.KILLFILEPATH.isSet():
            inputObject.setKILL_FILE(self.container.inputData.KILLFILEPATH.__str__())
        else:
            inputObject.setKILL_FILE(os.path.join(self.getWorkDirectory(),'INTERRUPT'))
        #Alternative space groups
        #print '\n\n\n****Dir of autoMR input object'
        #print [word+'\n' for word in dir(inputObject) if 'sg' in word.lower()]
        inp = self.container.inputData
        if inp.SGALT_SELECT.isSet():
            inputObject.setSGAL_SELE(str(inp.SGALT_SELECT))
            #print 'Setting SGAL_SELE to ',str(inp.SGALT_SELECT)
            if inp.SGALT_SELECT.__str__() == 'LIST' and inp.SGALT_TEST.isSet():
                for sgAltTest in inp.SGALT_TEST:
                    inputObject.addSGAL_TEST(sgAltTest.__str__())

        #Now run the main calculation....do something to catch the stdout from the
        #underlying C++ calls
        inputObject.setMUTE(False)
        self.prepareCaptureCPlusPlusStdoutToLog()
        inputObject.setKEYW(True)
        try:
            self.resultObject = phaser.runMR_AUTO(inputObject, outputObject)
        except RuntimeError as e:
            self.finishCaptureCPlusPlusStdout()
            self.appendErrorReport(105, str(e))
            return CPluginScript.FAILED


        self.finishCaptureCPlusPlusStdout()
        # Write phaser's internal log buffer for the main MR_AUTO phase (append mode)
        with open(self.makeFileName('LOG'), 'a') as logfile:
            logfile.write(self.resultObject.logfile())

        if not self.resultObject.Success():
            self.appendErrorReport(105, self.resultObject.ErrorName() + '-' + self.resultObject.ErrorMessage())
            return CPluginScript.FAILED

        self.analyseResults(self.resultObject)
        return CPluginScript.SUCCEEDED

    def processInputFiles(self):
        from ccp4i2.core import CCP4XtalData

        # Changed Mtz merging to included phases. Due to issues with makeHkln() (column names), I used cad to manually merge files.
        cnMtz = ['F_SIGF']
        if self.container.inputData.F_OR_I.isSet() and self.container.inputData.F_OR_I.__str__() == 'I':
            self.hklin,error = self.makeHklin([['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_IMEAN]])
        else:
            self.hklin,error = self.makeHklin([['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]])

        if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            for report in error._reports:
                if report['code'] == 32:
                    report['details'] = 'Observed data has no F/SIGF columns, required by Phaser. Check file import.'
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def addSearches(self, inputObject):
        inputData = self.container.inputData
        try:
            for i in range(len(inputData.ENSEMBLES)):
                if int(inputData.ENSEMBLES[i].number)>0:
                    inputObject.addSEAR_ENSE_NUM(str(inputData.ENSEMBLES[i].label), int(inputData.ENSEMBLES[i].number))
        except:
            self.appendErrorReport(202)
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    # process one or more output files
    # also writes the XML file, previously done by postProcess()
    def processOutputFiles(self):
        import sys

        import phaser
        resultObject = self.resultObject
        num_sol = len(resultObject.getPdbFiles())
        for i in range(1,num_sol+1):
            xyzout = os.path.join(self.getWorkDirectory(), "PHASER."+str(i)+".pdb")
            if os.path.exists(xyzout):
                self.container.outputData.XYZOUT.append(self.container.outputData.XYZOUT.makeItem())
                self.container.outputData.XYZOUT[-1].setFullPath(xyzout)
                self.container.outputData.XYZOUT[-1].annotation.set('Positioned coordinates for solution '+str(i))
            else:
                self.appendErrorReport(201,xyzout)
                return CPluginScript.FAILED

            hklout = os.path.join(self.getWorkDirectory(), "PHASER."+str(i)+".mtz")
            if os.path.exists(hklout):
                self.container.outputData.HKLOUT.append(self.container.outputData.HKLOUT.makeItem())
                self.container.outputData.HKLOUT[-1].setFullPath(hklout)
            else:
                self.appendErrorReport(201,hklout)
                return CPluginScript.FAILED

        from ccp4i2.core import CCP4XtalData
        self.splitHkloutList(miniMtzsOut=['MAPOUT','DIFMAPOUT','PHASEOUT'],programColumnNames=['FWT,PHWT','DELFWT,PHDELWT','PHIC,FOM'],outputBaseName=['MAPOUT','DIFMAPOUT','PHASEOUT'],outputContentFlags=[1,1,CCP4XtalData.CPhsDataFile.CONTENT_FLAG_PHIFOM],infileList=self.container.outputData.HKLOUT)

        for indx in range(len(self.container.outputData.MAPOUT)):
            self.container.outputData.MAPOUT[indx].annotation.set('Map for solution '+str(indx+1))
            self.container.outputData.MAPOUT[indx].contentFlag.set(1)
            self.container.outputData.MAPOUT[indx].subType.set(1)
            self.container.outputData.DIFMAPOUT[indx].annotation.set('Difference map for solution '+str(indx+1))
            self.container.outputData.DIFMAPOUT[indx].contentFlag.set(1)
            self.container.outputData.DIFMAPOUT[indx].subType.set(2)
            self.container.outputData.PHASEOUT[indx].annotation.set('Calculated phases for solution '+str(indx+1))

        solutions = resultObject.getDotSol()
        if len(solutions) > 0:
            picklePath = str(self.container.outputData.SOLOUT.fullPath)
            with open(picklePath,'wb') as pickleFile:
                try:
                    pickle.dump(solutions, pickleFile)
                except:
                    raise
                    print('Unable to Pickle solutions')
                self.container.outputData.SOLOUT.annotation.set('Solutions from Phaser')

        #Remove warnings and replace with ones parsed from the resultObject
        if len(self.xmlroot.xpath('PhaserWarnings')) > 0:
            phaser_warnings = [wrng for wrng in resultObject.warnings()]
            for warningsElement in self.xmlroot.xpath('PhaserWarnings')[0]:
                if warningsElement.text not in phaser_warnings:
                    advisoriesElement = etree.SubElement(self.xmlroot,'PhaserAdvisories')
                    advisoryElement = etree.SubElement(advisoriesElement,'Advisory')
                    advisoryElement.text = warningsElement.text
            for warningsElement in self.xmlroot.xpath('PhaserWarnings')[0]:
                warningsElement.getparent().remove(warningsElement)
            for warning in phaser_warnings:
                warningsElement = etree.SubElement(self.xmlroot,'PhaserWarnings')
                warningElement = etree.SubElement(warningsElement,'Warning')
                warningElement.text = warning
      
        #Remove old digested summaries and add new ones parsed from the result summary block
        for summaryNode in self.xmlroot.xpath('Summary'):
            summaryNode.getparent().remove(summaryNode)
        summary_buffer = '***'
        for text in resultObject.summary().split('\n'):
            if text.startswith("**********") and not summary_buffer.strip().endswith("***"):
                summaryNode = etree.SubElement(self.xmlroot,'Summary')
                summaryNode.text = summary_buffer
                summary_buffer = ""
            summary_buffer += (text+'\n')
        summaryNode = etree.SubElement(self.xmlroot,'Summary')
        summaryNode.text = summary_buffer
        
        self.flushXML(self.xmlroot)
        return CPluginScript.SUCCEEDED


    def analyseResults(self, results):
        import phaser      
        solutionsNode = etree.SubElement(self.xmlroot,'PhaserMrSolutions')
        
        if not results.foundSolutions():
            node=self.subElementWithNameAndText(solutionsNode,'solutionsFound','False')
        else:
            node=self.subElementWithNameAndText(solutionsNode,'solutionsFound','True')

        solutionListNode = etree.SubElement(solutionsNode,'Solutions')

        isol = 0
        # available items are in phaser/source/phaser/include/mr_set.h
        # example of usage in phaser/source/phaser/phaser/test_reporter.py

        for solution in results.getDotSol():
            if isol == 0:
                for oldNode in self.xmlroot.xpath('//PhaserCurrentBestSolution'):
                    oldNode.getparent().remove(oldNode)
                bcsNode = etree.SubElement(self.xmlroot,'PhaserCurrentBestSolution')
                solutionNode = etree.SubElement(bcsNode,'Solution')
                for nd in solution.KNOWN:
                    componentNode = etree.SubElement(solutionNode,'Component')
                    componentNameNode = etree.SubElement(componentNode,'Name')
                    componentNameNode.text = nd.MODLID
                spaceGroupNode = etree.SubElement(solutionNode,'spaceGroup')

                spaceGroupNode.text = solution.getSpaceGroupName()
                annotationNode = etree.SubElement(solutionNode,'Annotation')
                annotationNode.text = solution.ANNOTATION
                phaser_MR.expandSolutionAnnotation(solutionNode)
            if isol == 0 and self.inputSpaceGroup != solution.getSpaceGroupName():
                self.container.outputData.dataReindexed.set(True)
                warningsElements = self.xmlroot.xpath('PhaserWarnings')
                if len(warningsElements) > 0: warningsElement = warningsElements[0]
                else: warningsElement = etree.SubElement(self.xmlroot,'PhaserWarnings')
                warningElement = etree.SubElement(warningsElement,'Warning')
                warningElement.text = 'Spacegroup of best solution (%s) does not match input data spacegroup (%s)' % (str(solution.getSpaceGroupName()), str(self.inputSpaceGroup))
            elif isol == 0:
                self.container.outputData.dataReindexed.set(False)

            isol += 1
            solutionNode = etree.SubElement(solutionListNode,'Solution')
            node = self.subElementWithNameAndText(solutionNode,'ISOL',str(isol))
            node = self.subElementWithNameAndText(solutionNode,'SPG',str(solution.getSpaceGroupName()))
            #Properties not carried over:'CELL', 'DRMS','MAPCOEFS', 'NEWVRMS', 'RLIST', 'VRMS', 'TMPLT', 'KNOWN'
            for property in ['EQUIV', 'HALL', 'KEEP', 'LLG', 'NUM', 'ORIG_LLG', 'ORIG_NUM', 'ORIG_R', 'PAK', 'R', 'TF', 'TFZ', 'TFZeq']:
                value = getattr(solution,property,None)
                if value is not None:
                    node = self.subElementWithNameAndText(solutionNode, property,str(value))
            for nd in solution.KNOWN:
                componentNode = etree.SubElement(solutionNode,'COMPONENT')
                node = self.subElementWithNameAndText(componentNode,'modlid',str(nd.MODLID))

        #print dir(self.results.getTemplatesForSolution(0))
        #print dir(self.results.getDotSol()[0])
        #print etree.tostring(self.xmlroot, pretty_print=True)
        #print results.getTopPdbFile()

        return CPluginScript.SUCCEEDED
    
    def subElementWithNameAndText(self, parentNode, name, text):
        newNode = etree.SubElement(parentNode, name)
        newNode.text = text
        return newNode

    def flushXML(self, xml):
        import sys

        from lxml import etree
        tmpFilename = self.makeFileName('PROGRAMXML')+'_tmp'
        finalFilename = self.makeFileName('PROGRAMXML')
        with open(tmpFilename,'w') as tmpFile:
            xmlText = etree.tostring(xml, pretty_print=True)
            CCP4Utils.writeXML(tmpFile,xmlText)
            #Here adapt the update frequency to depend on the size of the current XML structure
            xmlUpdateDelay = max(5, int(len(xmlText)/100000))
            self.callbackObject.minimumDelayInSeconds = xmlUpdateDelay
        self.renameFile(tmpFilename, finalFilename)

    def prepareCaptureCPlusPlusStdoutToLog(self):
        # This suggested by Stack Overflow
        # http://stackoverflow.com/questions/24277488/in-python-how-to-capture-the-stdout-from-a-c-shared-library-to-a-variable
        # Highly modified to simply push the output onto the end of "log" file
        # Create pipe and dup2() the write end of it on top of stdout, saving a copy
        # of the old stdout
        self.stdout_fileno = sys.stdout.fileno()
        self.stdout_save = os.dup(self.stdout_fileno)
        
        self.logFile = open(self.makeFileName('LOG'),'a')
        os.dup2(self.logFile.fileno(), self.stdout_fileno)
        #os.close(self.logFile.fileno())
        #os.close(self.stdout_pipe[1])

    def finishCaptureCPlusPlusStdout(self):
        os.dup2(self.stdout_save, self.stdout_fileno)
        os.close(self.stdout_save)
        jobInfo = CCP4Modules.PROJECTSMANAGER().db().getJobInfo(jobId=self.jobId)
        if "jobtitle" in jobInfo and jobInfo["jobtitle"]:
            self.logFile.write(str(jobInfo["jobtitle"])+"\n")
        while "parentjobid" in jobInfo and jobInfo["parentjobid"]:
            jobInfo = CCP4Modules.PROJECTSMANAGER().db().getJobInfo(jobId=jobInfo["parentjobid"])
            if "jobtitle" in jobInfo and jobInfo["jobtitle"]:
                self.logFile.write(str(jobInfo["jobtitle"])+"\n")

        self.logFile.close()
