"""
Copyright (C) 2014 Newcastle University
"""

import os
import xml.etree.ElementTree as ET

from ....core import CCP4ErrorHandling
from ....core.CCP4PluginScript import CPluginScript


class pointless_reindexToMatch(CPluginScript):
    TASKNAME = 'pointless_reindexToMatch'
    TASKCOMMAND = 'pointless'
    TASKVERSION= 0.0
    ASYNCHRONOUS = False

    ERROR_CODES = { 201 : {'description' : 'Failed to make input files' }, 202 : {'description' : 'Failed to make output files' }, 203 : {'description' : 'Failed to identify a best reindexing' }}

    def formatCellLength(self, p):
        return "%7.1f" % float(p)

    def formatCellAngle(self, p):
        return "%7.1f" % (float(p) if float(p) > 3.141592653589793 else float(p) * 57.29577951308233)

    def shortformatCell(self, cell):
        s = ""
        s += self.formatCellLength(cell.a).strip()+','
        s += self.formatCellLength(cell.b).strip()+','
        s += self.formatCellLength(cell.c).strip()+','
        s += self.formatCellAngle(cell.alpha).strip()+','
        s += self.formatCellAngle(cell.beta).strip()+','
        s += self.formatCellAngle(cell.gamma).strip()
        return s

    def makeCommandAndScript(self):
        par = self.container.controlParameters
        self.mergedFilename = os.path.join(self.workDirectory,'MergedToReindex.mtz')
        # XMLOUT on command line so that syntax errors go into it
        self.appendCommandLine(['XMLOUT',str( self.makeFileName( 'PROGRAMXML' ) )])
        if str(self.container.controlParameters.REFERENCE) == 'HKLIN_FOBS_REF':
            self.appendCommandLine(['HKLREF',str(self.container.inputData.HKLIN_FOBS_REF.fullPath)])
        elif str(self.container.controlParameters.REFERENCE) == 'HKLIN_FC_REF':
            self.appendCommandLine(['HKLREF',str(self.container.inputData.HKLIN_FC_REF.fullPath)])
        elif str(self.container.controlParameters.REFERENCE) == 'HKLIN_FMAP_REF':
            self.appendCommandLine(['HKLREF',str(self.container.inputData.HKLIN_FMAP_REF.fullPath)])
        elif str(self.container.controlParameters.REFERENCE) == 'XYZIN_REF':
            self.appendCommandLine(['XYZIN',str(self.container.inputData.XYZIN_REF.fullPath)])
        elif str(self.container.controlParameters.REFERENCE) == 'SPECIFY':
            if self.container.controlParameters.CHOOSE_SPACEGROUP.isSet():
                if self.container.controlParameters.CHOOSE_SPACEGROUP.isSet():
                    self.appendCommandScript('spacegroup '+\
                                             str(self.container.controlParameters.CHOOSE_SPACEGROUP))
            if self.container.controlParameters.USE_REINDEX:
                self.appendCommandScript("reindex %s, %s, %s" % (par.REINDEX_OPERATOR.h,par.REINDEX_OPERATOR.k,par.REINDEX_OPERATOR.l))
        elif str(self.container.controlParameters.REFERENCE) == 'EXPAND':
            self.appendCommandScript('expand')
        elif self.container.controlParameters.LATTICE_CENTERING.isSet():
            lattype = str(self.container.controlParameters.LATTICE_CENTERING)
            if lattype != 'P':
                self.appendCommandScript('lattice '+lattype)
                
        mergedFilename = os.path.join(self.workDirectory,'MergedToReindex.mtz')
        self.appendCommandLine(['HKLIN',mergedFilename])

        if str(self.container.controlParameters.REFERENCE) != 'ANALYSE':
            reindexedFilename = os.path.join(self.workDirectory,'Reindexed.mtz')
            self.appendCommandLine(['HKLOUT',reindexedFilename])

        self.appendCommandScript('end')

        return CPluginScript.SUCCEEDED

    def processInputFiles(self):
        try:
            self.container.inputData.F_SIGF.loadFile()
            colList = self.container.inputData.F_SIGF.fileContent.getListOfColumns()
            print('\n\n\n\n F_SIGF',self.container.inputData.F_SIGF)
            print('\n\n\n\n fileContent',self.container.inputData.F_SIGF.fileContent)
            print('\n\n\n\n colList',colList)
            fSigFLabelList = [str(column.get('columnLabel')) for column in colList]
            self.fSigFLabelsString = ','.join(fSigFLabelList)
            fileList = []
            fileList += [(str(self.container.inputData.F_SIGF.fullPath), self.fSigFLabelsString, self.fSigFLabelsString)]
            if self.container.inputData.FREERFLAG.isSet():
                colList = self.container.inputData.FREERFLAG.fileContent.getListOfColumns()
                freeRFlagLabelList = [str(column.columnLabel) for column in colList]
                self.freeRFlagLabelsString = ','.join(freeRFlagLabelList)
                fileList += [(str(self.container.inputData.FREERFLAG.fullPath), self.freeRFlagLabelsString, self.freeRFlagLabelsString)]
            mergedFilename = os.path.join(self.workDirectory,'MergedToReindex.mtz')
            rv = self.joinMtz(mergedFilename, fileList)
            return rv
        except:
            self.appendErrorReport(201,'Pointless_reindexToMatch: processInputFiles: exception in process')
            return CPluginScript.FAILED

    def postProcessCheck(self, processId):
        status, exitStatus, exitCode = super(pointless_reindexToMatch,self).postProcessCheck(processId)
        if (exitStatus != CPluginScript.SUCCEEDED) or (exitCode != 0):
            print("postProcessCheck FAIL")
            return CPluginScript.UNSATISFACTORY, exitStatus, exitCode
        return status, exitStatus, exitCode

    def processOutputFiles(self):
        print('#PRM processOutputFiles')
        try:
            taskoption = str(self.container.controlParameters.REFERENCE)
            if taskoption != 'ANALYSE':
                # deal with output reflection files, unless ANALYSE
                self.processHKLOUT()
            # Process XML 
            try:
                with open(self.makeFileName('PROGRAMXML'),'r') as unfixedXMLFile:
                    text = unfixedXMLFile.read()
                try:
                    rootNode = ET.fromstring(text)
                except:
                    #MN Pointless's XML can be corrupt...I've spotted this in instances where the XYZIN_REF has been set, and
                    #the corresponding XML closing tag is split across two lines.
                    #Also, the attributes listed within markup may not be space delimited
                    #Here is a kludged fix
                    fixedText = ''
                    inMarkup = False
                    inQuote = False
                    for iChar in range(len(text)-1):
                        character = text[iChar:iChar+1]
                        if character == '<': inMarkup = True
                        elif character == '>': inMarkup = False
                        if inMarkup and character == '"':
                            if inQuote: character += ' '
                            inQuote = not inQuote
                        if (inMarkup and character != '\n') or not inMarkup: fixedText += character
                    with open(self.makeFileName('PROGRAMXML'),'w') as fixedXMLFile:
                        fixedXMLFile.write(fixedText)
                    rootNode = ET.fromstring(fixedText)
                #print( '#PRM rootNode',rootNode)
                
                if taskoption != 'ANALYSE' and taskoption != 'LATTICE' \
                       and taskoption != 'EXPAND':
                    bestReindexNodes = rootNode.xpath('//BestReindex')
                    scoreCountNodes = rootNode.xpath('//ScoreCount')
                    copyMessageNodes = rootNode.xpath('//CopyMessage')
                    if len(bestReindexNodes) == 0 and len(scoreCountNodes) == 0 and len(copyMessageNodes) == 0:
                        self.container.outputData.BestReindexIdentified = False
                        self.appendErrorReport(203, 'No reindex found')
                        print("reindex fail")
                        return CPluginScript.FAILED
                    else:
                        self.container.outputData.BestReindexIdentified = True
            except:
                print('Unable to complete check for best reindex found')
                self.appendErrorReport(203, 'Unable to complete check for best reindex found')
                return CPluginScript.FAILED
                    
        except:
            self.appendErrorReport(202,'Pointless_reindexToMatch: processOutputFiles exception in process')
            return CPluginScript.FAILED

        print("processOutput success")
        return CPluginScript.SUCCEEDED
    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
    def processHKLOUT(self):
        outputFilesList = ['F_SIGF_OUT']
        outputColumnsList  = [self.fSigFLabelsString]
        outputContentsList = [1]
        if self.container.inputData.FREERFLAG.isSet():
            outputFilesList += ['FREERFLAG_OUT']
            outputColumnsList  += [self.freeRFlagLabelsString]
            outputContentsList += [1]
        reindexedFilename = os.path.join(self.workDirectory,'Reindexed.mtz')
        # set contentFlag here so that splitmtz (as called from splitHklout) will know what columns to expect
        self.container.outputData.F_SIGF_OUT.contentFlag.set(self.container.inputData.F_SIGF.contentFlag)
        print("** splitHklout ", outputFilesList, outputColumnsList, reindexedFilename)
        print("CF", self.container.outputData.F_SIGF_OUT.contentFlag)
        print("FSO", self.container.outputData.F_SIGF_OUT)


        error = self.splitHklout(outputFilesList,outputColumnsList,infile = reindexedFilename)
        if error.maxSeverity()>CCP4ErrorHandling.Severity.WARNING:
            self.appendErrorReport(202,'Pointless_reindexToMatch: error in splitting file')
            return CPluginScript.FAILED
        else:
            # update the annotations
            self.container.outputData.F_SIGF_OUT.contentFlag.set(self.container.inputData.F_SIGF.contentFlag)
            # print "FSF", self.container.outputData.F_SIGF_OUT.fileContent
                
            #MN Did filecontent change from having spaceGroupName attribute to spaceGroup ????
            try: sgname = self.container.outputData.F_SIGF_OUT.fileContent.spaceGroupName
            except: sgname = str(self.container.outputData.F_SIGF_OUT.fileContent.spaceGroup)
            
            highres = "%7.2f" % self.container.outputData.F_SIGF_OUT.fileContent.resolutionRange.high
            title ='SG:'+str(sgname).strip()+';Resolution:'+highres.strip()+\
                    ";Cell:"+self.shortformatCell(self.container.outputData.F_SIGF_OUT.fileContent.cell)
            if self.container.controlParameters.REINDEX_OPERATOR.isSet():
                par = self.container.controlParameters
                reindexop = "Reindexed:"+("[%s, %s, %s]" % (par.REINDEX_OPERATOR.h,par.REINDEX_OPERATOR.k,par.REINDEX_OPERATOR.l)).strip()
                title = reindexop + ";" + title
            else:
                title = "New" + title
                
            self.container.outputData.F_SIGF_OUT.annotation = \
              str(self.container.inputData.F_SIGF.annotation) + " " + title

            self.container.outputData.F_SIGF_OUT.subType = self.container.inputData.F_SIGF.subType
            if self.container.inputData.FREERFLAG.isSet():
                self.container.outputData.FREERFLAG_OUT.annotation = title
