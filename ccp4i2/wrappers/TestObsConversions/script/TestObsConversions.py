"""
Copyright (C) 2015 Newcastle University
"""

import os
import xml.etree.ElementTree as ET

from ....core import CCP4ErrorHandling
from ....core import CCP4Utils
from ....core.CCP4PluginScript import CPluginScript


class TestObsConversions(CPluginScript):
    TASKTITLE = 'TestObsConversions'
    TASKNAME = 'TestObsConversions'
    TASKVERSION= 0.0
    RUNEXTERNALPROCESS= False
    PERFORMANCECLASS = 'CTestObsConversionsPerformance'
    ERROR_CODES = { 201 : {'description' : 'Failed to turn the provided data object into the requested input representation' },}
    ERROR_CODES = { 202 : {'description' : 'Failed to turn the intermediate data object into the requested output representation' },}

    def __init__(self,parent=None,name=None,workDirectory=None,dummy=False,taskName=None,**kw):
        super(TestObsConversions,self).__init__(parent, name, workDirectory, dummy, taskName, **kw)
        self.xmlroot = ET.Element('TestObsConversions')

    def startProcess(self, command):
        CCP4Utils.writeXml(self.xmlroot, self.makeFileName('PROGRAMXML'))
        return CPluginScript.SUCCEEDED

    def processInputFiles(self):

        if int(self.container.controlParameters.INPUT_REPRESENTATION) == 1:
            inputType = 'F_SIGF_AS_IPAIR'
        elif int(self.container.controlParameters.INPUT_REPRESENTATION) == 2:
            inputType = 'F_SIGF_AS_FPAIR'
        elif int(self.container.controlParameters.INPUT_REPRESENTATION) == 3:
            inputType = 'F_SIGF_AS_IMEAN'
        
        fileNode = ET.SubElement(self.xmlroot,'File')
        typeNode = ET.SubElement(fileNode,'Role')
        typeNode.text = 'Input'
        contentFlagNode = ET.SubElement(fileNode,'ContentFlag')
        contentFlagNode.text = str(getattr(self.container.inputData,inputType).contentFlag)
        columnsNode = ET.SubElement(fileNode,'Columns')
        fileToInterrogate = getattr(self.container.inputData,inputType)
        columnsInFile = fileToInterrogate.getFileContent().getListOfColumns()
        columnsNode.text = str([column.columnLabel.__str__() for column in columnsInFile])
        pathNode = ET.SubElement(fileNode,'Path')
        pathNode.text = str(getattr(self.container.inputData,inputType))
        
        self.hklin,error = self.makeHklin([[inputType, int(self.container.controlParameters.INPUT_REPRESENTATION)]])
        if error.maxSeverity()>CCP4ErrorHandling.Severity.WARNING:
            self.appendErrorReport(201, inputType + " conversion to type number " + str(self.container.controlParameters.INPUT_REPRESENTATION))
            return CPluginScript.FAILED
        intermediateFilePath = os.path.join(self.getWorkDirectory(),'intermediate.mtz')
        os.rename(self.hklin, intermediateFilePath)

        self.container.inputData.F_SIGF_INTERMEDIATE.setFullPath(intermediateFilePath)
        self.container.inputData.F_SIGF_INTERMEDIATE.setContentFlag(reset=True)
        
        inputType = 'F_SIGF_INTERMEDIATE'
        fileNode = ET.SubElement(self.xmlroot,'File')
        typeNode = ET.SubElement(fileNode,'Role')
        typeNode.text = 'Intermediate'
        contentFlagNode = ET.SubElement(fileNode,'ContentFlag')
        contentFlagNode.text = str(getattr(self.container.inputData,inputType).contentFlag)
        columnsNode = ET.SubElement(fileNode,'Columns')
        fileToInterrogate = getattr(self.container.inputData,inputType)
        columnsInFile = fileToInterrogate.getFileContent().getListOfColumns()
        columnsNode.text = str([column.columnLabel.__str__() for column in columnsInFile])
        pathNode = ET.SubElement(fileNode,'Path')
        pathNode.text = str(getattr(self.container.inputData,inputType))

        self.hklin,columns,error = self.makeHklin0([['F_SIGF_INTERMEDIATE', int(self.container.controlParameters.OUTPUT_REPRESENTATION)]])
        if error.maxSeverity()>CCP4ErrorHandling.Severity.WARNING:
            self.appendErrorReport(202, 'F_SIGF_INTERMEDIATE' + " conversion to type number " + str(self.container.controlParameters.OUTPUT_REPRESENTATION))
            return CPluginScript.FAILED

        outputFiles = ['F_SIGF_FINAL']
        outputColumns = [columns]
        error = self.splitHklout(outputFiles,outputColumns,infile=self.hklin)
        if error.maxSeverity()>CCP4ErrorHandling.Severity.WARNING:
            return CPluginScript.FAILED
        self.container.outputData.F_SIGF_FINAL.setContentFlag(reset=True)
        
        inputType = 'F_SIGF_FINAL'
        fileNode = ET.SubElement(self.xmlroot,'File')
        typeNode = ET.SubElement(fileNode,'Role')
        typeNode.text = 'Final'
        contentFlagNode = ET.SubElement(fileNode,'ContentFlag')
        contentFlagNode.text = str(getattr(self.container.outputData,inputType).contentFlag)
        columnsNode = ET.SubElement(fileNode,'Columns')
        fileToInterrogate = getattr(self.container.outputData,inputType)
        columnsInFile = fileToInterrogate.getFileContent().getListOfColumns()
        columnsNode.text = str([column.columnLabel.__str__() for column in columnsInFile])
        self.container.outputData.PERFORMANCEINDICATOR.columnLabelsString = columnsNode.text
        pathNode = ET.SubElement(fileNode,'Path')
        pathNode.text = str(getattr(self.container.outputData,inputType))

        return CPluginScript.SUCCEEDED
