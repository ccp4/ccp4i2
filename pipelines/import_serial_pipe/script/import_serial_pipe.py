import os
import re
import shutil
import subprocess
import glob
import json
from lxml import etree
from baselayer import QtCore

from core.CCP4PluginScript import CPluginScript
from core import CCP4XtalData
from core import CCP4ErrorHandling
from core import CCP4Utils
from core import CCP4Modules
# try:
#     from baselayer.QtCore import Slot
# except:
#     from PyQt4.QtCore import pyqtSlot as Slot

class import_serial_pipe(CPluginScript):
    TASKMODULE = 'data_entry'         # GIU menu location
    TASKTITLE = 'Import Serial Pipeline'       # Short title for GUI
    TASKNAME = 'import_serial_pipe'        # Task name - same as class name
    TASKVERSION = 1.1                 # plugin version
    COMTEMPLATE = None                # The program com file template
    COMTEMPLATEFILE = None            # Name of file containing com file template
    # PERFORMANCECLASS = ''
    ASYNCHRONOUS = False
    MAINTAINER = 'martin.maly@soton.ac.uk'
    # PURGESEARCHLIST = [ [ 'HKLIN*.mtz' , 1 ],
    #                     ['aimless_pipe%*/HKLOUT*.mtz', 1],
    #                     [ 'hklout.mtz' , 5 ],    
    #                     [ 'HKLOUT.mtz' , 5 ]
    #                   ]
    # ERROR_CODES = { 101 : {'description' : 'Blank for now, may need this ',
    #                        'severity':CCP4ErrorHandling.SEVERITY_ERROR } }

    def __init__(self, *args, **kwargs):
        # self.seqin = None
        self.hklin = None
        self.hklin1 = None
        self.hklin2 = None
        self.stream = None
        CPluginScript.__init__(self, *args, **kwargs)

    def process(self):
        self.importSerialProcess = self.makePluginObject('import_serial')
        self.importSerialProcess.container.inputParameters = self.container.inputParameters
        self.importSerialProcess.container.inputData.copyData(otherContainer=self.container.inputData)
        self.importSerialProcess.container.outputData.copyData(otherContainer=self.container.outputData)
        #self.importSerialProcess.container.outputData.copyData(otherContainer=self.container.outputData,dataList=['HKLOUT','OBSOUT'])
        # XML output 'program.xml' is produced by the command line application
        self.xmlout = self.makeFileName('PROGRAMXML')
        # rootNode = etree.Element("import_serial")
        # Save xml
        #xmlfile = open(self.xmlout, 'wb')
        #xmlString= etree.tostring(rootNode, pretty_print=True)
        #xmlfile.write(xmlString)
        #xmlfile.close()
        ###self.connectSignal(self.importSerialProcess, 'finished', self.process2)
        print("import_serial_pipe: import_serial starting")
        status = self.importSerialProcess.process()
        print("import_serial_pipe: import_serial ended")

        # Run aimless for a report on data quality
        self.aimlessPipe = self.makePluginObject('aimless_pipe', pluginTitle='DR run for data analysis')
        mergedList = self.aimlessPipe.container.inputData.UNMERGEDFILES
        if len(mergedList)==0: mergedList.addItem()
        # Always do analysis on the file which is saved as pipeline output
        mergedList[0].file.set(self.importSerialProcess.container.outputData.HKLOUT.__str__())
        ##xname = self.filteredName(str(self.container.inputData.CRYSTALNAME), 'X')
        ##dname = self.filteredName(str(self.container.inputData.DATASETNAME), 'D')
        ##unmergedList[0].crystalName.set('CRYSTAL')
        ##unmergedList[0].dataset.set('DATASET')
        ##unmergedList[0].cell.set(self.container.inputParameters.CELL)
        # parameters for Pointless
        self.aimlessPipe.container.controlParameters.MODE = 'CHOOSE'
        self.aimlessPipe.container.controlParameters.CHOOSE_MODE = 'SPACEGROUP'
        self.aimlessPipe.container.controlParameters.CHOOSE_SPACEGROUP = \
            str(self.container.inputParameters.SPACEGROUP).strip()
        # parameters for Aimless
        self.aimlessPipe.container.controlParameters.SCALING_PROTOCOL = 'CONSTANT'
        self.aimlessPipe.container.controlParameters.ONLYMERGE = True
        self.aimlessPipe.container.controlParameters.ANALYSIS_MODE = True
        self.aimlessPipe.container.controlParameters.OUTPUT_UNMERGED = False
        self.aimlessPipe.container.controlParameters.SDCORRECTION_OVERRIDE = True
        self.aimlessPipe.container.controlParameters.SDCORRECTION_REFINE = False
        self.aimlessPipe.container.controlParameters.SDCORRECTION_SET = True
        self.aimlessPipe.container.controlParameters.SDCORRECTION_SDFAC = 1.0
        self.aimlessPipe.container.controlParameters.SDCORRECTION_SDB = 0.0
        self.aimlessPipe.container.controlParameters.SDCORRECTION_SDADD = 0.0

        ####self.importXML = etree.Element('IMPORT_LOG')  # information about the import step
        #  +1 if intensity, -1 if amplitude, 0 if unknown
        # self.isintensity = 0
        ####self.makeReportXML(self.importXML)  # add initial stuff for XML into self.importXML
        ####print("aim 3")
        ####self.outputLogXML(self.importXML)  # send self.importXML to program.xml

        ####tempXML = self.importXML
        ####self.addElement(tempXML, "DRPIPE_RUNNING", "True") 
        ####self.outputLogXML(tempXML)  # SEND self.importXML to program.xml
        #  Start data reduction
        print("import_serial_pipe: import_merged aimless_pipe starting")
        # self.connectSignal(self.aimlessPipe, 'finished', self.nearlyDone)
        self.aimlessPipe.process()
        print("import_serial_pipe: import_merged aimless_pipe end")
        # return CPluginScript.SUCCEEDED
        # def processOutputFiles(self):
        # shutil.copyfile(str(self.aimlessPipe.container.outputData.HKLOUT[0]), str(self.container.outputData.HKLOUT))
        shutil.copyfile(str(self.aimlessPipe.container.outputData.IMEANOUT[0]), str(self.container.outputData.HKLOUT))
        self.container.outputData.HKLOUT.set(self.aimlessPipe.container.outputData.HKLOUT)
        self.container.outputData.HKLOUT.setAnnotation("Merged intensities")
        # self.container.outputData.HKLOUT.contentFlag = 3 # CCP4XtalData.CObsDataFile.CONTENT_FLAG_IMEAN
        self.container.outputData.HKLOUT.setContentFlag(reset=True) # this will set contentFlag to 3 (CONTENT_FLAG_IMEAN)

        # Save xml
        #xmlfile = open(self.xmlout, 'wb')
        #xmlString= etree.tostring(root, pretty_print=True)
        #xmlfile.write(xmlString)
        #xmlfile.close()
        self.xmlout = self.makeFileName('PROGRAMXML')
        importSerialProcessXMLpath = self.importSerialProcess.makeFileName('PROGRAMXML')
        #importSerialProcessEtree = etree.parse(importSerialProcessXMLpath)
        #importSerialProcessRoot = importSerialProcessEtree.getroot()
        #self.outputLogXML(self.xmlout, importSerialProcessRoot)  # and output it
        importSerialPipelineXMLpath = self.makeFileName('PROGRAMXML')
        shutil.copyfile(str(importSerialProcessXMLpath), str(importSerialPipelineXMLpath))
        self.reportStatus(CPluginScript.SUCCEEDED)
        # return CPluginScript.SUCCEEDED


        #def outputLogXML(self, x1XML, x2XML=None):
        #'output x1XML and optionally x2XML to program.xml'
        ##print "outputLogXML", x1XML
        #rootXML = etree.Element('import_serial') # Global XML for everything
        #rootXML.append(x1XML)
        #if x2XML is not None:
        #    rootXML.append(x2XML)
        #with open(self.makeFileName('PROGRAMXML'), "w") as outputXML:
        #    CCP4Utils.writeXML(outputXML, etree.tostring(rootXML, pretty_print=True))
