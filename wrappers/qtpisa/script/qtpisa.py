from __future__ import print_function


from core.CCP4PluginScript import CPluginScript
from PySide2 import QtCore
import os,re,time,sys
#from lxml import etree
from xml.etree import ElementTree as ET
from core import CCP4Utils

class qtpisa(CPluginScript):
    
    TASKMODULE = 'validation'            # Where this plugin will appear on the gui
    TASKTITLE = 'Analyze bilogical units with qtpisa'     # A short title for gui menu
    TASKNAME = 'qtpisa'                  # Task name - should be same as class name
    TASKCOMMAND = 'qtpisa'                          # The command to run the executable
    TASKVERSION= 0.1                                # Version of this plugin
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 9999999.9
    MAINTAINER = 'eugene.krissinel@stfc.ac.uk'

    ERROR_CODES = {  200 : { 'description' : 'QtPisa exited with error status' }, 201 : { 'description' : 'Failed in harvest operation' },202 : { 'description' : 'Failed in processOutputFiles' }}

    def makeCommandAndScript(self):
        from core import CCP4Utils
        self.dropDir = os.path.join(self.workDirectory,'QTPISA_FILE_DROP')
        if not os.path.exists(self.dropDir):
          try:
            os.mkdir(self.dropDir)
          except:
            self.dropDir = self.workDirectory
            print('Could not make dropDir reset to',self.dropDir)

        clArgs = [ '-i2dir',self.dropDir,str(self.container.inputData.XYZIN) ]

        self.appendCommandLine(clArgs)
        # Use Qt class to watch the drop directory
        self.fileSystemWatcher = QtCore.QFileSystemWatcher(parent=self)
        self.fileSystemWatcher.addPath(self.dropDir)
        self.fileSystemWatcher.directoryChanged.connect(self.handleFileDrop)

        return CPluginScript.SUCCEEDED


    def numberOfOutputFiles(self):
        import glob
        outList = glob.glob(os.path.normpath(os.path.join(self.dropDir,'*.xml')))
        #print 'numberOfOutputFiles outList',os.path.join(self.dropDir,'*.xml'),outList
        #print 'numberOfOutputFiles xmlList',glob.glob(os.path.normpath(os.path.join(self.workDirectory,'*.xml')))
        return len(outList)

    @QtCore.Slot(str)
    def handleFileDrop(self,directory):
        import time,glob
        print('qtpisa',time.time())
        print('qtpisa',glob.glob(os.path.join(self.workDirectory,'*.*')))
        #print 'handleFileDrop',directory
        #Note that I don't copy the file to the appropriate xyzout filename here, since the file may not yet
        #be closed and/or flushed
  
    def processOutputFiles(self):
        try:
            # First up import PDB files that have been output
            
            import os, glob, shutil
            globPath = os.path.normpath(os.path.join(self.dropDir,'*.xml'))
            outList = glob.glob(globPath)
            

            report_xml = ET.Element("qtpisa_report")

            iFile = 1
            xyzoutList = self.container.outputData.XYZOUT
            for outputXML in outList:
                fpath,fname = os.path.split(outputXML)
                xyzoutList.append(xyzoutList.makeItem())

                parser = ET.XMLParser()
                f = open(outputXML)
                s = f.read()
                f.close()
                tree = ET.fromstring(s, parser)

                pisa_type = tree.findall(".//type")[0].text
                pisa_no = tree.findall(".//ser_no")[0].text

                pisa_file = tree.findall(".//file")[0].text

                outputFilePath = os.path.normpath(os.path.join(self.workDirectory,'XYZOUT_'+pisa_no+'-'+pisa_type+'-coordinates.pdb'))

                shutil.copyfile(pisa_file, outputFilePath)

                xyzoutList[-1].setFullPath(outputFilePath)
                xyzoutList[-1].annotation = "QtPISA "+pisa_type+" number "+pisa_no
                xyzoutList[-1].subType = 1

                iFile += 1
                report_xml.append(tree)

            # Create a trivial xml output file
            self.xmlroot = ET.Element('qtpisa')
            e = ET.Element('number_output_files')
            e.text = str(self.numberOfOutputFiles())
            self.xmlroot.append(e)
            self.xmlroot.append(report_xml)
            
            #Separate out here activity to attempt merge into project dictionary....this seems flakey,
            #but is needed for ongoing work, so I am making it give an report a warning in case of failure, rather than
            #offer the sad face of doom
        except:
            exc_type, exc_value,exc_tb = sys.exc_info()[:3]
            sys.stderr.write(str(exc_type)+'\n')
            sys.stderr.write(str(exc_value)+'\n')

            self.appendErrorReport(202,'Data harvesting failed')
            
        with open(self.makeFileName("PROGRAMXML"),"w") as programXMLFile:
            CCP4Utils.writeXML(programXMLFile,ET.tostring(self.xmlroot))

        if ( len(outList) ) > 0:
          return CPluginScript.SUCCEEDED
        else:
          return CPluginScript.MARK_TO_DELETE

    def addReportWarning(self, text):
        warningsNode = None
        warningsNodes = self.xmlroot.findall('.//Warnings')
        if len(warningsNodes) == 0: warningsNode = ET.SubElement(self.xmlroot, 'Warnings')
        else: warningsNode = warningsNodes[0]
        warningNode = ET.SubElement(warningsNode,'Warning')
        warningNode.text = text

