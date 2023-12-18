from __future__ import print_function


from core.CCP4PluginScript import CPluginScript
from PySide2 import QtCore
import os,re,time,sys
from core import CCP4Utils
from xml.etree import ElementTree as ET

class pdbview_edit(CPluginScript):
    
    TASKMODULE = 'model_data_utility'            # Where this plugin will appear on the gui
    TASKTITLE = 'Edit PDB/CIF files by hand with the PdbView program'     # A short title for gui menu
    TASKNAME = 'pdbview_edit'                  # Task name - should be same as class name
    TASKCOMMAND = 'ccp4-python'                          # The command to run the executable
    TASKVERSION= 0.1                                # Version of this plugin
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 9999999.9
    MAINTAINER = 'stuart.mcnicholas@york.ac.uk'

    ERROR_CODES = {  200 : { 'description' : 'Coordinate editor exited with error status' }, 201 : { 'description' : 'Failed in harvest operation' },202 : { 'description' : 'Failed in processOutputFiles' }}

    def makeCommandAndScript(self):
        from core import CCP4Utils
        self.dropDir = os.path.join(self.workDirectory,'CCP4MG_FILE_DROP')
        if not os.path.exists(self.dropDir):
          try:
            os.mkdir(self.dropDir)
          except:
            self.dropDir = self.workDirectory
            print('Could not make dropDir reset to',self.dropDir)

        # Make a script file with additional menu options to save to i2
        self.mgStatusPath = os.path.normpath(os.path.join(self.workDirectory,'script.mgpic.xml'))
        if sys.platform == 'win32':
          self.mgStatusPath = re.sub(r'\\\\',r'\\',self.mgStatusPath)
        # Declare script text then re.sub in the variables
        if sys.platform == "win32":
          i2dir = CCP4Utils.getCCP4I2Dir().replace('\\','/')
        else:
          i2dir = CCP4Utils.getCCP4I2Dir()

        clArgs = [ ]

        origInterface_fname = os.path.join(os.path.dirname(__file__),"PdbView-i2.py")
        if not os.path.exists(origInterface_fname):
          origInterface_fname = os.path.join(os.path.dirname(__file__),"PdbView-i2.pyc")
        if not os.path.exists(origInterface_fname):
          origInterface_fname = os.path.join(os.path.dirname(__file__),"PdbView-i2.pyo")
        origInterface = open(origInterface_fname)
        origInterface_s = origInterface.read()
        origInterface.close()

        #interface_s = origInterface_s.replace('XXXXX_WORK_DIR_XXXXX',self.workDirectory)
        interface_s = origInterface_s
        interface_fname = os.path.join(self.workDirectory,"PdbView-i2.py")
        interface = open(interface_fname,"w+")
        interface.write(interface_s)
        interface.close()
        
        clArgs.append(interface_fname)
        clArgs.append(self.workDirectory)

        if self.container.inputData.XYZIN_LIST.isSet():
            try:
                iFile = 1
                for XYZIN in self.container.inputData.XYZIN_LIST:
                    if os.path.isfile(XYZIN.__str__()):
                        clArgs.append(XYZIN.__str__())
                    else:
                        print('pdbview_edit.makeCommandAndScript XYZIN does not exist:',XYZIN.__str__())
                    iFile += 1
            except:
                #an issue with the existence of files
                pass

        self.appendCommandLine(clArgs)
        # Use Qt class to watch the drop directory
        self.fileSystemWatcher = QtCore.QFileSystemWatcher(parent=self)
        self.fileSystemWatcher.addPath(self.dropDir)
        self.fileSystemWatcher.directoryChanged.connect(self.handleFileDrop)

        return CPluginScript.SUCCEEDED


    def numberOfOutputFiles(self):
        import glob
        outList = glob.glob(os.path.normpath(os.path.join(self.dropDir,'output*.pdb')))
        #print 'numberOfOutputFiles outList',os.path.join(self.dropDir,'output*.pdb'),outList
        #print 'numberOfOutputFiles xmlList',glob.glob(os.path.normpath(os.path.join(self.workDirectory,'*.xml')))
        maxIndx = 0
        for f in outList:
           fpath,fname = os.path.split(f)
           #print 'numberOfOutputFiles  fpath,fname', fpath,fname
           maxIndx =  max(maxIndx,int(fname[6:-4]))
        return maxIndx

    @QtCore.Slot(str)
    def handleFileDrop(self,directory):
        import time,glob
        print('pdbview_edit',time.time())
        print('pdbview_edit',glob.glob(os.path.join(self.workDirectory,'*.*')))
        #print 'handleFileDrop',directory
        #Note that I don't copy the file to the appropriate xyzout filename here, since the file may not yet
        #be closed and/or flushed

        
    def processOutputFiles(self):
        try:
            # First up import PDB files that have been output
            
            import os, glob, shutil
            globPath = os.path.normpath(os.path.join(self.dropDir,'output*.pdb'))
            outList = glob.glob(globPath)
            
            xyzoutList = self.container.outputData.XYZOUT
            for outputPDB in outList:
                fpath,fname = os.path.split(outputPDB)
                iFile = int(fname[6:-4])
                xyzoutList.append(xyzoutList.makeItem())
                outputFilePath = os.path.normpath(os.path.join(self.workDirectory,'XYZOUT_'+str(iFile)+'-coordinates.pdb'))
                shutil.copyfile(outputPDB, outputFilePath)
                xyzoutList[-1].setFullPath(outputFilePath)
                xyzoutList[-1].annotation = "Coordinate editor output file number"+str(iFile)
                xyzoutList[-1].subType = 1

            # Create a trivial xml output file
            self.xmlroot = ET.Element('pdbview_edit')
            e = ET.Element('number_output_files')
            e.text = str(self.numberOfOutputFiles())
            self.xmlroot.append(e)
            
            #Separate out here activity to attempt merge into project dictionary....this seems flakey,
            #but is needed for ongoing work, so I am making it give an report a warning in case of failure, rather than
            #offer the sad face of doom
        except:
            exc_type, exc_value,exc_tb = sys.exc_info()[:3]
            sys.stderr.write(str(exc_type)+'\n')
            sys.stderr.write(str(exc_value)+'\n')

            self.appendErrorReport(202,'Data harvesting failed')
            
        CCP4Utils.saveEtreeToFile(self.xmlroot,self.makeFileName('PROGRAMXML'))
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

