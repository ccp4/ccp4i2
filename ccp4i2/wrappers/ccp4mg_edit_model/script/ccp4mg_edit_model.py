from stat import S_ISREG, ST_CTIME, ST_MODE
import glob
import os
import re
import shutil
import sys
import time

from lxml import etree
from PySide2 import QtCore

from ....core import CCP4Utils
from ....core.CCP4MgImports import PhmmerReportNoGui
from ....core.CCP4PluginScript import CPluginScript


class ccp4mg_edit_model(CPluginScript):
    
    TASKTITLE = 'Interactive model preparation - CCP4mg and MrBUMP'     # A short title for gui menu
    TASKNAME = 'ccp4mg_edit_model'                  # Task name - should be same as class name
    TASKCOMMAND = 'ccp4mg'                          # The command to run the executable
    TASKVERSION= 0.1                                # Version of this plugin
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 9999999.9
    MAINTAINER = 'stuart.mcnicholas@york.ac.uk'

    ERROR_CODES = {  200 : { 'description' : 'CCP4MG exited with error status' }, 201 : { 'description' : 'Failed in harvest operation' },202 : { 'description' : 'Failed in processOutputFiles' }}

    def makeCommandAndScript(self):
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


        status_xml = ""

        NSMAP = {'xsi':"http://www.w3.org/2001/XMLSchema-instance"}
        NS = NSMAP['xsi']
        location_attribute = '{%s}noNamespaceSchemaLocation' % NS
        tree = etree.Element("CCP4MG_Status",nsmap = NSMAP,attrib={location_attribute: 'http://www.ysbl.york.ac.uk/~mcnicholas/schema/CCP4MGApplicationOutput.xsd'})

        status_xml += etree.tostring(tree,encoding='utf-8', pretty_print=True).decode("utf-8")

        print("Writing",self.mgStatusPath)
        print(status_xml)

        with open(self.mgStatusPath, 'wb') as xmlFile:
             xmlFile.write(bytes(status_xml,"utf-8"))
        
        if sys.platform == 'win32':
            clArgs = [ '-noredirect','-norestore',self.mgStatusPath ]
        else:
            clArgs = [ '-norestore',self.mgStatusPath ]

        origInterface_fname = os.path.join(CCP4Utils.getCCP4I2Dir(),"wrappers","ccp4mg_edit_model","script","ccp4i2CCP4MGInterface.py")
        print(origInterface_fname)
        origInterface = open(origInterface_fname)
        origInterface_s = origInterface.read()
        origInterface.close()

        interface_s = origInterface_s
        interface_fname = os.path.join(self.workDirectory,"ccp4i2CCP4MGInterface.py")
        interface = open(interface_fname,"w+")
        interface.write(interface_s)
        interface.close()
        
        clArgs.extend(['-scr',interface_fname])

        asuin = self.container.inputData.ASUIN
        if asuin.isSet():
          for idx in range(len(asuin.fileContent.seqList)):
            if asuin.isSelected(asuin.fileContent.seqList[idx]):
              seqFile = os.path.join(self.workDirectory,'SEQIN_'+str(asuin.fileContent.seqList[idx].name)+'.fasta')
              asuin.writeFasta(seqFile,idx)
              clArgs.append(seqFile)

        if self.container.inputData.PHMMERCUTOFF:
               clArgs.extend(["-scriptArg","mrBumpCutoff="+self.container.inputData.PHMMERCUTOFF.__str__()])

        if self.container.inputData.MRMAX:
            clArgs.extend(["-scriptArg","mrBumpMRNUM="+self.container.inputData.MRMAX.__str__()])
        if self.container.inputData.PDBLOCAL.isSet():
            clArgs.extend(["-scriptArg","mrBumpUsePDBLocal="+self.container.inputData.PDBLOCAL.__str__()])
        if self.container.inputData.HHPREDIN.isSet():
            clArgs.extend(["-scriptArg","mrBumpHHRfile="+self.container.inputData.HHPREDIN.__str__()])

        if self.container.inputData.SEARCH_PDB:
            if self.container.inputData.REDUNDANCYLEVEL:
                if str(self.container.inputData.REDUNDANCYLEVEL) == '110':
                   clArgs.extend(["-scriptArg","mrBumpSim=All"])
                elif str(self.container.inputData.REDUNDANCYLEVEL) == '100':
                   clArgs.extend(["-scriptArg","mrBumpSim=100"])
                elif str(self.container.inputData.REDUNDANCYLEVEL) == '95':
                   clArgs.extend(["-scriptArg","mrBumpSim=95"])
                elif str(self.container.inputData.REDUNDANCYLEVEL) == '90':
                   clArgs.extend(["-scriptArg","mrBumpSim=90"])
                elif str(self.container.inputData.REDUNDANCYLEVEL) == '70':
                   clArgs.extend(["-scriptArg","mrBumpSim=70"])
                elif str(self.container.inputData.REDUNDANCYLEVEL) == '50':
                   clArgs.extend(["-scriptArg","mrBumpSim=50"])
        else:
            clArgs.extend(["-scriptArg","mrBumpSim=0"])

        if self.container.inputData.SEARCH_AFDB:
            if self.container.inputData.AFDBLEVEL is not None:
                if str(self.container.inputData.AFDBLEVEL) == '90':
                   clArgs.extend(["-scriptArg","mrBumpAFLEVEL=90"])
                elif str(self.container.inputData.AFDBLEVEL) == '80':
                   clArgs.extend(["-scriptArg","mrBumpAFLEVEL=80"])
                elif str(self.container.inputData.AFDBLEVEL) == '70':
                   clArgs.extend(["-scriptArg","mrBumpAFLEVEL=70"])
                elif str(self.container.inputData.AFDBLEVEL) == '60':
                   clArgs.extend(["-scriptArg","mrBumpAFLEVEL=60"])
                elif str(self.container.inputData.AFDBLEVEL) == '50':
                   clArgs.extend(["-scriptArg","mrBumpAFLEVEL=50"])
                elif str(self.container.inputData.AFDBLEVEL) == '40':
                   clArgs.extend(["-scriptArg","mrBumpAFLEVEL=40"])
                elif str(self.container.inputData.AFDBLEVEL) == '30':
                   clArgs.extend(["-scriptArg","mrBumpAFLEVEL=30"])
                elif str(self.container.inputData.AFDBLEVEL) == '20':
                   clArgs.extend(["-scriptArg","mrBumpAFLEVEL=20"])
                elif str(self.container.inputData.AFDBLEVEL) == '10':
                   clArgs.extend(["-scriptArg","mrBumpAFLEVEL=10"])
                elif str(self.container.inputData.AFDBLEVEL) == '0':
                   clArgs.extend(["-scriptArg","mrBumpAFLEVEL=0"])
    
        #TODO - Add options for this in future.
        #clArgs.extend(["-scriptArg","mrBumpMRNUM=10"])

        self.appendCommandLine(clArgs)

        # Use Qt class to watch the drop directory
        self.fileSystemWatcher = QtCore.QFileSystemWatcher(parent=self)
        self.fileSystemWatcher.addPath(self.dropDir)
        self.fileSystemWatcher.directoryChanged.connect(self.handleFileDrop)

        return CPluginScript.SUCCEEDED


    def numberOfOutputFiles(self):
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
        print('ccp4mg_edit_model',time.time())
        print('ccp4mg_edit_model',glob.glob(os.path.join(self.workDirectory,'*.*')))
        #print 'handleFileDrop',directory
        #Note that I don't copy the file to the appropriate xyzout filename here, since the file may not yet
        #be closed and/or flushed

        
    def processOutputFiles(self):
        try:
            # First up import PDB files that have been output

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
                xyzoutList[-1].annotation.set("MrBUMP-CCP4MG output file number "+str(iFile))
                xyzoutList[-1].subType.set(1)

            # Create an xml output file and harvest the MrBUMP log files.
            self.xmlroot = etree.Element('ccp4mg_edit_model')
            e = etree.Element('number_output_files')
            e.text = str(self.numberOfOutputFiles())
            self.xmlroot.append(e)
            logFName = os.path.join(self.workDirectory,"log.txt")
            if os.path.exists(logFName):
                with open(logFName) as f:
                    logLines = f.readlines()
                try:
                    mrBumpDir = None
                    for logL in logLines:
                        if logL.startswith("MRBUMP_DIRECTORY"):
                            mrBumpDir = logL.split()[1].strip()
                    if mrBumpDir is not None:
                        edir = etree.Element('MGMRBUMPDIR')
                        edir.text = str(mrBumpDir)
                        self.xmlroot.append(edir)
                        logs = glob.glob(os.path.normpath(os.path.join(mrBumpDir,'logs','*.log')))
                        logs.extend(glob.glob(os.path.normpath(os.path.join(mrBumpDir,'logs','*.txt'))))
                        # Sort the logs by creation date. I guess XML keeps order? Is this guaranteed?
                        logs = ((os.stat(path), path) for path in logs)
                        logs = ((stat[ST_CTIME], path) for stat, path in logs if S_ISREG(stat[ST_MODE]))
                        
                        for cdate, log in sorted(logs):
                            elog = etree.Element(os.path.basename(log))
                            elog.text = str(log)
                            self.xmlroot.append(elog)
                            shutil.copyfile(log, os.path.join(self.workDirectory,"mrbump_"+os.path.basename(log)))
                    if mrBumpDir is not None:
                        try:
                            win = PhmmerReportNoGui()
                            win.setResultsDir(mrBumpDir)
                            svg = win.svg(500,short=True)
                            print(svg)
                            self.xmlroot.append(etree.fromstring(svg))
                        except:
                            exc_type, exc_value,exc_tb = sys.exc_info()[:3]
                            sys.stderr.write(str(exc_type)+'\n')
                            sys.stderr.write(str(exc_value)+'\n')
                except:
                    print("No MrBump directory reported in log file.")
            
        except:
            exc_type, exc_value,exc_tb = sys.exc_info()[:3]
            sys.stderr.write(str(exc_type)+'\n')
            sys.stderr.write(str(exc_value)+'\n')

            self.appendErrorReport(202,'Data harvesting failed')
            
        CCP4Utils.saveEtreeToFile(self.xmlroot,self.makeFileName('PROGRAMXML'))
        """
        if ( len(outList) ) > 0:
          return CPluginScript.SUCCEEDED
        else:
          return CPluginScript.MARK_TO_DELETE
        """
        return CPluginScript.SUCCEEDED

    def addReportWarning(self, text):
        warningsNode = None
        warningsNodes = self.xmlroot.xpath('//Warnings')
        if len(warningsNodes) == 0: warningsNode = etree.SubElement(self.xmlroot, 'Warnings')
        else: warningsNode = warningsNodes[0]
        warningNode = etree.SubElement(warningsNode,'Warning')
        warningNode.text = text

