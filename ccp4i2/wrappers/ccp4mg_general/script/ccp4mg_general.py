import glob
import os
import re
import shutil
import sys
import time
import xml.etree.ElementTree as ET

from PySide2 import QtCore

from ....core import CCP4Utils
from ....core.CCP4PluginScript import CPluginScript


class ccp4mg_general(CPluginScript):
    
    TASKMODULE = 'model_building'            # Where this plugin will appear on the gui
    TASKTITLE = 'Molecular graphics visualization and figure creation - CCP4MG'     # A short title for gui menu
    TASKNAME = 'ccp4mg_general'                  # Task name - should be same as class name
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
        tree = ET.Element("CCP4MG_Status",nsmap = NSMAP,attrib={location_attribute: 'http://www.ysbl.york.ac.uk/~mcnicholas/schema/CCP4MGApplicationOutput.xsd'})

        # FIXME - We should allow list of dicts.
        DICT = None
        if self.container.inputData.DICT.isSet():
            DICT = self.container.inputData.DICT
            
        if self.container.inputData.XYZIN_LIST.isSet():
            if len(self.container.inputData.XYZIN_LIST)>0:
               try:
                   View = ET.Element('View')
                   scale_auto = ET.Element('scale_auto')
                   scale_auto.text = "true"
                   View.append(scale_auto)
                   centre_molData = ET.Element('centre_MolData')
                   XYZIN = self.container.inputData.XYZIN_LIST[0]
                   centre_molData.text = os.path.splitext(os.path.basename(XYZIN.__str__()))[0]
                   View.append(centre_molData)
                   orientation_auto = ET.Element('orientation_auto')
                   molData = ET.Element('molData')
                   molData.text = os.path.splitext(os.path.basename(XYZIN.__str__()))[0]
                   orientation_auto.append(molData)
                   selection = ET.Element('selection')
                   selection.text = "all"
                   orientation_auto.append(selection)
                   View.append(orientation_auto)
                   tree.append(View)
               except:
                   #an issue with the existence of files
                   pass
            try:
                iFile = 1
                for XYZIN in self.container.inputData.XYZIN_LIST:
                    if os.path.isfile(XYZIN.__str__()):
                        molData = ET.Element('MolData')
                        name = ET.Element('name')
                        name.text = os.path.splitext(os.path.basename(XYZIN.__str__()))[0]
                        molData.append(name)
                        filename = ET.Element('filename')
                        filetype = ET.Element('filetype')
                        filetype.text = "FULLPATH"
                        filename.append(filetype)
                        shortPath = ET.Element('shortPath')
                        shortPath.text = os.path.basename(XYZIN.__str__())
                        filename.append(shortPath)
                        fullPath = ET.Element('fullPath')
                        fullPath.text = XYZIN.__str__()
                        filename.append(fullPath)
                        molData.append(filename)
                        if DICT is not None:
                            customResCIFFiles = ET.Element('customResCIFFiles')
                            cifmonomer = ET.Element('cifmonomer')
                            cifname = ET.Element('name')
                            ciffilename = ET.Element('filename')
                            dict_lines = []
                            dict_name = "DRG"
                            with open(DICT.__str__(), 'rb') as f:
                               dict_lines = f.readlines()
                            for dict_l in dict_lines:
                                if dict_l.startswith("data_comp_") and dict_l.strip() != "data_comp_list":
                                    dict_name = dict_l[len("data_comp_"):].strip()
                                    print("Set dictionary name to ",dict_name)
                                    break
                            cifname.text = dict_name
                            ciffilename.text = DICT.__str__()
                            cifmonomer.append(cifname)
                            cifmonomer.append(ciffilename)
                            customResCIFFiles.append(cifmonomer)
                            molData.append(customResCIFFiles)
                        molDisp = ET.Element("MolDisp")
                        style = ET.Element("style")
                        style.text = "BONDS"
                        molDisp.append(style)
                        molData.append(molDisp)
                        tree.append(molData)
                    else:
                        print('ccp4mg_general.makeCommandAndScript XYZIN does not exist:',XYZIN.__str__())
                    iFile += 1
            except:
                #an issue with the existence of files
                exc_type, exc_value,exc_tb = sys.exc_info()[:3]
                sys.stdout.write(str(exc_type)+'\n')
                sys.stdout.write(str(exc_value)+'\n')
                pass

        if self.container.inputData.FPHIIN_LIST.isSet():
            try:
                iFile = 1
                for FPHIIN in self.container.inputData.FPHIIN_LIST:
                    print(' reading file number ' + str ( iFile )) 
                    if os.path.isfile(FPHIIN.__str__()):
                        mapData = ET.Element('MapData')
                        name = ET.Element('name')
                        name.text = os.path.splitext(os.path.basename(FPHIIN.__str__()))[0]
                        mapData.append(name)
                        filename = ET.Element('filename')
                        filetype = ET.Element('filetype')
                        filetype.text = "FULLPATH"
                        filename.append(filetype)
                        shortPath = ET.Element('shortPath')
                        shortPath.text = os.path.basename(FPHIIN.__str__())
                        filename.append(shortPath)
                        fullPath = ET.Element('fullPath')
                        fullPath.text = FPHIIN.__str__()
                        filename.append(fullPath)
                        mapData.append(filename)
                        f = ET.Element('f')
                        f.text = "F"
                        mapData.append(f)
                        phi = ET.Element('phi')
                        phi.text = "PHI"
                        mapData.append(phi)
                        mapDisp = ET.Element("MapDisp")
                        style = ET.Element("style")
                        style.text = "surface_style_chickenwire"
                        mapDisp.append(style)
                        surface_style = ET.Element("surface_style")
                        surface_style.text = "surface_style_chickenwire"
                        mapDisp.append(surface_style)
                        colour = ET.Element("colour")
                        colour.text = "blue"
                        mapDisp.append(colour)
                        radius = ET.Element("radius")
                        radius.text = "10.0"
                        mapDisp.append(radius)
                        contour_level = ET.Element("contour_level")
                        contour_level.text = "1.5"
                        mapDisp.append(contour_level)
                        mapData.append(mapDisp)
                        tree.append(mapData)
                    else:
                        print('ccp4mg_general.makeCommandAndScript FPHIIN does not exist:',FPHIIN.__str__())
                    iFile += 1
            except:
                print(' Exception ')
                #an issue with the existence of files
                pass

        if self.container.inputData.DELFPHIIN_LIST.isSet():
            try:
                iFile = 1
                for DELFPHIIN in self.container.inputData.DELFPHIIN_LIST:
                    print(' reading diff file number ' + str ( iFile )) 
                    if os.path.isfile(DELFPHIIN.__str__()):
                        mapData = ET.Element('MapData')
                        name = ET.Element('name')
                        name.text = os.path.splitext(os.path.basename(DELFPHIIN.__str__()))[0]
                        mapData.append(name)
                        filename = ET.Element('filename')
                        filetype = ET.Element('filetype')
                        filetype.text = "FULLPATH"
                        filename.append(filetype)
                        shortPath = ET.Element('shortPath')
                        shortPath.text = os.path.basename(DELFPHIIN.__str__())
                        filename.append(shortPath)
                        fullPath = ET.Element('fullPath')
                        fullPath.text = DELFPHIIN.__str__()
                        filename.append(fullPath)
                        mapData.append(filename)
                        f = ET.Element('f')
                        f.text = "F"
                        mapData.append(f)
                        phi = ET.Element('phi')
                        phi.text = "PHI"
                        mapData.append(phi)
                        difference = ET.Element('difference')
                        difference.text = "1"
                        mapData.append(difference)
                        mapDisp = ET.Element("MapDisp")
                        style = ET.Element("style")
                        style.text = "surface_style_chickenwire"
                        mapDisp.append(style)
                        surface_style = ET.Element("surface_style")
                        surface_style.text = "surface_style_chickenwire"
                        mapDisp.append(surface_style)
                        colour = ET.Element("colour")
                        colour.text = "green"
                        mapDisp.append(colour)
                        radius = ET.Element("radius")
                        radius.text = "10.0"
                        mapDisp.append(radius)
                        contour_level = ET.Element("contour_level")
                        contour_level.text = "3.0"
                        mapDisp.append(contour_level)
                        mapData.append(mapDisp)
                        mapDisp = ET.Element("MapDisp")
                        style = ET.Element("style")
                        style.text = "surface_style_chickenwire"
                        mapDisp.append(style)
                        surface_style = ET.Element("surface_style")
                        surface_style.text = "surface_style_chickenwire"
                        mapDisp.append(surface_style)
                        colour = ET.Element("colour")
                        colour.text = "red"
                        mapDisp.append(colour)
                        radius = ET.Element("radius")
                        radius.text = "10.0"
                        mapDisp.append(radius)
                        contour_level = ET.Element("contour_level")
                        contour_level.text = "-3.0"
                        mapDisp.append(contour_level)
                        mapData.append(mapDisp)
                        tree.append(mapData)
                    else:
                        print('ccp4mg_general.makeCommandAndScript DELFPHIIN does not exist:',DELFPHIIN.__str__())
                    iFile += 1
            except:
                print(' Exception ')
                #an issue with the existence of files
                pass

        ET.indent(tree)
        status_xml += ET.tostring(tree)

        print("Writing",self.mgStatusPath)
        print(status_xml)

        with open(self.mgStatusPath, 'wb') as xmlFile:
             xmlFile.write(bytes(status_xml,"utf-8"))
        
        clArgs = [ '-norestore',self.mgStatusPath ]

        origInterface_fname = os.path.join(CCP4Utils.getCCP4I2Dir(),"wrappers","ccp4mg_general","script","ccp4i2CCP4MGInterface.py")
        origInterface = open(origInterface_fname)
        origInterface_s = origInterface.read()
        origInterface.close()

        interface_s = origInterface_s
        interface_fname = os.path.join(self.workDirectory,"ccp4i2CCP4MGInterface.py")
        interface = open(interface_fname,"w+")
        interface.write(interface_s)
        interface.close()
        
        clArgs.extend(['-scr',interface_fname])

        if self.container.inputData.SEQUENCE_LIST.isSet():
            for sequenceFile in self.container.inputData.SEQUENCE_LIST:
               clArgs.append(sequenceFile)

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
        print('ccp4mg_general',time.time())
        print('ccp4mg_general',glob.glob(os.path.join(self.workDirectory,'*.*')))
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
                xyzoutList[-1].annotation = "CCP4MG coordinate output file number"+str(iFile)
                xyzoutList[-1].subType = 1

            """
            # Comment out this for now, may be used in future ...
            # I believe that these will need processing into mini mtz files somehow.... So this is likely to be very different
            globPath = os.path.normpath(os.path.join(self.dropDir,'mtzout*.mtz'))
            outList = glob.glob(globPath)
            xyzoutList = self.container.outputData.FPHIOUT
            for outputMTZ in outList:
                fpath,fname = os.path.split(outputMTZ)
                iFile = int(fname[6:-4])
                xyzoutList.append(xyzoutList.makeItem())
                outputFilePath = os.path.normpath(os.path.join(self.workDirectory,'MTZOUT_'+str(iFile)+'-density.mtz'))
                shutil.copyfile(outputMTZ, outputFilePath)
                xyzoutList[-1].setFullPath(outputFilePath)
                xyzoutList[-1].annotation = "CCP4MG density output file number"+str(iFile)
                xyzoutList[-1].subType = 1

            globPath = os.path.normpath(os.path.join(self.dropDir,'delmtz*.mtz'))
            outList = glob.glob(globPath)
            xyzoutList = self.container.outputData.DELFPHIOUT
            for outputMTZ in outList:
                fpath,fname = os.path.split(outputMTZ)
                iFile = int(fname[6:-4])
                xyzoutList.append(xyzoutList.makeItem())
                outputFilePath = os.path.normpath(os.path.join(self.workDirectory,'DELMTZ_'+str(iFile)+'-density.mtz'))
                shutil.copyfile(outputMTZ, outputFilePath)
                xyzoutList[-1].setFullPath(outputFilePath)
                xyzoutList[-1].annotation = "CCP4MG difference density output file number"+str(iFile)
                xyzoutList[-1].subType = 1
            """

            # Create a trivial xml output file
            self.xmlroot = ET.Element('ccp4mg_general')
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
        warningsNodes = self.xmlroot.xpath('//Warnings')
        if len(warningsNodes) == 0: warningsNode = ET.SubElement(self.xmlroot, 'Warnings')
        else: warningsNode = warningsNodes[0]
        warningNode = ET.SubElement(warningsNode,'Warning')
        warningNode.text = text

