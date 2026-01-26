import os
import re
import sys

from lxml import etree

from ccp4i2.core import CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript


class ccp4mg_edit_nomrbump(CPluginScript):
    
    TASKMODULE = 'bioinformatics'            # Where this plugin will appear on the gui
    TASKTITLE = 'Interactive selection of MR model components - CCP4mg'     # A short title for gui menu
    TASKNAME = 'ccp4mg_edit_nomrbump'                  # Task name - should be same as class name
    TASKCOMMAND = 'ccp4mg'                          # The command to run the executable
    ASYNCHRONOUS = True
    MAINTAINER = 'stuart.mcnicholas@york.ac.uk'

    ERROR_CODES = {  200 : { 'description' : 'CCP4MG exited with error status' }, 201 : { 'description' : 'Failed in harvest operation' },202 : { 'description' : 'Failed in processOutputFiles' }}

    def makeCommandAndScript(self):
        from ccp4i2.core import CCP4Utils
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

        if self.container.inputData.XYZIN_LIST.isSet():
            if len(self.container.inputData.XYZIN_LIST)>0:
               try:
                   View = etree.Element('View')
                   scale_auto = etree.Element('scale_auto')
                   scale_auto.text = "true"
                   View.append(scale_auto)
                   centre_molData = etree.Element('centre_MolData')
                   XYZIN = self.container.inputData.XYZIN_LIST[0]
                   centre_molData.text = os.path.splitext(os.path.basename(XYZIN.__str__()))[0]
                   View.append(centre_molData)
                   orientation_auto = etree.Element('orientation_auto')
                   molData = etree.Element('molData')
                   molData.text = os.path.splitext(os.path.basename(XYZIN.__str__()))[0]
                   orientation_auto.append(molData)
                   selection = etree.Element('selection')
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
                        molData = etree.Element('MolData')
                        name = etree.Element('name')
                        name.text = os.path.splitext(os.path.basename(XYZIN.__str__()))[0]
                        molData.append(name)
                        filename = etree.Element('filename')
                        filetype = etree.Element('filetype')
                        filetype.text = "FULLPATH"
                        filename.append(filetype)
                        shortPath = etree.Element('shortPath')
                        shortPath.text = os.path.basename(XYZIN.__str__())
                        filename.append(shortPath)
                        fullPath = etree.Element('fullPath')
                        fullPath.text = XYZIN.__str__()
                        filename.append(fullPath)
                        molData.append(filename)
                        molDisp = etree.Element("MolDisp")
                        style = etree.Element("style")
                        style.text = "BONDS"
                        molDisp.append(style)
                        molData.append(molDisp)
                        tree.append(molData)
                    else:
                        print('ccp4mg_general.makeCommandAndScript XYZIN does not exist:',XYZIN.__str__())
                    iFile += 1
            except:
                #an issue with the existence of files
                pass

        status_xml += etree.tostring(tree,encoding='utf-8', pretty_print=True).decode("utf-8")

        print("Writing",self.mgStatusPath)
        print(status_xml)

        with open(self.mgStatusPath, 'wb') as xmlFile:
             xmlFile.write(bytes(status_xml,"utf-8"))

        
        clArgs = [ '-norestore',self.mgStatusPath ]

        origInterface_fname = os.path.join(CCP4Utils.getCCP4I2Dir(),"wrappers","ccp4mg_edit_nomrbump","script","ccp4i2CCP4MGInterfaceNoMRBUMP.py")
        print(origInterface_fname)
        origInterface = open(origInterface_fname)
        origInterface_s = origInterface.read()
        origInterface.close()

        interface_s = origInterface_s
        interface_fname = os.path.join(self.workDirectory,"ccp4i2CCP4MGInterfaceNoMRBUMP.py")
        interface = open(interface_fname,"w+")
        interface.write(interface_s)
        interface.close()
        
        clArgs.extend(['-scr',interface_fname])

        self.appendCommandLine(clArgs)

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

    def processOutputFiles(self):
        try:
            # First up import PDB files that have been output
            
            import glob
            import os
            import shutil
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
                xyzoutList[-1].annotation = "CCP4MG output file number"+str(iFile)
                xyzoutList[-1].subType = 1

            # Create a trivial xml output file
            from lxml import etree
            self.xmlroot = etree.Element('ccp4mg_edit_nomrbump')
            e = etree.Element('number_output_files')
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
        from lxml import etree
        warningsNode = None
        warningsNodes = self.xmlroot.xpath('//Warnings')
        if len(warningsNodes) == 0: warningsNode = etree.SubElement(self.xmlroot, 'Warnings')
        else: warningsNode = warningsNodes[0]
        warningNode = etree.SubElement(warningsNode,'Warning')
        warningNode.text = text

