from __future__ import print_function

#from lxml import etree
from xml.etree import ElementTree as ET

from core.CCP4PluginScript import CPluginScript
from PySide2 import QtCore
from core import CCP4Utils

class coot_script_lines(CPluginScript):
    
    TASKMODULE = 'model_building'                               # Where this plugin will appear on the gui
    TASKTITLE = 'Execute arbitrary script code within Coot'     # A short title for gui menu
    TASKNAME = 'coot_script_lines'                                  # Task name - should be same as class name
    TASKCOMMAND = 'coot'                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    WHATNEXT = ['prosmart_refmac']
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 9999999.9
    
    '''
        def __init__(self,parent=None,name=None,workDirectory=''):
        CPluginScript. __init__(self,parent=parent,name=name)
        '''
    
    def makeCommandAndScript(self):
        import os
        
        self.dropDir = os.path.join(self.workDirectory,'COOT_FILE_DROP')
        if not os.path.exists(self.dropDir):
            try:
                os.mkdir(self.dropDir)
            except:
                self.dropDir = self.workDirectory
    
        cootScriptPath = os.path.join(self.workDirectory,'script.py')
        
        self.appendCommandLine(['--no-graphics', '--script',cootScriptPath])

        cootScript = open(cootScriptPath,"w")
        cootScript.write('import coot\n')
        cootScript.write('import os\n')
        
        i = 1
        for XYZIN in self.container.inputData.XYZIN:
          if XYZIN.exists():
            cootScript.write ("MolHandle_"+str(i)+"=coot.read_pdb(r'"+str(XYZIN.fullPath)+"')\n\n")
            i += 1
          else: pass
            #print '\n\n ** Non-file :[' + str(i)+ ']'+str(XYZIN.fullPath)
                
        i = 1
        for FPHIIN in self.container.inputData.FPHIIN:
          if FPHIIN.exists():
            cootScript.write ("MapHandle_"+str(i)+"=coot.make_and_draw_map(r'" + str(FPHIIN.fullPath)+"', 'F', 'PHI', 'PHI', 0, 0)\n\n")
            i += 1
          else: pass
            #print '\n\n ** Non-file :[' +str(i)+ ']'+ str(FPHIIN.fullPath)
            
        i = 1
        for DELFPHIIN in self.container.inputData.DELFPHIIN:
          if DELFPHIIN.exists():
            cootScript.write ("DifmapHandle_"+str(i)+"=coot.make_and_draw_map(r'" + str(DELFPHIIN.fullPath)+"', 'F', 'PHI', 'PHI', 0, 1)\n\n")
            i += 1
          else: pass
            #print '\n\n ** Non-file :[' +str(i)+ ']'+str(DELFPHIIN.fullPath)
          
        if self.container.inputData.DICT.isSet():
            cootScript.write(f"coot.read_cif_dictionary('{self.container.inputData.DICT.fullPath.__str__()}')\n")

        cootScript.write ("dropDir='"+self.dropDir+"'\n\n")

        if self.container.controlParameters.SCRIPT.isSet():
            #Here try to except and exit if script is faulty
            scriptLines = self.container.controlParameters.SCRIPT.__str__().split('\n')
            cootScript.write ('try:\n')
            for scriptLine in scriptLines:
                cootScript.write(f'    {scriptLine}\n')
            cootScript.write('except Exception as err:\n')
            cootScript.write('    import traceback\n')
            cootScript.write('    traceback.print_exc()\n')
            cootScript.write('    coot.coot_real_exit(0)\n')
            cootScript.write ('\n')
          
        cootScript.write("coot.coot_real_exit(0)\n")
        cootScript.close()
        
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        print('#coot_script_lines.processOutputFiles')
        #First up check for exit status of the program
        from core.CCP4Modules import PROCESSMANAGER
        exitStatus = 0
        try:
            exitStatus = PROCESSMANAGER().getJobData(pid=self.getProcessId(), attribute='exitStatus')
            if exitStatus != CPluginScript.SUCCEEDED:
                return exitStatus
        except:
            print('Unable to recover exitStatus')
            return CPluginScript.FAILED


        iPDBOut = 0
        try:
            import os, glob, shutil
            outList = glob.glob(os.path.join(self.dropDir,'*.pdb'))
            xyzoutList = self.container.outputData.XYZOUT
            print(outList)
            for outputPDB in outList:
                fpath,fname = os.path.split(outputPDB)
                xyzoutList.append(xyzoutList.makeItem())
                outputFilePath = os.path.join(self.workDirectory,'XYZOUT_'+str(iPDBOut)+'-coordinates.pdb')
                shutil.copyfile(outputPDB, outputFilePath)
                xyzoutList[-1].setFullPath(outputFilePath)
                xyzoutList[-1].annotation=fname
                iPDBOut += 1
        except:
            return CPluginScript.FAILED

        try:
            import os, glob, shutil
            outList = glob.glob(os.path.join(self.dropDir,'*.cif'))
            xyzoutList = self.container.outputData.XYZOUT
            print(outList)
            for outputPDB in outList:
                fpath,fname = os.path.split(outputPDB)
                xyzoutList.append(xyzoutList.makeItem())
                outputFilePath = os.path.join(self.workDirectory,'XYZOUT_'+str(iPDBOut)+'-coordinates.cif')
                shutil.copyfile(outputPDB, outputFilePath)
                xyzoutList[-1].setFullPath(outputFilePath)
                xyzoutList[-1].annotation=fname
                xyzoutList[-1].contentFlag=2
                iPDBOut += 1
        except:
            return CPluginScript.FAILED

        iCIFOut = 0
        try:
            outList = glob.glob(os.path.join(self.dropDir,'coot-ccp4', 'prodrg-out*.cif'))
            dictoutList = self.container.outputData.DICTOUT
            for outputCIF in outList:
                fpath,fname = os.path.split(outputCIF)
                dictoutList.append(dictoutList.makeItem())
                outputFilePath = os.path.join(self.workDirectory,'COOT_'+str(iCIFOut)+'-prodrg.cif')
                shutil.copyfile(outputCIF, outputFilePath)
                dictoutList[-1].setFullPath(outputFilePath)
                dictoutList[-1].annotation=fname
                iCIFOut += 1
        except:
            return CPluginScript.FAILED

        # Create a trivial xml output file
        from core import CCP4File
        self.xmlroot = ET.Element('coot_script_lines')
        e = ET.Element('number_output_pdbs')
        e.text = str(iPDBOut)
        self.xmlroot.append(e)
        e = ET.Element('number_output_cifs')
        e.text = str(iCIFOut)
        self.xmlroot.append(e)
        tableElement = self.logToXML()
        if tableElement is not None: self.xmlroot.append(tableElement)
        
        try:
            aFile=open( self.makeFileName('PROGRAMXML'),'w')
            ET.indent(self.xmlroot)
            CCP4Utils.writeXML(aFile,ET.tostring(self.xmlroot))
            aFile.close()
        except:
            print('Oops')
            return CPluginScript.FAILED

        return CPluginScript.SUCCEEDED

    def logToXML(self):
        
        tableelement = ET.Element('Table', title='Per residue statistics')
        
        pairs = [('Col_0','N'),('Col_1','StartBonds'),('Col_2','FinalBonds')]
        for pair in pairs:
            headerElement = ET.SubElement(tableelement,'Header',label=pair[1],identifier=pair[0])
        
        cootlines = open(self.makeFileName('LOG')).readlines()
        iRow = 1
        currentChangeList = []
        for line in cootlines:
            if line.startswith('bonds:'):
                #Possibilities are that this is the first or second bonds reported, i.e. this is an initial or final value
                if len(currentChangeList) == 0:
                    currentChangeList.append(line.split()[1])
                else:
                    currentChangeList.append(line.split()[1])
                    rowElement = ET.SubElement(tableelement,"row")
                    rowNoElement = ET.SubElement(rowElement,'Col_0')
                    rowNoElement.text = str(iRow)
                    startBondsElement = ET.SubElement(rowElement,'Col_1')
                    startBondsElement.text = str(currentChangeList[0])
                    finalBondsElement = ET.SubElement(rowElement,'Col_2')
                    finalBondsElement.text = str(currentChangeList[1])
                    currentChangeList = []
                    iRow += 1
    
        if iRow == 1: return None
    
        graphElement = ET.SubElement(tableelement,"Graph", title = "By residue bonds")
        graphColumnElement = ET.SubElement(graphElement,"Column", label='N', positionInList=str(0))
        graphColumnElement = ET.SubElement(graphElement,"Column", label='StartBonds', positionInList=str(1))
        graphColumnElement = ET.SubElement(graphElement,"Column", label='FinalBonds', positionInList=str(2))
        
        return tableelement

#====================================================================================================
# PLUGIN TESTS
# See Python documentation on unittest module

import unittest

class testcoot_script_lines(unittest.TestCase):
    
    def setUp(self):
        # make all background jobs wait for completion
        # this is essential for unittest to work
        from core.CCP4Modules import QTAPPLICATION,PROCESSMANAGER
        self.app = QTAPPLICATION()
        PROCESSMANAGER().setWaitForFinished(10000)
    
    def tearDown(self):
        from core.CCP4Modules import PROCESSMANAGER
        PROCESSMANAGER().setWaitForFinished(-1)
    
    def test_1(self):
        from core.CCP4Modules import QTAPPLICATION
        wrapper = coot_script_lines(parent=QTAPPLICATION(),name='coot_script_lines_test1')
        wrapper.container.loadDataFromXml()


def TESTSUITE():
    suite = unittest.TestLoader().loadTestsFromTestCase(testcoot_script_lines)
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)
