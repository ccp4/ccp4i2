import ctypes
import os
import sys

from lxml import etree
import chapi
import gemmi

from core.CCP4PluginScript import CPluginScript
from core import CCP4File

class coot_fit_residues(CPluginScript):
    TASKTITLE = 'Coot fit residues'     # A short title for gui menu
    TASKNAME = 'coot_fit_residues'                                  # Task name - should be same as class name
    TASKVERSION= 0.0                                     # Version of this plugin
    WHATNEXT = ['prosmart_refmac']
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 9999999.9

    def process(self, command=None, handler=None, **kw):
        out_path = os.path.normpath(str(self.container.outputData.XYZOUT))
        if sys.platform.startswith("win"):
            out_path = out_path.replace("\\","\\\\")
        mc = chapi.molecules_container_py(True)
        mc.set_use_gemmi(True)
        imol = mc.read_pdb(str(self.container.inputData.XYZIN.fullPath))
        gemmiStructure = gemmi.read_structure(str(self.container.inputData.XYZIN.fullPath))
        imap = mc.read_mtz(str(self.container.inputData.FPHIIN.fullPath),"F","PHI","",False,False)
        mc.set_imol_refinement_map(imap)
        mc.fill_partial_residues(imol)
        for model in gemmiStructure:
            for chain in model:
                firstResidue, lastResidue = chain[0].seqid.num,chain[-1].seqid.num
                mc.refine_residue_range(imol,chain.name,firstResidue, lastResidue,4000)
        mc.write_coordinates(imol,out_path)
        libc = ctypes.CDLL(None)
        if sys.platform == "darwin":
            c_stdout = ctypes.c_void_p.in_dll(libc, '__stdoutp')
            libc.fflush(c_stdout)
        elif sys.platform.startwith("linux"):
            c_stdout = ctypes.c_void_p.in_dll(libc, 'stdout')
            libc.fflush(c_stdout)
        # Create a trivial xml output file
        root = etree.Element("coot_fit_residues")
        self.container.outputData.XYZOUT.subType = 1
        xml_file = CCP4File.CXmlDataFile(fullPath=self.makeFileName("PROGRAMXML"))
        xml_file.saveFile(root)
        self.reportStatus(CPluginScript.SUCCEEDED)
        return CPluginScript.SUCCEEDED

    def handleLogChanged(self, logFilename, inHandleFinish=None):
        if inHandleFinish is None:
            inHandleFinish=False
        # Create a trivial xml output file
        
        self.tableelement.clear()
        tableelement=self.tableelement
        
        pairs = [('Col_0','N'),('Col_1','StartBonds'),('Col_2','FinalBonds')]
        for pair in pairs:
            headerElement = etree.SubElement(tableelement,'Header',label=pair[1],identifier=pair[0])
        
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
                    rowElement = etree.SubElement(tableelement,"row")
                    rowNoElement = etree.SubElement(rowElement,'Col_0')
                    rowNoElement.text = str(iRow)
                    startBondsElement = etree.SubElement(rowElement,'Col_1')
                    startBondsElement.text = str(currentChangeList[0])
                    finalBondsElement = etree.SubElement(rowElement,'Col_2')
                    finalBondsElement.text = str(currentChangeList[1])
                    currentChangeList = []
                    iRow += 1
        
        graphElement = etree.SubElement(tableelement,"Graph", title = "By residue bonds")
        graphColumnElement = etree.SubElement(graphElement,"Column", label='N', positionInList=str(0))
        graphColumnElement = etree.SubElement(graphElement,"Column", label='StartBonds', positionInList=str(1))
        graphColumnElement = etree.SubElement(graphElement,"Column", label='FinalBonds', positionInList=str(2))
        
        if iRow%20 == 0 or inHandleFinish:
            from core import CCP4File
            f = CCP4File.CXmlDataFile(fullPath=self.makeFileName('PROGRAMXML'))
            newXml = etree.tostring(self.xmlroot,pretty_print=True)
            
            if len(newXml) > self.xmlLength:
                # Save the xml if it has grown
                f.saveFile(self.xmlroot)
                self.xmlLength = len(newXml)
    
    def processOutputFiles(self):
        from core.CCP4PluginScript import CPluginScript
        import os
        status = CPluginScript.FAILED
        if os.path.exists(self.container.outputData.XYZOUT.__str__()): status = CPluginScript.SUCCEEDED
        
        self.container.outputData.XYZOUT.subType = 1
        
        # Create a trivial xml output file
        self.handleLogChanged(self.makeFileName('LOG'), inHandleFinish=True)
        self.reportStatus(finishStatus=status)

#====================================================================================================
# PLUGIN TESTS
# See Python documentation on unittest module

import unittest

class test_coot_fit_residues(unittest.TestCase):
    
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
        wrapper = coot_fit_residues(parent=QTAPPLICATION(),name='coot_fit_residues_test1')
        wrapper.container.loadDataFromXml()


def TESTSUITE():
    suite = unittest.TestLoader().loadTestsFromTestCase(test_coot_fit_residues)
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)

