from lxml import etree

from ccp4i2.core.CCP4PluginScript import CPluginScript

class coot_fit_residues(CPluginScript):
    
    TASKTITLE = 'Coot fit residues'     # A short title for gui menu
    TASKNAME = 'coot_fit_residues'                                  # Task name - should be same as class name
    TASKCOMMAND = 'coot'                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    WHATNEXT = ['prosmart_refmac']
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 9999999.9
    
    '''
        def __init__(self,parent=None,name=None,workDirectory=''):
        CPluginScript. __init__(self,parent=parent,name=name)
        '''

    def processInputFiles(self):
        #Create root for Output XML
        self.xmlroot = etree.Element('Coot_fit_residues')
        self.tableelement = etree.SubElement(self.xmlroot, 'Table', title='Per residue statistics')
        self.xmlLength = 0
        # watch the log file
        logFilename = self.makeFileName('LOG')
        from ccp4i2.core import CCP4Utils
        CCP4Utils.saveFile(logFilename,'')
        self.watchFile(logFilename,self.handleLogChanged)
    
    def makeCommandAndScript(self):
        import os
        cootScriptPath = os.path.join(self.workDirectory,'script.py')
        self.appendCommandLine(['--no-state-script','--no-graphics','--python','--pdb',self.container.inputData.XYZIN.fullPath,'--script',cootScriptPath])
        
        cootScript = open(cootScriptPath,"w")
        cootScript.write("make_and_draw_map(r'" + str(self.container.inputData.FPHIIN.fullPath)+"', 'F', 'PHI', 'PHI', 0, 0)\n")
        refineCommand = "fit_protein(0)\n"
        cootScript.write(refineCommand)
        cootScript.write("write_pdb_file(0,r'"+str(self.container.outputData.XYZOUT.fullPath)+"')\n")
        cootScript.write("coot_real_exit(0)\n")
        cootScript.close()
        
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
            from ccp4i2.core import CCP4File
            f = CCP4File.CXmlDataFile(fullPath=self.makeFileName('PROGRAMXML'))
            newXml = etree.tostring(self.xmlroot,pretty_print=True)
            
            if len(newXml) > self.xmlLength:
                # Save the xml if it has grown
                f.saveFile(self.xmlroot)
                self.xmlLength = len(newXml)
    
    def processOutputFiles(self):
        from ccp4i2.core.CCP4PluginScript import CPluginScript
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
        from ccp4i2.core.CCP4Modules import QTAPPLICATION,PROCESSMANAGER
        self.app = QTAPPLICATION()
        PROCESSMANAGER().setWaitForFinished(10000)
    
    def tearDown(self):
        from ccp4i2.core.CCP4Modules import PROCESSMANAGER
        PROCESSMANAGER().setWaitForFinished(-1)
    
    def test_1(self):
        from ccp4i2.core.CCP4Modules import QTAPPLICATION
        wrapper = coot_fit_residues(parent=QTAPPLICATION(),name='coot_fit_residues_test1')
        wrapper.container.loadDataFromXml()


def TESTSUITE():
    suite = unittest.TestLoader().loadTestsFromTestCase(test_coot_fit_residues)
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)

