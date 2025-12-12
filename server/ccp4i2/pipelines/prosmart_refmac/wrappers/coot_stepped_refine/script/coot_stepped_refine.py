from __future__ import print_function

from lxml import etree

from ccp4i2.core.CCP4PluginScript import CPluginScript

class coot_stepped_refine(CPluginScript):
    
    #TASKMODULE = 'model_building'                               # Where this plugin will appear on the gui
    TASKTITLE = 'Coot stepped refinement'     # A short title for gui menu
    TASKNAME = 'coot_stepped_refine'                                  # Task name - should be same as class name
    TASKCOMMAND = 'coot'                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    WHATNEXT = ['prosmart_refmac']
    TIMEOUT_PERIOD = 9999999.9
    

    def processInputFiles(self):
      #Create root for Output XML
      self.xmlroot = etree.Element('Coot_stepped_refine')
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
        refineCommand = "stepped_refine_protein(0)\n"
        if self.container.controlParameters.USERAMA.isSet():
            if self.container.controlParameters.USERAMA:
                refineCommand = "stepped_refine_protein_for_rama(0)\n"
        cootScript.write(refineCommand)
        cootScript.write("write_pdb_file(0,r'"+str(self.container.outputData.XYZOUT.fullPath)+"')\n")
        cootScript.write("coot_real_exit(0)\n")
        cootScript.close()
        
        return CPluginScript.SUCCEEDED
    
    
    def handleLogChanged(self, logFilename, inHandleFinish=None):
        if inHandleFinish is None:
            inHandleFinish=False
        # Create a trivial xml output file
        
        print('coot_stepped_refine.handleLogChanged',inHandleFinish)

        if getattr(self,'xmlroot',None) is None:
          self.xmlroot = etree.Element('Coot_stepped_refine')
          self.tableelement = etree.SubElement(self.xmlroot, 'Table', title='Per residue statistics')
          self.xmlLength = 0
        
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
