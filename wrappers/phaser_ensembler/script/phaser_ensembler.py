
from core.CCP4PluginScript import CPluginScript
from core import CCP4Utils
import base64

class phaser_ensembler(CPluginScript):
    TASKNAME = 'phaser_ensembler'                                  # Task name - should be same as class name
    TASKCOMMAND = 'phaser.ensembler'                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    WHATNEXT = ['prosmart_refmac']
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 9999999.9


    def makeCommandAndScript(self):
        for iCoordSet, xyzin in enumerate(self.container.inputData.XYZIN_LIST):
            if not xyzin.isSelectionSet():
                self.appendCommandLine(['input.model='+xyzin.fullPath.__str__()])
            else:
                import os
                xyzin.selection.text='{'+xyzin.selection.__str__()+'} and {(ALA,CYS,ASP,GLU,PHE,GLY,HIS,ILE,LYS,LEU,MET,ASN,PRO,GLN,ARG,SER,THR,VAL,TRP,TYR)}'
                inputCoordPath = os.path.normpath(os.path.join(self.getWorkDirectory(),'selected_'+iCoordSet.__str__()+'.pdb'))
                xyzin.getSelectedAtomsPdbFile(inputCoordPath)
                self.appendCommandLine(['input.model='+inputCoordPath])
    
        if self.container.inputData.ALIGNIN.isSet():
            self.appendCommandLine(['input.alignment='+self.container.inputData.ALIGNIN.fullPath.__str__()])

        self.appendCommandLine(['output.location='+self.getWorkDirectory()])

        ignoreLevels = ['input']
        ignoreParameters = ['output__gui_output_dir','output__location','output__job_title','output__root','configuration__superposition__atoms']
        
        for attr in self.container.keywords.dataOrder():
            if attr not in ignoreParameters:
                parameter = self.container.keywords.get(attr)
                if parameter.isSet():
                    levels = attr.split('__')
                    if levels[0] not in ignoreLevels:
                        commandLineKey = ".".join(levels)
                        keyValuePair = commandLineKey+"="+parameter.__str__()
                        self.appendCommandLine([keyValuePair])
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        import os
        fileIfMadePath = os.path.join(self.getWorkDirectory(),'ensemble_merged.pdb')
        if os.path.isfile(fileIfMadePath):
            remarkedFilePath = os.path.join(self.getWorkDirectory(),'ensemble_remarked.pdb')
            if self.container.inputData.OVERRIDEID.isSet():
                with open(remarkedFilePath,'w') as remarkedFile:
                    from core.CCP4ModelData import CPdbData
                    ensembledUnremarked = CPdbData()
                    ensembledUnremarked.loadFile(fileIfMadePath)
                    mmdbManager = ensembledUnremarked.mmdbManager
                    for i in range(mmdbManager.GetNumberOfModels()):
                        remarkedFile.write('REMARK PHASER ENSEMBLE MODEL{0:2d} ID {1:4.1f}\n'.format(i+1, float(self.container.inputData.OVERRIDEID)))
                    remarkedFile.write(open(fileIfMadePath).read())
        
                self.container.outputData.XYZOUT.setFullPath(remarkedFilePath)
            else:
                self.container.outputData.XYZOUT.setFullPath(fileIfMadePath)
            
            
            self.container.outputData.XYZOUT.annotation = 'Merged ensemble'
            
            from lxml import etree
            logText = open(self.makeFileName('LOG'),"r").read()
            rootNode = etree.Element('PHASER_ENSEMBLER')
            logNode = etree.SubElement(rootNode,'LOGTEXT')
            #logNode.text = etree.CDATA(logText)
            logNode.text = base64.b64encode(logText)
            with open (self.makeFileName('PROGRAMXML'),'w') as outputFile:
                CCP4Utils.writeXML(outputFile,etree.tostring(rootNode))
            return CPluginScript.SUCCEEDED
        else:
            return CPluginScript.FAILED

