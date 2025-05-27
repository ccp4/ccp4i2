import os
import sys

from lxml import etree

from ......core.CCP4PluginScript import CPluginScript


class hklin2cif(CPluginScript):
    DESCRIPTION = 'Convert hklin file (+/- scaled unmerged data) creates as aprt of adding_stats_to_mmcif pipeline'
    TASKNAME = 'hklin2cif'                                  # Task name - should be same as class name
    TASKVERSION= 0.0                                     # Version of this plugin
    ASYNCHRONOUS = False
    TIMEOUT_PERIOD =3.
    TASKCOMMAND='gemmi'

    ERROR_CODES = {  200 : { 'description' : 'Failed to add item to mol list' },201 : { 'description' : 'Failed to setFullPath' },}
    
    def __init__(self,*args,**kws):
        innate=CPluginScript.__init__(self, *args, **kws)
        self.xmlroot = ET.Element('Hklin2cif')

    def makeCommandAndScript(self):
        self.outputCifPath = os.path.normpath(os.path.join(self.getWorkDirectory(),'Reflections.cif'))
        self.appendCommandLine(['mtz2cif']) 
        if self.container.inputData.SCALEDUNMERGED.isSet():
            self.appendCommandLine(['--depo'])  
        self.appendCommandLine(os.path.join(self.workDirectory, 'GEMMIFIED.mtz'))
        if self.container.inputData.SCALEDUNMERGED.isSet():
            self.appendCommandLine(str(self.container.inputData.SCALEDUNMERGED))
        self.appendCommandLine([self.outputCifPath])
        
        return CPluginScript.SUCCEEDED

    def processInputFiles(self):
        self.container.inputData.HKLIN.loadFile()
        columns = self.container.inputData.HKLIN.fileContent.listOfColumns
        columnLabels = [column.columnLabel.__str__() for column in columns]
        COLUMNS2GEMMI = {
            'FREERFLAG_FREER':'FREER', 
            'F_SIGF_Iplus':' I(+)', 
            'F_SIGF_SIGIplus':'SIGI(+)', 
            'F_SIGF_Iminus':'I(-)', 
            'F_SIGF_SIGIminus':'SIGI(-)', 
            'F_SIGF_Fplus':' F(+)', 
            'F_SIGF_SIGFplus':'SIGF(+)', 
            'F_SIGF_Fminus':'F(-)', 
            'F_SIGF_SIGFminus':'SIGF(-)', 
            'F_SIGF_Fmean':'F', 
            'F_SIGF_SIGFmean':'SIGF', 
            'F_SIGF_F':'F', 
            'F_SIGF_SIGF':'SIGF', 
            'F_SIGF_I':'I', 
            'F_SIGF_SIGI':'SIGI', 
            'F_SIGF_Imean':'I', 
            'F_SIGF_SIGImean':'SIGI', 
            'FPHIOUT_F':'FWT', 
            'FPHIOUT_PHI':'PHWT', 
            'DIFFPHIOUT_F':'DELFWT', 
            'DIFFPHIOUT_PHI':'PHDELWT'
        }
        inLineArgs = ['LABIN', 'FILE', '1']
        outLineArgs = ['LABOUT', 'FILE', '1']
        I=1
        for columnLabel in columnLabels:
            if columnLabel in COLUMNS2GEMMI:
                inLineArgs.append(f'E{I}={columnLabel}')
                outLineArgs.append(f'E{I}={COLUMNS2GEMMI[columnLabel]}')
                I+=1

        self.container.inputData.HKLIN.runCad(
            os.path.join(self.workDirectory, 'GEMMIFIED.mtz'), comLines = [' '.join(inLineArgs), ' '.join(outLineArgs)])
        sys.stdout.flush()
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        self.container.outputData.CIFFILE.setFullPath(self.outputCifPath)
        self.container.outputData.CIFFILE.annotation = 'Cif file of observations used in refinement'
        return CPluginScript.SUCCEEDED
