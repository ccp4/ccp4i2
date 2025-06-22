import os
import shutil
import xml.etree.ElementTree as ET

from ....core import CCP4ErrorHandling
from ....core.CCP4Modules import PROCESSMANAGER
from ....core.CCP4PluginScript import CPluginScript


class lorestr_i2(CPluginScript):
    TASKNAME = 'lorestr_i2'   # Task name - should be same as class name and match pluginTitle in the .def.xml file
    TASKVERSION= 0.1               # Version of this plugin
    MAINTAINER = 'kovoleg'
    ERROR_CODES = { 201 : {'description' : 'Failed to analyse output files' },
                    202 : {'description' : 'Failed applying selection ot PDB file' }
                    }
    PURGESEARCHLIST = [ [ 'hklin.mtz' , 0 ],
                       ['log_mtzjoin.txt', 0]
                       ]
    TASKCOMMAND="lorestr" #_dummy
    PERFORMANCECLASS = 'CRefinementPerformance'
    
    def __init__(self, *args, **kws):
        super().__init__(*args, **kws)
        self._readyReadStandardOutputHandler = self.handleReadyReadStandardOutput
        self.xmlroot = ET.Element('lorestr_i2')
        self.xmlLength = 0

    def processInputFiles(self):
        from ....core import CCP4XtalData
        error = None
        self.hklin = None
        dataObjects = []
        #print '\n\n\n***contentFlag',self.container.inputData.F_SIGF.contentFlag
        #Append Observation with representation dependent on whether we are detwining on Is or not

        obsTypeRoot = 'CONTENT_FLAG_F'
        obsPairOrMean = 'MEAN'
        obsType = getattr(CCP4XtalData.CObsDataFile, obsTypeRoot+obsPairOrMean)
        dataObjects += [['F_SIGF',obsType]] # ,obsType

        #Apply coordinate selection if set
        self.inputCoordPath = os.path.normpath(self.container.inputData.XYZIN.fullPath.__str__())
        if self.container.inputData.XYZIN.isSelectionSet():
            self.inputCoordPath = os.path.normpath(os.path.join(self.getWorkDirectory(),'selected.pdb'))
            self.container.inputData.XYZIN.getSelectedAtomsPdbFile(self.inputCoordPath)

        #Include FreeRflag if called for
        if self.container.inputData.FREERFLAG.isSet():
            dataObjects += ['FREERFLAG']
        self.hklin,error = self.makeHklin(dataObjects)
        if error.maxSeverity()>CCP4ErrorHandling.Severity.WARNING:
            return CPluginScript.FAILED
        else:
            return CPluginScript.SUCCEEDED

    def handleReadyReadStandardOutput(self):
        print('Runtime call handleReadyReadStandardOutput\n\n')
        if not hasattr(self,'logFileHandle'): self.logFileHandle = open(self.makeFileName('LOG'),'w')
        if not hasattr(self,'logFileBuffer'): self.logFileBuffer = ''
        pid = self.getProcessId()
        qprocess = PROCESSMANAGER().getJobData(pid,attribute='qprocess')
        availableStdout = qprocess.readAllStandardOutput()
        self.logFileHandle.write(availableStdout)
        self.logFileHandle.flush()

        self.xmlroot = ET.Element("lorestr_i2")
        self.flushXML()

    def flushXML(self):
        newXml = ET.tostring(self.xmlroot)
        if len(newXml)>self.xmlLength:
            self.xmlLength = len(newXml)
            with open (self.makeFileName('PROGRAMXML')+'_tmp','w') as programXmlFile:
                programXmlFile.write(newXml)
            shutil.move(self.makeFileName('PROGRAMXML')+'_tmp', self.makeFileName('PROGRAMXML'))


    def makeCommandAndScript(self, **kw):
        self.hklout = os.path.join(self.workDirectory,"hklout.mtz")

        # ccp4-python -m lorestr.main
        # -f 1.mtz
        # -p1 1.pdb
        # -xyzout 1_lorestr1.pdb
        # -hklout 1_lorestr1.mtz
        # -o ./8_lorestr
        # -nc 3
        # -minres 3.5
        # -nh 5
        # -auto
        # -save_space

        # Input
        self.appendCommandLine(['-p1', self.inputCoordPath])
        self.appendCommandLine(['-f', self.hklin])
        if self.container.inputData.TLSIN.isSet():
            self.appendCommandLine(['-tls', self.container.inputData.TLSIN.fullPath])
        if self.container.inputData.DICT.isSet():
            self.appendCommandLine(['-libin', self.container.inputData.DICT.fullPath])
        p2 = ['-p2']
        for iCoordSet, xyzRef in enumerate(self.container.inputData.REFERENCE_LIST):
            p2.append(str(xyzRef.fullPath))
        if len(p2) > 1:
            self.appendCommandLine(p2)

        # Output
        self.appendCommandLine(['-xyzout', self.container.outputData.XYZOUT.fullPath])
        self.appendCommandLine(['-hklout', self.hklout])
        self.appendCommandLine(['-o', self.workDirectory + '/lorestrOutput'])

        # Parameters
        self.appendCommandLine(['-xml', str(self.makeFileName('PROGRAMXML'))]) # ccp4i2 compatibility

        if self.container.controlParameters.AUTO == 'all':
            self.appendCommandLine(['-auto'])
        elif self.container.controlParameters.AUTO == 'pdb':
            self.appendCommandLine(['-auto', 'pdb'])
        elif self.container.controlParameters.AUTO == 'af':
            self.appendCommandLine(['-auto', 'af'])

        if self.container.controlParameters.MR:
            self.appendCommandLine(['-mr'])

        if self.container.controlParameters.SS:
            self.appendCommandLine(['-save_space'])

        if self.container.controlParameters.DNA:
            self.appendCommandLine(['-dna'])

        if self.container.controlParameters.OVB:
            self.appendCommandLine(['-overalb'])

        if self.container.controlParameters.CPU:
            self.appendCommandLine(['-cpu', str(self.container.controlParameters.CPU)])

        if self.container.controlParameters.MINRES:
            self.appendCommandLine(['-minres', str(self.container.controlParameters.MINRES)])

        if self.container.controlParameters.NC:
            self.appendCommandLine(['-nc', str(self.container.controlParameters.NC)])

        if self.container.controlParameters.NH:
            self.appendCommandLine(['-nh', str(self.container.controlParameters.NH)])

        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        self.container.outputData.XYZOUT.annotation = 'Refined model'
        self.container.outputData.FPHIOUT.annotation = 'Weighted map after refinement'
        self.container.outputData.DIFFPHIOUT.annotation = 'Weighted difference map after refinement'
        self.container.outputData.FPHIOUT.subType = 1
        self.container.outputData.DIFFPHIOUT.subType = 2

        outputFiles = ['FPHIOUT','DIFFPHIOUT']
        outputColumns = ['FWT,PHWT','DELFWT,PHDELWT']

        error = self.splitHklout(outputFiles,outputColumns)

        xmlRoot = ET.parse(self.makeFileName('PROGRAMXML')).getroot()
        bestProtocolNumber = int(xmlRoot.find("Protocols/BestProtocol").text)
        bestRestraintsFileName = str(xmlRoot.find("Protocols/BestProtocolRestraintsFileName").text) # Contains absolute path
        rFact = xmlRoot.find("Protocols/P%d/Rfact" % bestProtocolNumber).text
        rFree = xmlRoot.find("Protocols/P%d/Rfree" % bestProtocolNumber).text

        self.container.outputData.PERFORMANCE.RFactor = float(rFact)
        self.container.outputData.PERFORMANCE.RFree = float(rFree)

        # Apparently CCP4i2 wants all output files to be strictly in self.getWorkDirectory(), not in the sub-directories. Copying then.
        exportRestraintsFileName = ''
        if os.path.exists(bestRestraintsFileName): # if the best protocol is jelly-body, then no restraints generated
            try:
                exportRestraintsFileName = os.path.join(self.getWorkDirectory(), os.path.basename(bestRestraintsFileName))
                shutil.copyfile(bestRestraintsFileName, exportRestraintsFileName)
            except:
                print('Can`t copy generated external restraints from %s to %s' % (bestRestraintsFileName, exportRestraintsFileName))

        self.container.outputData.EXTERNAL_RESTRAINTS_FILE.setFullPath(exportRestraintsFileName)
        self.container.outputData.EXTERNAL_RESTRAINTS_FILE.annotation.set('External restraints generated by the best refinement protocol')

        # CCP4i2 validation
        try:
           self.validate = self.makePluginObject('validate_protein')
           self.validate.container.inputData.XYZIN_1 = self.container.outputData.XYZOUT
           self.validate.container.inputData.XYZIN_2 = self.container.outputData.XYZOUT
           self.validate.container.inputData.F_SIGF_1 = self.container.inputData.F_SIGF
           self.validate.container.inputData.F_SIGF_2 = self.container.inputData.F_SIGF
           self.validate.container.outputData.COOTSCRIPTOUT = self.container.outputData.COOTSCRIPTOUT

           self.validate.container.controlParameters.DO_BFACT = True
           self.validate.container.controlParameters.DO_RAMA = True
           self.validate.container.controlParameters.DO_MOLPROBITY = True

           self.validate.doAsync = False
           self.validate.waitForFinished = -1
           self.validate.process()


           print('############################')
           validateXMLPath = self.validate.makeFileName('PROGRAMXML')
           print(validateXMLPath)

           validateXML = ET.parse(validateXMLPath).getroot()
           print(validateXML)

           self.xmlroot = ET.parse(self.makeFileName('PROGRAMXML')).getroot()
           print(self.xmlroot)

           xml_validation = ET.SubElement(self.xmlroot,"Validation")
           print(xml_validation)

           xml_validation.append(validateXML.xpath("//Validate_geometry_CCP4i2/B_averages")[0])

           print(validateXML.xpath("//Validate_geometry_CCP4i2/B_averages")[0])

           xml_validation.append(validateXML.xpath("//Validate_geometry_CCP4i2/Ramachandran_maps")[0])

           print(validateXML.xpath("//Validate_geometry_CCP4i2/Ramachandran_maps")[0])

           xml_validation.append(validateXML.xpath("//Validate_geometry_CCP4i2/Molprobity")[0])

           print(validateXML.xpath("//Validate_geometry_CCP4i2/Molprobity")[0])

           self.flushXML()
        except Exception as err:
           print("...Failed validation run after refinement", err)

        # End of CCP4i2 validation

        if error.maxSeverity()>CCP4ErrorHandling.Severity.WARNING:
            return CPluginScript.FAILED

        return CPluginScript.SUCCEEDED
