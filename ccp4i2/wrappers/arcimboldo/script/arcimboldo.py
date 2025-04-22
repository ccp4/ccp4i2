"""
Copyright (C) 2014 The University of York
"""

## @package arcimboldo
# This script runs all three versions of arcimboldo

from distutils.dir_util import copy_tree
import os

from lxml import etree

from ....core import CCP4Utils
from ....core import CCP4XtalData
from ....core.CCP4Modules import PROCESSMANAGER
from ....core.CCP4PluginScript import CPluginScript


ccp4_home = os.environ.get ( "CCP4", "not_set" )

class arcimboldo(CPluginScript):
    TASKTITLE = 'Arcimboldo'
    TASKNAME = 'arcimboldo'
    TASKCOMMAND = 'ARCIMBOLDO_LITE'
    TASKVERSION= 0.1
    WHATNEXT = [ 'prosmart_refmac' ]
    MAINTAINER = 'jtvcri@ibmb.csic.es'
    PURGESEARCHLIST = [['hklin*.mtz', 0]]
    ASYNCHRONOUS = True

    def genHKL(self, hklin):
        binf = os.path.normpath(os.path.join( CCP4Utils.getCCP4Dir().__str__(), 'bin', 'mtz2hkl' ))
        arglist = ['-f', hklin.__str__()]
        pid = PROCESSMANAGER().startProcess(binf, arglist)
        return PROCESSMANAGER().getJobData(pid, 'exitCode')

    def generateBor (self, hklin, columns):
        inputData = self.container.inputData
        advancedData = self.container.advancedData
        controlParameters = self.container.controlParameters
        guiAdmin = self.container.guiAdmin
        f_bor = open(os.path.join(self.getWorkDirectory(),'setup.bor'),'w')
        f_bor.write('[CONNECTION]\n')
        f_bor.write('distribute_computing = %s\n' % (controlParameters.ARCIMBOLDO_RUN))
        if controlParameters.ARCIMBOLDO_RUN != 'multiprocessing':
            if controlParameters.RUN_MODE == 'SYSTEM':
                f_bor.write('setup_bor_path =  %s/lib/python3.7/site-packages/arcimboldo/setup.bor\n' % (ccp4_home))
            else:
                config_file = os.path.join(inputData.CONFIG_FILE.relPath.__str__(),inputData.CONFIG_FILE.baseName.__str__())
                f_bor.write('setup_bor_path = %s\n' % (config_file))
        f_bor.write('[GENERAL]\n')
        f_bor.write('working_directory = %s\n' % (self.getWorkDirectory()))
        f_bor.write('mtz_path = %s\n' % (hklin))
        base = os.path.splitext(hklin)[0]
        f_bor.write('hkl_path = %s.hkl\n' % (base))
        if controlParameters.ARCIMBOLDO_OPTIONS == 'LITE':
            f_bor.write('[ARCIMBOLDO]\n')
            f_bor.write('rmsd = %6.2f\n' % (controlParameters.LITE_RMSD))
            if controlParameters.LITE_MODELS == 'HELIX':
                f_bor.write('fragment_to_search = %d\n' % (controlParameters.N_FRAGMENTS)) 
                f_bor.write('helix_length = %d\n' % (controlParameters.HELIX_LENGTH)) 
            elif controlParameters.LITE_MODELS == 'CUSTOM':
                f_bor.write('fragment_to_search = %d\n' % (controlParameters.N_FRAGMENTS)) 
                f_bor.write('model_file = %d\n' % (inputData.PDB_LITE)) 
            elif controlParameters.LITE_MODELS == 'HELICES':
                frag_count = 0
                i = 0
                f_bor.write('fragment_to_search = %d\n' % (len(inputData.LITE_HELICES_LIST))) 
                while i < len(inputData.LITE_HELICES_LIST):
                    f_bor.write('helix_length_%d = %d\n' % (i+1,inputData.LITE_HELICES_LIST[i])) 
                    i += 1
            elif controlParameters.LITE_MODELS == 'CUSTOMS':
                frag_count = 0
                i = 0
                f_bor.write('fragment_to_search = %d\n' % (len(inputData.LITE_CUSTOMS_LIST))) 
                while i < len(inputData.LITE_CUSTOMS_LIST):
                    f_bor.write('model_file_%d = %s\n' % (i+1,inputData.LITE_CUSTOMS_LIST[i])) 
                    i += 1
            
            if controlParameters.LITE_PARTIAL:
                f_bor.write('fixed_models_directory = %s\n' % (inputData.LITE_FIXED))

        elif controlParameters.ARCIMBOLDO_OPTIONS == 'BORGES':
            f_bor.write('[ARCIMBOLDO-BORGES]\n')
            if controlParameters.BORGES_LIBRARY == 'CUSTOM':
                f_bor.write('library_path = %s\n' % (inputData.BORGES_CUSTOM))
            else:
                ccp4_master_home = os.environ.get ( 'CCP4_MASTER', 'not_set' )
                lib_path = os.path.join(ccp4_master_home,'BORGES_LIBS',controlParameters.BORGES_LIBRARY.__str__())
                f_bor.write('library_path = %s\n' % (lib_path))
            if controlParameters.BORGES_GYRE:
                f_bor.write('rotation_model_refinement = %s\n' % (controlParameters.BORGES_GYRE_T))
            if controlParameters.BORGES_GIMBLE:
                f_bor.write('gimble = %s\n' % (controlParameters.BORGES_GIMBLE_T))
            if controlParameters.BORGES_MULTICOPY:
                f_bor.write('multicopy = %s\n' % (controlParameters.BORGES_MULTICOPY_T))                    
        else:
            f_bor.write('[ARCIMBOLDO-SHREDDER]\n')
            f_bor.write('model_file = %s\n' % (inputData.PDB_SHREDDER))
            if controlParameters.SHREDDER_RMSD:
                f_bor.write('rmsd_shredder = %6.2f\n' % (controlParameters.SHREDDER_RMSD_T))
            if controlParameters.SHREDDER_CONVERT:
                f_bor.write('trim_to_polyala = %s\n' % (controlParameters.SHREDDER_CONVERT_T))
            if controlParameters.SHREDDER_MAKE:
                f_bor.write('bfacnorm = %s\n' % (controlParameters.SHREDDER_MAKE_T))
            f_bor.write('shred_method = %s\n' % (controlParameters.SHREDDER_OPTIONS))
            if controlParameters.SHREDDER_OPTIONS == 'spherical':
                if(controlParameters.FRAGMENT_SIZE and controlParameters.FRAGMENT_SIZE_T != None):
                    fragment_size = str(controlParameters.FRAGMENT_SIZE_T)
                else:
                    fragment_size = "default"
                if(controlParameters.SHREDDER_COIL):
                    f_bor.write('sphere_definition = %s 1 %s 7 4 0.45 0.3\n' % (fragment_size, controlParameters.SHREDDER_COIL_T))
                elif fragment_size != "default":
                    f_bor.write('sphere_definition = %s 1 remove_coil 7 4 0.45 0.3\n' % (fragment_size))
    
                if(controlParameters.SHREDDER_GYRE):
                    f_bor.write('rotation_model_refinement = %s\n' % (controlParameters.SHREDDER_GYRE_T))
                if(controlParameters.SHREDDER_GIMBLE):
                    f_bor.write('gimble = %s\n' % (controlParameters.SHREDDER_GIMBLE_T))
                if(controlParameters.SHREDDER_LLG):
                    f_bor.write('occ = %s\n' % (controlParameters.SHREDDER_LLG_T))
                if(controlParameters.SHREDDER_COMBINE):
                    f_bor.write('alixe = %s\n' % (controlParameters.SHREDDER_COMBINE_T))
                if(controlParameters.SHREDDER_MULTICOPY):
                    f_bor.write('multicopy = %s\n' % (controlParameters.SHREDDER_MULTICOPY_T))
                f_bor.write('predicted_model = %s\n' % (controlParameters.SHREDDER_PREDICTED))
        #Common configuration values
        f_bor.write('name_job = arcimboldo\n')
        f_bor.write('coiled_coil = %s\n' % (controlParameters.COIL_COILED))
        f_bor.write('molecular_weight = %10.2f\n' % (controlParameters.MOLECULAR_WEIGHT))
        f_bor.write('number_of_component = %d\n' % (controlParameters.N_COMPONENTS))
        if(controlParameters.TNCS):
            f_bor.write('TNCS = %s\n' % (controlParameters.TNCS_T))
        if advancedData.SHELXE_LINE.isSet():
            f_bor.write('shelxe_line = %s\n' % (advancedData.SHELXE_LINE))
        if advancedData.KEYWORDS.isSet():
            f_bor.write('%s\n' % (advancedData.KEYWORDS))
        if inputData.F_SIGF.contentFlag == 2 or inputData.F_SIGF.contentFlag == 4:
            f_bor.write('f_label = %s\n' % (columns[0]))
            f_bor.write('sigf_label = %s\n' % (columns[1]))
        elif inputData.F_SIGF.contentFlag == 1 or inputData.F_SIGF.contentFlag == 3:
            f_bor.write('i_label = %s\n' % (columns[0]))
            f_bor.write('sigi_label = %s\n' % (columns[1]))
        f_bor.write('[LOCAL]\n')
        f_bor.write('path_local_phaser = %s/bin/phaser\n' % (ccp4_home))
        f_bor.write('path_local_shelxe = %s/bin/shelxe\n' % (ccp4_home))
        f_bor.close()

    def generateProgram(self):
        guiAdmin = self.container.guiAdmin
        nameJob = str(guiAdmin.jobTitle)
        pathHtml = str(os.path.join(self.getWorkDirectory(),'arcimboldo.html'))
        pathXml = str(os.path.join(self.getWorkDirectory(),'arcimboldo.xml'))
        self.programXml = etree.Element('arcimboldo')
        element = etree.SubElement(self.programXml, 'nameJob')
        element.text = nameJob
        element = etree.SubElement(self.programXml, 'pathHtml')
        element.text = pathHtml
        element = etree.SubElement(self.programXml, 'pathXml')
        element.text = pathXml
        with open(self.makeFileName('PROGRAMXML'), 'w+') as xml:
            xml.write(etree.tostring(self.programXml, encoding='unicode', pretty_print=True))
        self.watchFile(pathXml, self.refreshXML)

    def processInputFiles ( self ):
        list_of_stuff = [ ]
        inputData = self.container.inputData
        developerOptions = self.container.developerOptions
        if inputData.F_SIGF.contentFlag == 1 or inputData.F_SIGF.contentFlag == 3:
            list_of_stuff.append ( [ 'F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_IMEAN ] )
        elif inputData.F_SIGF.contentFlag == 2 or inputData.F_SIGF.contentFlag == 4:
            list_of_stuff.append ( [ 'F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN ] )  
        self.hklin, self.columns, error = self.makeHklin0 ( list_of_stuff )
        self.columns = self.columns.split(',')
        exitCode = self.genHKL(self.hklin)
        if exitCode != 0:
            return CPluginScript.FAILED
        if developerOptions.DEVELOPER_MODE == 'EXISTING':
            copy_tree(str(developerOptions.EXISTING_FOLDER), str(self.getWorkDirectory()), update=1)
        self.generateBor(self.hklin, self.columns)
        self.generateProgram()
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self):
        controlParameters = self.container.controlParameters
        if controlParameters.ARCIMBOLDO_OPTIONS == 'LITE':
            self.TASKCOMMAND = 'ARCIMBOLDO_LITE'
        elif controlParameters.ARCIMBOLDO_OPTIONS == 'BORGES':
            self.TASKCOMMAND = 'ARCIMBOLDO_BORGES'
        else:
            self.TASKCOMMAND = 'ARCIMBOLDO_SHREDDER'
        if self.container.developerOptions.DEVELOPER_MODE != 'BOR':
            self.appendCommandLine([os.path.join(self.getWorkDirectory(),'setup.bor')])
        return CPluginScript.SUCCEEDED

    def process(self):
        CPluginScript.process(self)

    def processOutputFiles(self):
        outputData = self.container.outputData
        pdbout = os.path.join(self.getWorkDirectory(), "best.pdb")
        if os.path.exists(pdbout):
            outputData.XYZOUT.append(outputData.XYZOUT.makeItem())
            outputData.XYZOUT[-1].setFullPath(pdbout)
            outputData.XYZOUT[-1].annotation = 'Best pdb solution'
    #      phsout = os.path.join(self.getWorkDirectory(), "best.phs")
    #      if os.path.exists(phsout):
    #         outputData.PHSOUT.append(outputData.PHSOUT.makeItem())
    #         outputData.PHSOUT[-1].setFullPath(phsout)
    #         outputData.PHSOUT[-1].annotation = 'Best phs solution'
        return CPluginScript.SUCCEEDED

    def refreshXML(self, filename):
        tmpFilename = self.makeFileName('PROGRAMXML') + '_tmp'
        with open(tmpFilename, 'w+') as xmlFile:
            xmlFile.write(etree.tostring(self.programXml, encoding='unicode', pretty_print=True))
        if os.path.exists(self.makeFileName('PROGRAMXML')):
            os.remove(self.makeFileName('PROGRAMXML'))
        os.rename(tmpFilename, self.makeFileName('PROGRAMXML'))
