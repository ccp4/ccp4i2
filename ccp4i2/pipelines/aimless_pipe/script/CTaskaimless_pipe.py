"""
Copyright (C) 2012 STFC
Martyn Winn August 2012 - aimless_pipe gui
Phil Evans 2014
"""

import gemmi
from PySide2 import QtCore, QtWidgets

from ....pipelines.import_merged.script.dybuttons import ChoiceButtons
from ....pipelines.import_merged.script.mmcifutils import CifBlockInfo
from ....qtgui import CCP4TaskWidget
from .xdstype import Xdstype


class CTaskaimless_pipe(CCP4TaskWidget.CTaskWidget):
    
    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'aimless_pipe'
    TASKVERSION = 0.0
    TASKMODULE='data_reduction'
    TASKTITLE='Data reduction - AIMLESS'
    SHORTTASKTITLE='Data reduction'
    TASKLABEL='aimles'
    WHATNEXT = ['phaser_pipeline','molrep_pipe','crank2','shelx']
    PROGRAMHELP = ['aimless','pointless','ctruncate']
    RANK=1
    DESCRIPTION = '''Scale and analyse unmerged data and suggest space group (Pointless, Aimless, Ctruncate, FreeRflag)'''

    # -------------------------------------------------------------
    def __init__(self,parent):
        CCP4TaskWidget.CTaskWidget.__init__(self,parent)
    
    # -------------------------------------------------------------
    def drawContents(self):

        self.setProgramHelpFile('aimless')
        self.container.guiParameters.SCALING_PROTOCOL_SET.set(False)
        self.hasfreer = False  # true if any accepted CIF blocks have FreeR data
        

        # - - - - - - - - - - - - - - - - - - - - - - - - - - -  Input data folder
        folder = self.openFolder(folderFunction='inputData',title='Input Data')
                
        self.createLine( ['widget', '-title','Select unmerged data files',
                          'UNMERGEDFILES'])
        self.container.inputData.UNMERGEDFILES.dataChanged.connect( self.getUnknownCell)
        self.openSubFrame(toggle=[ 'SHOW_MMCIF_BLOCKS', 'open', [ True ] ] )
        self.setMMCIFframe()
        self.setCIFblocklist()
        self.closeSubFrame()

        s = 'NOTE: input file has multiple datasets, only one name reported in the list above'
        self.createLine(['subtitle',s],
                        toggle=['NFILES_WITH_MULTIPLE_DATASETS', 'closed', [0]])
        
        line = self.createLine( ['label','Resolution range (\xc5)',
                                 'widget', 'RESOLUTION_RANGE',
                                 'advice',' Maximum resolution in files ',
                                 'label','range to be determined'],
                                toggleFunction=[self.getMaximumResolution,['UNMERGEDFILES']] )
        # get label widget, note itemAt counts from 0
        self.maximum_resolution_label = line.layout().itemAt(3).widget()
        self.createLine(['widget','POINTLESS_USE_RESOLUTION_RANGE',
                         'label',
                         ' use explicit resolution range in symmetry determination as well as in scaling'])

        # Automatic resolution cutoff
        line_cutoff = self.createLine(['widget','AUTOCUTOFF',
                                'label',
                                ' automatically cut resolution range based on a first Aimless run'])
        #        self.createLine(['label', 'Minimum CC(1/2) for resolution estimation',
        #                        'widget','CCMINIMUM']
        #                        ,appendLine=line_cutoff,
        #                         toggle=['AUTOCUTOFF', 'closed', [False]])
        self.createLine(['label', 'Minimum information threshold',
                        'widget','INFOCONTENTTHRESHOLD',
                        'label', ' bits/reflection']
                        ,appendLine=line_cutoff,
                         toggle=['AUTOCUTOFF', 'closed', [False]])
                                      
        self.createLine( ['label','Options for symmetry determination', 'widget', 'MODE'] )
        self.container.controlParameters.MODE.dataChanged.connect( self.handleChangeMode)

        # CHOOSE mode options
        self.createLine( ['advice', 'Options for choice of space group or Laue group:'],
                         toggle=['MODE','open',['CHOOSE']]  )
        line_choose = self.createLine( ['widget', 'CHOOSE_MODE'],
                                       toggle=['MODE','open',['CHOOSE']]  )
        self.createLine( ['widget', 'CHOOSE_SOLUTION_NO'], appendLine=line_choose,
                        toggle=['CHOOSE_MODE', 'open', ['SOLUTION_NO']] )
        self.createLine( ['widget', 'CHOOSE_SPACEGROUP'], appendLine=line_choose,
                        toggle=['CHOOSE_MODE', 'open', ['SPACEGROUP','REINDEX_SPACE']] )
        self.createLine( ['Use reindex operator ', \
                          'widget', 'REINDEX_OPERATOR'],
                         toggle=['MODE', 'open', ['CHOOSE']] )

        #self.createLine( ['subtitle','<i>Optional input data</i>',
        #                  'Reference data and existing FreeR sets, if any'])


        # Reference data for MODE == MATCH: reflection data ...
        self.createLine( ['subtitle',
                          'Optional reference data to resolve indexing ambiguity and space group'+\
                          ' (also later option to scale against this reference set)',
                          'In cases where there are indexing ambiguities, these can be resolved '+\
                          'by providing as reference a previously processed data set or '+\
                          'coordinates. The space group in the reference will be assumed to be correct'])

        self.createLine( ['label','   ', 'widget', 'REFERENCE_FOR_AIMLESS',
                          'label', ' use reference data in analysis against Batch after scaling'],
                         toggle=['REFERENCE_DATASET_SET','closed',[False]])

        line_hklref = self.createLine( ['tip',
                          'The reference may be a merged or unmerged reflection list, of coordinates for Fcalc^2',
                                        'label','   Reference data are ','widget',  'REFERENCE_DATASET'])
        self.createLine( ['advice', ' and is optionally defined in next line'],
                          appendLine=line_hklref,
                          toggleFunction=[ self.refdataNotNeeded,
                                           [ 'MODE', 'HKLIN_REF', 'XYZIN_REF']] )
        self.createLine( ['advice', ' and MUST be defined in next line'],
                          appendLine=line_hklref,
                          toggleFunction=[ self.refdataNeeded,
                                           [ 'MODE', 'HKLIN_REF', 'XYZIN_REF']] )

        # reflection list ...
        self.createLine( ['widget','-browseDb', True,'HKLIN_REF'],
                        toggleFunction=[ self.openHKLIN_REF, [ 'MODE', 'REFERENCE_DATASET' ] ] )
        self.container.inputData.HKLIN_REF.dataChanged.connect( self.handleSelectHklinRef)
        # ... or coordinates
        self.createLine( ['widget', 'XYZIN_REF'],
                        toggleFunction=[ self.openXYZIN_REF, [ 'MODE', 'REFERENCE_DATASET' ] ] )
        self.container.inputData.XYZIN_REF.dataChanged.connect( self.handleSelectXyzinRef)

        # existing FreeR set
        self.createLine( ['subtitle',
                          'Optional existing FreeR set, define to copy or extend if necessary',
                          'If there are existing FreeR sets, you should choose a compatible one, '+\
                          'ie with same cell and point group'])
        
        self.createLine( ['widget', 'FREERFLAG'] )
        self.container.inputData.FREERFLAG.dataChanged.connect( self.handleSelectFreeRflag)
        
        self.createLine(['label', 'Fraction of reflections in generated freeR set',
                         'widget','FREER_FRACTION'],
                        toggle=['COMPLETE', 'closed', [True]] )
        self.createLine([ 'advice', 'Default fraction is 0.05.', 'advice', 'Potential twinning operations will be taken into account'],
                        toggle=['COMPLETE', 'closed', [True]] )

        self.createLine(['widget','CUTRESOLUTION',
                         'label',
                         'Cut resolution of FreeR set if necessary to match the data'],
                        toggle=['COMPLETE', 'closed', [False]] )


        # - - - - - - - - - - - - - - - - - - - - - - - - - - -  Important options folder
        self.openFolder(title='Important Options',drawFolder=self.drawImportant)

        

        # - - - - - - - - - - - - - - - - - - - - - - - - - - -  Additional options folder
        self.openFolder(title='Additional Options',drawFolder=self.drawAdditional)

    def drawImportant(self):

        self.createLine( ['subtitle',
                          'Options for symmetry determination in Pointless:   -------',
                       'Less common non-default options for crystal symmetry determination '])

        self.createLine( ['label','   Maximum resolution for scoring set by CC(1/2) in P1 > ',
                          'widget', 'CCHALFLIMIT','label',' [usual method, default 0.6]'])
        self.createLine( ['label','   Maximum resolution for scoring set by I/sigma(I) > ',
                          'widget', 'ISIGLIMIT','label',' [fall-back method, default 6]'])

        self.createLine( ['label','   Tolerance for comparing lattices (degrees or equivalent on lengths) ',
                          'widget', 'TOLERANCE','label',' [default 2.0]'])
        self.createLine(['label',' '])

        self.createLine( ['subtitle','Options for scaling and merging in Aimless:  -------',
                         'Less common non-default options in scaling and merging'])
        
        # Analysis, CC(1/2) & I/sigI thresholds
        self.analysis()

        # intensity (profile/summation) and partials
        self.intensitiesAndPartials()
        # parameters for SDcorrection
        self.SDcorrection()
        # scaling
        self.scalingDetails()
        unmergedline = self.createLine(['widget','OUTPUT_UNMERGED',
                         'advice','output unmerged data as well as merged'])
        self.createLine(['widget', 'ORIGINAL_HKL',
                         'advice','output original hkl instead of reduced hkl'],
                        appendLine=unmergedline,
                        toggle=['OUTPUT_UNMERGED','closed',[False]])

        self.createLine(['label',' '])

        # Parameters for Free R
        self.createLine( ['subtitle','Options for FreeR set extension:  -------'])
        self.openSubFrame(frame=True, toggle=['COMPLETE', 'closed', [False]])

        self.createLine( ['advice','If you are extending an existing FreeR set, it must match the observed data in unit cell and Laue group'])
        self.createLine( ['advice','The unit cells should match to the lower resolution of the two datasets'])
        self.createLine([
            'tip','DANGEROUS: only sensible if the unit cells are very similar',
            'widget', 'OVERRIDE_CELL_DIFFERENCE',
                         'label', 'allow existing freeR set to have different unit cells'])
        self.createLine(['advice',
          '<span style="color: DarkOrange;font-weight: bold;">Be sure you know what you are doing: the cells must be very similar even if outside the test limits</span>'])
        self.closeSubFrame()


        # - - - - - - - - - - - - - - - - - - - - - - - - - - -  folder

    def drawAdditional(self):

        self.createLine(['label',' '])
        self.createLine( ['subtitle',
                          'Options for symmetry determination in Pointless:   -------',
                          'Common non-default options for crystal symmetry determination '])
        self.createLine( ['subtitle','   Choice of cell setting conventions:',
                          'The default option is the IUCr standard, but you may prefer '+\
                          '"reference" settings, eg P2<sub>1</sub>2<sub>1</sub>2 rather than '+\
                          'P22<sub>1</sub>2<sub>1</sub>, and also not to use I2 instead of C2'])
        self.createLine( ['widget', 'SET_SETTING'] )

        self.createLine(['label',' '])

        # Aimless things
        self.createLine(['subtitle',
                          'Options for scaling in Aimless:   -------',
                         'Common non-default options in scaling and merging'])
        self.scaleProtocol()
        self.createLine(['widget','REFINE_REFERENCE',
                         'advice','determine scales relative to reference dataset (not normally recommended)'],
                        toggle=['REFERENCE_DATASET_SET','closed',[False]])
        self.container.controlParameters.REFINE_REFERENCE.dataChanged.connect\
           (self.handleSelectAimlessRefineReference)
        self.rejectOutliers()
        self.createLine(['label',' '])

        self.createLine(['subtitle',
                          'Options for both Pointless and Aimless:   -------',
                         'Run definitions by batch range apply to both steps'])
        self.createLine(['subtitle',
                         '   Override automatic definition of runs to mark discontinuities in data',
              'Batch numbers after the 1st file may be incremented by a multiple of 1000, so this must be taken into account here'])

        self.createLine(['tip','by file mode applies only to Pointless',
                         'label','Run selection options',
                         'widget','-guiMode','radio', 'RUN_MODE'])
        self.createLine(['widget',
                         '-title',
                         'Batch ranges to define runs (after any renumbering in Pointless), with optional high resolution limit',
                         'RUN_BATCHLIST'],
                        toggle=['RUN_MODE', 'open', ['BYRANGE']] )
        self.createLine(['label',' '])
        self.createLine(['label',' '])

####
        self.createLine( ['subtitle',
                          'Expert options, not for normal use:   -------',
                          '... use only if you know what you are doing',
                          'widget','EXPERT_OPTIONS',
                          'label','show expert options'])

        # . . . . . {
        self.openSubFrame(frame=True, toggle=['EXPERT_OPTIONS','open',[True]])
        self.latticeCentreOptions()

        self.createLine(['subtitle',
                         'Option to allow non-chiral space groups, not for biological macromolecules'])
        self.createLine(['widget','ALLOW_NONCHIRAL',
                         'label','allow non-chiral space groups'],
                        toggle=['EXPERT_OPTIONS', 'open'] )

        # Thresholds on CC(1/2) and multiplicity, testing for complet disastar
        self.createLine( ['subtitle',
          'Processing will be abandoned after Aimless scaling if CC(1/2) in the inner shell is very small, ie the data are too weak to be useful'])
        line = self.createLine(['label',
                                'Abandon job if CC(1/2) is less than ',
                                'widget','CCHALFDISASTER'],
                               toggle=['EXPERT_OPTIONS', 'open'] )
        self.createLine(['label',' and multiplicity is more then ',
                                'widget','MULTIPLICITYDISASTER'],
                               appendLine=line,
                               toggle=['EXPERT_OPTIONS', 'open'] )
        self.closeSubFrame()
        # . . . . . }
    # -------------------------------------------------------------
    def latticeCentreOptions(self):
        #  LATTICE option
        self.createLine(['subtitle',
                         'Option to remove centred lattice absences'])
        self.createLine(['widget','REMOVE_LATTICE_CENTERING',
                         'label','Remove centred lattice absences'])
        self.createLine(['label','Desired lattice centering type',
                         'widget','LATTICE_CENTERING'],
                        toggle=['REMOVE_LATTICE_CENTERING','open'])
        self.createLine(['advice',
          'This option should ONLY be used if you are sure that '+\
                         'the wrong cell was used in integration:'+\
                         ' probably better to redo the integration<br/>'+\
         'Note that not all centred lattices are consisent with all Bravais lattices,'+\
                         ' check the result carefully<br/>'+\
                         'Lattice type "P" is ignored'],
                        toggle=['REMOVE_LATTICE_CENTERING','open'])

    # -------------------------------------------------------------
    def getMaximumResolution( self ) :
        highRes = 1000000.0
        for i in range(len(self.container.inputData.UNMERGEDFILES)):
            if self.container.inputData.UNMERGEDFILES[i].file.fileContent.highRes.isSet():
                highRes = min(highRes, self.container.inputData.UNMERGEDFILES[i].file.fileContent.highRes.__float__())

        if  highRes > 999999.0:
            highRes = 0.0
        if self.maximum_resolution_label: self.maximum_resolution_label.setText("%5.2f\xc5" % highRes)
        return True

    # -------------------------------------------------------------
    @QtCore.Slot(bool)
    def getUnknownCell( self, force = False ) :
        # returns true is any files are of types which are missing the cell,
        # ie SCA unmerged or ShelX
        print("\n=====getUnknownCell\n")

        # 1) if any file is of type SCA unmerged or ShelX, set cell from GUI or from known cell
        # 2) if any file is of type SCA or ShelX, set wavelength from GUI or from known wavelength
        # 3) if any file is of type SCA or ShelX, or MERGED, set ONLYMERGE, else set DEFAULT
        #  Not fool-proof if there are multiple files, but these will be picked up by Pointless

        #        for i in range(len(self.container.inputData.UNMERGEDFILES)):
        #            # get dataset name etc, but don't signal, to avoid recursion
        #            self.container.inputData.UNMERGEDFILES[i].loadDatasetName(signal=False)
        
        anyunknowncell = False
        goodcell = None
        anyunknownwavelength = False
        goodwavelength = None
        anyunscaleablefile = False

        forceCell = False
        forceWavelength = False
        self.container.inputData.HKLIN_IS_MERGED.set(False)

        print("Number of files:", len(self.container.inputData.UNMERGEDFILES))
        for i in range(len(self.container.inputData.UNMERGEDFILES)):
            print("File:", self.container.inputData.UNMERGEDFILES[i].file)
            print("***getUnknownCell fileContent", i, \
                  self.container.inputData.UNMERGEDFILES[i].file.fileContent)
            print("**getUnknownCell self.container.inputData.UNMERGEDFILES[i].file",
                  self.container.inputData.UNMERGEDFILES[i].file)
            print("*getUnknownCell self.container.inputData.UNMERGEDFILES[i]",
                  self.container.inputData.UNMERGEDFILES[i])

            #self.container.inputData.UNMERGEDFILES[i].file.loadFile()
            #print("^^ filecontent",
            #     self.container.inputData.UNMERGEDFILES[i].file.fileContent)

            fformat = self.container.inputData.UNMERGEDFILES[i].file.fileContent.format
            xdsformat = False
            if fformat == 'mtz':
                self.container.inputData.HKLIN_FORMAT.set('MTZ')
            elif fformat == 'mmcif':
                self.container.inputData.HKLIN_FORMAT.set('MMCIF')
            elif fformat == 'sca':
                self.container.inputData.HKLIN_FORMAT.set('SCA')
            elif fformat == 'xds':
                self.container.inputData.HKLIN_FORMAT.set('XDS')
                xdsformat = True
            else:
                self.container.inputData.HKLIN_FORMAT.set('OTHER')
            if fformat == 'mmcif':
                if i == 0:
                    # 1st file is MMCIF
                    self.container.guiParameters.SHOW_MMCIF_BLOCKS.set(True)
                    #  index to mmcif file (always == 0)
                    self.container.inputData.MMCIF_FILE_INDEX.set(i)
                    self.openMmcifFile()
                    self.processMmcifFile(False)
                    anyunscaleablefile = True
                    if len(self.container.guiParameters.MMCIF_BLOCKNAMES) == 0:
                        message = 'No importable unmerged data blocks found'
                        QtWidgets.QMessageBox.warning(self,'Importing mmCIF file', message)
                        return
                    # default blockname is first one
                    blockname = self.container.guiParameters.MMCIF_BLOCKNAMES[0]
                    self.selectedBlock = blockname
                    self.extractMmcifInfo(blockname)
                    continue
                else:
                    if self.container.guiParameters.SHOW_MMCIF_BLOCKS:
                        message = "You can only have one mmCIF file"
                    else:
                        message = "mmCIF files cannot be mixed with other file types"
                    QtWidgets.QMessageBox.warning(self,'Importing mmCIF file', message)
                    return
            else:
                if i > 0 and self.container.guiParameters.SHOW_MMCIF_BLOCKS:
                    message = "mmCIF files cannot be mixed with other file types"
                    QtWidgets.QMessageBox.warning(self,'Importing mmCIF file', message)
                    return
                self.container.guiParameters.SHOW_MMCIF_BLOCKS.set(False)

            # not mmCIF
            print("shouldn't be mmcif", self.container.guiParameters.SHOW_MMCIF_BLOCKS)
            knowncell = self.container.inputData.UNMERGEDFILES[i].file.fileContent.knowncell
            # cell from file or GUI
            filecell = self.container.inputData.UNMERGEDFILES[i].cell
            contentcell = self.container.inputData.UNMERGEDFILES[i].file.fileContent.cell
            print("Filecell, knowncell: ", filecell, knowncell)
            print("contentcell", contentcell, contentcell.isSet())
            if knowncell:
                # cell from file or GUI
                if contentcell.isSet():
                    goodcell = contentcell
                    forceCell = True
                else:
                    goodcell = filecell
            else:
                anyunknowncell = True

            if xdsformat:
                # There are 4 types of XDS format with different information
                #  about wavelengths and scaleability
                xtype = Xdstype(str(self.container.inputData.UNMERGEDFILES[i].file))
                xtype.dump()
                knownwavelength = self.container.inputData.UNMERGEDFILES[i].file.fileContent.knownwavelength
                wavelength = self.container.inputData.UNMERGEDFILES[i].file.fileContent.wavelength
                if not xtype.scaleable:
                    anyunscaleablefile = True
            else:
                knownwavelength = self.container.inputData.UNMERGEDFILES[i].file.fileContent.knownwavelength
                wavelength = self.container.inputData.UNMERGEDFILES[i].file.fileContent.wavelength
                
                fileformat = self.container.inputData.UNMERGEDFILES[i].file.fileContent.format
                mergedfile = self.container.inputData.UNMERGEDFILES[i].file.fileContent.merged
                print("getUnknownCell fileformat, mergedfile", fileformat, mergedfile)
                # print("getUnknownCell filecontent",
                # self.container.inputData.UNMERGEDFILES[i].file.fileContent)
                if mergedfile == 'merged':
                    self.container.inputData.HKLIN_IS_MERGED.set(True)
                if (not knowncell) or (fileformat == 'sca') or (mergedfile == 'merged'):
                    print("getUnknownCell notscaleable")
                    anyunscaleablefile = True

            print('knownwavelength, wavelength', knownwavelength, wavelength)
            if knownwavelength:
                # wavelength from file or GUI
                goodwavelength = wavelength
                forceWavelength = True
            else:
                anyunknownwavelength = True


        # end loop files
                
        # Set cell if needed
        #print("getUnknownCell setting things", anyunscaleablefile)
        #print("auc, auw", anyunknowncell, anyunknownwavelength)
        #print("self.container.guiParameters.SCALING_PROTOCOL_SET", \
         #     self.container.guiParameters.SCALING_PROTOCOL_SET)
        if anyunknowncell:
            if goodcell is not None:
                self.container.controlParameters.CELL.set(goodcell)
            else:
                self.container.controlParameters.CELL.set(filecell)
        elif (force or forceCell):
            self.container.controlParameters.CELL.unSet()

        # Set wavelength if needed
        if anyunknownwavelength:
            if goodwavelength is not None:
                self.container.controlParameters.WAVELENGTH.set(goodwavelength)
            else:
                self.container.controlParameters.WAVELENGTH.set(wavelength)
        elif (force or forceWavelength):
             self.container.controlParameters.WAVELENGTH.unSet()

        if anyunscaleablefile:
            #print("**anyunscaleablefile")
            if len(self.container.inputData.UNMERGEDFILES) > 1:
                # Multiple input files, default CONSTANT
                self.container.controlParameters.ONLYMERGE.set(False)
                self.container.controlParameters.SCALING_PROTOCOL.set('CONSTANT')
                self.container.guiParameters.SCALING_PROTOCOL_SET.set(False)
            else:
                self.container.controlParameters.ONLYMERGE.set(True)
                self.container.controlParameters.SCALING_PROTOCOL.set('ONLYMERGE')
                self.container.guiParameters.SCALING_PROTOCOL_SET.set(False)
        else:
            if not self.container.guiParameters.SCALING_PROTOCOL_SET:
                #print("SCALING_PROTOCOL_SET false")
                self.container.controlParameters.ONLYMERGE.set(False)
                self.container.controlParameters.SCALING_PROTOCOL.set('DEFAULT')

        if self.container.inputData.HKLIN_IS_MERGED:
            # switch off optimisation of SD corrections
            self.container.controlParameters.SDCORRECTION_OVERRIDE.set(True)
            self.container.controlParameters.SDCORRECTION_REFINE.set(False)
                
        nfmult = 0
        for i in range(len(self.container.inputData.UNMERGEDFILES)):
            print(">> getUnknownCell  filecontent",
                self.container.inputData.UNMERGEDFILES[i].file.fileContent)
            ndatasets = self.container.inputData.UNMERGEDFILES[i].file.fileContent.numberofdatasets
            print("ndatasets", ndatasets, type(ndatasets))
            if ndatasets:
                ndatasets = int(ndatasets)
            else:
                ndatasets = 1
            print("ndatasets", ndatasets, type(ndatasets))
            self.container.inputData.NUMBER_OF_DATASETS.append(ndatasets)
            if ndatasets > 1:
                xdname = 'Multiple'
                self.container.inputData.UNMERGEDFILES[i].file.fileContent.crystalName.set(xdname)
                self.container.inputData.UNMERGEDFILES[i].file.fileContent.datasetName.set(xdname)
                nfmult += 1

            self.container.guiParameters.NFILES_WITH_MULTIPLE_DATASETS.set(nfmult)
            print(">> getUnknownCell end filecontent",
                self.container.inputData.UNMERGEDFILES[i].file.fileContent)

        print("**returning")
        return (not anyunknowncell)
 
    # -------------------------------------------------------------
    def scaleProtocol(self):
        #  define Aimless scaling protocol and parameters, for all runs (for now)

        line = self.createLine(['label','Scale','widget','SCALING_PROTOCOL'])
        self.container.controlParameters.SCALING_PROTOCOL.dataChanged.connect( self.handleChangeProtocol)

        self.createLine(['label',' with relative B-factor','widget','BFACTOR_SCALE'],
                        appendLine=line,
                        toggleFunction=[self.scaleProtocolBfactorOn,
                                        ['SCALING_PROTOCOL', 'BFACTOR_SCALE']])
        self.createLine(['label',' with no relative B-factor','widget','BFACTOR_SCALE'],
                        appendLine=line,
                        toggleFunction=[self.scaleProtocolBfactorOff,
                                        ['SCALING_PROTOCOL', 'BFACTOR_SCALE']])
        # . . . {
        self.openSubFrame(toggle=['SCALING_PROTOCOL','open',['ROTATION','SECONDARY']])
        line = self.createLine(['label','Define scale ranges along rotation axis by     ',
                                'widget','SCALES_ROTATION_TYPE'])
        self.createLine(['widget','SCALES_ROTATION_SPACING','label',' degrees'],
                        appendLine=line,
                        toggle=['SCALES_ROTATION_TYPE','open',['SPACING']])
        self.createLine(['widget','SCALES_ROTATION_NBINS'], appendLine=line,
                        toggle=['SCALES_ROTATION_TYPE','open',['NBINS']])
        line = self.createLine(['label','Define B-factor ranges along rotation axis by',
                                'widget','SCALES_BROTATION_TYPE'])
        self.createLine(['widget','SCALES_BROTATION_SPACING','label',' degrees'],
                        appendLine=line,
                        toggle=['SCALES_BROTATION_TYPE','open',['SPACING']])
        self.createLine(['widget','SCALES_BROTATION_NBINS'], appendLine=line,
                        toggle=['SCALES_BROTATION_TYPE','open',['NBINS']])
        
        self.closeSubFrame()
        # . . . }
        self.createLine(['label',
                         'Maximum order of spherical harmonics for secondary beam correction (eg 4 or 6) ',
                         'widget','SCALES_SECONDARY_NSPHHARMONICS'],
                        toggle=['SCALING_PROTOCOL','open',['SECONDARY']])
        # . . . {
        self.openSubFrame(toggle=['SCALING_PROTOCOL','open',['SECONDARY']])
        self.createLine(['label','Tile scaling for CCD detectors',
                         'widget','SCALES_TILETYPE'])
        self.createLine(['label','       number of tiles on X', 'widget', 'SCALES_NTILEX',
                         'label',' and on Y', 'widget', 'SCALES_NTILEY'],
                        toggle=['SCALES_TILETYPE','closed',['DEFAULT', 'NONE']])
        self.closeSubFrame()
        # . . . }
        
        self.closeSubFrame()
        # . . . . . }

    # -------------------------------------------------------------
    def rejectOutliers(self):
        # parameters for outlier rejection
        self.createLine(['widget','OUTLIER_OVERRIDE',
            'advice','override default parameters for outlier rejection'])

        # . . . . . {
        self.openSubFrame(frame=True, toggle=['OUTLIER_OVERRIDE','open',[True]])

        self.createLine(['label','Reject outliers if > ',
                                'widget','OUTLIER_SDMAX',
                                'label',' from mean, or > ',
                                'widget','OUTLIER_SDMAX2',
                                'label',' if 2 observations'])
        self.createLine(['label','Reject outliers between I+ and I- if > ', 
                         'widget','OUTLIER_SDMAXALL','label',', ',
                         'widget','OUTLIER_SDMAXALL_ADJUST',
                         'label',' increase for large anomalous differences'])
        line = self.createLine(['label','         ','widget','OUTLIER_COMBINE'])
        self.createLine(['label', ' compare outliers across all datasets'], appendLine=line,
                        toggle=['OUTLIER_COMBINE','open'])
        self.createLine(['label', ' compare outliers within each datasets'], appendLine=line,
                        toggle=['OUTLIER_COMBINE','closed'])

        self.createLine(['label', 'Set maximum E to reject unreasonably large intensities',
                        'widget','OUTLIER_EMAX'])

        self.closeSubFrame()
        # . . . . . }

    # -------------------------------------------------------------
    def SDcorrection(self):
        # parameters for SD correction
        self.createLine(['widget','SDCORRECTION_OVERRIDE',
                         'advice','override default parameters for SD correction'])

        # . . . . . {
        self.openSubFrame(frame=True, toggle=['SDCORRECTION_OVERRIDE','open',[True]])
        line = self.createLine(['widget','SDCORRECTION_REFINE'])
        self.createLine(['label','refine SDcorrection parameters'],appendLine=line,
                        toggle=['SDCORRECTION_REFINE','open',[True]])
        self.createLine(['label','do not refine SDcorrection parameters'],appendLine=line,
                        toggle=['SDCORRECTION_REFINE','open',[False]])

        self.createLine(['widget','SDCORRECTION_SET',
                         'label','set SD correction parameters: SdFac',
                         'widget','SDCORRECTION_SDFAC',
                         'label','   SdB','widget','SDCORRECTION_SDB',
                         'label','   SdAdd','widget','SDCORRECTION_SDADD'],
                        toggle=['SDCORRECTION_REFINE','open',[False]])
        self.container.controlParameters.SDCORRECTION_SDFAC.dataChanged.connect(self.handleSetSDcorrection)
        self.container.controlParameters.SDCORRECTION_SDB.dataChanged.connect(self.handleSetSDcorrection)
        self.container.controlParameters.SDCORRECTION_SDADD.dataChanged.connect(self.handleSetSDcorrection)

        
        line = self.createLine(['widget','SDCORRECTION_OPTIONS',
                                'label','to be determined'],
                               toggleFunction=[self.getSDcorrectionlabel,
                                               ['SDCORRECTION_REFINE','SDCORRECTION_OPTIONS']] )
        # get label widget, note itemAt counts from 0
        self.sdcorrection_label = line.layout().itemAt(1).widget()

        line = self.createLine(['widget','SDCORRECTION_FIXSDB','label','fix sdB parameters'],
                        toggle=['SDCORRECTION_REFINE','open',[True]])
        self.createLine(['spacing',150,'label','SD to tie sdB to zero',
                         'widget','SDCORRECTION_TIESDB_SD'],appendLine=line)

        line = self.createLine(['advice','Similarity restraint SDs: '],
                        toggleFunction=[self.SDcorrectionSimilar,
                                        ['SDCORRECTION_REFINE','SDCORRECTION_OPTIONS']])
        self.createLine(['label','sdFac','widget','SDCORRECTION_SIMILARITY_SDFAC'],
                        appendLine=line)
        self.createLine(['label',' sdB','widget','SDCORRECTION_SIMILARITY_SDB'],
                        appendLine=line)
        self.createLine(['label',' sdAdd','widget','SDCORRECTION_SIMILARITY_SDADD'],
                        appendLine=line)
        self.closeSubFrame()
        # . . . . . }

    # -------------------------------------------------------------
    @QtCore.Slot()
    def handleSetSDcorrection(self):
        self.container.controlParameters.SDCORRECTION_SET.set(True)

    # -------------------------------------------------------------
    def SDcorrectionSimilar(self):
        if self.container.controlParameters.SDCORRECTION_REFINE:
            if self.container.controlParameters.SDCORRECTION_OPTIONS == 'SIMILAR':
                return True
        return False
        
    # -------------------------------------------------------------
    def getSDcorrectionlabel(self):

        s = ""
        if self.container.controlParameters.SDCORRECTION_REFINE:
            if self.container.controlParameters.SDCORRECTION_OPTIONS == 'INDIVIDUAL':
                s = " different parameters for each run"
            elif self.container.controlParameters.SDCORRECTION_OPTIONS == 'SAME':
                s = " same parameters for each run"
            elif self.container.controlParameters.SDCORRECTION_OPTIONS == 'SIMILAR':
                s = " similar parameters for each run"
            self.sdcorrection_label.setText(s)
            return True
        else:
            self.sdcorrection_label.setText(s)
            return False

    # -------------------------------------------------------------
    def analysis(self):
        # parameters for analysis
        self.createLine(['widget','ANALYSIS_OVERRIDE',
                         'advice',
       'override default parameters for estimation of maximum resolution'])

        # . . . . . {
        self.openSubFrame(frame=True, toggle=['ANALYSIS_OVERRIDE','open',[True]])

        self.createLine(['label',
                         'Set minimum information content (average bits/reflection)',
                        'widget','INFOCONTENTTHRESHOLD'])

        self.createLine(['label', 'Set minimum CC(1/2)',
                        'widget','CCMINIMUM'])

        self.createLine(['label', 'Set minimum anomalous CC(1/2)',
                        'widget','CCANOMMINIMUM'])

        self.createLine(['label', 'Set minimum MnI/sigI',
                        'widget','ISIGMINIMUM'])
    
        self.closeSubFrame()
        # . . . . . }

    # -------------------------------------------------------------
    def intensitiesAndPartials(self):
        # parameters for intensity selection
        self.createLine(['widget','INTENSITIES_OVERRIDE',
                         'advice',
       'override default parameters for selection of intensities and treatment of partials'])

        # . . . . . {
        self.openSubFrame(frame=True, toggle=['INTENSITIES_OVERRIDE','open',[True]])
        line = self.createLine(['label','Use','widget','INTENSITIES_OPTIONS'])
        self.createLine(['label','(profile-fitted for weak intensities, summation for strong)'],
                        appendLine=line,
                        toggle=['INTENSITIES_OPTIONS','open',['COMBINE']])
        self.createLine(['widget','PARTIALS_TEST',
                         'label','Only accept partials with total fraction between',
                         'widget','PARTIALS_FRACLOW','label',' and',
                         'widget','PARTIALS_FRACHIGH'])
        self.container.controlParameters.PARTIALS_FRACLOW.dataChanged.connect(self.handlePartialstest)
        self.container.controlParameters.PARTIALS_FRACHIGH.dataChanged.connect(self.handlePartialstest)

        line = self.createLine(['widget','PARTIALS_CHECK'])
        self.createLine(['label','Scale partials outside rejection range'],
                        appendLine=line,
                        toggle=['PARTIALS_CHECK','open',[False]])
        self.createLine(['label','Scale partials in range',
                         'widget','PARTIALS_SCALE_MIN',
                         'label','to lower acceptance limit'],
                        appendLine=line,
                        toggle=['PARTIALS_CHECK','open',[True]])
        self.createLine(['widget','ACCEPT_OVERLOADS',
                         'label',' accept overloaded observations'])
        self.createLine(['widget','ACCEPT_EDGES',
                         'label',' accept observations on edge of tile or detector'])
        self.createLine(['widget','ACCEPT_XDS_MISFITS',
                         'label',' accept observations flagged by XDS as outliers (MISFITS)'])
    
        self.closeSubFrame()
        # . . . . . }

    # -------------------------------------------------------------
    @QtCore.Slot()
    def handlePartialstest(self):
        self.container.controlParameters.PARTIALS_TEST.set(True)

    # -------------------------------------------------------------
    def scalingDetails(self):
        # parameters for scaling
        self.createLine(['widget','SCALING_DETAILS',
                         'advice','override default parameters for scaling details'])

        # . . . . . {
        self.openSubFrame(frame=True, toggle=['SCALING_DETAILS','open',[True]])
        self.createLine(['widget','CYCLES_FLAG',
                         'label','Refine scale factors for',
                         'widget','CYCLES_N','label','cycles'])
        self.container.controlParameters.CYCLES_N.dataChanged.connect(self.handleNcycles)
        line = self.createLine(['widget','PARALLEL',
                                'label','use multiple processors'])
        self.container.controlParameters.PARALLEL.dataChanged.connect(self.handleParallel)
        self.createLine(['widget','PARALLEL_MODE'],
                         appendLine=line,
                         toggle=['PARALLEL','open',[True]])
        self.createLine(['label', 'determined from number of reflections'],
                         appendLine=line,
                         toggle=['PARALLEL_MODE','open',['AUTO']])
        self.createLine(['widget','NPROC','label','number of processors'],
                         appendLine=line,
                         toggle=['PARALLEL_MODE','open',['NUMBER']])
        self.createLine(['widget','FRACPROC','label','fraction of processors'],
                         appendLine=line,
                         toggle=['PARALLEL_MODE','open',['FRACTION']])
        self.createLine(['widget','SELECT1','widget','SELECT_IOVSDMIN',
                         'label','minimum I/sd for 1st round scaling   ',
                         'widget','SELECT2','widget','SELECT_EMIN',
                         'label','minimum E for 2nd round scaling'])
        self.container.controlParameters.SELECT_IOVSDMIN.dataChanged.connect(self.handleSelect1)
        self.container.controlParameters.SELECT_EMIN.dataChanged.connect(self.handleSelect2)
        self.createLine(['widget','TIE_ROTATION','label',
                         'Restrain neighbouring scale factors on rotation axis with SD',
                         'widget','TIE_ROTATION_SD'])
        self.container.controlParameters.TIE_ROTATION_SD.dataChanged.connect(self.handleTie_rotation_sd)
        self.createLine(['widget','TIE_SURFACE','label',
                         'Restrain surface parameters to a sphere with SD',
                         'widget','TIE_SURFACE_SD'])
        self.container.controlParameters.TIE_SURFACE_SD.dataChanged.connect(self.handleTie_surface_sd)
        self.createLine(['widget','TIE_BFACTOR','label',
                         'Restrain neighbouring B-factors on rotation axis with SD',
                         'widget','TIE_BFACTOR_SD'])
        self.container.controlParameters.TIE_BFACTOR_SD.dataChanged.connect(self.handleTie_bfactor_sd)
        self.createLine(['widget','TIE_BZERO','label',
                         'Restrain B-factors to zero with SD',
                         'widget','TIE_BZERO_SD'])
        self.container.controlParameters.TIE_BZERO_SD.dataChanged.connect(self.handleTie_bzero_sd)
        
        self.closeSubFrame()
        # . . . . . }

    # -------------------------------------------------------------
    @QtCore.Slot()
    def handleNcycles(self):
        self.container.controlParameters.CYCLES_FLAG.set(True)

    # -------------------------------------------------------------
    @QtCore.Slot()
    def handleParallel(self):
        if self.container.controlParameters.PARALLEL == False:
            self.container.controlParameters.PARALLEL_MODE.set('OFF')
        else:
            self.container.controlParameters.PARALLEL_MODE.set('AUTO')
                       
    # -------------------------------------------------------------
    @QtCore.Slot()
    def handleSelect1(self):
        self.container.controlParameters.SELECT1.set(True)

    # -------------------------------------------------------------
    @QtCore.Slot()
    def handleSelect2(self):
        self.container.controlParameters.SELECT2.set(True)

    # -------------------------------------------------------------
    @QtCore.Slot()
    def handleTie_rotation_sd(self):
        self.container.controlParameters.TIE_ROTATION.set(True)

    # -------------------------------------------------------------
    @QtCore.Slot()
    def handleTie_surface_sd(self):
        self.container.controlParameters.TIE_SURFACE.set(True)

    # -------------------------------------------------------------
    @QtCore.Slot()
    def handleTie_bfactor_sd(self):
        self.container.controlParameters.TIE_BFACTOR.set(True)

    # -------------------------------------------------------------
    @QtCore.Slot()
    def handleTie_bzero_sd(self):
        self.container.controlParameters.TIE_BZERO.set(True)

    # -------------------------------------------------------------
    def openHKLIN_REF( self ) :
        
        par = self.container.controlParameters
        if par.REFERENCE_DATASET == 'HKL':
            return True
        else:
            return False

    # -------------------------------------------------------------
    def openXYZIN_REF( self ) :

        par = self.container.controlParameters
        if par.REFERENCE_DATASET == 'XYZ':
            return True
        else:
            return False

    # -------------------------------------------------------------
    @QtCore.Slot()
    def handleSelectAimlessRefineReference(self):
        # from REFINE_REFERENCE
        if self.container.guiParameters.REFERENCE_DATASET_SET:
            if self.container.controlParameters.REFINE_REFERENCE:
                # Refine to reference selected, 
                self.container.controlParameters.REFERENCE_FOR_AIMLESS.set(True)
        else:
            self.container.controlParameters.REFINE_REFERENCE.set(False)
            self.container.controlParameters.REFERENCE_FOR_AIMLESS.set(False)

    # -------------------------------------------------------------
    @QtCore.Slot()
    def handleSelectHklinRef(self):

     if self.container.inputData.HKLIN_REF.isSet():
         self.container.controlParameters.MODE.set('MATCH')
         # set MODE = COMBINE if hklin is merged
         print(">> handleSelectHklinRef ", self.container.inputData.HKLIN_IS_MERGED)
         if self.container.inputData.HKLIN_IS_MERGED:
             self.container.controlParameters.MODE.set('LAUE')
             print("MODE = LAUE", self.container.inputData.HKLIN_REF)
         self.container.guiParameters.REFERENCE_DATASET_SET.set(True)
     else:
         # Don't forget users can unselect something!
         self.container.guiParameters.REFERENCE_DATASET_SET.set(False)
         self.container.controlParameters.REFERENCE_FOR_AIMLESS.set(False)
         self.container.controlParameters.REFINE_REFERENCE.set(False)
         if self.container.controlParameters.MODE == 'MATCH':
             self.container.controlParameters.MODE.set('LAUE')

    # -------------------------------------------------------------
    @QtCore.Slot()
    def handleSelectXyzinRef(self):

     if self.container.inputData.XYZIN_REF.isSet():
         self.container.controlParameters.MODE.set('MATCH')
         self.container.guiParameters.REFERENCE_DATASET_SET.set(True)
     else:
         # Don't forget users can unselect something!
         self.container.guiParameters.REFERENCE_DATASET_SET.set(False)
         if self.container.controlParameters.MODE == 'MATCH':
             self.container.controlParameters.MODE.set('LAUE')

    # -------------------------------------------------------------
    @QtCore.Slot()
    def handleChangeMode(self):

        if self.container.controlParameters.MODE != 'MATCH' and \
               not self.container.inputData.HKLIN_IS_MERGED:
            self.container.inputData.HKLIN_REF.unSet()
            self.container.inputData.XYZIN_REF.unSet()


    # -------------------------------------------------------------
    @QtCore.Slot()
    def handleSelectFreeRflag(self):

        if self.container.inputData.FREERFLAG.isSet():
            self.container.controlParameters.COMPLETE.set(True)
        else:
            self.container.controlParameters.COMPLETE.set(False)

    # -------------------------------------------------------------
    @QtCore.Slot()
    def handleChangeProtocol(self):
        #print("ChangeProtocol", self.container.controlParameters.SCALING_PROTOCOL)
        #print("Protocol", self.container.guiParameters.SCALING_PROTOCOL_SET)
        self.container.guiParameters.SCALING_PROTOCOL_SET.set(True)
        if self.container.controlParameters.SCALING_PROTOCOL != 'ONLYMERGE':
            self.container.controlParameters.ONLYMERGE.set(False)

    # -------------------------------------------------------------
    def refdataNotNeeded(self):
#        if self.container.inputData.HKLIN_REF.isSet() or \
#           self.container.inputData.XYZIN_REF.isSet():
#            return True

        if self.container.controlParameters.MODE != 'MATCH':
            return True

        return False
        
    # -------------------------------------------------------------
    def refdataNeeded(self):
        
        if self.container.controlParameters.MODE == 'MATCH':
            return True

        return False

    # -------------------------------------------------------------
    def scaleProtocolVariable(self):
        # True if scaling protocol not DEFAULT or ONLYMERGE
        if ((self.container.controlParameters.SCALING_PROTOCOL != 'DEFAULT') and
            (self.container.controlParameters.SCALING_PROTOCOL != 'ONLYMERGE')):
            return True
        return False

    # -------------------------------------------------------------
    def scaleProtocolBfactorOn(self):
        if (self.scaleProtocolVariable()
            and self.container.controlParameters.BFACTOR_SCALE):
            return True
        return False

    # -------------------------------------------------------------
    def scaleProtocolBfactorOff(self):
        if (self.scaleProtocolVariable()
            and (not self.container.controlParameters.BFACTOR_SCALE)):
            return True
        return False

    # -------------------------------------------------------------
    def openMmcifFile(self):
        i = 0
        print("openMmcifFile", self.container.inputData.UNMERGEDFILES[i].file)
        filename = str(self.container.inputData.UNMERGEDFILES[i].file.fullPath)
        self.mmcif = gemmi.cif.read(filename)
        self.rblocks = gemmi.as_refln_blocks(self.mmcif)
        if len(self.rblocks) == 0 or not self.rblocks[0]:
            # not a reflection mmcif file
            mess = QtWidgets.QMessageBox.warning(self,'Importing unmerged file',
                   'This does not seem to be a reflection mmcif file: may be coordinates?')
            return

        self.cifblockinfo = []
        for rb in self.rblocks:
            self.cifblockinfo.append(CifBlockInfo(rb))

        print("end openMmcifFile: nblocks = ", len(self.cifblockinfo))

    # -------------------------------------------------------------
    def processMmcifFile(self, acceptMerged=True):
        # acceptMerged True for merged blocks, False for unmerged
        print("** processMmcifFile", acceptMerged)

        nblock = len(self.rblocks)

        # for each accepted reflection block
        #  only list "accepted" blocks which have I or F data
        ids = []          #  block names
        detailslist = []  #  _diffrn.details
        infolist = []     #  info, hkl list type  
        columnlist = []   #  column labels
        #  information about blocks not accepted for input
        #    eg merged data,
        #       or phases etc, not Is or Fs    FIXME
        otherlist = []
        idxblkinfo = []  # index into main list for accepted blocks
        otheridxblkinfo = [] # ... and for unaccepted ones
        self.hasfreer = False  # true if any accepted blocks have FreeR data
        for idx, cifinfo in enumerate(self.cifblockinfo):
            if (acceptMerged and cifinfo.merged_diffn_data()) or \
               ((not acceptMerged) and (not cifinfo.ismerged())):
                # merged and merged or unmerged and unmerged
                #print("cifinfo", idx, cifinfo)
                idxblkinfo.append(idx)
                ids.append(cifinfo.bname)
                details = cifinfo.details   # _diffrn.details if present
                if details is None: details = ''
                detailslist.append(details)
                info = cifinfo.info   # formatted column info
                info += "\n    hkl list type: "+cifinfo.hklcheckformat
                freerwarning = self.formatFreeRinfo(cifinfo)
                if freerwarning is not None:
                    info += freerwarning
                if info is None: info = ''
                infolist.append(info)
                if cifinfo.haveFreeR():
                    self.hasfreer = True
                columns = cifinfo.columnnames
                if columns is None: columns = ''
                columnlist.append(columns)
            else:
                # accumulate strings identifying rejected blocks
                otherlist.append(self.otherinfo(cifinfo, acceptMerged))
                otheridxblkinfo.append(idx)

        #print("processMmcifFile ids", ids, "\n detailslist", detailslist,
        #     "\ninfolist", infolist, ",  idxblkinfo", idxblkinfo)
        #print("Lengths:", len(ids), len(detailslist), len(infolist))

        self.container.guiParameters.MMCIF_INDICES.set(idxblkinfo)
        self.container.guiParameters.MMCIF_BLOCKNAMES.set(ids)
        self.container.guiParameters.MMCIF_BLOCK_DETAILS.set(detailslist)
        self.container.guiParameters.MMCIF_BLOCK_INFO.set(infolist)
        self.setColumnNames(columnlist)
        # non-accepted blocks
        self.container.guiParameters.MMCIF_BLOCK_OTHER.set(otherlist)
        
        mbnl = list(self.container.guiParameters.MMCIF_BLOCKNAMES)
        for b in self.container.guiParameters.MMCIF_BLOCKNAMES:
            s = str(b)
            self.extractMmcifInfo(s)
            self.setCIFblocklist()

    # -------------------------------------------------------------
    def otherinfo(self, cifinfo, acceptMerged=True):
        # Information about an mmcif not suitable for merged or unmerged input
        s = ">>> mmCIF block " + cifinfo.bname + ": "
        details = cifinfo.details   # _diffrn.details if present
        if not cifinfo.ismerged() and acceptMerged:
            s += "Unmerged data"
            if details is not None:
                s += ": "+details
        else:
            info = cifinfo.info   # formatted column info
            info += "\n   hkl list type: "+cifinfo.hklcheckformat
            s += info
            
            if details is not None:
                s += ": "+details

            #  cifinfo.columnnames is a dictionary
            columnlist = self.formatColumnNames(cifinfo.columnnames)
            s += '\n   '+columnlist
            #print("otherinfo columnnames", columnlist,"\n>",s)

        freerwarning = self.formatFreeRinfo(cifinfo)
        if freerwarning is not None:
            s += freerwarning

        return s  

    # -------------------------------------------------------------
    def setColumnNames(self, columnlist):
        # columnlist is a dictionary of columns names indexed by content
        print("setColumnNames columnlist", columnlist)
        # Possible content types, in order of priority
        contenttypes = ['unmerged I', 'I+- anomalous', 'F+- anomalous', 'Imean', 'Fmean']
        columnlists = []
        for clist in columnlist:
            collist = None
            for ctype in contenttypes:
                if ctype in clist:
                    collist = clist[ctype]
                    fclmns = self.formatColumnlist(collist)
                    columnlists.append(fclmns)
                    break

        #print("columnlists", columnlists)
        self.container.guiParameters.MMCIF_BLOCK_COLUMNS.set(columnlists)
    # -------------------------------------------------------------
    def formatColumnNames(self, columnlist):
        # columnlist is a dictionary of columns names indexed by content
        # Return list of column names, line-warpped as needed
        
        # Possible content types, in order of priority
        s = ''
        maxlinelength = 90  # characters
        nkeys = len(columnlist)
        i = 0
        for key in columnlist:
            collist = columnlist[key]
            fclmns = self.formatColumnlist(collist)
            if i < nkeys:
                # not last
                if len(s)+len(fclmns) > maxlinelength:
                    s += "\n   "
            s += str(key) + ': [' + fclmns + '], '
            i += 1
        return s[:-2]

    # -------------------------------------------------------------
    def formatColumnlist(self, collist):
        # collist is a list of strings, return string
        s = ''
        if len(collist) == 0: return s
        for label in collist:
            s += label + ", "
        return s[:-2]
    # -------------------------------------------------------------
    def formatFreeRinfo(self, cifinfo):
        # return string with warnings about invalid FreeR, = None if no FreeR
        messagelist = cifinfo.freerWarning()
        if messagelist is None: return None
        message = ''
        for m in messagelist:
            message += '\n  '+m
        return message
    # -------------------------------------------------------------
    def extractMmcifInfo(self, blockname=None):
        print("**extractMmcifInfo", blockname)

        nblocks = len(self.rblocks)
        naccepted = len(self.container.guiParameters.MMCIF_BLOCK_DETAILS)

        cifinfo = None
        jblock = 0
        if blockname is None or nblocks == 1:
            # take the first one
            rblock = self.rblocks[0]
            cifinfo = self.cifblockinfo[0]
            blockname = cifinfo.bname
        else:
            #print("len(self.cifblockinfo)", len(self.cifblockinfo))
            # Find block in list of accepted blocks
            for i in range(naccepted):
                j = self.container.guiParameters.MMCIF_INDICES[i]
                if self.cifblockinfo[j].bname == blockname:
                    cifinfo = self.cifblockinfo[j]
                    jblock = i  # selected block
                    break

        if cifinfo is None:
            print("***Failed") # shouldn't happen111111111111111111111111111

        #print("extractMmcifInfo jblock, len", jblock,
        #     len(self.container.guiParameters.MMCIF_BLOCK_DETAILS))

        # print("*Rcell:",cifinfo.cell)
        # print("RSG:", cifinfo.spacegroup_name)
        # print("*Rwvl:",cifinfo.wavelength)

        #rc = cifinfo.cell
        #self.container.inputData.SPACEGROUPCELL.cell.set(rc)

        #self.container.inputData.SPACEGROUPCELL.spaceGroup.set(cifinfo.spacegroup_name)
        #self.container.inputData.WAVELENGTH.set(cifinfo.wavelength)

        #self.container.inputData.DATASETNAME.set(blockname)
        #self.container.inputData.CRYSTALNAME.set(blockname)

        fidx = self.container.inputData.MMCIF_FILE_INDEX   # index to input file list, = 0

        # set things into file list
        #   disconnect dataChanged connection to avoid recursion
        self.container.inputData.UNMERGEDFILES.dataChanged.disconnect( self.getUnknownCell)
        self.container.inputData.UNMERGEDFILES[fidx].file.fileContent.cell.set(cifinfo.cell)
        self.container.inputData.UNMERGEDFILES[fidx].file.fileContent.spaceGroup.set(cifinfo.spacegroup_name)
        self.container.inputData.UNMERGEDFILES[fidx].file.fileContent.wavelength.set(cifinfo.wavelength)

        self.container.inputData.UNMERGEDFILES[fidx].crystalName.set(blockname)
        self.container.inputData.UNMERGEDFILES[fidx].dataset.set(blockname)

        self.getWidget('UNMERGEDFILES').validate()
        self.container.inputData.UNMERGEDFILES.dataChanged.connect( self.getUnknownCell)
        #self.getWidget('CRYSTALNAME').validate()

        #print("extractMmcifInfo lengths",
        #     len(self.container.guiParameters.MMCIF_BLOCK_INFO[jblock]),
        #     len(self.container.guiParameters.MMCIF_BLOCK_COLUMNS[jblock]))
        self.container.inputData.MMCIF_SELECTED_BLOCK.set(blockname)
        self.container.inputData.MMCIF_SELECTED_DETAILS.set(\
            self.container.guiParameters.MMCIF_BLOCK_DETAILS[jblock])
        self.container.inputData.MMCIF_SELECTED_INFO.set(\
          self.container.guiParameters.MMCIF_BLOCK_INFO[jblock])
        if len(self.container.guiParameters.MMCIF_BLOCK_COLUMNS) > 0:
            self.container.inputData.MMCIF_SELECTED_COLUMNS.set(\
                self.container.guiParameters.MMCIF_BLOCK_COLUMNS[jblock])

        #  +1 if intensity, -1 if amplitude, 0 if unknown
        isintensity = -1
        if 'I' in str(self.container.inputData.MMCIF_SELECTED_INFO):
            # only intensity contents contains the letter 'I'
            isintensity = +1
        else:
            isintensity = -1

        self.container.inputData.MMCIF_SELECTED_ISINTENSITY.set(isintensity)

        merged = cifinfo.ismerged()
        if merged:
            mess = 'Block '+blockname+\
                   ' apppears to be an merged data block - please use the Import Merged task to import it'
            mess = QtWidgets.QMessageBox.warning(self,
                                           'Importing unmerged file', mess)
            return
        
        # Check reflection label sets
        if not cifinfo.OK:
            # Probably not a reflection block, eg maybe coordinates
            mess = QtWidgets.QMessageBox.warning(self,'Importing unmerged file',
     'This does not seem to be a reflection mmcif file: may be coordinates?')
            return
      
        message = ""
        if nblocks == 1:
            message = "This file"
        else:
            message = "This data block"
          
        status, mess = cifinfo.columnsOK()
        message += mess

        if status < 0:
            mess = QtWidgets.QMessageBox.warning(self,
                                               self.windowTitle(), message)
      
    # -------------------------------------------------------------
    def setMMCIFframe(self):
        # just make an empty area for later, if needed
        # and if HKLIN is MMCIF, create file data objects
        print("markCIFframe")
        if self.widget.subFrame is not None:
            self.cifpane = self.widget.subFrame.layout()
        else:
            self.cifpane = self.widget.currentFolderLayout

        self.cifbuttons = ChoiceButtons()
        self.cifpane.addWidget(self.cifbuttons)
        self.selectedBlock = ''
        self.cifbuttons.clickedSignal.connect(self.cifblockClicked)

        if self.container.inputData.HKLIN_FORMAT == 'MMCIF':
            # make self.mmcif & .rblocks
            self.openMmcifFile()

    # -------------------------------------------------------------
    @QtCore.Slot(str)
    def cifblockClicked(self):
        s = self.cifbuttons.selected
        # print("Clicked", s)
        self.selectedBlock = s
        self.extractMmcifInfo(s)
      
    # -------------------------------------------------------------
    def setCIFblocklist(self):
        if self.container.inputData.HKLIN_FORMAT != 'MMCIF':
            return

        infolist = ['']
        if self.container.guiParameters.MMCIF_BLOCK_INFO:
            infolist = self.strlist(self.container.guiParameters.MMCIF_BLOCK_INFO)
            infolist = ['Column content type: '+x for x in infolist]

        # print("Info", self.container.guiParameters.MMCIF_BLOCK_INFO)
        
        nblocks = len(self.container.guiParameters.MMCIF_BLOCKNAMES)
        if nblocks > 1:
            title = 'Select which mmCIF unmerged reflection block to use (first is default)'
        elif nblocks == 1:
            title = 'mmCIF file contains one appropriate unmerged reflection block'
        else:
            title = 'NB mmCIF file contains NO appropriate unmerged reflection blocks'

        subtitle = None
        if self.hasfreer:
            subtitle = "NOTE FreeR data will NOT be inported from unmerged blocks"

        self.cifbuttons.setChoices(title,
            self.strlist(self.container.guiParameters.MMCIF_BLOCKNAMES),
            tags=self.strlist(self.container.guiParameters.MMCIF_BLOCK_DETAILS),
            notes=infolist,subtitle=subtitle)

        if self.container.guiParameters.MMCIF_BLOCK_OTHER.isSet():
            if len(self.container.guiParameters.MMCIF_BLOCK_OTHER) > 0:
                  self.cifbuttons.addOtherText(\
                 "Reflection blocks not containing importable unmerged reflection data",
                      self.strlist(self.container.guiParameters.MMCIF_BLOCK_OTHER))

    # -------------------------------------------------------------
    def strlist(self, cstringlist):
        # string list from CString list
        sl = []
        for i in range(len(cstringlist)):
            sl.append(str(cstringlist[i]))
        return sl
