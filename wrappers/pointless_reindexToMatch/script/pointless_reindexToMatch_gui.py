from __future__ import print_function
"""
    pipelines/pohaser_pipeline/wrappers/pointless_reindexToMatch/script/pointless_reindexToMatch_gui.py
    Copyright (C) 2014 Newcastle University
    Author: Martin Noble
    
    """

from qtgui.CCP4TaskWidget import CTaskWidget
from ccp4i2.baselayer import QtCore

#-------------------------------------------------------------------
class pointless_reindexToMatch_gui(CTaskWidget):
    #-------------------------------------------------------------------
    
    # Subclass CTaskWidget to give specific task window
    TASKMODULE=['expt_data_utility']
    TASKTITLE='Reindex reflections or change spacegroup'
    DESCRIPTION = 'Reindex: match to reference data/coordinates; change space group; analyse symmetry; or expand to P1 (Pointless)'
    TASKNAME = 'pointless_reindexToMatch'
    TASKVERSION = 0.1
    RANK = 2

    ERROR_CODES = { 301 : { 'description' : 'You must specify a space group or reindex operator' } }

    
    def drawContents(self):
        self.openFolder(folderFunction='inputData', title='Input')
        self.createLine(['subtitle','Reflection objects to manipulate'])
        self.openSubFrame(frame=True)

        self.createLine(['widget','-browseDb',True,'F_SIGF'])
        self.container.inputData.F_SIGF.dataChanged.connect(self.getSpacegroup)

        self.createLine(['widget','-browseDb',True,'FREERFLAG'])
        self.closeSubFrame()
        
        self.createLine(['subtitle','Spacegroup and indexing'])
        self.createLine(['spacing',5,
                         'label','Define new indexing and spacegroup using',
                         'widget','REFERENCE'],
                        toggle=['REFERENCE','closed',['ANALYSE','EXPAND']])
        self.createLine(['spacing',5,
                         'label','Analyse data symmetry',
                         'widget','REFERENCE'],
                        toggle=['REFERENCE','open',['ANALYSE']])

        self.createLine(['spacing',5,
                         'label','Expand to space group P1',
                         'widget','REFERENCE'],
                        toggle=['REFERENCE','open',['EXPAND']])

        self.createLine(['widget','HKLIN_FOBS_REF'],toggle=['REFERENCE','open',['HKLIN_FOBS_REF']])
        self.createLine(['widget','HKLIN_FC_REF'],toggle=['REFERENCE','open',['HKLIN_FC_REF']])
        self.createLine(['widget','HKLIN_FMAP_REF'],toggle=['REFERENCE','open',['HKLIN_FMAP_REF']])
        self.createLine(['widget','XYZIN_REF'],toggle=['REFERENCE','open',['XYZIN_REF']])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - -  Explicit reindex
        self.openSubFrame(toggle=['REFERENCE','open',['SPECIFY']])
        self.createLine(['subtitle',
                         "Note that reindexing a merged file is only valid within the same point group"])
        self.mtzSpgp = str(self.container.inputData.F_SIGF.fileContent.spaceGroup)
        #print "mtzSpgp", self.mtzSpgp
        line = self.createLine(['label','Space group from input data','label',self.mtzSpgp])
        # get label widget, note itemAt counts from 0
        self.inputSpaceGroupLabel = line.layout().itemAt(1).widget()

        self.createLine(['widget','CHOOSE_SPACEGROUP'])
        self.createLine(['widget','REINDEX_OPERATOR'])
        self.createLine(['widget','USE_REINDEX','label', 'use reindex operator'])

        self.createLine(['advice',
                         'You can give either a spacegroup, or a reindex operator (eg k,l,h), or both'])
        self.createLine(['advice',
                         'If only one of these is given, Pointless is usually able to generate the other automatically, but you should check'])
        self.createLine(['advice',
                         'The validity and consistency will be checked by Pointless: note warnings'])
        self.closeSubFrame()
        self.container.controlParameters.REINDEX_OPERATOR.dataChanged.connect(self.handleReindexOperator)
        # - - - - - - - - - - - - - - - - - - - - - - - - -   Explicit reindex

        # - - - - - - - - - - - - - - - - - - - - - -  Remove lattice centering
        self.openSubFrame(toggle=['REFERENCE','open',['LATTICE']])
        self.createLine(['subtitle',
                         "Remove centred lattice absences: BEWARE dangerous"])
        self.inputSpaceGroupLabel = line.layout().itemAt(1).widget()

        self.createLine(['label','Desired lattice centering type',
                         'widget','LATTICE_CENTERING'])
        self.createLine(['advice',
          'This option should ONLY be used if you are sure that the wrong cell was used in integration'])
        self.createLine(['advice',
         'Note that not all centred lattices are consisent with all Bravais lattices,'+\
                         ' check the result carefully'])
        self.createLine(['advice','Lattice type "P" is ignored'])

        self.closeSubFrame()
        # - - - - - - - - - - - - - - - - - - - - - -  Remove lattice centering

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    def isValid(self):
        invalidElements = super(pointless_reindexToMatch_gui, self).isValid()
        #Length of "solution elements" list can take any value if SOLIN is not set
        if not str(self.container.controlParameters.REFERENCE) == 'HKLIN_FOBS_REF':
            if self.container.inputData.HKLIN_FOBS_REF in invalidElements:
                invalidElements.remove(self.container.inputData.HKLIN_FOBS_REF)
                self.container.inputData.HKLIN_FOBS_REF.unSet()
        if not str(self.container.controlParameters.REFERENCE) == 'HKLIN_FC_REF':
            if self.container.inputData.HKLIN_FC_REF in invalidElements:
                invalidElements.remove(self.container.inputData.HKLIN_FC_REF)
                self.container.inputData.HKLIN_FC_REF.unSet()
        if not str(self.container.controlParameters.REFERENCE) == 'HKLIN_FMAP_REF':
            if self.container.inputData.HKLIN_FMAP_REF in invalidElements:
                invalidElements.remove(self.container.inputData.HKLIN_FMAP_REF)
                self.container.inputData.HKLIN_FMAP_REF.unSet()
        if not str(self.container.controlParameters.REFERENCE) == 'XYZIN_REF':
            if self.container.inputData.XYZIN_REF in invalidElements:
                invalidElements.remove(self.container.inputData.XYZIN_REF)
                self.container.inputData.XYZIN_REF.unSet()
        return invalidElements

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    @QtCore.Slot()
    def getSpacegroup(self):
        self.container.inputData.F_SIGF.loadFile()
        SG = self.container.inputData.F_SIGF.fileContent.spaceGroup
        #print "getSpacegroup", SG, type(SG)
        self.mtzSpgp = str(SG)
        self.inputSpaceGroupLabel.setText(self.mtzSpgp)
        
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    @QtCore.Slot()
    def handleReindexOperator(self):
        self.container.controlParameters.USE_REINDEX = True

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # from dummy job example
    def taskValidity(self):
        # based on dummy job example
        from ccp4i2.core import CCP4ErrorHandling
        rv = CCP4ErrorHandling.CErrorReport()
        if str(self.container.controlParameters.REFERENCE) == 'SPECIFY':
            # for explicit reindex, must define either Spacegroup or Reindex or both
            if (not self.container.controlParameters.CHOOSE_SPACEGROUP.isSet()) and \
               (not self.container.controlParameters.USE_REINDEX):
                # Fail
                ''' print('PRM .taskValidity fail')
                message = 'For explicit reindexing, you must set a space group or a reindex operator'
                rv.append(self.__class__,101,details=message,stack=False)
                '''
        return rv
