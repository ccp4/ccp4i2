"""
    pipelines/SubstituteLigand/SubstituteLigand_gui.py: CCP4 GUI Project
    Copyright (C) 2015 Newcastle University
    """

"""
    Maritn Noble
    """

from PySide2 import QtGui, QtWidgets,QtCore
from qtgui.CCP4TaskWidget import CTaskWidget
from qtgui import CCP4Widgets

class SubstituteLigand_gui(CTaskWidget):
    
    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'SubstituteLigand'
    TASKVERSION = 0.0
    TASKMODULE=['ligands','bigpipes']
    TASKTITLE='Automated solution of isomorphous ligand complex'
    SHORTTASKTITLE='Isomorphous ligand solution'
    WHATNEXT = ['prosmart_refmac','coot_rebuild']
    DESCRIPTION = '''A ligand workflow, starting from merged or unmerged reflections, SMILES, and an isomorphous parent structure'''
    RANK=1

    # -------------------------------------------------------------
    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)
    
    # -------------------------------------------------------------
    def drawContents(self):
        folder = self.openFolder(folderFunction='inputData',title='Input Data')
        self.createLine( [ 'subtitle', 'Ligand geometry' ])
        self.createLine(['label','Format in which geometry will be specified: ','stretch','widget','LIGANDAS'])
        self.container.controlParameters.OBSAS.dataChanged.connect(self.LIGANDASChanged)
        
        self.openSubFrame(frame=True,toggle=['LIGANDAS','open',['MOL']])
        self.createLine(['widget','MOLIN'])
        self.closeSubFrame()

        self.openSubFrame(frame=True,toggle=['LIGANDAS','open',['DICT']])
        self.createLine(['widget','-browseDb','True','DICTIN'])
        self.closeSubFrame()

        self.openSubFrame(frame=True,toggle=['LIGANDAS','open',['SMILES']])
        self.createLine ( [ 'widget', '-guiMode', 'multiLine', 'SMILESIN' ] )
        self.closeSubFrame()

        self.openSubFrame(frame=True)
        self.createLine(['label','For rigid body refinement use','widget','PIPELINE'])
        self.closeSubFrame()
        
        self.createLine( [ 'subtitle', 'Starting PDB <i>In the appropriate spacegroup</i>.' ])
        self.createLine(['widget','-title','Parent structure','XYZIN'])
        self.getWidget('XYZIN').showAtomSelection()
        
        
        self.createLine( [ 'subtitle', 'Reflection data' ])
        self.createLine(['label','Format in which reflections will be specified:  ','stretch','widget','OBSAS'])
        self.container.controlParameters.OBSAS.dataChanged.connect(self.OBSASChanged)
        
        self.openSubFrame(frame=True,toggle=['OBSAS','open',['UNMERGED']])
        self.createLine(['widget','UNMERGEDFILES'])
        self.closeSubFrame()

        self.openSubFrame(frame=True,toggle=['OBSAS','open',['MERGED']])
        self.createLine(['widget','-browseDb','True','F_SIGF_IN'])
        self.closeSubFrame()

        self.createLine( [ 'subtitle', 'Free-R flags' ])
        self.openSubFrame(frame=True)
        self.createLine(['widget','-browseDb','True','FREERFLAG_IN'])
        self.closeSubFrame()

        self.OBSASChanged()

    @QtCore.Slot()
    def OBSASChanged(self):
        if self.container.controlParameters.OBSAS.__str__() == 'MERGED':
            self.container.inputData.F_SIGF_IN.setQualifiers({'allowUndefined' : False } )
            self.container.inputData.FREERFLAG_IN.setQualifiers({'allowUndefined' : True } )
        else:
            self.container.inputData.F_SIGF_IN.setQualifiers({'allowUndefined' : True } )
            self.container.inputData.FREERFLAG_IN.setQualifiers({'allowUndefined' : True } )

    @QtCore.Slot()
    def LIGANDASChanged(self):
        if self.container.controlParameters.LIGANDAS.__str__() == 'MOL':
            self.container.inputData.MOLIN.setQualifiers({'allowUndefined' : False } )
            self.container.inputData.DICTIN.setQualifiers({'allowUndefined' : True } )
            self.container.inputData.SMILESIN.setQualifiers({'minLength' : 0 } )
        elif self.container.controlParameters.LIGANDAS.__str__() == 'DICT':
            self.container.inputData.MOLIN.setQualifiers({'allowUndefined' : True } )
            self.container.inputData.DICTIN.setQualifiers({'allowUndefined' : False } )
            self.container.inputData.SMILESIN.setQualifiers({'minLength' : 0 } )
        if self.container.controlParameters.LIGANDAS.__str__() == 'SMILES':
            self.container.inputData.MOLIN.setQualifiers({'allowUndefined' : True } )
            self.container.inputData.DICTIN.setQualifiers({'allowUndefined' : True } )
            self.container.inputData.SMILESIN.setQualifiers({'minLength' : 1 } )

    def isValid(self):
        invalidElements = CTaskWidget.isValid(self)
        copiedInvalidElements = []
        for invalidElement in invalidElements:
            doCopy = True
            if self.container.controlParameters.OBSAS.__str__() == 'MERGED':
                for toRemove in  [self.container.inputData.UNMERGEDFILES, self.container.inputData.UNMERGEDFILES[0], self.container.inputData.UNMERGEDFILES[0].file, self.container.inputData.UNMERGEDFILES[0].crystalName,self.container.inputData.UNMERGEDFILES[0].dataset,self.container.inputData.UNMERGEDFILES[0].cell,self.container.inputData.UNMERGEDFILES[0].cell.a,self.container.inputData.UNMERGEDFILES[0].cell.b,self.container.inputData.UNMERGEDFILES[0].cell.c]:
                    if invalidElement is toRemove:
                        doCopy = False
            if doCopy:
                copiedInvalidElements.append(invalidElement)
        return copiedInvalidElements
