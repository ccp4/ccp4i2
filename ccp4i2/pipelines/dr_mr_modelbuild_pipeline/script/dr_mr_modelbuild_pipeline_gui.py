"""
Copyright (C) 2011 University of York
Stuart McNicholas 2020-2021 adapted from Andrey Lebedev September 2011 - molrep_mr gui
"""

from PySide2 import QtCore

from qtgui import CCP4TaskWidget


class Cdr_mr_modelbuild_pipeline(CCP4TaskWidget.CTaskWidget):

  TASKNAME = 'dr_mr_modelbuild_pipeline'
  TASKVERSION = 0.0
  TASKMODULE= ['model_building','bigpipes']
  WHATNEXT = ['prosmart_refmac','coot_rebuild']
  TASKTITLE = 'Data Reduction, MR, Model build pipeline'
  SHORTTASKTITLE = 'Data Rection, MR, model building'
  DESCRIPTION='Data Reduction, MR, Model build pipeline'
  MGDISPLAYFILES = ['XYZOUT']  
  RANK=1
  
  def __init__(self,parent):
    CCP4TaskWidget.CTaskWidget.__init__(self,parent)

  def drawContents(self):

    self.setProgramHelpFile('molrep_mr')

#-  --------------------          --------------------          --------------------

    folder = self.openFolder(folderFunction='inputData',title='Input Data and Protocol')


    self.openSubFrame()
    self.createLine( [ 'subtitle', 'Input data' ])
    self.createLine(['label',"Input data type (run import merged, use merged from 'import merged' job or unmerged: ",'stretch','widget','MERGED_OR_UNMERGED'])
    self.container.controlParameters.MERGED_OR_UNMERGED.dataChanged.connect(self.DataTypeChanged)
    self.createLine( ['widget', '-title','Select unmerged data files', 'UNMERGEDFILES'], toggle=['MERGED_OR_UNMERGED','open',[ 'UNMERGED' ]]  )
    self.createLine( [ 'widget', '-browseDb', True, 'F_SIGF_IN'], toggle=['MERGED_OR_UNMERGED','open',[ 'MERGED' ]] )
    self.createLine( [ 'widget', '-browseDb', True, 'FREER_IN'], toggle=['MERGED_OR_UNMERGED','open',[ 'MERGED' ]] )
    self.createLine( [ 'widget', '-browseDb', True, 'HKLIN'], toggle=['MERGED_OR_UNMERGED','open',[ 'MERGED_F' ]] )
    
    self.closeSubFrame()

    self.openSubFrame()
    indent = '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'
    self.createLine( [ 'subtitle', 'Search model' ] )
    self.createLine( [ 'widget', 'XYZINORMRBUMP' ] )
    self.createLine( [ 'widget', 'XYZIN' ] , toggle=['XYZINORMRBUMP','open',['XYZINPUT']])
    self.openSubFrame(frame=[False],toggle=['XYZINORMRBUMP','open',['MRBUMP']])
    self.createLine( [ 'widget', 'SEARCH_PDB', 'label', 'Search PDB for for possible MR search models' ])
    self.createLine( [ 'advice', indent+'Non-redundancy level for homologue search:'], toggle=['SEARCH_PDB','open',[ True ]])
    self.createLine( [ 'label', indent, 'widget', 'REDUNDANCYLEVEL' ], toggle=['SEARCH_PDB','open',[ True ]])
    self.createLine( [ 'widget', 'SEARCH_AFDB', 'label', 'Search EBI-AFDB for possible MR search models' ])
    self.createLine( [ 'label', indent, 'advice', indent+'EBI-AFDB pLDDT residue score cut-off:' ], toggle=['SEARCH_AFDB','open',[ True ]])
    self.createLine( [ 'label', indent, 'widget', 'AFDBLEVEL' ], toggle=['SEARCH_AFDB','open',[ True ]])
    self.createLine( [ 'advice', 'Maximum no. of search models to create:'] , toggle=['XYZINORMRBUMP','open',['MRBUMP']] )
    self.createLine( [ 'widget', 'MRMAX' ] , toggle=['XYZINORMRBUMP','open',['MRBUMP']])
    self.closeSubFrame()
    self.getWidget('XYZIN').showAtomSelection()
    self.createLine( [ 'subtitle', 'Sequence of target model' ] )
    self.createLine( [ 'widget', 'ASUIN' ] )
    self.createLine( [ 'label', 'The number of monomers to search for', 'widget', 'NMON' ] )
    print('molrep_mr setWhatsThis',self.getWidget('NMON'))
    self.getWidget('NMON').setWhatsThis("The number of monomers in the asymmetric unit - recommended that you choose 'auto' and let the program decide.")
    self.closeSubFrame()

    """
    self.openSubFrame( toggle = [ 'PERFORM','open', [ 'pat' ]] )
    self.createLine( [ 'subtitle', 'Fixed Model' ] )
    self.createLine( [ 'widget', 'XYZIN_FIX' ] )
    self.closeSubFrame()
    """

    self.openSubFrame()
    self.createLine( [ 'subtitle', 'Options' ] )
    self.createLine(['widget', 'AUTOCUTOFF', 'label', 'Run Aimless twice, first to find resolution limit'])
    self.createLine(['widget', 'RUNACORN', 'label', 'Run phase refinement with acorn before model building'])
    self.createLine(['label', 'Run', 'widget', 'REFMAC_NCYC', 'label', 'cycles of restrained refinement after MR'])
    self.createLine(['label', 'Run', 'widget', 'BUCC_NCYC', 'label', 'model building pipeline iterations'])
    self.closeSubFrame()

    self.createLine( [ 'subtitle', 'Ligand geometry' ])
    self.createLine(['label','Format in which geometry will be specified: ','stretch','widget','LIGANDAS'])
    self.container.controlParameters.LIGANDAS.dataChanged.connect(self.LIGANDASChanged)
        
    self.openSubFrame(frame=True,toggle=['LIGANDAS','open',['MOL']])
    self.createLine(['widget','MOLIN'])
    self.closeSubFrame()

    self.openSubFrame(frame=True,toggle=['LIGANDAS','open',['DICT']])
    self.createLine(['widget','-browseDb','True','DICTIN'])
    self.closeSubFrame()

    self.openSubFrame(frame=True,toggle=['LIGANDAS','open',['SMILES']])
    self.createLine ( [ 'widget', '-guiMode', 'multiLine', 'SMILESIN' ] )
    self.closeSubFrame()

    advanced_folder = self.openFolder(folderFunction='inputData',title='Advanced options', drawFolder=self.drawAdvanced)

  def drawAdvanced(self):
    self.openSubFrame()
    self.createLine( [ 'subtitle', 'Advanced model building pipeline options' ] )
    self.createLine( ['label', 'Use', 'widget', 'BUCCANEER_OR_MODELCRAFT', 'label', 'model building pipleline (Buccaneer or Modelcraft)'] )
    self.createLine( [ 'advice', 'The modelcraft pipeline is now the default option and recommended in most cases'] )
    self.closeSubFrame()

  @QtCore.Slot()
  def DataTypeChanged(self):
      if self.container.controlParameters.MERGED_OR_UNMERGED.__str__() == 'MERGED':
            self.container.inputData.UNMERGEDFILES.setQualifiers({'allowUndefined' : True } )
            self.container.inputData.UNMERGEDFILES.setQualifiers({'mustExist' : False } )
            self.container.inputData.FREER_IN.setQualifiers({'allowUndefined' : False } )
            self.container.inputData.FREER_IN.setQualifiers({'mustExist' : True } )
            self.container.inputData.F_SIGF_IN.setQualifiers({'allowUndefined' : False } )
            self.container.inputData.F_SIGF_IN.setQualifiers({'mustExist' : True } )
            self.container.inputData.HKLIN.setQualifiers({'allowUndefined' : True } )
            self.container.inputData.HKLIN.setQualifiers({'mustExist' : False } )
      elif self.container.controlParameters.MERGED_OR_UNMERGED.__str__() == 'MERGED_F':
            self.container.inputData.UNMERGEDFILES.setQualifiers({'allowUndefined' : True } )
            self.container.inputData.UNMERGEDFILES.setQualifiers({'mustExist' : False } )
            self.container.inputData.FREER_IN.setQualifiers({'allowUndefined' : True } )
            self.container.inputData.FREER_IN.setQualifiers({'mustExist' : False } )
            self.container.inputData.F_SIGF_IN.setQualifiers({'allowUndefined' : True } )
            self.container.inputData.F_SIGF_IN.setQualifiers({'mustExist' : False } )
            self.container.inputData.HKLIN.setQualifiers({'allowUndefined' : False } )
            self.container.inputData.HKLIN.setQualifiers({'mustExist' : True } )
      else:
            self.container.inputData.UNMERGEDFILES.setQualifiers({'allowUndefined' : False } )
            self.container.inputData.UNMERGEDFILES.setQualifiers({'mustExist' : True } )
            self.container.inputData.FREER_IN.setQualifiers({'allowUndefined' : True } )
            self.container.inputData.FREER_IN.setQualifiers({'mustExist' : False } )
            self.container.inputData.F_SIGF_IN.setQualifiers({'allowUndefined' : True } )
            self.container.inputData.F_SIGF_IN.setQualifiers({'mustExist' : False } )
            self.container.inputData.HKLIN.setQualifiers({'allowUndefined' : True } )
            self.container.inputData.HKLIN.setQualifiers({'mustExist' : False } )
      unmerged_widget = self.getWidget("UNMERGEDFILES")
      freer_widget = self.getWidget("FREER_IN")
      f_sigf_widget = self.getWidget("F_SIGF_IN")
      hklin_widget = self.getWidget("HKLIN")
      unmerged_widget.validate()
      freer_widget.validate()
      f_sigf_widget.validate()
      hklin_widget.validate()

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
        self.DataTypeChanged()
        invalidElements = super(Cdr_mr_modelbuild_pipeline, self).isValid()
        print("##################################################")
        print("original",invalidElements)
        print(self.container.controlParameters.MERGED_OR_UNMERGED.__str__())
        if self.container.inputData.XYZINORMRBUMP.__str__() == 'MRBUMP' or self.container.inputData.XYZINORMRBUMP.__str__() == 'MRPARSE':
            invalidElements = [x for x in invalidElements if x.objectName() != "XYZIN"]
        if self.container.controlParameters.MERGED_OR_UNMERGED.__str__() == 'MERGED':
            invalidElements = [x for x in invalidElements if x.objectName() != "HKLIN"]
            invalidElements = [x for x in invalidElements if x.objectName() != "UNMERGEDFILES"]
        elif self.container.controlParameters.MERGED_OR_UNMERGED.__str__() == 'MERGED_F':
            invalidElements = [x for x in invalidElements if x.objectName() != "FREER_IN"]
            invalidElements = [x for x in invalidElements if x.objectName() != "F_SIGF_IN"]
            invalidElements = [x for x in invalidElements if x.objectName() != "UNMERGEDFILES"]
        else:
            invalidElements = [x for x in invalidElements if x.objectName() != "FREER_IN"]
            invalidElements = [x for x in invalidElements if x.objectName() != "F_SIGF_IN"]
            invalidElements = [x for x in invalidElements if x.objectName() != "HKLIN"]
        print("final",invalidElements)
        print("##################################################")
        return invalidElements
