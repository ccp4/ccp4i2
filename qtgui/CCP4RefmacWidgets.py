
from PySide2 import QtGui, QtWidgets
from PySide2 import QtCore
from qtgui import CCP4Widgets
from core import CCP4RefmacData

# -----------------------------------------------------------------------------------

class CRefmacRigidGroupSegmentView(CCP4Widgets.CComplexLineWidget):
  MODEL_CLASS = CCP4RefmacData.CRefmacRigidGroupSegment

  def __init__(self,parent=None,model=None,qualifiers={}):
    qualis = {'vboxLayout':True,'iconButton':False}
    qualis.update(qualifiers)
    CCP4Widgets.CComplexLineWidget.__init__(self,parent,qualis)
    if model is not None: self.setModel(model)
    self.widgets = {}
    self.widgets['chain_id'] = CCP4Widgets.CStringView(self)
    self.widgets['residue_1'] = CCP4Widgets.CIntView(self)
    self.widgets['residue_2'] = CCP4Widgets.CIntView(self)
    line = QtWidgets.QHBoxLayout()
    line.addWidget(QtWidgets.QLabel('Chain',self))
    line.addWidget(self.widgets['chain_id'])
    line.addWidget(QtWidgets.QLabel('from residue',self))
    line.addWidget(self.widgets['residue_1'])
    line.addWidget(QtWidgets.QLabel('to residue',self))
    line.addWidget(self.widgets['residue_2'])
    self.layout().addLayout(line)

  def setModel(self,model):
    CCP4Widgets.CComplexLineWidget.setModel(self,model)
    if model is not None and self.editable:
      for item in ['chain_id','residue_1','residue_2']:
        model.get(item).dataChanged.connect(self.validate)


class CRefmacRigidGroupHeaderView(CCP4Widgets.CComplexLineWidget):
  def __init__(self,parent=None,model=None,qualifiers={}):
    CCP4Widgets.CComplexLineWidget.__init__(self,parent,qualifiers=qualifiers)
    self.widgets['rigid_group_id'] = CCP4Widgets.CStringView(self, qualifiers = { 'editable' : self.editable } )
    self.layout().addWidget(QtWidgets.QLabel('Add Domain Definition',self))
    self.layout().addWidget( self.widgets['rigid_group_id'])
    self.setModel(model)


class CRefmacRigidGroupListView(CCP4Widgets.CTreeView):
  MODEL_CLASS = CCP4RefmacData.CRefmacRigidGroupList
  def __init__(self,parent=None,model=None,qualifiers={}):
    qualis = { 'hierarchy' : [ {'name':'self','label':'Domain Definition','editorClass':CRefmacRigidGroupHeaderView, 'list':
                                   [{'name':'segmentList','label':'Residue Range','editorClass':CRefmacRigidGroupSegmentView }]
                             } ],
               'columnHeaders':['Chain','from residue','to residue'],
               }
    qualis.update(qualifiers)
    CCP4Widgets.CTreeView.__init__(self,parent,model=model,qualifiers=qualis)

# -----------------------------------------------------------------------------------

class CRefmacAnomalousAtomView(CCP4Widgets.CComplexLineWidget):

   MODEL_CLASS = CCP4RefmacData.CRefmacAnomalousAtom

   def __init__(self,parent=None,model=None,qualifiers={}):
      qualis = qualifiers
      CCP4Widgets.CComplexLineWidget.__init__(self,parent,qualis)

      self.widgets = {}
      if model is not None:
        self.widgets['atomType'] = CCP4Widgets.CStringView(self,model=model.atomType)
        self.widgets['Fp'] = CCP4Widgets.CFloatView(self,model=model.Fp)
        self.widgets['Fpp'] = CCP4Widgets.CFloatView(self,model=model.Fpp)
      else:
        self.widgets['atomType'] = CCP4Widgets.CStringView(self)
        self.widgets['Fp'] = CCP4Widgets.CFloatView(self)
        self.widgets['Fpp'] = CCP4Widgets.CFloatView(self)
        

      line = self.layout()
#     line.removeWidget(self.iconButton)
      line.addWidget(QtWidgets.QLabel('Anomalous element   ',self))
      line.addWidget(self.widgets['atomType'])
      line.addWidget(QtWidgets.QLabel("       f'   ",self))
      line.addWidget(self.widgets['Fp'])
      line.addWidget(QtWidgets.QLabel("       f''   ",self))
      line.addWidget(self.widgets['Fpp'])

      self.setModel(model)

