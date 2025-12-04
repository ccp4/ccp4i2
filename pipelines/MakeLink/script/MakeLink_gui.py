"""
    MakeLink_gui.py: CCP4 GUI Project
    
    This library is free software: you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    version 3, modified in accordance with the provisions of the
    license to address the requirements of UK law.
    
    You should have received a copy of the modified GNU Lesser General
    Public License along with this library.  If not, copies may be
    downloaded from http://www.ccp4.ac.uk/ccp4license.php
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
    """

from qtgui.CCP4TaskWidget import CTaskWidget
from baselayer import QtCore
import sys, os
from core import CCP4Utils
import gemmi

#-------------------------------------------------------------------
class MakeLink_gui(CTaskWidget):
    #-------------------------------------------------------------------

    # Subclass CTaskWidget to give specific task window
    TASKNAME = 'MakeLink' # this has to match the pluginName given in the corresponding .def.xml
    TASKVERSION = 0.1
    TASKMODULE = [ 'ligands' ] #Section in the task list where this task will be listed e.g. 'refinement','model_building' for full list see MODULE_ORDER in core/CCP4TaskManager.py
    SHORTTASKTITLE='Make Covalent Link - AceDRG'
    TASKTITLE='Make Covalent Link - AceDRG'
    DESCRIPTION = '''Generate a link dictionary to describe a covalent bond between two monomers, allowing modification of the linked monomers'''
    MGDISPLAYFILES = ['XYZIN']
    WHATNEXT = ['coot_rebuild','prosmart_refmac']

    def __init__(self,parent):
        CTaskWidget.__init__(self,parent)

    def ToggleDict1(self):
       if str(self.container.inputData.MON_1_TYPE) == 'CIF':
          self.container.inputData.DICT_1.setQualifiers({'allowUndefined':False})
          self.container.inputData.RES_NAME_1_CIF.setQualifiers({'allowUndefined':False})
          self.container.inputData.RES_NAME_1_TLC.setQualifiers({'allowUndefined':True})
          self.container.inputData.ATOM_NAME_1_CIF.setQualifiers({'allowUndefined':False})
          self.container.inputData.ATOM_NAME_1_TLC.setQualifiers({'allowUndefined':True})
          self.validate()
          return True
       self.container.inputData.DICT_1.setQualifiers({'allowUndefined':True})
       self.container.inputData.RES_NAME_1_CIF.setQualifiers({'allowUndefined':True})
       self.container.inputData.RES_NAME_1_TLC.setQualifiers({'allowUndefined':False})
       self.container.inputData.ATOM_NAME_1_CIF.setQualifiers({'allowUndefined':True})
       self.container.inputData.ATOM_NAME_1_TLC.setQualifiers({'allowUndefined':False})
       self.validate()
       return False

    def ToggleDict2(self):
       if str(self.container.inputData.MON_2_TYPE) == 'CIF':
          self.container.inputData.DICT_2.setQualifiers({'allowUndefined':False})
          self.container.inputData.RES_NAME_2_CIF.setQualifiers({'allowUndefined':False})
          self.container.inputData.RES_NAME_2_TLC.setQualifiers({'allowUndefined':True})
          self.container.inputData.ATOM_NAME_2_CIF.setQualifiers({'allowUndefined':False})
          self.container.inputData.ATOM_NAME_2_TLC.setQualifiers({'allowUndefined':True})
          self.validate()
          return True
       self.container.inputData.DICT_2.setQualifiers({'allowUndefined':True})
       self.container.inputData.RES_NAME_2_CIF.setQualifiers({'allowUndefined':True})
       self.container.inputData.RES_NAME_2_TLC.setQualifiers({'allowUndefined':False})
       self.container.inputData.ATOM_NAME_2_CIF.setQualifiers({'allowUndefined':True})
       self.container.inputData.ATOM_NAME_2_TLC.setQualifiers({'allowUndefined':False})
       self.validate()
       return False

    def ToggleTLC1(self):
       if str(self.container.inputData.MON_1_TYPE) == 'TLC':
          return True
       return False

    def ToggleTLC2(self):
       if str(self.container.inputData.MON_2_TYPE) == 'TLC':
          return True
       return False

    def ToggleCIF1(self):
       if str(self.container.inputData.MON_1_TYPE) == 'CIF':
          return True
       return False

    def ToggleCIF2(self):
       if str(self.container.inputData.MON_2_TYPE) == 'CIF':
          return True
       return False

    def ToggleDelete1(self):
       if self.container.inputData.TOGGLE_DELETE_1:
          self.container.inputData.DELETE_1_LIST.setQualifiers({'allowUndefined':False})
          self.DeleteListChanged1()
          self.validate()
          return True
       self.container.inputData.DELETE_1_LIST.setQualifiers({'allowUndefined':True})
       self.DeleteListChanged1()
       self.validate()
       return False

    def ToggleDelete2(self):
       if self.container.inputData.TOGGLE_DELETE_2:
          self.container.inputData.DELETE_2_LIST.setQualifiers({'allowUndefined':False})
          self.DeleteListChanged2()
          self.validate()
          return True
       self.container.inputData.DELETE_2_LIST.setQualifiers({'allowUndefined':True})
       self.DeleteListChanged2()
       self.validate()
       return False

    def ToggleChange1(self):
       if self.container.inputData.TOGGLE_CHANGE_1:
          self.container.inputData.CHANGE_BOND_1_LIST.setQualifiers({'allowUndefined':False})
          self.container.inputData.CHANGE_BOND_1_TYPE.setQualifiers({'allowUndefined':False})
          self.validate()
          return True
       self.container.inputData.CHANGE_BOND_1_LIST.setQualifiers({'allowUndefined':True})
       self.container.inputData.CHANGE_BOND_1_TYPE.setQualifiers({'allowUndefined':True})
       self.validate()
       return False

    def ToggleChange2(self):
       if self.container.inputData.TOGGLE_CHANGE_2:
          self.container.inputData.CHANGE_BOND_2_LIST.setQualifiers({'allowUndefined':False})
          self.container.inputData.CHANGE_BOND_2_TYPE.setQualifiers({'allowUndefined':False})
          self.validate()
          return True
       self.container.inputData.CHANGE_BOND_2_LIST.setQualifiers({'allowUndefined':True})
       self.container.inputData.CHANGE_BOND_2_TYPE.setQualifiers({'allowUndefined':True})
       self.validate()
       return False

    def ToggleCharge1(self):
       if self.container.inputData.TOGGLE_CHARGE_1:
          self.container.inputData.CHARGE_1_LIST.setQualifiers({'allowUndefined':False})
          self.container.inputData.CHARGE_1_VALUE.setQualifiers({'allowUndefined':False})
          self.ChargeListChanged1()
          self.validate()
          return True
       self.container.inputData.CHARGE_1_LIST.setQualifiers({'allowUndefined':True})
       self.container.inputData.CHARGE_1_VALUE.setQualifiers({'allowUndefined':True})
       self.ChargeListChanged1()
       self.validate()
       return False

    def ToggleCharge2(self):
       if self.container.inputData.TOGGLE_CHARGE_2:
          self.container.inputData.CHARGE_2_LIST.setQualifiers({'allowUndefined':False})
          self.container.inputData.CHARGE_2_VALUE.setQualifiers({'allowUndefined':False})
          self.ChargeListChanged2()
          self.validate()
          return True
       self.container.inputData.CHARGE_2_LIST.setQualifiers({'allowUndefined':True})
       self.container.inputData.CHARGE_2_VALUE.setQualifiers({'allowUndefined':True})
       self.ChargeListChanged2()
       self.validate()
       return False

#    def ToggleAutoLink(self):
#       if self.container.controlParameters.TOGGLE_AUTOLINK:
#          return True
#       return False

    def ToggleLinkMode(self):
       if self.container.controlParameters.TOGGLE_LINK:
          if self.container.inputData.XYZIN.isSet():
             return True
       return False

    def ToggleAutoLinkMode(self):
       if self.container.controlParameters.TOGGLE_LINK:
          if self.container.inputData.XYZIN.isSet():
             if str(self.container.controlParameters.LINK_MODE) == 'AUTO':
                return True
       return False

    def ToggleManualLinkMode(self):
       if self.container.controlParameters.TOGGLE_LINK:
          if self.container.inputData.XYZIN.isSet():
             if str(self.container.controlParameters.LINK_MODE) == 'MANUAL':
                return True
       return False
    
    def getAtomsFromCif(self,cif_path,comp_id):
       atom_list = []
       try:
         if not os.path.exists(cif_path):
            print("Cannot find file: "+cif_path)
            return atom_list
         input_cif = gemmi.cif.read_file(cif_path)
         block = input_cif.find_block("comp_{}".format(comp_id))
         if block:
            for (atom_id,type) in block.find(['_chem_comp_atom.atom_id','_chem_comp_atom.type_symbol']):
               if type != 'H':
                  atom_list.append(str(atom_id))
       except Exception as e:
         print("Error reading atoms from CIF: %s" % e)
       if len(atom_list)>0:
         atom_list = [atom.strip('\"') for atom in atom_list]
       return atom_list
    
    def getBondsFromCif(self,cif_path,comp_id):
       bond_list = []
       try:
         if not os.path.exists(cif_path):
            print("Cannot find file: "+cif_path)
            return bond_list
         input_cif = gemmi.cif.read_file(cif_path)
         block = input_cif.find_block("comp_{}".format(comp_id))
         if block:
            atom_list = []
            for (atom_id,type) in block.find(['_chem_comp_atom.atom_id','_chem_comp_atom.type_symbol']):
               if type != 'H':
                  atom_list.append(str(atom_id))
            for (atom1,atom2,type) in block.find(['_chem_comp_bond.atom_id_1','_chem_comp_bond.atom_id_2','_chem_comp_bond.type']):
               if all(x in atom_list for x in [atom1,atom2]):
                  bond_list.append([str(atom1).strip('\"'),str(atom2).strip('\"'),str(type).lower()])
       except Exception as e:
         print("Error reading bonds from CIF: %s" % e)
       return bond_list

    @QtCore.Slot()
    def updateAtomList1TLC(self):
       # Called when residue name is changed
       if self.container.inputData.MON_1_TYPE == 'TLC':
          # Get info from monomer library
          if self.container.inputData.RES_NAME_1_TLC.isSet():
             comp_id = self.container.inputData.RES_NAME_1_TLC.__str__()
             if len(comp_id)>0:
                cif_path = os.path.join(CCP4Utils.getCCP4Dir().__str__(), 'lib', 'data', 'monomers', comp_id[0].lower(), comp_id+'.cif')
                if len(cif_path)>0:
                   atom_list = self.getAtomsFromCif(cif_path,comp_id)
                   if len(atom_list)>0:
                      self.container.inputData.ATOM_NAME_1_TLC.setQualifiers({'enumerators':atom_list,'menuText':atom_list})
                      if self.container.inputData.ATOM_NAME_1.__str__() in self.container.inputData.ATOM_NAME_1_TLC.qualifiers('enumerators')[:]:
                         # If job is cloned or reloaded, reselect correct element of combobox
                         self.container.inputData.ATOM_NAME_1_TLC.set(self.container.inputData.ATOM_NAME_1.__str__())
                      else:
                         self.container.inputData.ATOM_NAME_1_TLC.unSet()
                      self.getWidget('ATOM_NAME_1_TLC').populateComboBox(self.container.inputData.ATOM_NAME_1_TLC)
                      self.getWidget('ATOM_NAME_1_TLC').updateViewFromModel()
                      self.updateChangeList1()
                      return
          self.container.inputData.ATOM_NAME_1_TLC.setQualifiers({'enumerators':[],'menuText':[]})
          self.container.inputData.ATOM_NAME_1_TLC.unSet()
          self.getWidget('ATOM_NAME_1_TLC').populateComboBox(self.container.inputData.ATOM_NAME_1_TLC)
          self.getWidget('ATOM_NAME_1_TLC').updateViewFromModel()
          self.updateChangeList1()
       return

    @QtCore.Slot()
    def updateAtomList2TLC(self):
       # Called when residue name is changed
       if self.container.inputData.MON_2_TYPE == 'TLC':
          # Get info from monomer library
          if self.container.inputData.RES_NAME_2_TLC.isSet():
             comp_id = self.container.inputData.RES_NAME_2_TLC.__str__()
             if len(comp_id)>0:
                cif_path = os.path.join(CCP4Utils.getCCP4Dir().__str__(), 'lib', 'data', 'monomers', comp_id[0].lower(), comp_id+'.cif')
                if len(cif_path)>0:
                   atom_list = self.getAtomsFromCif(cif_path,comp_id)
                   if len(atom_list)>0:
                      self.container.inputData.ATOM_NAME_2_TLC.setQualifiers({'enumerators':atom_list,'menuText':atom_list})
                      if self.container.inputData.ATOM_NAME_2.__str__() in self.container.inputData.ATOM_NAME_2_TLC.qualifiers('enumerators')[:]:
                         # If job is cloned or reloaded, reselect correct element of combobox
                         self.container.inputData.ATOM_NAME_2_TLC.set(self.container.inputData.ATOM_NAME_2.__str__())
                      else:
                         self.container.inputData.ATOM_NAME_2_TLC.unSet()
                      self.getWidget('ATOM_NAME_2_TLC').populateComboBox(self.container.inputData.ATOM_NAME_2_TLC)
                      self.getWidget('ATOM_NAME_2_TLC').updateViewFromModel()
                      self.updateChangeList2()
                      return
          self.container.inputData.ATOM_NAME_2_TLC.setQualifiers({'enumerators':[],'menuText':[]})
          self.container.inputData.ATOM_NAME_2_TLC.unSet()
          self.getWidget('ATOM_NAME_2_TLC').populateComboBox(self.container.inputData.ATOM_NAME_2_TLC)
          self.getWidget('ATOM_NAME_2_TLC').updateViewFromModel()
          self.updateChangeList2()
       return

    @QtCore.Slot()
    def updateAtomList1CIF(self):
       # Called when residue name is changed
       if self.container.inputData.MON_1_TYPE == 'CIF':
          # Get info from CIF file
          if self.container.inputData.DICT_1.isSet() and self.container.inputData.RES_NAME_1_CIF.isSet():
             comp_id = self.container.inputData.RES_NAME_1_CIF.__str__()
             cif_path = self.container.inputData.DICT_1.fullPath.__str__()
             if len(comp_id)>0 and len(cif_path)>0:
                atom_list = self.getAtomsFromCif(cif_path,comp_id)
                if len(atom_list)>0:
                   self.container.inputData.ATOM_NAME_1_CIF.setQualifiers({'enumerators':atom_list,'menuText':atom_list})
                   if self.container.inputData.ATOM_NAME_1.__str__() in self.container.inputData.ATOM_NAME_1_CIF.qualifiers('enumerators')[:]:
                      # If job is cloned or reloaded, reselect correct element of combobox
                      self.container.inputData.ATOM_NAME_1_CIF.set(self.container.inputData.ATOM_NAME_1.__str__())
                   else:
                      self.container.inputData.ATOM_NAME_1_CIF.unSet()
                   self.getWidget('ATOM_NAME_1_CIF').populateComboBox(self.container.inputData.ATOM_NAME_1_CIF)
                   self.getWidget('ATOM_NAME_1_CIF').updateViewFromModel()
                   self.updateChangeList1()
                   return
          self.container.inputData.ATOM_NAME_1_CIF.setQualifiers({'enumerators':[],'menuText':[]})
          self.container.inputData.ATOM_NAME_1_CIF.unSet()
          self.getWidget('ATOM_NAME_1_CIF').populateComboBox(self.container.inputData.ATOM_NAME_1_CIF)
          self.getWidget('ATOM_NAME_1_CIF').updateViewFromModel()
          self.updateChangeList1()
       return
    
    @QtCore.Slot()
    def updateAtomList2CIF(self):
       # Called when residue name is changed
       if self.container.inputData.MON_2_TYPE == 'CIF':
          # Get info from CIF file
          if self.container.inputData.DICT_2.isSet() and self.container.inputData.RES_NAME_2_CIF.isSet():
             comp_id = self.container.inputData.RES_NAME_2_CIF.__str__()
             cif_path = self.container.inputData.DICT_2.fullPath.__str__()
             if len(comp_id)>0 and len(cif_path)>0:
                atom_list = self.getAtomsFromCif(cif_path,comp_id)
                if len(atom_list)>0:
                   self.container.inputData.ATOM_NAME_2_CIF.setQualifiers({'enumerators':atom_list,'menuText':atom_list})
                   if self.container.inputData.ATOM_NAME_2.__str__() in self.container.inputData.ATOM_NAME_2_CIF.qualifiers('enumerators')[:]:
                      # If job is cloned or reloaded, reselect correct element of combobox
                      self.container.inputData.ATOM_NAME_2_CIF.set(self.container.inputData.ATOM_NAME_2.__str__())
                   else:
                      self.container.inputData.ATOM_NAME_2_CIF.unSet()
                   self.getWidget('ATOM_NAME_2_CIF').populateComboBox(self.container.inputData.ATOM_NAME_2_CIF)
                   self.getWidget('ATOM_NAME_2_CIF').updateViewFromModel()
                   self.updateChangeList2()
                   return
          self.container.inputData.ATOM_NAME_2_CIF.setQualifiers({'enumerators':[],'menuText':[]})
          self.container.inputData.ATOM_NAME_2_CIF.unSet()
          self.getWidget('ATOM_NAME_2_CIF').populateComboBox(self.container.inputData.ATOM_NAME_2_CIF)
          self.getWidget('ATOM_NAME_2_CIF').updateViewFromModel()
          self.updateChangeList2()
       return
    
    def updateDeleteList1(self,atom_list):
       # Called after the atom list changes
       self.container.inputData.DELETE_1_LIST.setQualifiers({'enumerators':atom_list,'menuText':atom_list})
       if self.container.inputData.DELETE_1.__str__() in self.container.inputData.DELETE_1_LIST.qualifiers('enumerators')[:]:
          self.container.inputData.DELETE_1_LIST.set(self.container.inputData.DELETE_1.__str__())
       else:
          self.container.inputData.DELETE_1_LIST.unSet()
          self.container.inputData.DELETE_1.unSet()
       self.getWidget('DELETE_1_LIST').populateComboBox(self.container.inputData.DELETE_1_LIST)
       self.getWidget('DELETE_1_LIST').updateViewFromModel()
       self.validate()
       return

    def updateDeleteList2(self,atom_list):
       # Called after the atom list changes
       self.container.inputData.DELETE_2_LIST.setQualifiers({'enumerators':atom_list,'menuText':atom_list})
       if self.container.inputData.DELETE_2.__str__() in self.container.inputData.DELETE_2_LIST.qualifiers('enumerators')[:]:
          self.container.inputData.DELETE_2_LIST.set(self.container.inputData.DELETE_2.__str__())
       else:
          self.container.inputData.DELETE_2_LIST.unSet()
          self.container.inputData.DELETE_2.unSet()
       self.getWidget('DELETE_2_LIST').populateComboBox(self.container.inputData.DELETE_2_LIST)
       self.getWidget('DELETE_2_LIST').updateViewFromModel()
       self.validate()
       return

    @QtCore.Slot()
    def DeleteListChanged1(self):
       # Called when the delete list is changed
       if self.container.inputData.DELETE_1_LIST.isSet():
          self.container.inputData.DELETE_1.set(self.container.inputData.DELETE_1_LIST.__str__())
       else:
          if self.container.inputData.DELETE_1.isSet():
             if self.container.inputData.DELETE_1.__str__() in self.container.inputData.DELETE_1_LIST.qualifiers('enumerators')[:]:
                self.container.inputData.DELETE_1_LIST.set(self.container.inputData.DELETE_1.__str__())
                self.getWidget('DELETE_1_LIST').populateComboBox(self.container.inputData.DELETE_1_LIST)
                self.getWidget('DELETE_1_LIST').updateViewFromModel()
                self.validate()
             else:
                self.container.inputData.DELETE_1.unSet()
       self.updateChangeList1()
       self.reloadType1()
       return

    @QtCore.Slot()
    def DeleteListChanged2(self):
       # Called when the delete list is changed
       if self.container.inputData.DELETE_2_LIST.isSet():
          self.container.inputData.DELETE_2.set(self.container.inputData.DELETE_2_LIST.__str__())
       else:
          if self.container.inputData.DELETE_2.isSet():
             if self.container.inputData.DELETE_2.__str__() in self.container.inputData.DELETE_2_LIST.qualifiers('enumerators')[:]:
                self.container.inputData.DELETE_2_LIST.set(self.container.inputData.DELETE_2.__str__())
                self.getWidget('DELETE_2_LIST').populateComboBox(self.container.inputData.DELETE_2_LIST)
                self.getWidget('DELETE_2_LIST').updateViewFromModel()
                self.validate()
             else:
                self.container.inputData.DELETE_2.unSet()
       self.updateChangeList2()
       self.reloadType2()
       return

    def updateChargeList1(self,atom_list):
       # Called after the atom list changes
       self.container.inputData.CHARGE_1_LIST.setQualifiers({'enumerators':atom_list,'menuText':atom_list})
       if self.container.inputData.CHARGE_1.__str__() in self.container.inputData.CHARGE_1_LIST.qualifiers('enumerators')[:]:
          self.container.inputData.CHARGE_1_LIST.set(self.container.inputData.CHARGE_1.__str__())
       else:
          self.container.inputData.CHARGE_1_LIST.unSet()
          self.container.inputData.CHARGE_1.unSet()
       self.getWidget('CHARGE_1_LIST').populateComboBox(self.container.inputData.CHARGE_1_LIST)
       self.getWidget('CHARGE_1_LIST').updateViewFromModel()
       self.validate()
       return
          
    def updateChargeList2(self,atom_list):
       # Called after the atom list changes
       self.container.inputData.CHARGE_2_LIST.setQualifiers({'enumerators':atom_list,'menuText':atom_list})
       if self.container.inputData.CHARGE_2.__str__() in self.container.inputData.CHARGE_2_LIST.qualifiers('enumerators')[:]:
          self.container.inputData.CHARGE_2_LIST.set(self.container.inputData.CHARGE_2.__str__())
       else:
          self.container.inputData.CHARGE_2_LIST.unSet()
          self.container.inputData.CHARGE_2.unSet()
       self.getWidget('CHARGE_2_LIST').populateComboBox(self.container.inputData.CHARGE_2_LIST)
       self.getWidget('CHARGE_2_LIST').updateViewFromModel()
       self.validate()
       return

    def ChargeListChanged1(self):
       # Called when the charge list is changed
       if self.container.inputData.CHARGE_1_LIST.isSet():
          self.container.inputData.CHARGE_1.set(self.container.inputData.CHARGE_1_LIST.__str__())
       else:
          if self.container.inputData.CHARGE_1.isSet():
             if self.container.inputData.CHARGE_1.__str__() in self.container.inputData.CHARGE_1_LIST.qualifiers('enumerators')[:]:
                self.container.inputData.CHARGE_1_LIST.set(self.container.inputData.CHARGE_1.__str__())
                self.getWidget('CHARGE_1_LIST').populateComboBox(self.container.inputData.CHARGE_1_LIST)
                self.getWidget('CHARGE_1_LIST').updateViewFromModel()
                self.validate()
             else:
                self.container.inputData.CHARGE_1.unSet()
       #self.reloadType1()
       return

    def ChargeListChanged2(self):
       # Called when the charge list is changed
       if self.container.inputData.CHARGE_2_LIST.isSet():
          self.container.inputData.CHARGE_2.set(self.container.inputData.CHARGE_2_LIST.__str__())
       else:
          if self.container.inputData.CHARGE_2.isSet():
             if self.container.inputData.CHARGE_2.__str__() in self.container.inputData.CHARGE_2_LIST.qualifiers('enumerators')[:]:
                self.container.inputData.CHARGE_2_LIST.set(self.container.inputData.CHARGE_2.__str__())
                self.getWidget('CHARGE_2_LIST').populateComboBox(self.container.inputData.CHARGE_2_LIST)
                self.getWidget('CHARGE_2_LIST').updateViewFromModel()
                self.validate()
             else:
                self.container.inputData.CHARGE_2.unSet()
       #self.reloadType2()
       return

    def updateChangeList1(self):
       # Called after updating the atom list enumerator list
       comp_id = ''
       cif_path = ''
       bond_list = ''
       if self.container.inputData.MON_1_TYPE == 'CIF':
          if self.container.inputData.DICT_1.isSet() and self.container.inputData.RES_NAME_1_CIF.isSet():
             comp_id = self.container.inputData.RES_NAME_1_CIF.__str__()
             cif_path = self.container.inputData.DICT_1.fullPath.__str__()
       else:
          if self.container.inputData.RES_NAME_1_TLC.isSet():
             comp_id = self.container.inputData.RES_NAME_1_TLC.__str__()
             if len(comp_id)>0:
                cif_path = os.path.join(CCP4Utils.getCCP4Dir().__str__(), 'lib', 'data', 'monomers', comp_id[0].lower(), comp_id+'.cif')
       if len(comp_id)>0 and len(cif_path)>0:
          bond_list = self.getBondsFromCif(cif_path,comp_id)
       self.bond_dict1 = {}
       for bond in bond_list:
          if self.container.inputData.TOGGLE_DELETE_1:
             if self.container.inputData.DELETE_1_LIST.isSet():
                if self.container.inputData.DELETE_1_LIST.__str__() in [bond[0],bond[1]]:
                   continue
          self.bond_dict1[bond[0]+" -- "+bond[1]] = bond[2]
#FIXME - I have to convert dict_keys to list because of all the isinstance stuff in CData ...
       self.container.inputData.CHANGE_BOND_1_LIST.setQualifiers({'enumerators':list(self.bond_dict1.keys()),'menuText':list(self.bond_dict1.keys())})
       self.container.inputData.CHANGE_BOND_1_LIST.unSet()
       self.getWidget('CHANGE_BOND_1_LIST').populateComboBox(self.container.inputData.CHANGE_BOND_1_LIST)
       self.getWidget('CHANGE_BOND_1_LIST').updateViewFromModel()
       self.updateChangeType1()
       return

    def updateChangeList2(self):
       # Called after updating the atom list enumerator list
       comp_id = ''
       cif_path = ''
       bond_list = ''
       if self.container.inputData.MON_2_TYPE == 'CIF':
          if self.container.inputData.DICT_2.isSet() and self.container.inputData.RES_NAME_2_CIF.isSet():
             comp_id = self.container.inputData.RES_NAME_2_CIF.__str__()
             cif_path = self.container.inputData.DICT_2.fullPath.__str__()
       else:
          if self.container.inputData.RES_NAME_2_TLC.isSet():
             comp_id = self.container.inputData.RES_NAME_2_TLC.__str__()
             if len(comp_id)>0:
                cif_path = os.path.join(CCP4Utils.getCCP4Dir().__str__(), 'lib', 'data', 'monomers', comp_id[0].lower(), comp_id+'.cif')
       if len(comp_id)>0 and len(cif_path)>0:
          bond_list = self.getBondsFromCif(cif_path,comp_id)
       self.bond_dict2 = {}
       for bond in bond_list:
          if self.container.inputData.TOGGLE_DELETE_2:
             if self.container.inputData.DELETE_2_LIST.isSet():
                if self.container.inputData.DELETE_2_LIST.__str__() in [bond[0],bond[1]]:
                   continue
          self.bond_dict2[bond[0]+" -- "+bond[1]] = bond[2]
#FIXME - I have to convert dict_keys to list because of all the isinstance stuff in CData ...
       self.container.inputData.CHANGE_BOND_2_LIST.setQualifiers({'enumerators':list(self.bond_dict2.keys()),'menuText':list(self.bond_dict2.keys())})
       self.container.inputData.CHANGE_BOND_2_LIST.unSet()
       self.getWidget('CHANGE_BOND_2_LIST').populateComboBox(self.container.inputData.CHANGE_BOND_2_LIST)
       self.getWidget('CHANGE_BOND_2_LIST').updateViewFromModel()
       self.updateChangeType2()
       return

    @QtCore.Slot()
    def updateChangeType1(self):
       if self.container.inputData.CHANGE_BOND_1_LIST.isSet():
          self.container.inputData.CHANGE_BOND_1.set(self.container.inputData.CHANGE_BOND_1_LIST.__str__())
          type = self.bond_dict1[self.container.inputData.CHANGE_BOND_1_LIST.__str__()]
          other_types = [x for x in ["single","double","triple"] if x != type]
          self.container.inputData.CHANGE_BOND_1_TYPE.setQualifiers({'enumerators':other_types,'menuText':other_types})
          self.container.inputData.CHANGE_BOND_1_TYPE.unSet()
       else:
          if self.container.inputData.CHANGE_BOND_1.isSet():
             if self.container.inputData.CHANGE_BOND_1.__str__() in self.container.inputData.CHANGE_BOND_1_LIST.qualifiers('enumerators')[:]:
                self.container.inputData.CHANGE_BOND_1_LIST.set(self.container.inputData.CHANGE_BOND_1.__str__())
                self.updateChangeType1()
                return
             self.container.inputData.CHANGE_BOND_1.unSet()
          self.container.inputData.CHANGE_BOND_1_TYPE.unSet()
          self.container.inputData.CHANGE_BOND_1_TYPE.setQualifiers({'enumerators':[],'menuText':[]})
          self.container.inputData.CHANGE_1_TYPE.unSet()
       self.getWidget('CHANGE_BOND_1_TYPE').populateComboBox(self.container.inputData.CHANGE_BOND_1_TYPE)
       self.getWidget('CHANGE_BOND_1_TYPE').updateViewFromModel()
       self.validate()
       return

    @QtCore.Slot()
    def updateChangeType2(self):
       if self.container.inputData.CHANGE_BOND_2_LIST.isSet():
          self.container.inputData.CHANGE_BOND_2.set(self.container.inputData.CHANGE_BOND_2_LIST.__str__())
          type = self.bond_dict2[self.container.inputData.CHANGE_BOND_2_LIST.__str__()]
          other_types = [x for x in ["single","double","triple"] if x != type]
          self.container.inputData.CHANGE_BOND_2_TYPE.setQualifiers({'enumerators':other_types,'menuText':other_types})
          self.container.inputData.CHANGE_BOND_2_TYPE.unSet()
       else:
          if self.container.inputData.CHANGE_BOND_2.isSet():
             if self.container.inputData.CHANGE_BOND_2.__str__() in self.container.inputData.CHANGE_BOND_2_LIST.qualifiers('enumerators')[:]:
                self.container.inputData.CHANGE_BOND_2_LIST.set(self.container.inputData.CHANGE_BOND_2.__str__())
                self.updateChangeType2()
                return
             self.container.inputData.CHANGE_BOND_2.unSet()
          self.container.inputData.CHANGE_BOND_2_TYPE.unSet()
          self.container.inputData.CHANGE_BOND_2_TYPE.setQualifiers({'enumerators':[],'menuText':[]})
          self.container.inputData.CHANGE_2_TYPE.unSet()
       self.getWidget('CHANGE_BOND_2_TYPE').populateComboBox(self.container.inputData.CHANGE_BOND_2_TYPE)
       self.getWidget('CHANGE_BOND_2_TYPE').updateViewFromModel()
       self.validate()
       return

    @QtCore.Slot()
    def updateType1(self):
       if self.container.inputData.CHANGE_BOND_1_TYPE.isSet():
          self.container.inputData.CHANGE_1_TYPE.set(self.container.inputData.CHANGE_BOND_1_TYPE.__str__())
       return

    @QtCore.Slot()
    def updateType2(self):
       if self.container.inputData.CHANGE_BOND_2_TYPE.isSet():
          self.container.inputData.CHANGE_2_TYPE.set(self.container.inputData.CHANGE_BOND_2_TYPE.__str__())
       return

    def reloadType1(self):
       # Reload bond type when refreshing the page, cloning the job, or changing the delete atom
       if not self.container.inputData.CHANGE_BOND_1_TYPE.isSet():
          if self.container.inputData.CHANGE_1_TYPE.isSet():
             if self.container.inputData.CHANGE_1_TYPE.__str__() in self.container.inputData.CHANGE_BOND_1_TYPE.qualifiers('enumerators')[:]:
                self.container.inputData.CHANGE_BOND_1_TYPE.set(self.container.inputData.CHANGE_1_TYPE)
       return

    def reloadType2(self):
       # Reload bond type when refreshing the page, cloning the job, or changing the delete atom
       if not self.container.inputData.CHANGE_BOND_2_TYPE.isSet():
          if self.container.inputData.CHANGE_2_TYPE.isSet():
             if self.container.inputData.CHANGE_2_TYPE.__str__() in self.container.inputData.CHANGE_BOND_2_TYPE.qualifiers('enumerators')[:]:
                self.container.inputData.CHANGE_BOND_2_TYPE.set(self.container.inputData.CHANGE_2_TYPE)
       return

    @QtCore.Slot()
    def AtomListChanged1TLC(self):
       # Called whenever the atom list changes
       if self.container.inputData.MON_1_TYPE == 'TLC':
          atom_list = []
          atom_list_full = []
          if self.container.inputData.ATOM_NAME_1_TLC.isSet():
             self.container.inputData.ATOM_NAME_1.set(self.container.inputData.ATOM_NAME_1_TLC.__str__())
             atom_list_full = self.container.inputData.ATOM_NAME_1_TLC.qualifiers('enumerators')[:]
             atom_list = self.container.inputData.ATOM_NAME_1_TLC.qualifiers('enumerators')[:]
             atom_list.remove(self.container.inputData.ATOM_NAME_1_TLC.__str__())
          else:
             self.container.inputData.ATOM_NAME_1.unSet()
          self.updateDeleteList1(atom_list)
          self.updateChargeList1(atom_list_full)
       return

    @QtCore.Slot()
    def AtomListChanged2TLC(self):
       # Called whenever the atom list changes
       if self.container.inputData.MON_2_TYPE == 'TLC':
          atom_list = []
          atom_list_full = []
          if self.container.inputData.ATOM_NAME_2_TLC.isSet():
             self.container.inputData.ATOM_NAME_2.set(self.container.inputData.ATOM_NAME_2_TLC.__str__())
             atom_list_full = self.container.inputData.ATOM_NAME_2_TLC.qualifiers('enumerators')[:]
             atom_list = self.container.inputData.ATOM_NAME_2_TLC.qualifiers('enumerators')[:]
             atom_list.remove(self.container.inputData.ATOM_NAME_2_TLC.__str__())
          else:
             self.container.inputData.ATOM_NAME_2.unSet()
          self.updateDeleteList2(atom_list)
          self.updateChargeList2(atom_list_full)
       return

    @QtCore.Slot()
    def AtomListChanged1CIF(self):
       # Called whenever the atom list changes
       if self.container.inputData.MON_1_TYPE == 'CIF':
          atom_list = []
          atom_list_full = []
          if self.container.inputData.ATOM_NAME_1_CIF.isSet():
             self.container.inputData.ATOM_NAME_1.set(self.container.inputData.ATOM_NAME_1_CIF.__str__())
             atom_list_full = self.container.inputData.ATOM_NAME_1_CIF.qualifiers('enumerators')[:]
             atom_list = self.container.inputData.ATOM_NAME_1_CIF.qualifiers('enumerators')[:]
             atom_list.remove(self.container.inputData.ATOM_NAME_1_CIF.__str__())
          else:
             self.container.inputData.ATOM_NAME_1.unSet()
          self.updateDeleteList1(atom_list)
          self.updateChargeList1(atom_list_full)
       return

    @QtCore.Slot()
    def AtomListChanged2CIF(self):
       # Called whenever the atom list changes
       if self.container.inputData.MON_2_TYPE == 'CIF':
          atom_list = []
          atom_list_full = []
          if self.container.inputData.ATOM_NAME_2_CIF.isSet():
             self.container.inputData.ATOM_NAME_2.set(self.container.inputData.ATOM_NAME_2_CIF.__str__())
             atom_list_full = self.container.inputData.ATOM_NAME_2_CIF.qualifiers('enumerators')[:]
             atom_list = self.container.inputData.ATOM_NAME_2_CIF.qualifiers('enumerators')[:]
             atom_list.remove(self.container.inputData.ATOM_NAME_2_CIF.__str__())
          else:
             self.container.inputData.ATOM_NAME_2.unSet()
          self.updateDeleteList2(atom_list)
          self.updateChargeList2(atom_list_full)
       return

    def getResiduesFromCif(self,cif_path):
       residue_list = []
       try:
         input_cif = gemmi.cif.read_file(cif_path)
         block = input_cif.find_block("comp_list")
         if block:
            for (id,three_letter_code) in block.find(['_chem_comp.id','_chem_comp.three_letter_code']):
               if id != three_letter_code:
                  raise Exception("Comp ID not equal to three letter code: {} {}".format(id,three_letter_code))
               residue_list.append(str(id))
       except Exception as e:
         print("Error reading residues from CIF: %s" % e)
       return residue_list

    @QtCore.Slot()
    def updateResidueList1(self):
       self.container.inputData.RES_NAME_1_CIF.setQualifiers({'enumerators':[],'menuText':[]})
       self.container.inputData.RES_NAME_1_CIF.unSet()
       self.getWidget('RES_NAME_1_CIF').populateComboBox(self.container.inputData.RES_NAME_1_CIF)
       self.getWidget('RES_NAME_1_CIF').updateViewFromModel()
       if self.container.inputData.MON_1_TYPE == 'CIF':
          if self.container.inputData.DICT_1.isSet():
             cif_path = self.container.inputData.DICT_1.fullPath.__str__()
             residue_list = self.getResiduesFromCif(cif_path)
             if len(residue_list)>0:
                self.container.inputData.RES_NAME_1_CIF.setQualifiers({'enumerators':residue_list,'menuText':residue_list})
                self.container.inputData.RES_NAME_1_CIF.set(next(iter(residue_list),None))
                self.getWidget('RES_NAME_1_CIF').populateComboBox(self.container.inputData.RES_NAME_1_CIF)
                self.getWidget('RES_NAME_1_CIF').updateViewFromModel()
          self.updateAtomList1CIF()
       return

    @QtCore.Slot()
    def updateResidueList2(self):
       self.container.inputData.RES_NAME_2_CIF.setQualifiers({'enumerators':[],'menuText':[]})
       self.container.inputData.RES_NAME_2_CIF.unSet()
       self.getWidget('RES_NAME_2_CIF').populateComboBox(self.container.inputData.RES_NAME_2_CIF)
       self.getWidget('RES_NAME_2_CIF').updateViewFromModel()
       if self.container.inputData.MON_2_TYPE == 'CIF':
          if self.container.inputData.DICT_2.isSet():
             cif_path = self.container.inputData.DICT_2.fullPath.__str__()
             residue_list = self.getResiduesFromCif(cif_path)
             if len(residue_list)>0:
                self.container.inputData.RES_NAME_2_CIF.setQualifiers({'enumerators':residue_list,'menuText':residue_list})
                self.container.inputData.RES_NAME_2_CIF.set(next(iter(residue_list),None))
                self.getWidget('RES_NAME_2_CIF').populateComboBox(self.container.inputData.RES_NAME_2_CIF)
                self.getWidget('RES_NAME_2_CIF').updateViewFromModel()
          self.updateAtomList2CIF()
       return

    @QtCore.Slot()
    def mon1TypeChanged(self):
       self.container.inputData.DICT_1.unSet()
       self.container.inputData.RES_NAME_1_CIF.unSet()
       self.container.inputData.RES_NAME_1_CIF.setQualifiers({'enumerators':[],'menuText':[]})
       self.getWidget('RES_NAME_1_CIF').populateComboBox(self.container.inputData.RES_NAME_1_CIF)
       self.getWidget('RES_NAME_1_CIF').updateViewFromModel()
       self.container.inputData.RES_NAME_1_TLC.unSet()
       self.container.inputData.ATOM_NAME_1_CIF.unSet()
       self.container.inputData.ATOM_NAME_1_CIF.setQualifiers({'enumerators':[],'menuText':[]})
       self.getWidget('ATOM_NAME_1_CIF').populateComboBox(self.container.inputData.ATOM_NAME_1_CIF)
       self.getWidget('ATOM_NAME_1_CIF').updateViewFromModel()
       self.container.inputData.ATOM_NAME_1_TLC.unSet()
       self.container.inputData.ATOM_NAME_1_TLC.setQualifiers({'enumerators':[],'menuText':[]})
       self.getWidget('ATOM_NAME_1_TLC').populateComboBox(self.container.inputData.ATOM_NAME_1_TLC)
       self.getWidget('ATOM_NAME_1_TLC').updateViewFromModel()
       self.container.inputData.ATOM_NAME_1.unSet()
       self.container.inputData.CHANGE_BOND_1_LIST.unSet()
       self.container.inputData.CHANGE_BOND_1_LIST.setQualifiers({'enumerators':[],'menuText':[]})
       self.getWidget('CHANGE_BOND_1_LIST').populateComboBox(self.container.inputData.CHANGE_BOND_1_LIST)
       self.getWidget('CHANGE_BOND_1_LIST').updateViewFromModel()
       self.container.inputData.CHANGE_BOND_1.unSet()
       self.container.inputData.CHANGE_BOND_1_TYPE.unSet()
       self.container.inputData.CHANGE_BOND_1_TYPE.setQualifiers({'enumerators':[],'menuText':[]})
       self.getWidget('CHANGE_BOND_1_TYPE').populateComboBox(self.container.inputData.CHANGE_BOND_1_TYPE)
       self.getWidget('CHANGE_BOND_1_TYPE').updateViewFromModel()
       self.container.inputData.CHANGE_1_TYPE.unSet()
       self.updateDeleteList1([])
       self.validate()
       return

    @QtCore.Slot()
    def mon2TypeChanged(self):
       self.container.inputData.DICT_2.unSet()
       self.container.inputData.RES_NAME_2_CIF.unSet()
       self.container.inputData.RES_NAME_2_CIF.setQualifiers({'enumerators':[],'menuText':[]})
       self.getWidget('RES_NAME_2_CIF').populateComboBox(self.container.inputData.RES_NAME_2_CIF)
       self.getWidget('RES_NAME_2_CIF').updateViewFromModel()
       self.container.inputData.RES_NAME_2_TLC.unSet()
       self.container.inputData.ATOM_NAME_2_CIF.unSet()
       self.container.inputData.ATOM_NAME_2_CIF.setQualifiers({'enumerators':[],'menuText':[]})
       self.getWidget('ATOM_NAME_2_CIF').populateComboBox(self.container.inputData.ATOM_NAME_2_CIF)
       self.getWidget('ATOM_NAME_2_CIF').updateViewFromModel()
       self.container.inputData.ATOM_NAME_2_TLC.unSet()
       self.container.inputData.ATOM_NAME_2_TLC.setQualifiers({'enumerators':[],'menuText':[]})
       self.getWidget('ATOM_NAME_2_TLC').populateComboBox(self.container.inputData.ATOM_NAME_2_TLC)
       self.getWidget('ATOM_NAME_2_TLC').updateViewFromModel()
       self.container.inputData.ATOM_NAME_2.unSet()
       self.container.inputData.CHANGE_BOND_2_LIST.unSet()
       self.container.inputData.CHANGE_BOND_2_LIST.setQualifiers({'enumerators':[],'menuText':[]})
       self.getWidget('CHANGE_BOND_2_LIST').populateComboBox(self.container.inputData.CHANGE_BOND_2_LIST)
       self.getWidget('CHANGE_BOND_2_LIST').updateViewFromModel()
       self.container.inputData.CHANGE_BOND_2.unSet()
       self.container.inputData.CHANGE_BOND_2_TYPE.unSet()
       self.container.inputData.CHANGE_BOND_2_TYPE.setQualifiers({'enumerators':[],'menuText':[]})
       self.getWidget('CHANGE_BOND_2_TYPE').populateComboBox(self.container.inputData.CHANGE_BOND_2_TYPE)
       self.getWidget('CHANGE_BOND_2_TYPE').updateViewFromModel()
       self.container.inputData.CHANGE_2_TYPE.unSet()
       self.updateDeleteList2([])
       self.validate()
       return

    @QtCore.Slot()
    def updateModelResLists(self):
      if self.container.inputData.XYZIN.isSet():
        path = self.container.inputData.XYZIN.__str__()
      else:
        return
      if self.container.inputData.MON_1_TYPE.__str__() == 'CIF':
        if self.container.inputData.RES_NAME_1_CIF.isSet():
          res1 = self.container.inputData.RES_NAME_1_CIF.__str__()
        else:
          return
      elif self.container.inputData.MON_1_TYPE.__str__() == 'TLC':
        if self.container.inputData.RES_NAME_1_TLC.isSet():
          res1 = self.container.inputData.RES_NAME_1_TLC.__str__()
        else:
          return
      if self.container.inputData.MON_2_TYPE.__str__() == 'CIF':
        if self.container.inputData.RES_NAME_2_CIF.isSet():
          res2 = self.container.inputData.RES_NAME_2_CIF.__str__()
        else:
          return
      elif self.container.inputData.MON_2_TYPE.__str__() == 'TLC':
        if self.container.inputData.RES_NAME_2_TLC.isSet():
          res2 = self.container.inputData.RES_NAME_2_TLC.__str__()
        else:
          return
      print("Reading residues from model: "+path)
      print("Searching for residues: "+res1+" and "+res2)
      res_list1 = []
      res_list2 = []
      try: # try to read CIF file
        doc = gemmi.cif.read(path)
      except:
        try: # try to read PDB file
          st = gemmi.read_structure(path)
          if st:
            for model in st:
              for chain in model:
                for res in chain:
                  if res.name == res1:
                    res_list1.append(res)
                  if res.name == res2:
                    res_list2.append(res)
        except Exception as e:
          print("Error: %s" % e)
          print("Cannot continue - program terminated.")
      self.container.controlParameters.MODEL_RES_LIST.unSet()
      self.container.controlParameters.MODEL_RES_LIST.setQualifiers({'enumerators':[],'menuText':[]})
      hit_list = []
      for hit1 in res_list1:
        for hit2 in res_list2:
          hit_list.append(str(hit1.name+"_"+str(hit1.seqid)+" -- "+hit2.name+"_"+str(hit2.seqid)))
      self.container.controlParameters.MODEL_RES_LIST.setQualifiers({'enumerators':hit_list,'menuText':hit_list})
      self.getWidget('MODEL_RES_LIST').populateComboBox(self.container.controlParameters.MODEL_RES_LIST)
      self.getWidget('MODEL_RES_LIST').updateViewFromModel()
      return

    @QtCore.Slot()
    def updateLinkMode(self):
      if self.container.controlParameters.LINK_MODE.__str__() == 'AUTO':
         self.container.controlParameters.MODEL_RES_LIST.setQualifiers({'allowUndefined':True})
      else:
         self.container.controlParameters.MODEL_RES_LIST.setQualifiers({'allowUndefined':False})
            
    def drawContents(self):
        self.openFolder(folderFunction='inputData',followFrom=False)
        
        ###
        self.createLine(['subtitle','First monomer to be linked'])
        self.openSubFrame(frame=True)
        self.createLine ( [ 'label','Get ligand description from','widget','MON_1_TYPE' ] )
        self.createLine ( [ 'widget','DICT_1' ], toggleFunction=[self.ToggleDict1,['MON_1_TYPE','DICT_1']])

        # Custom CIF:
        self.createLine ( [ 'label','Residue name','widget','RES_NAME_1_CIF','label','linking atom','widget','ATOM_NAME_1_CIF' ], toggleFunction=[self.ToggleCIF1,['MON_1_TYPE']] )
        # Three letter code:
        self.createLine ( [ 'label','Residue name','widget','RES_NAME_1_TLC','label','linking atom','widget','ATOM_NAME_1_TLC' ], toggleFunction=[self.ToggleTLC1,['MON_1_TYPE']] )

        toggledelete1 = self.createLine ( [ 'widget','TOGGLE_DELETE_1','label','Delete atom' ] )
        self.createLine ( [ 'widget','DELETE_1_LIST' ], toggleFunction=[self.ToggleDelete1,['TOGGLE_DELETE_1','DELETE_1_LIST']], appendLine=toggledelete1 )

        togglechange1 = self.createLine ( [ 'widget','TOGGLE_CHANGE_1','label','Change order of bond' ] )
        self.createLine ( [ 'widget','CHANGE_BOND_1_LIST','label','to','widget','CHANGE_BOND_1_TYPE' ], toggleFunction=[self.ToggleChange1,['TOGGLE_CHANGE_1','CHANGE_BOND_1_LIST','CHANGE_BOND_1_TYPE']], appendLine=togglechange1 )
        
        togglecharge1 = self.createLine ( [ 'widget','TOGGLE_CHARGE_1','label','Change the formal charge of an atom' ] )
        self.createLine ( [ 'widget','CHARGE_1_LIST','label','to','widget','CHARGE_1_VALUE' ], toggleFunction=[self.ToggleCharge1,['TOGGLE_CHARGE_1','CHARGE_1_LIST']], appendLine=togglecharge1 )
        
        self.closeSubFrame()
        
        ###
        self.createLine(['subtitle','Second monomer to be linked'])
        self.openSubFrame(frame=True)
        self.createLine ( [ 'label','Get ligand description from','widget','MON_2_TYPE' ] )
        self.createLine ( [ 'widget','DICT_2' ], toggleFunction=[self.ToggleDict2,['MON_2_TYPE','DICT_2']])

        # Custom CIF:
        self.createLine ( [ 'label','Residue name','widget','RES_NAME_2_CIF','label','linking atom','widget','ATOM_NAME_2_CIF' ], toggleFunction=[self.ToggleCIF2,['MON_2_TYPE']] )
        # Three letter code:
        self.createLine ( [ 'label','Residue name','widget','RES_NAME_2_TLC','label','linking atom','widget','ATOM_NAME_2_TLC' ], toggleFunction=[self.ToggleTLC2,['MON_2_TYPE']] )

        toggledelete2 = self.createLine ( [ 'widget','TOGGLE_DELETE_2','label','Delete atom' ] )
        self.createLine ( [ 'widget','DELETE_2_LIST' ], toggleFunction=[self.ToggleDelete2,['TOGGLE_DELETE_2','DELETE_2_LIST']], appendLine=toggledelete2 )

        togglechange2 = self.createLine ( [ 'widget','TOGGLE_CHANGE_2','label','Change order of bond' ] )
        self.createLine ( [ 'widget','CHANGE_BOND_2_LIST','label','to','widget','CHANGE_BOND_2_TYPE' ], toggleFunction=[self.ToggleChange2,['TOGGLE_CHANGE_2','CHANGE_BOND_2_LIST','CHANGE_BOND_2_TYPE']], appendLine=togglechange2 )
        
        togglecharge2 = self.createLine ( [ 'widget','TOGGLE_CHARGE_2','label','Change the formal charge of an atom' ] )
        self.createLine ( [ 'widget','CHARGE_2_LIST','label','to','widget','CHARGE_2_VALUE' ], toggleFunction=[self.ToggleCharge2,['TOGGLE_CHARGE_2','CHARGE_2_LIST']], appendLine=togglecharge2 )
        
        self.closeSubFrame()
        
        ###
        self.createLine ( [ 'label','Order of the bond between linked atoms','widget','BOND_ORDER' ] )

        ###
        self.createLine(['widget','TOGGLE_LINK','subtitle','Apply links to model (optional)'])
        self.openSubFrame(frame=True, toggle = ['TOGGLE_LINK', 'open', [ True ] ])
        #self.createLine ( [ 'widget','TOGGLE_LINK','label','Apply links to model' ] )
        self.createLine ( [ 'widget','XYZIN' ], toggle = ['TOGGLE_LINK', 'open', [ True ] ] )
        toggleautolink = self.createLine ( [ 'label', 'Apply links', 'widget','LINK_MODE' ], toggleFunction=[self.ToggleLinkMode,['TOGGLE_LINK','XYZIN']] )
        self.createLine ( [ 'label','within','widget','LINK_DISTANCE','label','times the dictionary value for this bond' ], toggleFunction=[self.ToggleAutoLinkMode,['TOGGLE_LINK','XYZIN','LINK_MODE']], appendLine=toggleautolink )
        self.createLine ( [ 'label','Create link between residues:','widget','MODEL_RES_LIST','label','Filter list by atom proximity: ', 'widget', 'MODEL_LINK_DISTANCE' ], toggleFunction=[self.ToggleManualLinkMode,['TOGGLE_LINK','XYZIN','LINK_MODE']] )

        #toggleautolink = self.createLine ( [ 'widget','TOGGLE_AUTOLINK','label','Create links between all matching atom-pairs' ] )
        #self.createLine ( [ 'label','within','widget','LINK_DISTANCE','label','times the dictionary value for this bond' ], toggleFunction=[self.ToggleAutoLink,['TOGGLE_AUTOLINK']], appendLine=toggleautolink )
        self.closeSubFrame()

        ###
        self.openFolder(folderFunction='controlParameters',title='Advanced')
        self.createLine( [ 'widget', '-guiMode','multiLine','EXTRA_ACEDRG_INSTRUCTIONS' ] )
        self.createLine( [ 'widget', '-guiMode','multiLine','EXTRA_ACEDRG_KEYWORDS' ] )

        # Watch for changes
        if self.isEditable():
           self.container.inputData.MON_1_TYPE.dataChanged.connect( self.mon1TypeChanged)
           self.container.inputData.DICT_1.dataChanged.connect( self.updateResidueList1)
           self.container.inputData.RES_NAME_1_CIF.dataChanged.connect( self.updateAtomList1CIF)
           self.container.inputData.RES_NAME_1_TLC.dataChanged.connect( self.updateAtomList1TLC)
           self.container.inputData.ATOM_NAME_1_CIF.dataChanged.connect( self.AtomListChanged1CIF)
           self.container.inputData.ATOM_NAME_1_TLC.dataChanged.connect( self.AtomListChanged1TLC)
           self.container.inputData.DELETE_1_LIST.dataChanged.connect( self.DeleteListChanged1)
           self.container.inputData.CHARGE_1_LIST.dataChanged.connect( self.ChargeListChanged1)
           self.container.inputData.CHANGE_BOND_1_LIST.dataChanged.connect( self.updateChangeType1)
           self.container.inputData.CHANGE_BOND_1_TYPE.dataChanged.connect( self.updateType1)
           ###
           self.container.inputData.MON_2_TYPE.dataChanged.connect( self.mon2TypeChanged)
           self.container.inputData.DICT_2.dataChanged.connect( self.updateResidueList2)
           self.container.inputData.RES_NAME_2_CIF.dataChanged.connect( self.updateAtomList2CIF)
           self.container.inputData.RES_NAME_2_TLC.dataChanged.connect( self.updateAtomList2TLC)
           self.container.inputData.ATOM_NAME_2_CIF.dataChanged.connect( self.AtomListChanged2CIF)
           self.container.inputData.ATOM_NAME_2_TLC.dataChanged.connect( self.AtomListChanged2TLC)
           self.container.inputData.DELETE_2_LIST.dataChanged.connect( self.DeleteListChanged2)
           self.container.inputData.CHARGE_2_LIST.dataChanged.connect( self.ChargeListChanged2)
           self.container.inputData.CHANGE_BOND_2_LIST.dataChanged.connect( self.updateChangeType2)
           self.container.inputData.CHANGE_BOND_2_TYPE.dataChanged.connect( self.updateType2)
           ###
           self.container.inputData.XYZIN.dataChanged.connect( self.updateModelResLists)
           self.container.controlParameters.LINK_MODE.dataChanged.connect( self.updateLinkMode)

        # Refresh combo boxes when reloading the page
        if self.container.inputData.MON_1_TYPE == 'CIF':
           self.updateResidueList1()
        else:
           self.updateAtomList1TLC()
        if self.container.inputData.MON_2_TYPE == 'CIF':
           self.updateResidueList2()
        else:
           self.updateAtomList2TLC()
        self.reloadType1()
        self.reloadType2()
        self.updateModelResLists()
        self.updateLinkMode()

