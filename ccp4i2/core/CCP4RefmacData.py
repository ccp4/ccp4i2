
from core import CCP4Data
from core import CCP4File

# -----------------------------------------------------------------------------------

class CRefmacRestraintsDataFile(CCP4File.CDataFile):
    QUALIFIERS = {'fileLabel' : 'restraints', 'mimeTypeName' : 'application/refmac-external-restraints',
                  'mimeTypeDescription' : 'Refmac external restraints', 'guiLabel' : 'Additional restraints',
                  'fileExtensions' : ['txt'], 'fileContentClassName' : NotImplemented}

class CRefmacRigidGroupSegment(CCP4Data.CData):
    CONTENTS = {'chain_id' : {'class' : CCP4Data.CString, 'qualifiers' : {'charWidth': 1, 'allowUndefined' : False, 'mustExist': True}},
                'residue_1' : {'class' : CCP4Data.CInt, 'qualifiers' : {'mustExist': True, 'allowUndefined': False}},
                'residue_2' : {'class' : CCP4Data.CInt, 'qualifiers' : {'mustExist': True, 'allowUndefined': False}}}
    CONTENTS_ORDER = ['chain_id', 'residue_1', 'residue_2']
    ERROR_CODES = {101 : {'description' : 'No sequence identity or structure RMS to target set'}}

    def getTableTextItems(self):
        items = []
        if self.chain_id.isSet():
            items.append(self.chain_id.__str__())
        else:
            items.append('--')
        for key in ['residue_1', 'residue_2']:
            if self.__dict__['_value'][key].isSet():
                items.append( self.__dict__['_value'][key].__str__())
            else:
                items.append('')
        return items

# def validity(self,arg):
#   err = CCP4ErrorHandling.CErrorReport()
#   err.extend(self.chain_id.validity(self.chain_id.get()))
#   if not self.residue_1.isSet() and not self.residue_2.isSet():
#     err.append(self.__class__,101,name=self.objectPath())
#   else:
#     err.extend(self.residue_1.validity(self.residue_1.get()))
#     err.extend(self.residue_2.validity(self.residue_2.get()))
#   print 'CRefmacRigidGroupSegment.validity',err.report()
#   return err


class CRefmacRigidGroupItem(CCP4Data.CData):
    CONTENTS = {'rigid_group_id' : {'class' : CCP4Data.CString},
                'segmentList' : {'class' : CCP4Data.CList, 'subItem' : {'class' : CRefmacRigidGroupSegment}, 'qualifiers' : {'listMinLength' : 1}}}

    def getTextItem(self):
        return str(self.rigid_group_id)


class CRefmacRigidGroupList(CCP4Data.CList):
    SUBITEM = {'class' : CRefmacRigidGroupItem}

    def addItem(self, value={}, index=-1):
        obj = CCP4Data.CList.addItem(self, value=NotImplemented, index=index)
        if index < 0:
            name = 'Rigid group ' + str(self.__len__())
        else:
            name = 'Rigid group ' + str(index + 1)
        while not self.isUniqueName(name):
            name = name + 'x'
        obj.rigid_group_id = name
        return obj

    def isUniqueName(self, name):
        for item in self.__dict__['_value']:
            if item.rigid_group_id == name:
                return False
        return True

# -----------------------------------------------------------------------------------

class CRefmacAnomalousAtom(CCP4Data.CData):

    CONTENTS = {'atomType': {'class': CCP4Data.CString, 'qualifiers': {'charWidth': 5, # 'mustExist': True, 'allowUndefined': False, 'onlyEnumerators': False, 'default': 'Se',
                                                                       'toolTip': 'Element name as in PDB file'}},
                'Fp': {'class': CCP4Data.CFloat, 'qualifiers' : {#'mustExist': True, 'allowUndefined': False,
                                                                 'toolTip' : "Form factor f' for element at given wavelength"}},
                'Fpp': {'class': CCP4Data.CFloat, 'qualifiers': {# 'mustExist': True,'allowUndefined': False,
                                                                 'toolTip' : "Form factor f'' for element at given wavelength"}}}

    def getTextItem(self):
        model = getattr(self, '_value', None)
        atomType = '  --  '
        item = model['atomType']
        if item.isSet() :
            value = getattr(item, '_value', None)
            if value :
                atomType = "%6s" %(value)
        Fp = '  --  '
        item = model['Fp']
        if item.isSet() :
            value = getattr(item, '_value', None)
            if value :
                Fp = "%6.3f" %(value)
        Fpp = '  --  '
        item = model['Fpp']
        if item.isSet() :
            value = getattr(item, '_value', None)
            if value :
                Fpp = "%6.3f" %(value)
        return "  %s:  %s  %s" %(atomType, Fp, Fpp)

