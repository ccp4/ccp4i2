from __future__ import print_function

"""
     CCP4XtalData.py: CCP4 GUI Project
     Copyright (C) 2010 University of York

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

"""
   Liz Potterton Jan 2010 - Created. Classes for CCP4 xtal data
                 Sep 2010 - Converted to 'generic' data style
"""

## @package CCP4XtalData (QtCore) Data objects for CCP4 crystallographic data
import os
import re
import glob
import types
import math
import sys

from PySide2 import QtCore

from core import CCP4Data
from core import CCP4File
from core import CCP4Utils
from core import CCP4Modules
from core import CCP4ModelData

from core.CCP4Config import XMLPARSER,QT
from core.CCP4ErrorHandling import *

if XMLPARSER() == 'lxml':
    from lxml import etree
else:
    from elementtree import ElementTree as etree

def SYMMETRYMANAGER():
    # Horrible mess if loadSymLib crashes!
    if CSymmetryManager.insts is None:
        CSymmetryManager.insts = CSymmetryManager()
        try:
            CSymmetryManager.insts.loadSymLib()
        except CException as e:
            pass
    return CSymmetryManager.insts


class CSymmetryManager:
    insts = None
    ERROR_CODES= {101 : 'CCP4 directory undefined - can not read symmetry library',
                  102 : 'Failed to find/open symmetry library',
                  103 : 'Error reading symmetry library' }

    def __init__(self):
        self.hmSpaceGroupList = []
        self.crystalSystems = ['triclinic','monoclinic','orthorhombic','tetragonal','trigonal','hexagonal','cubic']
        self.chiralSpaceGroups = {'triclinic' : ['P 1'],
                                  'monoclinic' : ['P 1 2 1','P 1 21 1','C 1 2 1','I 1 2 1'],
                                  'orthorhombic' : ['P 2 2 2','P 2 2 21','P 2 21 2', 'P 21 2 2',
                                                    'P 21 21 2','P 21 2 21','P 2 21 21','P 21 21 21','C 2 2 21',
                                                    'C 2 2 2','F 2 2 2','I 2 2 2','I 21 21 21'],
                                  'tetragonal' : ['P 4','P 41','P 42','P 43','I 4','I 41','P 4 2 2',
                                                  'P 4 21 2','P 41 2 2','P 41 21 2','P 42 2 2','P 42 21 2',
                                                  'P 43 2 2','P 43 21 2','I 4 2 2','I 41 2 2'],
                                  'trigonal' : ['P 3', 'P 31','P 32','H 3',
                                                'P 3 1 2','P 3 2 1','P 31 1 2','P 31 2 1',
                                                'P 32 1 2', 'P 32 2 1','H 3 2',
                                                'R 3 :H', 'R 3 :R','R 3 2 :H', 'R 3 2 :R'],
                                  'hexagonal' : ['P 6', 'P 61', 'P 65','P 62','P 64', 'P 63','P 6 2 2',
                                                 'P 61 2 2','P 65 2 2','P 62 2 2','P 64 2 2','P 63 2 2'],
                                  'cubic'  : ['P 2 3','F 2 3','I 2 3','P 21 3','I 21 3',
                                              'P 4 3 2', 'P 42 3 2','F 4 3 2',
                                              'F 41 3 2','I 4 3 2','P 43 3 2','P 41 3 2','I 41 3 2']}
        self.laueGroups = [['P 21 3', 'P 2 3'],
                           ['P 41 3 2', 'P 43 3 2', 'P 42 3 2', 'P 4 3 2'],
                           ['I 21 3', 'I 2 3'], ['I 4 3 2', 'I 41 3 2'],
                           ['F 41 3 2', 'F 4 3 2'],
                           ['P 31', 'P 3', 'P 32'],
                           ['P 31 1 2', 'P 3 1 2', 'P 32 1 2'],
                           ['P 31 2 1', 'P 3 2 1', 'P 32 2 1'],
                           ['P 61', 'P 6', 'P 65', 'P 62', 'P 64', 'P 63'],
                           ['P 61 2 2', 'P 2 2 21', 'P 1 c 1', 'P 2 2 21', 'C 1 c 1', 'P 62 2 2', 'P 64 2 2', 'P 63 2 2'],
                           ['P 41', 'P 4', 'P 42', 'P 43'],
                           ['P 41 21 2', 'P 4 2 2', 'P 4 21 2', 'P 41 2 2', 'P 43 2 2', 'P 42 2 2', 'P 42 21 2', 'P 43 21 2', 1094],
                           ['I 41', 'I 4'], ['I 41 2 2', 'I 4 2 2'],
                           ['P 21 21 21', 'P 2 2 2', 'P 2 2 21', 'P 21 21 2', 'P 21 2 2', 'P 2 21 2', 'P 21 2 21', 'P 2 21 21'],
                           ['C 2 2 21', 'C 2 2 2', 1020, 1021],
                           ['I 21 21 21', 'I 2 2 2', 1023],
                           ['P 1 21 1', 'P 1 2 1', 'P 1 1 2', 'P 1 1 21']]

    def convertLaue(self):
        out = []
        for lG in self.laueGroups:
            out.append([])
            for sG in lG:
                if self.ccp4NumberList.count(sG):
                    out[-1].append(self.hmSpaceGroupList[self.ccp4NumberList.index(sG)])
                    if out[-1][-1] == '' :
                        out[-1][-1] = sG
                else:
                    out[-1].append( sG)
        print(out)

    def loadSymLib(self, fileName=None):
        if fileName is None:
            path = CCP4Utils.getCCP4Dir()
            if path is None or len(path) == 0:
                raise CException(self.__class__, 101, name=self.objectPath())
            fileName = os.path.normpath(os.path.join(path, 'lib', 'data', 'syminfo.lib'))
        try:
            text = CCP4Utils.readFile(fileName)
        except CException:
            raise CException(self.__class__, 102, fileName, name=self.objectPath())
        lineList = text.split('\n')
        il = -1
        self.hmSpaceGroupList = []
        self.oldSpaceGroupList = []
        self.pointGroupList = []
        self.ccp4NumberList = []
        xHM = None
        old = []
        pgrp = None
        while il < len(lineList):
            il = il + 1
            if lineList[il][0:16] == 'begin_spacegroup':
                il = il + 1
                while il < len(lineList):
                    if lineList[il][0:6] == 'symbol':
                        if lineList[il].count('xHM'):
                            s = re.search(r'(.*?)\'(.*?)\'(.*)',lineList[il])
                            if s is not None:
                                xHM = s.groups()[1]
                        elif lineList[il].count('pgrp'):
                            s = re.search(r'(.*?)\'(.*?)\'(.*)',lineList[il])
                            if s is not None:
                                pgrp = s.groups()[1]
                        elif lineList[il].count('old'):
                            s = re.search(r'(.*?)\'(.*?)\'(.*)',lineList[il])
                            while s is not None:
                                o = s.groups()[1]
                                if len(o)>0:
                                    old.append(o)
                                s = re.search(r'(.*?)\'(.*?)\'(.*)',s.groups()[2])
                        elif lineList[il].count('ccp4'):
                            ccp4 = int(lineList[il].split()[-1])
                    elif lineList[il][0:14] == 'end_spacegroup':
                        if xHM is not None:
                            self.hmSpaceGroupList.append(xHM)
                            self.oldSpaceGroupList.append(old)
                            self.pointGroupList.append(pgrp)
                            self.ccp4NumberList.append(ccp4)
                        xHM = None
                        old = []
                    il = il + 1

    def isChiral(self, spaceGroup):
        for key, value in list(self.chiralSpaceGroups.items()):
            if value.count(spaceGroup) > 0:
                return True
        return False

    def crystalSystem(self, spaceGroup):
        for key, value in list(self.chiralSpaceGroups.items()):
            if value.count(spaceGroup) > 0:
                return key
        return None

    def nonEnantiogenicList(self):
        l = []
        for system in self.crystalSystems:
            l.extend(self.chiralSpaceGroups[system]);
        # for system, nameList in self.chiralSpaceGroups.items(): l.extend(nameList)
        return l

    def spaceGroupValidity(self, spaceGroup=None):
        if spaceGroup is None or len(spaceGroup)==0:
            return (2, None)
        up = spaceGroup.upper()
        if up != spaceGroup:
            code = 100
            spaceGroup = up
        else:
            code = 0
        if len(self.hmSpaceGroupList) == 0:
            # No data from syminfo.lib so can only check if it is chiral
            if self.isChiral(spaceGroup):
                return (0 + code, spaceGroup)
            else:
                return (6 + code, spaceGroup)
        if self.hmSpaceGroupList.count(spaceGroup):
            ii = self.hmSpaceGroupList.index(spaceGroup)
            if self.isChiral(spaceGroup):
                return (0 + code, spaceGroup)
            else:
                return (1 + code, spaceGroup)
        # Is it an old space group?
        for ii in range(len(self.oldSpaceGroupList)):
            if self.oldSpaceGroupList[ii].count(spaceGroup):
                return (3 + code, self.hmSpaceGroupList[ii])
        # Is it a cut-down version of hm, match without spaces
        hits = []
        for item in self.hmSpaceGroupList:
            if re.match(spaceGroup + ' ', item):
                hits.append(item)
            elif re.match(spaceGroup + ' ', re.sub(' ', '', item)):
                hits.append(item)
        if len(hits) == 1:
            return (4 + code, hits[0])
        elif len(hits) > 1:
            return (5 + code, hits)
        # Is it a cut-down version of old
        hits = []
        for item in self.oldSpaceGroupList:
            if len(item) > 0:
                if re.match(spaceGroup, item[0]):
                    hits.append(item)
                elif re.match(spaceGroup,re.sub(' ', '', item[0])):
                    hits.append(item)
        if len(hits) == 1:
            return (3+code, hits[0])
        elif len(hits) > 1:
            return (5+code, hits)
        # Is it wrong spaces
        noSpace = re.sub(' ', '', spaceGroup)
        for item in self.hmSpaceGroupList:
            if noSpace == re.sub(' ', '', item):
                hits.append(item)
        if len(hits) > 0:
            return (7 + code, hits)
        return (2 + code, spaceGroup)

    def spaceGroupCompleter(self,spaceGroup=None):
        up = spaceGroup.upper().strip()
        hits = []
        for system, groups in list(self.chiralSpaceGroups.items()):
            for gp in groups:
                if re.match(spaceGroup, gp):
                    hits.append(gp)
        return hits

    def ccp4toHm(self,num):
        if self.ccp4NumberList.count(num):
            return self.hmSpaceGroupList[self.ccp4NumberList.index(num)]
        else:
            return None

    def hmtoCcp4(self,name):
        if name is None or name is NotImplemented:
            return None
        if self.hmSpaceGroupList.count(name):
            return self.ccp4NumberList[self.hmSpaceGroupList.index(name)]
        else:
            return None

class CSpaceGroup(CCP4Data.CString):
    '''A string holding the space group'''
    QUALIFIERS = {'allowUndefined' : True, 'toolTip' : 'Hermann-Mauguin space group name',
                  'helpFile' : 'crystal_data#space_group'}
    ERROR_CODES = {101 : {'description' : 'Invalid space group'},
                   102 : {'description' : 'Space group is not chiral', 'severity' : SEVERITY_WARNING},
                   103 : {'description' : 'Space group is not Hermann-Mauguin standard'},
                   104 : {'description' : 'Space group is not a chiral Hermann-Mauguin standard. Full syminfo.lib information not loaded.'},
                   105 : {'description' : 'Space group is not Hermann-Mauguin standard - has wrong number of spaces?'},
                   106 : {'description' : 'Space group is undefined', 'severity' : SEVERITY_UNDEFINED},
                   107 : {'description' : 'Space group is undefined'},
                   108 : {'description' : 'Space group is incomplete', 'severity' : SEVERITY_WARNING}}

    def validity(self,arg):
        ''' Need to check is valid space group name '''
        err= CErrorReport()
        if arg is None or len(str(arg)) == 0:
            if self.qualifiers('allowUndefined'):
                err.append(self.__class__, 106, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            else:
                err.append(self.__class__, 107, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            return err
        arg = str(arg)
        status, correctedSG = SYMMETRYMANAGER().spaceGroupValidity(str(arg))
        if correctedSG != arg:
            details = 'Could be: '
            if isinstance(correctedSG, list):
                for item in correctedSG:
                    details = details + ' ' + str(item) + ','
                details = details[0:-1]
            else:
                details = details + ' ' + str(correctedSG)
        else:
            details = arg
        status = status%100  # +100 indicates lowercase corrected
        if status == 1:
            err.append(self.__class__, 102, details, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        elif status == 2:
            err.append(self.__class__, 101, details, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        elif 3 <= status <= 4:
            err.append(self.__class__, 103, details, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        elif status == 5:
            err.append(self.__class__, 108, details, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        elif status == 6:
            err.append(self.__class__, 104, details, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        return err

    def set(self, value=None, checkValidity=True):
        if isinstance(value, int):
            CCP4Data.CString.set(self, SYMMETRYMANAGER().ccp4toHm(value), checkValidity=checkValidity)
        else:
            CCP4Data.CString.set(self, value, checkValidity=checkValidity)

    def fix(self, arg):
        if arg is None or len(arg) == 0:
            return arg
        status, correctedSG = SYMMETRYMANAGER().spaceGroupValidity(arg)
        if isinstance(correctedSG, list):
            if len(correctedSG) == 1:
                return correctedSG[0]
            else:
                return arg
        else:
            return correctedSG

    def number(self):
        return SYMMETRYMANAGER().hmtoCcp4(self.__dict__['_value'])


# What extra functionality is needed here?
class CReindexOperator(CCP4Data.CData):

    CONTENTS = {'h' : {'class' : CCP4Data.CString, 'qualifiers' : { 'default' : 'h'}},
                'k' : {'class' : CCP4Data.CString, 'qualifiers' : { 'default' : 'k'}},
                'l' : {'class' : CCP4Data.CString, 'qualifiers' : { 'default' : 'l'}}}
    CONTENTS_ORDER = ['h','k','l']
    ERROR_CODES = {201 : {'description' : 'Operator has bad syntax (needs three comma-separated fields)'},
                   202 : {'description' : 'Operator contains invalid characters'},
                   203 : {'description' : 'Operator is not set'}}

    def validity(self,arg):
        err = CErrorReport()
        if arg is None:
            err.append(self.__class__, 106, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            return err
        for key in  ['h','k','l']:
            op = self.__dict__['_value'][key].__str__().strip()
            if len(op) < 1:
                err.append(self.__class__, 203, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            s = re.search('[^0-9,\,,\-,\+,\/,h,k,l]', op)
            if s is not None:
                err.append(self.__class__, 202, name=self.objectPath())
        return err

    def isSet(self, allowUndefined=False, allowDefault=True, allSet=True):
        for key in ['h', 'k', 'l']:
            if self.__dict__['_value'][key].__str__() != key:
                return True
        return False


class CCellLength(CCP4Data.CFloat):
    '''A cell length'''
    QUALIFIERS = {'min' : 0.0, 'default' : None, 'allowUndefined' : False, 'toolTip' : 'Cell length in A' }

    def setNm(self,value):
        if value is None:
            self.set(None)
        else:
            self.set(value*10.0)

    def getNm(self):
        if self._value is None:
            return None
        else:
            return self._value / 10.0

    PROPERTIES = {'nm' : {'fget' : getNm, 'fset' : setNm}}


class CCellAngle(CCP4Data.CFloat):
    '''A cell angle'''

    QUALIFIERS = {'min' : 0.0, 'max' : 180.0, 'default' : None, 'allowUndefined' : True, 'toolTip' : 'Cell angle in degrees'}

    def getRadians(self):
        if self._value is None:
            return 0.0
        else:
            return self._value * (math.pi / 180.0)

    def setRadians(self,value):
        if value is None:
            self.set(value)
        else:
            self.set(value * (180.0 / math.pi))

    PROPERTIES = {'rad' : {'fget' : getRadians , 'fset' : setRadians}}


class CCell(CCP4Data.CData):
    '''A unit cell'''
    CONTENTS = {'a' : { 'class' : CCellLength, 'qualifiers' : {'toolTip' : 'Cell length a in A', 'guiLabel' : 'a'}},
                'b' : { 'class' : CCellLength, 'qualifiers' : {'toolTip' : 'Cell length b in A', 'guiLabel' : 'b'}},
                'c' : { 'class' : CCellLength, 'qualifiers' : {'toolTip' : 'Cell length c in A', 'guiLabel' : 'c'}},
                'alpha' : {'class' : CCellAngle, 'qualifiers' : {'toolTip' : 'Cell angle alpha in degrees', 'guiLabel' : 'alpha'}},
                'beta' : {'class' : CCellAngle, 'qualifiers' : {'toolTip' : 'Cell angle beta in degrees', 'guiLabel' : 'beta' }},
                'gamma' : {'class' : CCellAngle,  'qualifiers' : {'toolTip' : 'Cell angle gamma in degrees', 'guiLabel' : 'gamma'}}}
    CONTENTS_ORDER = ['a' , 'b', 'c', 'alpha', 'beta', 'gamma']
    QUALIFIERS = {'toolTip' : 'Cell lengths and angles', 'helpFile' : 'crystal_data#cell'}

    def fix(self, arg):
        if arg['a'] is not None and arg['b'] is not None and arg['c'] is not None:
            for item in ['alpha', 'beta', 'gamma']:
                if arg[item] is None:
                    arg[item] = 90.0
                elif arg[item] < 0.0:
                    arg[item] = 0.0
                elif arg[item] > 180.0:
                    arg[item] = 180.0
                elif arg[item] < 3.2:
                    arg[item] = round(arg[item] * (180.0 / math.pi), 4)
                else:
                    arg[item] = round(arg[item], 4)
            for item in ['a', 'b', 'c']:
                arg[item] = round(arg[item], 4)
            return arg
        else:
            return {'a' :None , 'b' : None , 'c' :None , 'alpha' :None , 'beta' : None, 'gamma' : None}

    def __cmp__(self,other):
        for item in self.CONTENTS_ORDER:
            c = self.__dict__['_value'][item].__cmp__(other.__dict__['_value'][item])
            if c != 0:
                return c
        return 0

    def isSet(self, allowUndefined=False, allowDefault=True, allSet=True):
        # Allow the cell angles to be unset
        return CCP4Data.CData.isSet(self, allSet=False)

    def set(self, value=None):
        if isinstance(value,list) and len(value) in [3,6]:
            d = {'a' : (value[0]), 'b' : (value[1]), 'c' : (value[2])}
            if len(value) == 6:
                d.update ({'alpha' : (value[3]), 'beta' : (value[4]), 'gamma' : (value[5])})
            else:
                d.update ({ 'alpha' : None, 'beta' :None, 'gamma' : None})
            value = d
        CCP4Data.CData.set(self,value)

    def guiLabel(self):
        if self.__dict__['_value']['alpha'].isSet():
            s = '%7.1f,%7.1f,%7.1f,%7.1f,%7.1f,%7.1f'%( self.__dict__['_value']['a'].__float__(), self.__dict__['_value']['b'].__float__(),
                                                        self.__dict__['_value']['c'].__float__(), self.__dict__['_value']['alpha'].__float__(),
                                                        self.__dict__['_value']['beta'].__float__(), self.__dict__['_value']['gamma'].__float__())
        elif self.__dict__['_value']['a'].isSet():
            s = '%7.1f,%7.1f,%7.1f,90,90,90'%(self.__dict__['_value']['a'], self.__dict__['_value']['b'], self.__dict__['_value']['c'])
        else:
            s = 'Unknown'
        s = re.sub(' ', '', s)
        return s

class CSpaceGroupCell(CCP4Data.CData):
    '''Cell space group and parameters'''
    CONTENTS = {'spaceGroup' : {'class' : CSpaceGroup, 'qualifiers' : {'guilabel' : 'space group'}},
                'cell' : {'class' : CCell, 'qualifiers' : { 'guilabel' : 'cell'}}}
    CONTENTS_ORDER = ['spaceGroup', 'cell']
    QUALIFIERS = {'toolTip' : 'Space group and cell length and angles',
                  'helpFile' : 'crystal_data#cell_space_group'}

    ERROR_CODES = {101 : {'description' : 'Cell lengths should NOT be identical'},
                   102 : {'description' : 'Cell angles should NOT be identical'},
                   103 : {'description' : 'Cell angle should be 90'},
                   104 : {'description' : 'Cell angle should NOT be 90'},
                   105 : {'description' : 'Cell lengths should be identical'},
                   106 : {'description' : 'Cell angle should be 120'},
                   107 : {'description' : 'Cell angle should be identical'}}

    def isNinety(self, arg):
        if arg is None or abs(arg - 90.0) < 0.0001:
            return True
        else:
            return False

    def validity(self, arg):
        ''' Needs checking of cell paramenters to be consistent with space group '''
        #
        # triclinic     a != b != c; alpha != beta != gamma  or anything
        # monoclinic    a != b != c; alpha= gamma = 90; beta != 90  a,b,c,beta unspecified
        # orthorhombic  a != b != c; alpha = beta = gamma = 90 a,b,c unspecified
        # tetragonal    a = b != c; alpha = beta = gamma = 90  c may be = a,b
        # rhombohedral  a = b = c;  alpha = beta = gamma != 90
        # hexagonal     a = b != c; alpha = beta = 90; gamma = 120 c may be = a,b
        # cubic         a = b = c;  alpha = beta = gamma = 90
        v = self.itemValidity()
        if len(v) > 0:
            print('CSpaceGroupCellitem.Validity', v.report())
            return v
        cell = arg.get('cell')
        if cell is None:
            return v
        for ang in ['alpha', 'beta', 'gamma']:
            if cell[ang] is None or cell[ang] == 0.0:
                cell[ang] = 90.0
            elif cell[ang] < 3.0:
                cell[ang] = cell[ang] * 180.0/math.pi
        xtlSys = SYMMETRYMANAGER().crystalSystem(arg.get('spaceGroup', None))
        if xtlSys is None:
            pass
        elif xtlSys == 'triclinic':
            pass  # allow anything
            """
            if cell['a'] == cell['b'] or cell['a'] == cell['c'] or cell['b'] == cell['c']:
                v.append(self.__class__, 101, 'all cell lengths in a triclinic space group',
                         name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            if cell['alpha'] == cell['beta'] or cell['alpha'] == cell['gamma'] or cell['beta'] == cell['gamma']:
                v.append(self.__class__, 102, 'all cell angles in a triclinic space group',
                         name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            """
        elif xtlSys == 'monoclinic':
            """
            if cell['a'] == cell['b'] or cell['a'] == cell['c'] or cell['b'] == cell['c']:
                v.append(self.__class__, 101, 'all cell lengths in a monoclinic space group',
                         name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            """
            if not self.isNinety(cell['alpha']) or not self.isNinety(cell['gamma']):
                v.append(self.__class__, 103, 'alpha and gamma in monoclinic space group',
                         name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            """
            if self.isNinety(cell['beta']):
                v.append(self.__class__, 104, 'beta in monoclinic space group',
                         name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            """
        elif xtlSys == 'orthorhombic':
            """
            if cell['a'] == cell['b'] or cell['a'] == cell['c'] or cell['b'] == cell['c']:
                v.append(self.__class__, 101, 'all cell lengths in an orthorhombic space group',
                         name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            """
            if not self.isNinety(cell['alpha']) or not self.isNinety(cell['beta']) or not self.isNinety(cell['gamma']):
                v.append(self.__class__, 103, 'all cell angles in an orthorhombic space group',
                         name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        elif xtlSys == 'tetragonal':
            if cell['a'] != cell['b']:
                v.append(self.__class__, 105, 'a and b in a tetragonal space group',
                         name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            """
            if cell['a'] == cell['c']:
                v.append(self.__class__, 101, 'a/b and c in a tetragonal space group',
                         name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            """
            if not self.isNinety(cell['alpha']) or not self.isNinety(cell['beta']) or not self.isNinety(cell['gamma']):
                v.append(self.__class__, 103, 'all cell angles in a tetragonal space group',
                         name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        elif xtlSys == 'rhombohedral':
            if cell['a'] != cell['b'] or cell['a'] != cell['c'] or cell['b'] != cell['c']:
                v.append(self.__class__, 105, 'All cell lengths in a rhombohedral space group',
                         name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            """
            if self.isNinety(cell['alpha']) or self.isNinety(cell['beta']) or self.isNinety(cell['gamma']):
                v.append(self.__class__, 104, 'All cell angles in a rhombohedral space group',
                         name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            """
            if cell['alpha'] != cell['beta'] or cell['alpha'] != cell['gamma']:
                v.append(self.__class__, 107, 'All cell angles in a rhombohedral space group',
                         name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        elif xtlSys == 'hexagonal':
            if cell['a'] != cell['b']:
                v.append(self.__class__, 105, 'a and b in a hexagonal space group',
                         name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            """
            if cell['a'] == cell['c'] or cell['b'] == cell['c']:
                v.append(self.__class__, 101, 'a/b and c in a hexagonal space group',
                         name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            """
            if not self.isNinety(cell['alpha']) or not self.isNinety(cell['beta']):
                v.append(self.__class__, 103, 'alpha and beta in a hexagonal space group',
                         name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            if cell['gamma'] != 120.:
                v.append(self.__class__, 106, 'gamma in a hexagonal space group',
                         name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        elif xtlSys == 'cubic':
            if cell['a'] != cell['b'] or cell['a'] != cell['c'] or cell['b'] != cell['c']:
                v.append(self.__class__, 105, 'All cell lengths in a cubic space group',
                         name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            if not self.isNinety(cell['alpha']) or not self.isNinety(cell['beta']) or not self.isNinety(cell['gamma']):
                v.append(self.__class__, 104, 'all cell angles in a cubic space group',
                         name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        return v


class CResolutionRange(CCP4Data.CData):
    CONTENTS = {'low' : {'class' : CCP4Data.CFloat, 'qualifiers' : {'min' : 0.0, 'allowUndefined' : True}},
                'high' : {'class' : CCP4Data.CFloat , 'qualifiers' : {'min' : 0.0, 'allowUndefined' : True}}}
    ERROR_CODES = {201 : {'description' : 'High/low resolution wrong way round?'}}

    def validity(self, arg):
        v = self.itemValidity(arg)
        if v.maxSeverity() > 0:
            return v
        if arg.get('high').__gt__(arg.get('low')):
            v.append(self.__class__, 201, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        return v


class CWavelength(CCP4Data.CFloat):
    '''Wavelength in Angstrom'''
    QUALIFIERS = {'min' : 0.0, 'toolTip' : 'Data collection wavelength in Angstrom' }

    def fix(self, arg = None):
        if arg is None or arg is NotImplemented:
            return arg
        arg = round(arg, 5)
        return arg

    def setNm(self, value):
        if value is None:
            self.set(None)
        else:
            self.set(value * 10.0)

    def getNm(self):
        if self._value is None:
            return None
        else:
            return self._value / 10.0

    PROPERTIES = {'nm' : {'fget' : getNm , 'fset' : setNm } }


class CAltSpaceGroup(CSpaceGroup):
    def guiLabel(self):
        if self.__dict__['_value'] is None:
            return '-'
        else:
            return str(self.__dict__['_value'])


class CAltSpaceGroupList(CCP4Data.CList):
    SUBITEM = { 'class' : CAltSpaceGroup }

    def validity(self,arg):
        try:
            mode = self.parent().get('SGALT_SELECT').__str__()
        except:
            mode = 'LIST'
        if mode != 'LIST':
            return CErrorReport()
        else:
            err = CCP4Data.CList.validity(self, arg)
            if len(arg) < 1:
                err.append(self.__class__, 101, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        return err

class CMapDataFile(CCP4File.CDataFile):
    '''A CCP4 Map file'''
    QUALIFIERS = {'mimeTypeName' : 'application/CCP4-map',
                  'mimeTypeDescription' : 'CCP4 Electron density map',
                  'fileExtensions' : ['map'],
                  'fileContentClassName' : None,
                  'guiLabel' : 'Electron Density Map',
                  'toolTip' : 'A map in CCP4 format',
                  'helpFile' : 'data_files#map_files' }

class CGenericReflDataFile(CCP4File.CDataFile):
    QUALIFIERS = {'guiLabel' : 'Reflection data',
                  'mimeTypeName' : "application/CCP4-generic-reflections",
                  'toolTip' : 'A reflection data file in MTZ or a non-CCP4 format',
                  'fileContentClassName' : 'CUnmergedDataContent',
                  'fileExtensions' :  ['mtz', 'hkl', 'HKL', 'sca', 'SCA', 'mmcif', 'cif', 'ent'],
                  'downloadModes' : ['ebiSFs'],
                  'helpFile' : 'import_merged#file_formats' }

    def __init__(self, value={}, qualifiers={}, parent=None, name=None, fullPath=None, keywords={}, **kw):
        CCP4File.CDataFile.__init__(self, value=value, qualifiers=qualifiers, parent=parent, name=name, fullPath=fullPath, keywords=keywords, **kw)

    def getFormat(self):
        if isinstance(self.__dict__['_fileContent'],CUnmergedDataContent):
            return str(self.__dict__['_fileContent'].format)
        elif  isinstance(self.__dict__['_fileContent'],CMtzData):
            return 'mtz'
        elif  isinstance(self.__dict__['_fileContent'],CMmcifReflData):
            return 'mmcif'

    def getMerged(self):
        if isinstance(self.__dict__['_fileContent'],CUnmergedDataContent):
            return self.__dict__['_fileContent'].merged
        elif isinstance(self.__dict__['_fileContent'],CMtzData):
            return 'merged'
        elif  isinstance(self.__dict__['_fileContent'],CMmcifReflData):
            return 'merged'

    def getFileContent(self):
        from core import CCP4DataManager
        contentClass = self.fileContentClass()
        if self.__dict__['_fileContent'] is None or self.__dict__['_fileContent'].__class__ != contentClass:
            self.__dict__['_fileContent'] = None
        rv = self.loadFile()
        if len(rv) > 0:
            print(rv.report())
        if isinstance(self.__dict__['_fileContent'],CUnmergedDataContent) and self.__dict__['_fileContent'].format == 'mtz' and self.__dict__['_fileContent'].merged == 'merged':
            cls = CCP4DataManager.DATAMANAGER().getClass('CMtzData')
            self.__dict__['_fileContent'] = cls()
            self.loadFile()
        return self.__dict__['_fileContent']

    def fileContentClass(self,className=None):
        from core import CCP4DataManager
        if className is None:
            if not self.exists() or (self.getExt() not in ['.mmcif','.cif','.ent']):
                className = 'CUnmergedDataContent'
            else:
                className = 'CMmcifReflData'
        cls = CCP4DataManager.DATAMANAGER().getClass(className)
        if cls is None:
            raise CException(self.__class__, 105, 'Contents class name:' + str(className), name=self.objectPath())
        return cls

    def importFileName(self, jobId=None, jobDirectory=None, ext=None):
        if ext is None:
            try:
                ext = self.getExt()
            except:
                pass
        return CCP4File.CDataFile.importFileName(self, jobId=jobId, jobDirectory=jobDirectory, ext=ext)

class CMmcifReflDataFile(CCP4File.CMmcifDataFile):
    '''A reflection file in mmCIF format'''
    QUALIFIERS = {'guiLabel' : 'mmCIF reflection data',
                  'mimeTypeName' : 'chemical/x-cif',
                  'toolTip' : 'A reflection file in mmCIF format',
                  'fileContentClassName' : 'CMmcifReflData',
                  'helpFile' : 'data_files#mmCIF'}

    def __init__(self, value={}, qualifiers={}, parent=None, name=None, fullPath=None, keywords={}, **kw):
        CCP4File.CDataFile.__init__(self, value=value, qualifiers=qualifiers, parent=parent, name=name,
                                    fullPath=fullPath, keywords=keywords, **kw)

    def updateData(self):
        self.loadFile()
        self.dataChanged.emit()

class CMmcifReflData(CCP4File.CMmcifData):
    dataLoaded = QtCore.Signal()
    '''Reflection data in mmCIF format'''

    CONTENTS = {'cell' : {'class' : CCell },
                'spaceGroup' : {'class' : CSpaceGroup },
                'wavelength' : {'class' : CWavelength },
                'haveFreeRColumn' : {'class' : CCP4Data.CBoolean},
                'haveFobsColumn' : {'class' : CCP4Data.CBoolean},
                'haveFpmObsColumn' : {'class' : CCP4Data.CBoolean},
                'haveIobsColumn' : {'class' : CCP4Data.CBoolean},
                'haveIpmObsColumn' : {'class' : CCP4Data.CBoolean}}
    ERROR_CODES = {101 : {'description' : 'Attempting to load mmCIF data from non-existant/broken file'},
                   102 : {'description' : 'Error reading interpreting line in cif file'}}

    def loadFile(self, fileName=None):
        print('CMmcifReflData.loadFile',fileName)
        if fileName is None or not os.path.exists(str(fileName)):
            self.unSet()
            self.__dict__['lastLoadedFile'] = None
            return
        fileName = str(fileName)
        # Beware lastLoadedFile might not have been unSet
        err = CErrorReport()
        mmcifLines = CCP4Utils.readFile(fileName).split("\n")[0:500]
        cell = {'a':None, 'b':None, 'c':None, 'alpha':None, 'beta':None, 'gamma':None}
        wavelength = None
        pyStrSpaceGroupName = ''
        iSpaceGroupNumber = -1
        haveFreeRColumn = False
        haveFobsColumn = False
        haveFpmObsColumn = False
        haveIobsColumn = False
        haveIpmObsColumn = False
        from core.CCP4Utils import safeFloat
        for j, pyStrLine in enumerate(mmcifLines):
            try:
                if "_symmetry.Int_Tables_number" in pyStrLine:
                    iSpaceGroupNumber = int(pyStrLine.split()[1])
                if "_symmetry.space_group_name_H-M" in pyStrLine:
                    pyStrSpaceGroupName = pyStrLine[pyStrLine.index("_symmetry.space_group_name_H-M")+31:]
                    pyStrSpaceGroupName = pyStrSpaceGroupName.strip().strip('"').strip("'").strip()
                if "_refln.status" in pyStrLine:
                    haveFreeRColumn = True
                if "_refln.F_meas" in pyStrLine or "_refln.F_meas_au" in pyStrLine:
                    haveFobsColumn = True
                if "_refln.pdbx_F_plus" in pyStrLine:
                    haveFpmObsColumn = True
                if "_refln.intensity_meas" in pyStrLine or "_refln.F_squared_meas" in pyStrLine:
                    haveIobsColumn = True
                if "_refln.pdbx_I_plus" in pyStrLine:
                    haveIpmObsColumn = True
                if "_diffrn_radiation_wavelength.wavelength" in pyStrLine:
                    print('wavelength',pyStrLine.split()[1])
                    wavelength = safeFloat(pyStrLine.split()[1])
                if "_cell." in pyStrLine:
                    for cif_item in ['length_a','length_b','length_c','angle_alpha','angle_beta','angle_gamma']:
                        if "_cell." + cif_item in pyStrLine:
                            value = safeFloat(pyStrLine.split()[1])
                            typ,item = cif_item.split('_')
                            if typ == 'angle' and value < 3.0:
                                cell[item] = value * (180.0/math.pi)
                            else:
                                cell[item] = value
            except:
                err.append(self.__class__, 102, str(pyStrLine))
        self.__dict__['_value']['cell'].set(cell)
        self.__dict__['_value']['wavelength'].set(wavelength)
        self.__dict__['_value']['spaceGroup'].set(pyStrSpaceGroupName)
        self.__dict__['_value']['haveFreeRColumn'].set(haveFreeRColumn)
        self.__dict__['_value']['haveFobsColumn'].set(haveFobsColumn)
        self.__dict__['_value']['haveFpmObsColumn'].set(haveFpmObsColumn)
        self.__dict__['_value']['haveIobsColumn'].set(haveIobsColumn)
        self.__dict__['_value']['haveIpmObsColumn'].set(haveIpmObsColumn)
        self.__dict__['lastLoadedFile'] = fileName
        self.dataLoaded.emit()
        return err


class CMtzDataFile(CCP4File.CDataFile):
    '''An MTZ experimental data file'''
    QUALIFIERS = {'mimeTypeName' : 'application/CCP4-mtz',
                  'mimeTypeDescription' : 'MTZ experimental data',
                  'fileExtensions' : ['mtz'],
                  'fileContentClassName' : 'CMtzData',
                  'guiLabel' : 'Experimental data',
                  'toolTip' : "Experimental data in CCP4's MTZ format",
                  'sameCrystalAs' : NotImplemented,
                  'sameCrystalLevel' : NotImplemented,
                  'helpFile' : 'data_files#MTZ' }

    QUALIFIERS_DEFINITION = {'sameCrystalAs' : {'type' : str, 'description' : 'Name of CMtzDataFile object that crystal parameters should match - probably the observed data' },
                             'sameCrystalLevel' : {'type' : int , 'description' : 'Rigour of same crystal test'}}
    ERROR_CODES = {151 : {'description' :'Failed converting MTZ file to alternative format'},
                   152 : {'description' :'Failed merging MTZ file - invalid input'},
                   153 : {'description' :'Failed merging MTZ files - error running cmtzjoin - see log'},
                   154 : {'description' :'Failed merging MTZ files - error running cad - see log'},
                   401 : {'description' :'MTZ file header data differs'},
                   402 : {'description' :'MTZ file columns differ'},
                   403 : {'description' :'Error trying to access number of reflections' , 'severity' : SEVERITY_WARNING},
                   404 : {'description' :'MTZ files have different number of reflections'},
                   405 : {'description' :'MTZ column mean value differs'},
                   406 : {'description' :'MTZ file header data differs - may be autogenerated names', 'severity' : SEVERITY_WARNING},
                   407 : {'description' : 'Error splitting MTZ file - failed creating input command to cmtzsplit'},
                   408 : {'description' : 'Error splitting MTZ file - output file missing'}}

    # Seems to be necessary to have this __init__ defined otherwise fails to load demo task widget. This needs looking into
    def __init__(self, value={}, qualifiers={}, parent=None, name=None, fullPath=None,**kw):
        qualis = {}
        qualis.update(qualifiers)
        qualis.update(kw)
        CCP4File.CDataFile.__init__(self, value=value, qualifiers=qualis, parent=parent, name=name, fullPath=fullPath, keywords=kw)

    def validity(self, arg={}):
        err = CCP4File.CDataFile.validity(self, arg)
        if err.maxSeverity() > SEVERITY_WARNING or not self.isSet():
            return err
        otherMtz = self.getDataByKey('sameCrystalAs')
        if otherMtz is None or not otherMtz.isSet():
            return err
        testLevel = self.qualifiers('sameCrystalLevel')
        if testLevel is not NotImplemented and testLevel is not None:
            err.extend(self.sameCrystal(otherMtz, testLevel))
        else:
            err.extend(self.sameCrystal(otherMtz))
        return err

    def sameCrystal(self, other=None, testLevel=None):
        if self._fileContent is None:
            self.loadFile()
        if other._fileContent is None:
            other.loadFile()
        return self.getFileContent().sameCrystal(other._fileContent,testLevel=testLevel)

    def runMtz2various(self, hklin=None, labin=[], hklout=None, output='SHELX', keywords={}, **args):
        ''' Run mtz2various with self as the input file
        keywords are the keyword and text value of mtz2various keyword input excluding OUTPUT,LABIN, and END
        Returns the hklout(str),err  (CErrorReport) '''
        labelMapping = {'F' : 'FP' , 'SIGF' : 'SIGFP' , 'I' : 'I' , 'SIGI' : 'I',
                        'Fplus' : 'F(+)', 'SIGFplus' : 'SIGF(+)', 'Fminus' : 'F(-)', 'SIGFminus' : 'SIGF(-)',
                        'Iplus' : 'I(+)', 'SIGIplus' : 'SIGI(+)', 'Iminus' : 'I(-)', 'SIGIminus' : 'SIGI(-)',
                        'FREER' : 'FREE' }
        kw = {}
        kw.update(keywords)
        kw.update(args)
        error = CErrorReport()
        cbin = os.path.normpath(os.path.join( CCP4Utils.getOSDir(), 'bin', 'mtz2various' ))
        if not os.path.exists(cbin):
            cbin =  os.path.normpath(os.path.join( CCP4Utils.getCCP4Dir(), 'bin', 'mtz2various' ))
        if hklin is None:
            hklin=self.__str__()
        arglist = ['HKLIN', hklin, 'HKLOUT', hklout]
        comText = 'OUTPUT ' + output + '\nLABIN'
        for fileLabel in labin:
            progLabel = labelMapping.get(fileLabel,None)
            if progLabel is not None:
                comText = comText + ' ' + progLabel + '=' + fileLabel
        comText = comText +'\n'
        for key,value in list(kw.items()):
            if key.lower() not in ['end','output','labin']:
                comText = comText + key + ' ' + value + '\n'
        comText = comText + 'END\n'
        inputFile = os.path.normpath(os.path.splitext(hklout)[0] + '_mtz2various.com')
        CCP4Utils.saveFile(inputFile,comText)
        logFile = os.path.normpath(os.path.splitext(hklout)[0] + '_mtz2various.log')
        pid = CCP4Modules.PROCESSMANAGER().startProcess(cbin, arglist, logFile=logFile, inputFile=inputFile)
        status = CCP4Modules.PROCESSMANAGER().getJobData(pid)
        exitCode = CCP4Modules.PROCESSMANAGER().getJobData(pid, 'exitCode')
        if status == 0 and os.path.exists(hklout):
            return hklout, error
        else:
            error.append(self.__class__, 151, self.__str__())
            return None, error

    def runMtzjoin(self, outfile, infiles):
        error = CErrorReport()
        logFile = os.path.normpath(os.path.splitext(outfile)[0] + '_cmtzjoin.log')
        cbin = CCP4Utils.getCCP4Exe('cmtzjoin')
        arglist = ['-mtzout', outfile]
        try:
            for name,cols in infiles:
                arglist.append('-mtzin')
                arglist.append(name)
                if len(cols) > 0:
                    arglist.append('-colout')
                    arglist.append(cols)
        except:
            error.append(self.__class__, 152, str(name) + ' ' + str(cols))
        pid = CCP4Modules.PROCESSMANAGER().startProcess(cbin, arglist, logFile=logFile)
        status = CCP4Modules.PROCESSMANAGER().getJobData(pid)
        exitCode = CCP4Modules.PROCESSMANAGER().getJobData(pid, 'exitCode')
        if status not in [0, 101] or not os.path.exists(outfile):
            error.append(self.__class__, 153, logFile)
            outfile = None
        return outfile, error

    def runCad(self, hklout=None, hklinList=[], comLines=[]):
        ''' Run cad
        hklout - output file name
        hklinList - list of mtz to be merged with self.__str__()
        comLines - list of strings to go into command file (END not necessary)
        If comLines is zero length then the compulsary LABIN input will be set to 'ALL' for all input files
        Returns the hklout(str),err  (CErrorReport) '''
        error = CErrorReport()
        cbin = os.path.normpath(os.path.join(CCP4Utils.getCCP4Dir(), 'bin', 'cad'))
        arglist = ['HKLOUT', hklout, 'HKLIN1', self.__str__()]
        nHklin = 1
        for hklin in hklinList:
            nHklin += 1
            arglist.extend([ 'HKLIN'+str(nHklin) , hklin ] )
        print('CMtzDataFile.runCad', arglist)
        comText = 'SYSAB_KEEP\n'
        labinLines = []
        for line in comLines:
            words = line.upper().split()
            if words[0] == 'LABIN':
                labinLines.append(int(words[2]))
        for n in range(1, nHklin+1):
            if n not in labinLines:
                comText = comText + 'LABIN FILE_NUMBER ' + str(n) + ' ALL\n'
        for line in comLines:
            if line.lower() not in ['end']:
                comText = comText + line + '\n'
        comText = comText + 'END\n'
        print('CMtzDataFile.runCad', comText)
        inputFile = os.path.normpath(os.path.splitext(hklout)[0] + '_cad.com')
        CCP4Utils.saveFile(inputFile, comText)
        logFile =  os.path.normpath(os.path.splitext(hklout)[0] + '_cad.log')
        pid = CCP4Modules.PROCESSMANAGER().startProcess(cbin, arglist, logFile=logFile, inputFile=inputFile)
        status = CCP4Modules.PROCESSMANAGER().getJobData(pid)
        exitCode = CCP4Modules.PROCESSMANAGER().getJobData(pid,'exitCode')
        if status == 0 and os.path.exists(hklout):
            return hklout, error
        else:
            error.append(self.__class__, 154, self.__str__())
            return None, error

    def assertSame(self, other, diagnostic=False, **kw):
        report = CCP4File.CDataFile.assertSame(self, other, diagnostic=diagnostic, testChecksum=False, **kw)
        if report.maxSeverity() > SEVERITY_WARNING or report.count(code=315) > 0:
            return report
        try:
            self.loadFile()
        except:
            report.append(self.__class__, 311, self.__str__(), name=self.objectPath(False))
            return report
        try:
            other.loadFile()
        except:
            report.append(self.__class__, 312, other.__str__(), name=self.objectPath(False))
            return report
        for item in ['cell', 'spaceGroup']:
            if self.fileContent.__getattr__(item) != other.fileContent.__getattr__(item):
                report.append(self.__class__, 401,item + ' : ' + str( self.fileContent.__getattr__(item) ) + ' : ' + str(other.fileContent.__getattr__(item)), stack=False, name=self.objectPath(False) )
        for item in ['low', 'high']:
            if self.fileContent.resolutionRange.__getattr__(item) != other.fileContent.resolutionRange.__getattr__(item):
                lerrStr = item + ' : ' + str(self.fileContent.resolutionRange.__getattr__(item)) + ' : ' \
                                       + str(other.fileContent.resolutionRange.__getattr__(item))
                report.append(self.__class__, 401, lerrStr, stack=False, name=self.objectPath(False))
        ok = 0
        if len(self.fileContent.datasets) != len(other.fileContent.datasets):
            ok = 2
        else:
            for idx in range(len(self.fileContent.datasets)):
                selfDName = str(self.fileContent.datasets[idx])
                otherDName = str(other.fileContent.datasets[idx])
                if selfDName != otherDName and selfDName:
                    if selfDName.count('unknown') and otherDName.count('unknown'):
                        if ok == 0:
                            ok = 1
                    else:
                        ok = 2
        if ok == 2:
            lerrStr = 'datasets ' + ' : ' + str( self.fileContent.datasets ) + ' : ' + str(other.fileContent.datasets)
            report.append(self.__class__, 401, lerrStr, stack=False, name=self.objectPath(False) )
        elif ok == 1:
            lerrStr = 'datasets ' + ' : ' + str( self.fileContent.datasets ) + ' : ' + str(other.fileContent.datasets)
            report.append(self.__class__, 406, lerrStr, stack=False, name=self.objectPath(False) )
        if len(self.fileContent.listOfColumns) != len(other.fileContent.listOfColumns):
            report.append(self.__class__, 402, stack=False, name=self.objectPath(False))
        # Test the number of reflections
        myRefnList = self.hklfileReflectionList()
        otherRefnList = other.hklfileReflectionList()
        if myRefnList is None or otherRefnList is None:
            report.append(self.__class__, 403, stack=False, name=self.objectPath(False))
        else:
            if myRefnList.NumberReflections() != otherRefnList.NumberReflections():
                lerrStr = str(myRefnList.NumberReflections() )+' : '+str (otherRefnList.NumberReflections() )
                report.append(self.__class__, 404, lerrStr, stack=False , name=self.objectPath(False))
            else:
                myRefnList.ReadData()
                otherRefnList.ReadData()
                myStats = myRefnList.GetColumnStatistics()
                otherStats = otherRefnList.GetColumnStatistics()
                for icol in range(min(myRefnList.NumberColumns(), otherRefnList.NumberColumns()) - 3):
                    if diagnostic:
                        print('Reflection data mean values for ',self.fileContent.listOfColumns[icol].columnLabel,' mine:',myStats[icol].Mean(),'other:',otherStats[icol].Mean())
                    if not abs(otherStats[icol].Mean() - myStats[icol].Mean()) < max(abs(otherStats[icol].Mean()/100.0), 0.0001):
                        report.append(self.__class__, 405, str(self.fileContent.listOfColumns[icol].columnLabel) + ' ' +  str( myStats[icol].Mean())+ ' : ' + str(otherStats[icol].Mean()), name=self.objectPath(False) )
        if diagnostic:
            print('CMtzDataFile.assertSame NumberReflections', myRefnList.NumberReflections(), otherRefnList.NumberReflections())
        if len(report) == 0:
            report.append(self.__class__, 300, name=self.objectPath(False), stack=False)
        #print 'CMtzDataFile.assertSame done',report.report()
        return report

    def hklfileReflectionList(self):
        if not self.exists():
            return None
        try:
            import ccp4mg
            import hklfile
            return hklfile.ReflectionList(self.__str__())
        except:
            return None

    def runMtzSplit(self, colinList, outfiles, coloutList=None, logFile=None ):
        print('CMtzDataFile',self, outfiles)
        ret = CErrorReport()
        if logFile is None:
            logFile = CCP4Utils.makeTmpFile(cdir=False)
        cbin = os.path.join(CCP4Utils.getOSDir(), 'bin', 'cmtzsplit')
        if not os.path.exists(cbin):
            cbin = os.path.join(CCP4Utils.getCCP4Dir(), 'bin', 'cmtzsplit')
        arglist = ['-mtzin', str(self)]
        for ii in range(min(len(colinList), len(outfiles))):
            try:
                name = outfiles[ii]
                colin = colinList[ii]
                arglist.append('-mtzout')
                arglist.append(name)
                arglist.append('-colin')
                arglist.append(colin)
                arglist.append('-colout')
                if coloutList is not None:
                    arglist.append( coloutList[ii])
                else:
                    signature = self.fileContent.getSignatureForColumns(colin.split(','))
                    print('CMtzDataFile.runMtzSplit', colin, signature)
                    for cls in (CObsDataFile, CPhsDataFile, CMapCoeffsDataFile, CFreeRDataFile):
                        contentFlag = 0
                        for correctSig in cls.QUALIFIERS['correctColumns']:
                            contentFlag += 1
                            if signature == correctSig:
                                colout = cls.CONTENT_SIGNATURE_LIST[contentFlag-1][0]
                                for col in cls.CONTENT_SIGNATURE_LIST[contentFlag-1][1:]:
                                    colout = colout + ',' + col
                                arglist.append(colout)
                                break
            except:
                ret.append(self.__class__, 407, name + ' ' + colin)
        pid = CCP4Modules.PROCESSMANAGER().startProcess(cbin, arglist, logFile=logFile)
        status = CCP4Modules.PROCESSMANAGER().getJobData(pid)
        exitCode = CCP4Modules.PROCESSMANAGER().getJobData(pid, 'exitCode')
        if status == 0:
            for outfile in outfiles:
                if not os.path.exists(outfile[0]):
                    ret.append(self.__class__,408,  outfile)
        return ret


class CCrystalName(CCP4Data.CString):
    QUALIFIERS = {'allowUndefined' : False, 'minLength' : 1, 'allowedChars' : 1, 'toolTip' : 'Unique identifier for crystal (one word)'}


class CDatasetName(CCP4Data.CString):
    QUALIFIERS  = {'allowUndefined' : False, 'allowedChars' : 1, 'minLength' : 1, 'toolTip' : 'Unique identifier for dataset (one word)' }


def getClipperCell(clipperCell):
    cell = {}
    for item,roundTo in [['a', 4], ['b', 4], ['c', 4], ['alpha', 5], ['beta', 5], ['gamma', 5]]:
        value = getattr(clipperCell, item)().__float__()
        if math.isnan(value) or value < 0.0001:
            cell[item] = None
        elif value < 3.2 and item in ['alpha', 'beta', 'gamma']:
            cell[item] = round(value * ( 180.0 / math.pi),roundTo)
        else:
            cell[item] = round(value, roundTo)
    return cell


class CUnmergedDataContent(CCP4File.CDataFileContent):
    CONTENTS = {'format' : {'class' : CCP4Data.CString,
                            'qualifiers' : {'onlyEnumerators':True , 'enumerators' : ['unk','mtz' ,'xds', 'sca', 'saint', 'shelx', 'mmcif'], 'default' : 'unk'}},
                'merged' : {'class' : CCP4Data.CString,
                            'qualifiers' : {'onlyEnumerators':True , 'enumerators' : ['unk','merged' ,'unmerged'], 'default' : 'unk'}},
                'crystalName' : {'class' : CCrystalName}, 'datasetName' : {'class' : CDatasetName}, 'cell' : {'class' : CCell},
                'spaceGroup' : {'class' : CSpaceGroup}, 'batchs' : {'class' : CCP4Data.CString}, 'lowRes' : {'class' : CCP4Data.CFloat},
                'highRes' : {'class' : CCP4Data.CFloat}, 'knowncell' : {'class' : CCP4Data.CBoolean},
                'knownwavelength' : {'class' : CCP4Data.CBoolean}, 'numberLattices' : {'class' : CCP4Data.CInt},
                'wavelength' : {'class' : CWavelength},
                'numberofdatasets':{'class' : CCP4Data.CInt}}

    def loadFile(self,fileName=None):
        self.unSet()
        fileName = str(fileName)
        # hklfile.ReflectionFileType NOTSET=0 ABSENT=1 UNKNOWN=2 MTZ = 3 SCA_MERGED=4 SCA_UNMERGED=5
        # SHELX=6 SAINT=7 XDS_INTEGRATE+8 XDS_ASCII=9
        if fileName is None or not os.path.exists(fileName):
            return
        self.__dict__['_value']['knowncell'].set(True)
        self.__dict__['_value']['knownwavelength'].set(True)
        # get the file format
        try:
            import ccp4mg
            import hklfile
        except Exception as e:
            print('FAILED IMPORTING HKLFILE')
            print(e)
        self.format.set('unk')
        if os.path.splitext(fileName)[1] in ['.mmcif', '.cif', '.ent']:
            self.format.set('mmcif')
        else:
            reflectionList = hklfile.ReflectionList(fileName)
            ftype = reflectionList.FileType()
            if ftype.FileType() in [hklfile.ReflectionFileType.ABSENT, hklfile.ReflectionFileType.UNKNOWN]:
                raise CException(self.__class__, 101, fileName, name=self.objectPath())
            self.merged.set('unk')
            # there may be settings of unk needed for some file types
            if (reflectionList.Merged()):
                self.merged.set('merged')
            else:
                self.merged.set('unmerged')
            if ftype.FileType() in [hklfile.ReflectionFileType.SCA_MERGED]:
                self.merged = 'merged'
            if ftype.FileType() in [hklfile.ReflectionFileType.SCA_UNMERGED]:
                self.merged.set('unmerged')
        self.__dict__['_value']['knowncell'].set(True)
        if ftype.FileType() in [hklfile.ReflectionFileType.SCA_UNMERGED, hklfile.ReflectionFileType.SHELX]:
            self.__dict__['_value']['knowncell'].set(False)
            self.__dict__['_value']['knownwavelength'].set(False)
        if ftype.FileType() in [hklfile.ReflectionFileType.SCA_MERGED, hklfile.ReflectionFileType.SHELX]:
            self.__dict__['_value']['knownwavelength'].set(False)
        if ftype.FileType() in [hklfile.ReflectionFileType.SCA_UNMERGED, hklfile.ReflectionFileType.SCA_MERGED]:
            self.format.set('sca')
        if ftype.FileType() == hklfile.ReflectionFileType.SHELX:
            self.format.set('shelx')
        if ftype.FileType() in [hklfile.ReflectionFileType.XDS_INTEGRATE, hklfile.ReflectionFileType.XDS_ASCII]:
            self.format.set('xds')
            # assuming xds in unmerged may be faulty
            self.merged.set('unmerged')
        if ftype.FileType() == hklfile.ReflectionFileType.SAINT:
            self.format.set('saint')
            self.merged.set('unmerged')
        if ftype.FileType() == hklfile.ReflectionFileType.MTZ:
            self.format.set('mtz')
        self.extractMtzData(ftype, reflectionList)

    def extractMtzData(self, ftype, reflectionList):
        #This is broken without swigged clipper
        listOfDatasets = []
        listOfXnames = []
        listOfWavelengths = []
        listOfCells = []
        xtalDataInf = reflectionList.GetDatasets()
        firstrealdataset = False
        self.wavelength = float(0.0)
        for n in range(len(xtalDataInf)):
            listOfDatasets.append(xtalDataInf[n].Dname().__str__())
            listOfXnames.append(xtalDataInf[n].Xname().__str__())
            listOfWavelengths.append(xtalDataInf[n].Wavelength().__float__())
            listOfCells.append(getClipperCell(xtalDataInf[n].Cell()))
            if not firstrealdataset and (xtalDataInf[n].Dname() != 'HKL_base'):
                firstrealdataset = True
                self.datasetName.set(listOfDatasets[n])
                self.crystalName.set(listOfXnames[n])
                self.wavelength.set(self.wavelength.fix(listOfWavelengths[n]))
        print('extractMtzData listOfCells', listOfCells, "Ndatasets", len(listOfDatasets))
        # We should probably store other dataset information, but for now
        # just recognising multiple datasets is useful  Phil Evans June 2023
        self.numberofdatasets.set(len(listOfDatasets))
        self.__dict__['_value']['cell'].set(getClipperCell(reflectionList.Cell()))
        try:
            hm = reflectionList.Spacegroup().symbol_hm()
            if not isinstance(hm, str):
                import ccp4mg
                import hklfile   # KJS : put this in here for now to fix line below
                hm = hklfile.ClipperStringAsString(hm)
            self.__dict__['_value']['spaceGroup'].set(hm)
        except:
            self.__dict__['_value']['spaceGroup'].set(reflectionList.CCP4_SpacegroupNumber().__int__())
        self.__dict__['_value']['batchs'].set(reflectionList.formatBatchNumberlist())
        self.__dict__['_value']['lowRes'].set(reflectionList.Resolution().ResLow())
        self.__dict__['_value']['highRes'].set(reflectionList.Resolution().ResHigh())
        self.__dict__['_value']['numberLattices'].set(reflectionList.NumberLattices())


class CUnmergedDataFile(CCP4File.CDataFile):
    ''' Handle MTZ, XDS and scalepack files. Allow wildcard filename'''
    QUALIFIERS = {'mimeTypeName' : 'application/CCP4-unmerged-experimental', 'mimeTypeDescription' : 'Unmerged experimental data',
                  'fileExtensions' : ['mtz', 'hkl', 'HKL', 'sca', 'SCA', 'ent', 'cif'],
                  'fileContentClassName' : 'CUnmergedDataContent',
                  'guiLabel' : 'Unmerged reflections', 'toolTip' : "Unmerged experimental data in any format",
                  'helpFile' : 'data_files#unmerged_data' }

    def __init__(self, value={}, qualifiers={}, parent=None, name=None, fullPath=None, **kw):
        qualis = {}
        qualis.update(qualifiers)
        qualis.update(kw)
        CCP4File.CDataFile.__init__(self, value=value, qualifiers=qualis, parent=parent, name=name, fullPath=fullPath, keywords=kw)
        self.baseName.setQualifier('allowedCharacters', '*?')
        self.relPath.setQualifier('allowedCharacters', '*?')

    def isWildCard(self):
        if self.baseName.isSet() and (str(self.baseName).count('*') + str(self.baseName).count('?')) > 0:
            return True
        else:
            return False

    def globFiles(self):
        allFiles = glob.glob(self.fullPath().__str__())
        return allFiles

    def fileFormat(self):
        if not self.baseName.isSet():
            return None
        base, ext = os.path.splitext( self.baseName.__str__())
        if ext == '.mtz':
            return 'mtz'


class CUnmergedDataFileList(CCP4Data.CList):
    SUBITEM = {'class' : CUnmergedDataFile}


class CImportUnmerged(CCP4Data.CData):
    CONTENTS = {'file' : {'class' :CUnmergedDataFile , 'qualifiers' : { 'allowUndefined' : False, 'mustExist': True, 'fromPreviousJob' : True }},
                'cell' : {'class' :CCell}, 'wavelength' : {'class' : CWavelength},
                'crystalName' : {'class': CCP4Data.CString , 'qualifiers' : {'allowUndefined' : False, 'minLength' : 1, 'guiLabel' : 'crystal name', 'allowedCharsCode' : 1}},
                'dataset' : {'class': CCP4Data.CString, 'qualifiers' : {'allowUndefined' : False, 'minLength' : 1, 'guiLabel' : 'dataset name','allowedCharsCode' : 1}},
                'excludeSelection' : {'class': CCP4Data.CRangeSelection, 'qualifiers' : {'allowUndefined' : True}},}
    CONTENTS_ORDER = ['file', 'crystalName', 'dataset', 'excludeSelection']
    QUALIFIERS = {'toolTip' : 'Unmerged experimental data file and name for crystal and dataset in the file', 'helpFile' : 'import_merged#file_formats',
                  'toolTip' : 'Imported data file, cell parameters and crystal/dataset identifiers'}   # KJS : Duplicate entry ?

    def __init__(self, **kw):
        CCP4Data.CData.__init__(self, **kw)
        # These are calling loadDatasetName to overwrite valid dataset/crystalNames
        # try just calling it from the gui

    def validity(self, args):
        keys = ['file','crystalName','dataset','excludeSelection']
        if self.file.isSet() and not self.file.fileContent.knowncell:
            keys.append('cell')
        if self.file.isSet() and not self.file.fileContent.knownwavelength:
            keys.append('wavelength')
        return self.itemValidity(args, keys=keys)

    def aimlessExcludeBatch(self):
        #Return a string of 'EXCLUDE BATCH ..' lines
        if not self.__dict__['_value']['excludeSelection'].isSet():
            return ''
        ret = ''
        batchList = self.excludeSelection.removeWhiteSpace(str(self.__dict__['_value']['excludeSelection'])).split(",")
        listout = ''
        for batch in batchList:
            if "-" in batch:
                ret +=  "EXCLUDE BATCH %s\n" % batch.replace("-", " TO ")
            else:
                listout += " " + batch
        if len(listout) > 0:
            ret +=  "EXCLUDE BATCH %s\n" % listout
        return ret

    def loadDatasetName(self,signal=True):
        disallowed =  ['dummy', 'new', 'none']
        self.file.loadFile()
        if self.file.fileContent.datasetName.isSet() and len(self.file.fileContent.datasetName) > 0:
            name = self.file.fileContent.datasetName.__str__()
            if name.lower() not in disallowed:
                self.dataset.set(name)
        if self.file.fileContent.crystalName.isSet() and len(self.file.fileContent.crystalName) > 0:
            name = self.file.fileContent.crystalName.__str__()
            if name.lower() not in disallowed:
                self.crystalName.set(name)
        if signal: self.dataChanged.emit()

    def getTextItem(self):
        return self.file.baseName.__str__()

    def getTableTextItems(self):
        return [self.file.baseName.__str__(), self.crystalName.__str__(), self.dataset.__str__(), self.excludeSelection.__str__()]

    # This is needed for the CImportUnmergedView to handle icon menu
    def exists(self):
        return self.__dict__['_value']['file'].exists()


class CImportUnmergedList(CCP4Data.CList):
    QUALIFIERS = {'listMinLength' : 1}
    SUBITEM = {'class' : CImportUnmerged}

    def getCrystalList(self):
        xtalList = []
        for item in self.__dict__['_value']:
            if item.crystalName.isSet():
                xtal = item.crystalName.__str__()
                if xtalList.count(xtal) == 0:
                    xtalList.append(xtal)
        return xtalList

    def saveToDb(self):
        fileObjList = []
        for item in self.__dict__['_value']:
            fileObjList.append(item.file)
        return fileObjList, None, {}


class CImageFile(CCP4File.CDataFile):
    QUALIFIERS = {'mimeTypeName' : 'application/CCP4-image', 'mimeTypeDescription' : 'Image file',
                  'fileExtensions' : ['img', 'cbf', 'mccd', 'mar1600', 'h5', 'nxs'],
                  'fileContentClassName' : None, 'guiLabel' : 'Image file', 'toolTip' : "First image file in a directory"}

    # Override to avoid copying image
    def importFile(self, **kw):
        return


class CImageFileList(CCP4Data.CList):
    SUBITEM = {'class' : CImageFile}

class CXia2ImageSelection(CCP4Data.CData):
    CONTENTS = {'imageFile': {'class':CImageFile},
                'imageStart': {'class': CCP4Data.CInt, 'qualifiers': {'allowUndefined': True, 'min' : 0}},
                'imageEnd': {'class': CCP4Data.CInt, 'qualifiers': {'allowUndefined': True, 'min' : 0}}}
    CONTENTS_ORDER = ['imageFile', 'imageStart', 'imageEnd']
    QUALIFIERS = {'toolTip' : 'select an image file and an optional range of files to define a dataset'}

    def __init__(self, **kw):
        CCP4Data.CData.__init__(self, **kw)

    def validity(self, args):
        keys = ['imageFile','imageStart','imageEnd']
        return self.itemValidity(args, keys=keys)

    def getTextItem(self):
        return str(self.imageFile.baseName)

    def getTableTextItems(self):
        return [str(self.imageFile.baseName), str(self.imageStart), str(self.imageEnd)]

class CXia2ImageSelectionList(CCP4Data.CList):
    SUBITEM = {'class' : CXia2ImageSelection}


class CUnmergedMtzDataFile(CMtzDataFile):
    QUALIFIERS = {'mimeTypeName' : 'application/CCP4-mtz-unmerged', 'mimeTypeDescription' : 'MTZ unmerged experimental data',
                  'fileExtensions' : ['mtz'], 'fileContentClassName' : None, 'guiLabel' : 'Unmerged MTZ reflections',
                  'toolTip' : "Unmerged experimental data in CCP4's MTZ format"}


class CColumnType(CCP4Data.CString):
    '''A list of recognised MTZ column types'''
    QUALIFIERS = {'enumerators' : ['H', 'J', 'F', 'D', 'Q', 'G', 'L', 'K', 'M', 'E', 'P', 'W', 'A', 'B', 'Y', 'I', 'R'],
                  'onlyEnumerators' : True, 'default' : 'F'}
    DESCRIPTION = {'H' : 'index h,k,l', 'J' : 'intensity', 'F' : 'structure amplitude, F', 'D' : 'anomalous difference',
                   'Q' : 'standard deviation of J,F,D or other (but see L and M below)',
                   'G' : 'structure amplitude associated with one member of an hkl -h-k-l pair, F(+) or F(-)',
                   'L' : 'standard deviation of a column of type G',
                   'K' : 'intensity associated with one member of an hkl -h-k-l pair, I(+) or I(-)', 'M' : 'standard deviation of a column of type K',
                   'E' : 'structure amplitude divided by symmetry factor ("epsilon"). Normally scaled as well to give normalised structure factor',
                   'P' : 'phase angle in degrees', 'W' : 'weight (of some sort)', 'A' : 'phase probability coefficients (Hendrickson/Lattman)',
                   'B' : 'BATCH number', 'Y' : 'M/ISYM, packed partial/reject flag and symmetry number',
                   'I' : 'any other integer', 'R' : 'any other real'}

    def __init__(self, value=NotImplemented, qualifiers={}, parent=None, name=None, **kw):
        qualis = {}
        qualis.update(qualifiers)
        qualis.update(kw)
        CCP4Data.CString.__init__(self, value=value, qualifiers=qualis, parent=parent, name=name)

    def getDescription(self):
        return CColumnType.DESCRIPTION[self.__dict__['_value']]

    PROPERTIES = {'description' : {'fget' : getDescription}}

class CColumnTypeList(CCP4Data.CList):
    '''A list of acceptable MTZ column types'''
    SUBITEM = {'class' : CColumnType}
    '''
    def __init__(self,value=NotImplemented,qualifiers={},parent=None,name=None,subItem={},**kw):
        qualis = {}
        qualis.update(qualifiers)
        qualis.update(kw)
        CCP4Data.CList.__init__(self,value=value,qualifiers=qualis,parent=parent,name=name,subItem=subItem)
    '''

    def setEtree(self, element, checkValidity=True):
        txt = str(element.text)
        cList = txt.split(',')
        for item in cList:
            self.addItem()
            self.__dict__['_value'][-1].set(item)

    def getEtree(self, tag=None):
        if tag is None:
            tag = self.objectName()
        if tag is None or len(tag) == 0:
            tag = self.className()
        element = etree.Element(tag)
        txt = ''
        for item in self.__dict__['_value']:
            txt = txt + item.get() + ','
        if len(txt) > 0:
            txt = txt[0:-1]
        element.text = txt
        return element


class CMtzColumn(CCP4Data.CData):
    '''An MTZ column with column label and column type'''
    CONTENTS = {'columnLabel' : {'class' : CCP4Data.COneWord, 'qualifiers' : {'allowUndefined' : True}}, 'columnType' : {'class' : CColumnType},
                'dataset' : {'class' : CCP4Data.COneWord}, 'groupIndex' : { 'class' : CCP4Data.CInt}}

    def guiLabel(self):
        if self.__dict__['_value']['dataset'].isSet() and self.__dict__['_value']['dataset'] != 'HKL_base':
            return str(self.__dict__['_value']['dataset']) + '/' + str(self.__dict__['_value']['columnLabel'])
        else:
            return str(self.__dict__['_value']['columnLabel'])


class CMtzColumnGroupType(CColumnType):
    pass


class CMtzColumnGroup(CCP4Data.CData):
    CONTENTS = {'groupType' : {'class' : CMtzColumnGroupType }, 'columns' : {'class' : CCP4Data.CList, 'subItem' : {'class' : CMtzColumn}}}

    def guiLabel(self):
        text = ''
        for col in self.__dict__['_value']['columns']:
            text = text + str(col.guiLabel()) + ' '
        return text

class CMtzDataset(CCP4Data.CData):
    CONTENTS = {'name' : {'class' : CCP4Data.CString}, 'columnGroups' : {'class' : CCP4Data.CList, 'subItem' : {'class' : CMtzColumnGroup}}}
    pass


class CMtzData(CCP4File.CDataFileContent):
    dataLoaded = QtCore.Signal()
    '''Some of the data contents of an MTZ file'''
    CONTENTS = {'cell' : {'class' : CCell}, 'spaceGroup' : {'class' : CSpaceGroup},    # HM space group name
                'resolutionRange' : {'class' : CResolutionRange}, 'listOfColumns' : {'class' : CCP4Data.CList, 'subItem' : {'class' : CMtzColumn}},
                'datasets' : {'class' : CCP4Data.CList, 'subItem' : {'class' : CCP4Data.CString}},
                'crystalNames' : {'class' : CCP4Data.CList, 'subItem' : {'class' : CCP4Data.CString}},
                'wavelengths' : {'class' : CCP4Data.CList, 'subItem' : {'class' : CWavelength}},
                'datasetCells' : {'class' : CCP4Data.CList, 'subItem' : {'class' : CCell}}, 'merged' : {'class' : CCP4Data.CBoolean}}

    ERROR_CODES = {101 : {'description' : 'Attempting to load MTZ data from non-existant/broken file'},
                   102 : {'description' : 'Error creating command file for mtzdump'},
                   103 : {'description' : 'No log file found from mtzdump'},
                   104 : {'description' : 'Error reading log file from mtzdump'},
                   105 : {'severity': SEVERITY_WARNING, 'description' : 'Different spacegroup'},
                   106 : {'severity': SEVERITY_WARNING, 'description' : 'Different cell parameter'},
                   107 : {'severity': SEVERITY_WARNING, 'description' : 'Different cell parameters'},
                   108 : {'severity': SEVERITY_ERROR, 'description' : 'Different Laue group'},
                   109 : {'severity': SEVERITY_ERROR, 'description' : 'Different point group'},
                   410 : {'description' :'Invalid CSeqDataFile passed to matthewCoeff'},
                   411 : {'description' :'Failed to run matthewCoeff'}}

    def loadFile(self, fileName=None):
        if fileName is None:
            return
        self.__dict__['lastLoadedFile'] = str(fileName)
        if not os.path.exists(self.__dict__['lastLoadedFile']):
            self.unSet()
            raise CException(self.__class__, 101, fileName, name=self.objectPath())
        try:
            import ccp4mg
            import hklfile   # KJS : Perhaps re-write this.
            useHklfile = True
        except:
            useHklfile = False
        if useHklfile:
            self.extractMtzData(self.__dict__['lastLoadedFile'])
            self.dataLoaded.emit()
        else:
            try:
                inputFile = os.path.normpath(os.path.join(CCP4Utils.getTMP(), 'runMtzdump.script'))
                if not os.path.exists(inputFile):
                    CCP4Utils.saveFile(fileName=inputFile, text='HEADER\nEND\n')
            except:
                raise CException(self.__class__, 102, inputFile, name=self.objectPath())
            logFile = self.logFileName()
            pid = CCP4Modules.PROCESSMANAGER().startProcess(command='mtzdump', args=['HKLIN', self.__dict__['lastLoadedFile']],
                                                            inputFile = inputFile, logFile = logFile, waitForFinished=1000)
            try:
                self.scrapeMtzdumpLog(pid)
            except CException as e:
                print(e.report())
            except Exception as e:
                print('ERROR extracting data from MTZ file: ', fileName)

    def logFileName(self):
        baseName = os.path.splitext(os.path.basename(self.__dict__['lastLoadedFile']))[0]
        return  os.path.normpath(os.path.join(CCP4Utils.getTMP(), baseName + '_mtzdump.log'))

    def scrapeMtzdumpLog(self, processId, **kw):
        logFile = self.logFileName()
        if not os.path.exists(logFile):
            raise CException(self.__class__, 103, logFile,name=self.objectPath())
        logText = CCP4Utils.readFile(logFile)
        if logText == '':
            raise CException(self.__class__, 104, logFile,name=self.objectPath())
        fileData = self.parseMtzdumpLog(logText)
        rv = CErrorReport()
        for item in self.dataOrder():
            try:
                self.__dict__['_value'][item].set(fileData[item])
            except CException as e:
                rv.append(e)
            except:
                pass
        self.dataLoaded.emit()
        return rv

    def parseMtzdumpLog(self, logText=''):
        from core.CCP4Utils import safeFloat
        # Extract data from log file. Code taken from EDNA example
        pyListLogLines = logText.split("\n")
        cell = []
        listOfColumns = []
        column_name_list = []
        column_type_list = []
        pyStrSpaceGroupName = ''
        iSpaceGroupNumber = -1
        lowerResolutionLimit = None
        upperResolutionLimit = None
        for j, pyStrLine in enumerate(pyListLogLines):
            if "* Dataset ID, project/crystal/dataset names, cell dimensions, wavelength:" in pyStrLine:
                cell = list(map(safeFloat, pyListLogLines[j + 5].split()))
            if " * Space group = " in pyStrLine:
                pyStrSpaceGroupName = pyStrLine.split("'")[1].strip()
                iSpaceGroupNumber = int(pyStrLine.replace("(", " ").replace(")", " ").split()[-1])
            if "*  Resolution Range" in pyStrLine:
                lowerResolutionLimit = float(((pyListLogLines[j + 2].split("(")[1]).split())[0])
                upperResolutionLimit = float(((pyListLogLines[j + 2].split("(")[1]).split())[2])
            if "* Column Labels" in pyStrLine:
                column_name_list = pyListLogLines[j + 2].split()
            if "* Column Types" in pyStrLine:
                column_type_list = pyListLogLines[j + 2].split()
        for j, column_name in enumerate(column_name_list):
            column_type = column_type_list[j]
            listOfColumns.append(CMtzColumn(columnLabel=column_name, columnType=column_type))
        rv = {}
        if len(cell) == 6:
            rv['cell'] = {'a' : cell[0], 'b' : cell[1], 'c' : cell[2], 'alpha' : cell[3], 'beta' : cell[4], 'gamma' : cell[5]}
        rv['spaceGroup'] = pyStrSpaceGroupName
        rv['resolutionRange'] = {'high' : lowerResolutionLimit, 'low' : upperResolutionLimit}
        rv['listOfColumns'] = listOfColumns[3:]
        return rv

    def getNColumns(self):
        return len(self.__dict__['_value']['listOfColumns'])

    def getColumn(self, guiLabel):
        if guiLabel is None:
            return None
        columnLabel = guiLabel.split('/')[-1]
        for column in self.__dict__['_value']['listOfColumns']:
            if column.columnLabel == columnLabel:
                return column
        return None

    def getColumnIndex(self, guiLabel):
        if guiLabel is None:
            return -1
        columnLabel = guiLabel.split('/')[-1]
        ii = -1
        for column in self.__dict__['_value']['listOfColumns']:
            ii = ii + 1
            if column.columnLabel == columnLabel:
                return ii
        return -1

    def getListOfWavelengths(self):
        #MN CHECK ME CHECKME expose list of wavelengths of datasets in an MTZ file
        return self.__dict__['_value']['wavelengths']

    def getListOfColumns(self,columnTypes=[]):
        if len(columnTypes) == 0:
            return self.__dict__['_value']['listOfColumns']
        else:
            listOfColumns = []
            for column in self.__dict__['_value']['listOfColumns']:
                if columnTypes.count(column.columnType):
                    listOfColumns.append(column)
            return listOfColumns

    def getColumnType(self, guiLabel=''):
        if guiLabel is None:
            return ''
        columnLabel = guiLabel.split('/')[-1]
        for column in self.__dict__['_value']['listOfColumns']:
            if column.colmnLabel == columnLabel:
                return column.columnType
        return ''

    def getPartnerColumn(self, firstColumn='', columnGroupItem=None):
        # NB returns a CMtzColumn
        partnerTypes = columnGroupItem.columnType
        if not columnGroupItem.partnerOffset.isSet():
            partnerOffset = 1
        else:
            partnerOffset = int(columnGroupItem.partnerOffset)
        columnIndex = self.getColumnIndex(firstColumn)
        if columnIndex < 0:
            return None
        columnIndex = columnIndex + partnerOffset
        if columnIndex >= len(self.__dict__['_value']['listOfColumns']):
            return None
        if partnerTypes.count(self.__dict__['_value']['listOfColumns'][columnIndex].columnType):
            return self.__dict__['_value']['listOfColumns'][columnIndex].guiLabel()
        else:
            return None

    def extractMtzData(self, fileName):
        try:
            import ccp4mg
            import hklfile
        except:
            print('FAILED IMPORTING HKLFILE')
        reflectionList = hklfile.ReflectionList(fileName)
        ftype = reflectionList.FileType()
        if ftype.FileType() in [hklfile.ReflectionFileType.ABSENT, hklfile.ReflectionFileType.UNKNOWN]:
            raise CException(self.__class__, 101, fileName, name=self.objectPath())
        #This is broken without swigged clipper
        xtalDataInf = reflectionList.GetDatasets()
        listOfDatasets = []
        listOfCrystals = []
        listOfWavelengths = []
        listOfCells = []
        for n in range(len(xtalDataInf)):
            listOfDatasets.append(xtalDataInf[n].Dname().__str__())
            listOfCrystals.append(xtalDataInf[n].Xname().__str__())
            listOfWavelengths.append(xtalDataInf[n].Wavelength().__str__())
            listOfCells.append(getClipperCell(xtalDataInf[n].Cell()))
        self.__dict__['_value']['datasets'].set(listOfDatasets)
        self.__dict__['_value']['crystalNames'].set(listOfCrystals)
        self.__dict__['_value']['datasetCells'].set(listOfCells)
        self.__dict__['_value']['wavelengths'].unSet()
        for item in listOfWavelengths:
            try:
                wl = float(item)
            except:
                wl = None
            #MN I sem to have an MTZ for which one of the wavelengths is -1.0, i.e. outside the allowed range
            #To avoid the exception being thrown, I will try except this
            try:
                self.__dict__['_value']['wavelengths'].append(wl)
            except CException as e:
                print(e)
        nCol = reflectionList.NumberColumns()
        colLabels = reflectionList.ColumnLabels()
        colTypes = reflectionList.ColumnTypes()
        groupIndex = reflectionList.ColumnGroupIndex()
        dataSetIndex = reflectionList.GetColumnDataset()
        listOfColumns = []
        for n in range(3, nCol):
            listOfColumns.append(CMtzColumn(columnLabel=colLabels[n], columnType=colTypes[n], groupIndex=groupIndex[n]))
        self.__dict__['_value']['listOfColumns'].set(listOfColumns)
        for n in range(nCol - 3):
            try:
                self.__dict__['_value']['listOfColumns'][n].dataset.set(self.__dict__['_value']['datasets'][dataSetIndex[n+3]])
            except:
                print('ERROR: extractMtzData listOfColumns', self.__dict__['_value']['listOfColumns'], 'nCol', nCol, 'n', n)
        self.__dict__['_value']['merged'].set(reflectionList.Merged())
        cell = reflectionList.Cell()
        for item in ['a', 'b', 'c', 'alpha', 'beta', 'gamma']:
            value = getattr(cell,item)().__float__()
            if not math.isnan ( value ) and value > 0.0001:
                if item in ['alpha', 'beta', 'gamma'] and value < 3.0:
                    getattr(self.__dict__['_value']['cell'],item).setRadians(value)
                else:
                    getattr(self.__dict__['_value']['cell'],item).set(value)
        try:
            hm = reflectionList.Spacegroup().symbol_hm()
            if not isinstance(hm,str):
                hm = hklfile.ClipperStringAsString(hm)
            self.__dict__['_value']['spaceGroup'].set(hm)
        except:
            self.__dict__['_value']['spaceGroup'].set(reflectionList.CCP4_SpacegroupNumber().__int__())
        self.__dict__['_value']['resolutionRange'].low.set(reflectionList.Resolution().ResLow().__str__())
        self.__dict__['_value']['resolutionRange'].high.set(reflectionList.Resolution().ResHigh().__str__())
        return listOfDatasets, listOfColumns

    def columnSignature(self):
        sig = ''
        for column in self.__dict__['_value']['listOfColumns']:
            sig += str(column.columnType)
        return sig

    def columnNames(self, indx=0, nColumns=None):
        nameList = []
        if nColumns is None:
            nColumns = len(self.__dict__['_value']['listOfColumns'])
        for ic in range(indx, indx + nColumns):
            nameList.append(str(self.__dict__['_value']['listOfColumns'][ic].columnLabel))
        return nameList

    def getSignatureForColumns(self, columnNameList):
        sig = ''
        nColumns = len(self.__dict__['_value']['listOfColumns'])
        for col in columnNameList:
            for ic in range(nColumns):
                if col == self.__dict__['_value']['listOfColumns'][ic].columnLabel:
                    sig += str(self.__dict__['_value']['listOfColumns'][ic].columnType)
                    break
        return sig

    def sameCrystal(self,other=None,testLevel=None):
        # testLevel = 0 no spacegroup test
        # testLevel = 1 same point group
        # testLevel = 2 same Laue group
        # testLevel = 3 same space group
        # testLevel = 4 same space group and cell
        if testLevel is None:
            testLevel = 4
        err= CErrorReport()
        if other is None:
            return err
        details = '(' + str(self.spaceGroup) + ') to ' + str(other.parent().qualifiers('guiLabel')) + ' (' + str(other.spaceGroup) + ')'
        path = str(self.parent().qualifiers('guiLabel'))
        if not self.__dict__['_value']['spaceGroup'].isSet() or not other.__dict__['_value']['spaceGroup'].isSet():
            return err
        if testLevel == 1:
            # test same point group
            same = False
            SG = self.__dict__['_value']['spaceGroup'].__str__()
            SGother = other.__dict__['_value']['spaceGroup'].__str__()
            if SG == SGother:
                same = True
            else:
                for name, sGList in list(SYMMETRYMANAGER().chiralSpaceGroups.items()):
                    if self.__dict__['_value']['spaceGroup'].__str__() in sGList:
                        if other.__dict__['_value']['spaceGroup'].__str__() in sGList:
                            same = True
                        break
                if not same:
                    err.append(self.__class__, 109, details, name=path, stack=False)
        elif testLevel == 2:
            # test same Laue group
            same=False
            for lG in SYMMETRYMANAGER().laueGroups:
                if  self.__dict__['_value']['spaceGroup'].__str__() in lG:
                    if  other.__dict__['_value']['spaceGroup'].__str__() in lG:
                        same = True
                    break
            if not same:
                err.append(self.__class__, 108, details, name=path, stack=False)
        elif testLevel in [3, 4]:
            # test same space group [ and cell]
            if self.__dict__['_value']['spaceGroup'] != other.__dict__['_value']['spaceGroup']:
                err.append(self.__class__, 105, details, name=path, stack=False)
            if testLevel == 4:
                for item in ['a', 'b', 'c', 'alpha', 'beta', 'gamma']:
                    dif = (float(self.__dict__['_value']['cell'].__dict__['_value'][item]) - float(other.__dict__['_value']['cell'].__dict__['_value'][item])) / float(self.__dict__['_value']['cell'].__dict__['_value'][item])
                    if dif < -1.0 or dif > 1.0:
                        err.append(self.__class__, 106, details, name=path, stack=False)
                        break
        return err

    def clipperSameCellCoor(self, other, tolerance=1.0):
        import clipper
        cell = other.mmdbManager.GetCell()
        otherA = float(cell[1])
        otherB = float(cell[2])
        otherC = float(cell[3])
        otherAlpha = float(cell[4])
        otherBeta = float(cell[5])
        otherGamma = float(cell[6])
        myCell = clipper.Cell(clipper.Cell_descr(float(self.cell.a), float(self.cell.b), float(self.cell.c),
                                                 float(self.cell.alpha), float(self.cell.beta), float(self.cell.gamma)))
        otherCell = clipper.Cell(clipper.Cell_descr(float(otherA), float(otherB), float(otherC),
                                                    float(otherAlpha), float(otherBeta), float(otherGamma)))
        recipdifference = self.clipperRecipCellDifference(myCell, otherCell)
        difference = self.clipperCellDifference(myCell, otherCell)
        equals = myCell.equals(otherCell, tolerance)  # Boolean
        """ Return dictionary of:
          'validity'           True if cells are simlar within resolution of tolerance
          'maximumResolution1' maximum allowed resolution in cell1
          'maximumResolution2' maximum allowed resolution in cell2
          'difference'  average cell difference in A
          'tolerance'   in A """
        result = {'validity': equals, 'maximumResolution1': recipdifference[0],
                  'maximumResolution2': recipdifference[1], 'difference': difference, 'tolerance': tolerance}
        return result

    def clipperSameCell(self, other, tolerance=None):
        import clipper
        if tolerance is None:
            myHigh = self.resolutionRange.high.get()
            otherHigh = other.resolutionRange.high.get()
            if myHigh is None:
                if otherHigh is None:
                    tolerance = 1.0
                else:
                    tolerance = 1.5 * otherHigh
            elif otherHigh is None:
                tolerance = 1.5 * myHigh
            else:
                tolerance = 1.5 * max(myHigh, otherHigh)
        myCell = clipper.Cell(clipper.Cell_descr(float(self.cell.a), float(self.cell.b), float(self.cell.c),
                                                 float(self.cell.alpha), float(self.cell.beta), float(self.cell.gamma)))
        otherCell = clipper.Cell(clipper.Cell_descr(float(other.cell.a), float(other.cell.b), float(other.cell.c),
                                                    float(other.cell.alpha), float(other.cell.beta), float(other.cell.gamma)))
        recipdifference = self.clipperRecipCellDifference(myCell, otherCell)
        difference = self.clipperCellDifference(myCell, otherCell)
        equals = myCell.equals(otherCell, tolerance)  # Boolean
        """ Return dictionary of:
          'validity'           True if cells are simlar within resolution of tolerance
          'maximumResolution1' maximum allowed resolution in cell1
          'maximumResolution2' maximum allowed resolution in cell2
          'difference'  average cell difference in A
          'tolerance'   in A """
        result = {'validity': equals, 'maximumResolution1': recipdifference[0],
                  'maximumResolution2': recipdifference[1], 'difference': difference, 'tolerance': tolerance}
        return result

    def clipperRecipCellDifference(self, cell1, cell2):
        fracmat1 = cell1.matrix_frac()
        fracmat2 = cell2.matrix_frac()
        s = 0.0
        for j in range(3):
            for i in range(3):
                s += (fracmat1[i][j] - fracmat2[i][j])**2
        s = math.sqrt(s)
        volume1 = cell1.volume()
        volume2 = cell2.volume()
        rv = s * math.pow(volume1, 0.66666666667), s * math.pow(volume2, 0.66666666667)
        return rv

    def clipperCellDifference(self, cell1, cell2):
        "Get resolution at which cell2 is different from cell1"
        orthmat1 = cell1.matrix_orth()
        orthmat2 = cell2.matrix_orth()
        s = 0.0
        for j in range(3):
            for i in range(3):
                s += (orthmat1[i][j] - orthmat2[i][j])**2
        return math.sqrt(s)

    def getColumnGroups(self):
        # Notes: This function builds a list containing mtz column info (for mini-mtz's)
        groupIndex = 0 # MapCoeffs broken KJS
        groupList = []
        # Sort listOfColumns into a list of grouped columns
        for ii in range(len(self.__dict__['_value']['listOfColumns'])):
            fileColumn = self.__dict__['_value']['listOfColumns'][ii]
            if fileColumn.groupIndex != groupIndex:
                groupList.append(CColumnGroup(dataset=str(fileColumn.dataset)))
                groupIndex = int(fileColumn.groupIndex)
            groupList[-1].columnList.append(CMtzColumn(columnLabel=fileColumn.columnLabel, columnType=fileColumn.columnType,
                                                       groupIndex=fileColumn.groupIndex, dataset=fileColumn.dataset))
        ## loop over and catch map-coefs.
        cols = self.__dict__['_value']['listOfColumns']
        for i, col1 in enumerate(cols[:-1]):
            col2 = cols[i+1]
            MapT = (col1.get('columnType') == 'F' and col2.get('columnType') == 'P')  # Map (FP): F-Phi / Phases (PW): Phi-FOM
            if MapT:
                groupList.append(CColumnGroup(dataset=str(col1.dataset)))
                groupList[-1].columnList.append(CMtzColumn(columnLabel=col1.columnLabel, columnType=col1.columnType,
                                                            groupIndex=col1.groupIndex, dataset=col1.dataset))
                groupList[-1].columnList.append(CMtzColumn(columnLabel=col2.columnLabel, columnType=col2.columnType,
                                                            groupIndex=col2.groupIndex, dataset=col2.dataset))
        # Assign the columnGroupType and contentFlag
        for group in groupList:
            signature = ''
            labels = []
            for col in group.columnList:
                signature = signature + str(col.columnType)
                labels.append(str(col.columnLabel))
            for cls, label in [[CObsDataFile, 'Obs'], [CPhsDataFile, 'Phs'], [CMapCoeffsDataFile, 'MapCoeffs'], [CFreeRDataFile, 'FreeR']]:
                if cls.QUALIFIERS['correctColumns'].count(signature):
                    #print 'correctColumns', cls.QUALIFIERS['correctColumns']
                    if label == 'FreeR':
                        # Is this really FreeR? check column label
                        if (len(labels) == 1) and ('free' in str(labels[0]).lower()):
                            print("Column recognised as FreeR, label:", str(labels[0]))
                            group.columnGroupType.set(label)
                            group.contentFlag.set(cls.QUALIFIERS['correctColumns'].index(signature) + 1)
                    else:
                        group.columnGroupType.set(label)
                        group.contentFlag.set(cls.QUALIFIERS['correctColumns'].index(signature) + 1)
            if not group.columnGroupType.isSet():
                if signature == 'FPW':
                    group.columnGroupType.set('Phs')
                    group.contentFlag.set(2)
                    group.columnList.pop(0)
                elif signature == 'PWF':
                    group.columnGroupType.set('Phs')
                    group.contentFlag.set(2)
                    group.columnList.pop(2)
                elif signature == 'FQDQ':
                    group.columnGroupType.set('Obs')
                    group.columnList.pop(3)
                    group.columnList.pop(2)
                    group.contentFlag.set(4)
        return groupList

    def matthewsCoeff(self, seqDataFile=None, nRes=None, molWt=None, polymerMode=""):
        if seqDataFile is not None:
            try:
                molWt = seqDataFile.fileContent.getAnalysis('molecularWeight')
            except:
                molWt = 0.0
        elif nRes is not None:
            #  Estimated residue wt as per ccp4 matthews_coeff documentation
            molWt = 112.5 * float(nRes)
        if molWt < 0.01:
            raise CException(self.__class__, 410, str(seqDataFile))
        # temporary log and xml files
        import tempfile
        f1 = tempfile.mkstemp()
        os.close(f1[0])
        f2 = tempfile.mkstemp()
        os.close(f2[0])
        comText = 'MOLWEIGHT ' + str(molWt) + '\nCELL'
        for p in ['a', 'b', 'c']:
            comText =  comText + ' ' + str(self.cell.get(p))
        for p in ['alpha','beta','gamma']:
            a = float(self.cell.get(p))
            if a < 3.0:
                a = a * 180. / math.pi
            comText = comText + ' ' + str(a)
        comText = comText + '\nSYMM ' + str(self.spaceGroup.number()) + '\n'
        comText = comText + 'XMLO\nAUTO\n'
        if polymerMode != "":
            comText =  comText + '\nMODE ' + polymerMode
        argList = ['XMLFILE' , f1[1]]
        pid = CCP4Modules.PROCESSMANAGER().startProcess('matthews_coef', argList, logFile=f2[1], inputText=comText)
        if not os.path.exists(f1[1]):
            raise CException(self.__class__,411,str(seqDataFile))
        # Get results from xml file
        rv = {'results' : []}
        xTree = CCP4Utils.openFileToEtree(fileName=f1[1])
        try:
            rv['cell_volume'] = float( xTree.xpath('cell')[0].get('volume'))
        except:
            pass
        xResultList = xTree.xpath('result')
        for xResult in xResultList:
            rv['results'].append({'nmol_in_asu' : int(xResult.get('nmol_in_asu'))})
            for item in ['matth_coef','percent_solvent','prob_matth']:
                rv['results'][-1][item] = float(xResult.get(item))
        return rv

    def getWavelength(self):
        for item in self.__dict__['_value']['wavelengths']:
            if item.isSet() and item > 0.01:
                return item
        return self.__dict__['_value']['wavelengths'][0]

    PROPERTIES  = {'wavelength' : {'fget' : getWavelength}}


class CColumnGroup(CCP4Data.CData):
    """Groups of columns in MTZ - probably from analysis by hklfile"""
    CONTENTS = {'columnGroupType' : {'class' : CCP4Data.COneWord,
                                     'qualifiers' : {'onlyEnumerators' : True,
                                                     'enumerators' : ['Obs', 'Phs', 'MapCoeffs', 'FreeR']}},
                'contentFlag' : {'class' : CCP4Data.CInt},
                'dataset' :  {'class' : CCP4Data.CString},
                'columnList' : {'class' : CCP4Data.CList, 'subItem' : {'class' : CMtzColumn}},
                'selected' : {'class' : CCP4Data.CBoolean}}

    def columnListStr(self, withTypes=True, splitter=','):
        columnList = self.__dict__['_value']['columnList']
        if len(columnList) > 0:
            label = str(columnList[0].columnLabel)
            if withTypes:
                label = label + '(' + str(columnList[0].columnType) + ')'
            for col in columnList[1:]:
                label = label + splitter + str(col.columnLabel)
                if withTypes:
                    label = label + '(' + str(col.columnType) + ')'
            return label
        else:
            return ''


class CColumnGroupList(CCP4Data.CList):
    SUBITEM = {'class' : CColumnGroup}


# This is a simple class to hold items of the requiredColumnsList qualifier for CMtzColumnGroup
# *****  Needs some validity checking **************
class CColumnGroupItem(CCP4Data.CData):
    """Definition of set of columns that form a 'group'"""
    # defaultList is string that is a comma separated list
    CONTENTS =  {'columnName' : {'class' : CCP4Data.COneWord},
                 'defaultList' : {'class' : CCP4Data.CString},
                 'columnType' : {'class' : CColumnTypeList},
                 'partnerTo' : {'class' : CCP4Data.COneWord},
                 'partnerOffset' : {'class' : CCP4Data.CInt}}

    ERROR_CODES = {1 : {'description' : 'Attempting to change immutable object'},
                   2 : {'description' : 'Attempting to access unknown attribute'}}

    def getEtree(self):
        ele = etree.Element('columnGroupItem')
        ele.set('id',str(self.__dict__['_value']['columnName']))
        ele.append(self.__dict__['_value']['columnType'].getEtree())
        if self.__dict__['_value']['partnerTo'].get() is not None:
            rcData = etree.Element('partnerTo')
            rcData.text = str(self.__dict__['_value']['partnerTo'])
            ele.append(rcData)
            if self.__dict__['_value']['partnerOffset'].get() is not None:
                rcOffset = etree.Element('partnerOffset')
                rcOffset.text = str(self.__dict__['_value']['partnerOffset'])
                ele.append(rcOffset)
        if self.__dict__['_value']['defaultList'].get() is not None:
            e = self.__dict__['_value']['defaultList'].getEtree()
            ele.append(e)
        return ele

    def setEtree(self, element, checkValidity=True):
        self.__dict__['_value']['columnName'].set(str(element.get('id')))
        for ele in element.iterchildren():
            name = str(ele.tag)
            if name in self.__dict__['_value']:
                self.__dict__['_value'][name].setEtree(ele)


class CProgramColumnGroup0(CCP4Data.CData):

    CONTENTS = {'columnGroup' : {'class' : CMtzColumnGroup},
                'datasetName' : {'class' : CCP4Data.CString}}
    QUALIFIERS = {'mustExist' : False, 'mtzFileKey' : '', 'groupTypes' : []}
    QUALIFIERS_ORDER = ['groupTypes','mtzFileKey','mustExist']
    QUALIFIERS_DEFINITION = {'groupTypes' : {'type' : list,
                                             'description' : 'Type of columnGroup required by program'},
                             'mtzFileKey' : {'type' :str,
                                             'description' : 'The key for a CMtxDataFile in the same CContainer'},
                             'mustExist' : {'type' : bool,
                                            'description' : 'Flag if the parameter must be set at run time'}}
    ERROR_CODES = {101 : {'description' : 'Column not in MTZ file'},
                   102 : {'description' : 'Column wrong type'},
                   103 : {'description' : 'MTZ file is not defined', 'severity' : SEVERITY_WARNING},
                   104 : {'description' : 'No column group selected'},
                   105 : {'description' : 'No column group selected', 'severity' : SEVERITY_WARNING}}

    def __getattr__(self,name):
        # This method necessary to get reference to columns (used most in scripts) to work
        # e.g. c.inputData.F_SIGF.FP
        if 'columnGroup' in self.__dict__['_qualifiers'] and 'columnGroup' in self.__dict__['_value']:
            ii = -1
            for column in self.__dict__['_qualifiers']['columnGroup']:
                ii = ii + 1
                if str(column.columnName) == name:
                    if len( self.__dict__['_value']['columnGroup'].columns) <= ii:
                        return None
                    else:
                        return self.__dict__['_value']['columnGroup'].columns[ii].columnLabel.__str__()
        return CCP4Data.CData.__getattr__(self, name)

    def __setattr__(self,name,value):
        if 'columnGroup' in self.__dict__['_qualifiers'] and 'columnGroup' in self.__dict__['_value']:
            ii = -1
            for column in self.__dict__['_qualifiers']['columnGroup']:
                ii = ii + 1
                if str(column.columnName) == name:
                    if len( self.__dict__['_value']['columnGroup'].columns)<=ii:
                        pass
                    else:
                        return self.__dict__['_value']['columnGroup'].columns[ii].set({'columnLabel' : value})
        return CCP4Data.CData.__setattr__(self,name,value)

    def getMtzData(self):
        mtz = self.getDataByKey('mtzFileKey')
        if mtz is None:
            return None
        return mtz.fileContent

    def setColumnGroup(self, datasetIndex=None, columnGroupIndex=None):
        mtzData = self.getMtzData()
        if mtzData is None:
            return
        try:
            cg = mtzData.datasets[datasetIndex].columnGroups[columnGroupIndex]
        except:
            print('ERROR attemting to get dataset', datasetIndex, 'and columngroup', columnGroupIndex, 'from mtzdata')
            return
        self.__dict__['_value']['columnGroup'].set(cg)
        self.__dict__['_value']['datasetName'].set(mtzData.datasets[datasetIndex].name.__str__())

    def validity(self, value ={}):
        v = CErrorReport()
        if len(value) == 0 or 'columnGroup' not in value or value['columnGroup'].isSet():
            if  not self.qualifiers('allowUndefined'):
                v.append(self.__class__,104, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            else:
                v.append(self.__class__,105, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            return v
        if value['columnGroup'].type not in self.qualifers('types'):
            v.append(self.__class__, 102, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        mtzData = self.getMtxData()
        if mtzData is None:
            v.append(self.__class__, 103, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        else:
            group = mtzData.fileContent.getColumnGroup(value)    # KJS : Changed mtz to mtzData
            if group is None:
                v.append(self.__class__,103, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        return v

    def setColumnGroupQualifier(self, columnGroup=[], parent=None, initialise=False):
        # Expecting columnGroup to be a list of dicts or CColumnGroupItem
        # Each dict should have attributes for CColumnGroupItem object
        if '_qualifiers' not in self.__dict__:
            self.__dict__['_qualifiers'] = {}
        if initialise or ('columnGroup' not in self.__dict__['_qualifiers']):
            self.__dict__['_qualifiers']['columnGroup'] = []
        indx = -1
        for item in columnGroup:
            indx = indx + 1
            if isinstance(item, CColumnGroupItem):
                itemObj = item
            else:
                itemObj = CColumnGroupItem(item, parent=parent)
            if indx > 0:
                if itemObj.partnerTo is None:
                    itemObj.partnerTo = str(self.__dict__['_qualifiers']['columnGroup'][0].columnName)
                if itemObj.partnerOffset is None:
                    itemObj.partnerOffset = indx
            self.__dict__['_qualifiers']['columnGroup'].append(itemObj)

    def qualifiers(self, name=None, default=True, custom=True, contentQualifiers=False):
        return CCP4Data.CData.qualifiers(self, name=name, default=default, custom=custom, contentQualifiers=False)

    def setQualifiers(self, qualifiers={}, **kw):
        qualis = {}
        qualis.update(qualifiers)
        if 'columnGroup' in qualis:
            self.setColumnGroupQualifier(qualis['columnGroup'])
            del qualis['columnGroup']
        CCP4Data.CData.setQualifiers(self, qualis)

    def setQualifiersEtree(self, element=None):
        # Custom parsing etree to extract the columnGroup info
        rv = CErrorReport()
        self.__dict__['_qualifiers']['columnGroup'] = []
        if element is not None:
            cGele = element.find('columnGroup')
            if cGele is not None:
                for ele in cGele.iterchildren():
                    try:
                        if str(ele.tag) == 'columnGroupItem':
                            item = CColumnGroupItem(parent=self)
                            item.setEtree(ele)
                            self.__dict__['_qualifiers']['columnGroup'].append(item)
                    except CException as e:
                        rv.extend(e)
                    except:
                        rv.append(self.__class__, 107, name=self.objectName())
                # Remove the columnGroup element from the etree
                # and parse the rest with the base class method
                element.remove(cGele)
            rv.extend(CCP4Data.CData.setQualifiersEtree(self, element))
            return rv

    def qualifiersEtree(self, customOnly=True, tag=None):
        error = CErrorReport()
        if tag is None:
            tag = self.objectName()
        if tag is None or len(tag) == 0:
            tag = self.__class__.__name__
            error.append(CCP4Data.CData, 14)
        root, err = CCP4Data.CData.qualifiersEtree(self, customOnly=customOnly, tag=tag)
        error.extend(err)
        # Remove the columnGroup element which needs doing differently
        try:
            ele = root.find('columnGroup')
            if ele is not None:
                root.remove(ele)
        except:
            pass
        ele = etree.Element('columnGroup')
        columnGroups = self.qualifiers('columnGroup')
        if columnGroups is not NotImplemented:
            for item in columnGroups:
                e = item.getEtree()
                ele.append(e)
        root.append(ele)
        return root,error

    def columnGroupNames(self):
        columnGroups = self.qualifiers('columnGroup')
        if columnGroups is NotImplemented:
            return []
        columnNameList = []
        for column in columnGroups:
            columnNameList.append(column.guiLabel())
        return columnNameList


class CProgramColumnGroup(CCP4Data.CData):
    '''A group of MTZ columns required for program input'''
    ERROR_CODES = {101 : {'description' : 'Column not in MTZ file'},
                   102 : {'description' : 'Column wrong type'},
                   103 : {'description' : 'Error setting columnGroup qualifier'},
                   104 : {'description' : 'Missing column selection'},
                   105 : {'description' : 'Specified column not found in MTZ file'},
                   106 : {'description' : 'Specified column has wrong type in MTZ file'},
                   107 : {'description' : 'Error reading columnGroup qualifier from XML file'},
                   108 : {'description' : 'No columnGroup qualifier'}}
    CONTENTS = {}
    QUALIFIERS = {'mustExist' : False, 'mtzFileKey' : '', # Name of the data item in parent dataContainer
                  'toolTipList' : [], 'default' : []} # A list of CColumnGroupItem
    QUALIFIERS_ORDER = ['mtzFileKey', 'mustExist', 'toolTipList', 'default']
    QUALIFIERS_DEFINITION = {'mtzFileKey'  : {'type' :str,
                                              'description' : 'The key for a CMtxDataFile in the same CContainer' },
                             'mustExist' : {'type' : bool,
                                            'description' : 'Flag if the parameter must be set at run time' },
                             'toolTipList' : {'type' : list,
                                              'description' : 'Tooltips for columns in group'},
                             'default' : {'type' : list, 'listItemType': str,
                                          'description' : 'Preferred values for column names'}}

    def build(self,qualifiers={}):
        self.__dict__['_value'] = {}
        for columnName in self.columnGroupNames():
            self.__dict__['_value'][columnName] = None

    def __init__(self,value={},qualifiers={},parent=None,name=None,**kw):
        qualis = {}
        qualis.update(qualifiers)
        qualis.update(kw)
        CCP4Data.CData.__init__(self,qualifiers=qualis,parent=parent,name=name)

    def contents(self):
        c = {}
        for columnName in self.columnGroupNames():
            c[columnName] = {'class' : CCP4Data.CString}
        return c

    def getColumnGroupQualifier(self):
        if 'columnGroup' not in self.__dict__['_qualifiers']:
            return []
        colGroupList = []
        columnGroups = self.__dict__['_qualifiers']['columnGroup']
        for colGroupObj in columnGroups:
            colGroup = colGroupObj.get()
            colGroupList.append(colGroup)
        return colGroupList

    def nColumns(self):
        return len(self.__dict__['_value'])

    def findColumnsInMtz(self,data={}):
        ''' Input is set of string values for columnLabels - Test if these columnLabels exist in the MTZ
        and return the CMtzColumns '''
        mtz = self.getDataByKey('mtzFileKey')
        if mtz is None or mtz.fileContent is None:
            return {}
        mtzColumns = {}
        for key,value in list(data.items()):
            mtzColumns[key] = mtz.fileContent.getColumn(value)
        return mtzColumns

    def validity(self,data = {}):
        if self.qualifiers('mustExist'):
            mtzColumns = self.findColumnsInMtz(data=data)
        else:
            mtzColumns = {}
        v = CErrorReport()
        unsetItems = 0
        columnGroups = self.qualifiers('columnGroup')
        if columnGroups is NotImplemented: return v
        for column in self.qualifiers('columnGroup'):
            columnName = str(column.columnName)
            if columnName not in data:
                if not self.qualifiers('allowUndefined'):
                    v.append(self.__class__,104,columnName,name=self.objectPath(),label=self.qualifiers('guiLabel'),stack=False)
            elif data[columnName] is None or  data[columnName] is NotImplemented or data[columnName] == '':
                if not self.qualifiers('allowUndefined'):
                    v.append(self.__class__,104,columnName,name=self.objectPath(),label=self.qualifiers('guiLabel'),stack=False)
            else:
                if self.qualifiers('mustExist'):
                    if columnName not in mtzColumns or mtzColumns[columnName] is None:
                        v.append( self.__class__,105,columnName,name=self.objectPath(),label=self.qualifiers('guiLabel'),stack=False)
                    elif column.columnType.count(str(mtzColumns[columnName].columnType))==0:
                        v.append( self.__class__,106,columnName,name=self.objectPath(),label=self.qualifiers('guiLabel'),stack=False)
        return v

    def fix(self,data={}):
        mtz = self.getDataByKey('mtzFileKey')
        if mtz is None:
            return -1
        nUnset = 0
        columnGroups = self.qualifiers('columnGroup')
        if columnGroups is NotImplemented:
            return {}
        for column in columnGroups:
            columnName = str(column.columnName)
            if data.get(columnName,None) is None:
                nUnset = nUnset + 1
                if column.partnerTo is not None and data.get(str(column.partnerTo),None) is not None:
                    partnerColumn = mtz.fileContent.getPartnerColumn(firstColumn=data.get(str(column.partnerTo)),columnGroupItem=column)
                    if  partnerColumn is not None:
                        data[columnName] = str(partnerColumn)
                        nUnset = nUnset - 1
        return data

    def __setattr__(self,name,value):
        if self.columnGroupNames().count(name):
            self.__dict__['_value'][name] = value
        else:
            CCP4Data.CData.__setattr__(self,name,value)

    def set(self,data={},fix=False,**kw):
        if isinstance(data,CProgramColumnGroup):
            data = data.get()
        data.update(kw)
        validity = self.validity(data)
        # If insufficient data try filling out from partner column information
        if validity.count(code=104) > 0 or fix:
            data = self.fix(data)
            validity = self.validity(data)
        if validity.maxSeverity()<SEVERITY_ERROR:
            for column in self.qualifiers('columnGroup'):
                columnName = str(column.columnName)
                if columnName in data:
                    if data[columnName] is None:
                        self.__dict__['_value'][columnName] = None
                    else:
                        self.__dict__['_value'][columnName] = data[columnName].split('/')[-1]
            self.dataChanged.emit()
        else:
            e = CException()
            e.extend(validity)
            raise e

    def unSet(self):
        for column in self.qualifiers('columnGroup'):
            self.__dict__['_value'][str(column.columnName)] = None

    def get(self,name='columnGroup'):
        if name == 'guiLabels':
            data = {}
            mtzColumns = self.findColumnsInMtz(data=self.__dict__['_value'])
            for key,columnObj in list(mtzColumns.items()):
                if columnObj is None:
                    data[key] = None
                else:
                    data[key] = columnObj.guiLabel()
            return data
        elif name == 'columnGroup':
            return self.__dict__['_value']
        else:
            return CCP4Data.CData.get(self,name)

    def setEtree(self,element=None,checkValidity=True):
        rv = CErrorReport()
        if element is not None:
            data = {}
            for columnName in self.columnGroupNames():
                ele = element.find(columnName)
                if ele is not None:
                    if ele.text is None:
                        data[columnName] = None
                    else:
                        data[columnName] = str(ele.text)
            rv.extend(self.set(data))
            return rv

    def getEtree(self):
        name = self.objectName()
        if name is None or len(name) == 0:
            name = self.className()
        element = etree.Element(name)
        for columnName in self.columnGroupNames():
            ele = etree.Element(columnName)
            if self.__dict__['_value'][columnName] is not None:
                ele.text = self.__dict__['_value'][columnName]
            else:
                ele.text = ''
            element.append(ele)
        return element

    def setColumnGroupQualifier(self, columnGroup=[], parent=None, initialise=False):
        # Expecting columnGroup to be a list of dicts or CColumnGroupItem
        # Each dict should have attributes for CColumnGroupItem object
        if '_qualifiers' not in self.__dict__:
            self.__dict__['_qualifiers'] = {}
        if initialise or ('columnGroup' not in self.__dict__['_qualifiers']):
            self.__dict__['_qualifiers']['columnGroup'] = []
        indx = -1
        for item in columnGroup:
            indx = indx + 1
            if isinstance(item,CColumnGroupItem):
                itemObj = item
            else:
                itemObj = CColumnGroupItem(item,parent=parent)
            if indx > 0:
                if itemObj.partnerTo is None:
                    itemObj.partnerTo = str(self.__dict__['_qualifiers']['columnGroup'][0].columnName)
                if itemObj.partnerOffset is None:
                    itemObj.partnerOffset = indx
            self.__dict__['_qualifiers']['columnGroup'].append(itemObj)
        self.build()

    def qualifiers(self, name=None, default=True, custom=True, contentQualifiers=False):
        return CCP4Data.CData.qualifiers(self, name=name, default=default, custom=custom, contentQualifiers=False)

    def setQualifiers(self,qualifiers={},**kw):
        qualis = {}
        qualis.update(qualifiers)
        if 'columnGroup' in qualis:
            self.setColumnGroupQualifier(qualis['columnGroup'])
            del qualis['columnGroup']
        CCP4Data.CData.setQualifiers(self,qualis)

    def setQualifiersEtree(self,element=None):
        # Custom parsing etree to extract the columnGroup info
        rv = CErrorReport()
        self.__dict__['_qualifiers']['columnGroup'] = []
        if element is not None:
            cGele = element.find('columnGroup')
            if cGele is not None:
                for ele in cGele.iterchildren():
                    try:
                        if str(ele.tag) == 'columnGroupItem':
                            item = CColumnGroupItem(parent=self)
                            item.setEtree(ele)
                            self.__dict__['_qualifiers']['columnGroup'].append(item)
                    except CException as e:
                        rv.extend(e)
                    except:
                        rv.append(self.__class__, 107, name=self.objectName())
                # Remove the columnGroup element from the etree and parse the rest with the base class method
                element.remove(cGele)
            rv.extend(CCP4Data.CData.setQualifiersEtree(self,element))
            self.build()
            return rv

    def qualifiersEtree(self,customOnly=True,tag=None):
        error = CErrorReport()
        if tag is None:
            tag = self.objectName()
        if tag is None or len(tag) == 0:
            tag = self.__class__.__name__
            error.append(CCP4Data.CData, 14)
        root,err = CCP4Data.CData.qualifiersEtree(self, customOnly=customOnly, tag=tag)
        error.extend(err)
        # Remove the columnGroup element which needs doing differently
        try:
            ele = root.find('columnGroup')
            if ele is not None:
                root.remove(ele)
        except:
            pass
        ele = etree.Element('columnGroup')
        columnGroups = self.qualifiers('columnGroup')
        if columnGroups is not NotImplemented:
            for item in columnGroups:
                e = item.getEtree()
                ele.append(e)
        root.append(ele)
        return root, error

    def columnGroupNames(self):
        columnGroups = self.qualifiers('columnGroup')
        if columnGroups is NotImplemented:
            return []
        columnNameList = []
        for column in columnGroups:
            columnNameList.append(str(column.columnName))
        return columnNameList

    def isSet(self, allowUndefined=False, allowDefault=True, allSet=True):
        for key in self.columnGroupNames():
            if self._value[key] is None:
                return False
        return True

class CFSigFColumnGroup(CProgramColumnGroup):
    QUALIFIERS = {'guiLabel' : 'Structure factor and sigma'}

    def build(self,qualifiers={}):
        # print("BUILD F-SIGF")
        self.__dict__['_qualifiers']['columnGroup'] = []
        item = CColumnGroupItem(parent = self, value = {'columnName' : 'F', 'columnType' : ['F']})
        self.__dict__['_qualifiers']['columnGroup'].append(item)
        item = CColumnGroupItem(parent = self, value = {'columnName' : 'SIGF', 'columnType' : ['Q'],
                                                        'partnerTo': 'F', 'partnerOffset' : 1})
        self.__dict__['_qualifiers']['columnGroup'].append(item)
        CProgramColumnGroup.build(self)

class CISigIColumnGroup(CProgramColumnGroup):
    QUALIFIERS = {'guiLabel' : 'Intensity and sigma'}

    def build(self,qualifiers={}):
        self.__dict__['_qualifiers']['columnGroup'] = []
        item = CColumnGroupItem(parent=self, value = {'columnName' : 'I', 'columnType' : ['J']})
        self.__dict__['_qualifiers']['columnGroup'].append(item)
        item = CColumnGroupItem(parent=self, value = {'columnName' : 'SIGI', 'columnType' : ['Q'],
                                                      'partnerTo': 'I', 'partnerOffset' : 1})
        self.__dict__['_qualifiers']['columnGroup'].append(item)
        CProgramColumnGroup.build(self)

class CFPairColumnGroup(CProgramColumnGroup):
    QUALIFIERS = { 'guiLabel' : 'Anomalous structure factors and sigma' }

    def build(self, qualifiers={}):
        self.__dict__['_qualifiers']['columnGroup'] = []
        item = CColumnGroupItem( parent=self, value ={ 'columnName' : 'Fplus',
                                                        'columnType' : [ 'G' ] })
        self.__dict__['_qualifiers']['columnGroup'].append(item)
        offset = 0
        for name,coltype in [['SIGFplus','L'],['Fminus','G'],['SIGFminus','L']]:
            offset += 1
            item =  CColumnGroupItem( parent = self, value = {'columnName' : name, 'columnType' : [coltype],
                                                              'partnerTo': 'Fplus', 'partnerOffset' : offset})
            self.__dict__['_qualifiers']['columnGroup'].append(item)
        CProgramColumnGroup.build(self)

class CIPairColumnGroup(CProgramColumnGroup):
    QUALIFIERS = {'guiLabel' : 'Anomalous intensities and sigma' }

    def build(self,qualifiers={}):
        self.__dict__['_qualifiers']['columnGroup'] = []
        item = CColumnGroupItem( parent=self, value ={ 'columnName' : 'Iplus',
                                                        'columnType' : [ 'K' ] })
        self.__dict__['_qualifiers']['columnGroup'].append(item)
        offset = 0
        for name,coltype in [['SIGIplus','M'],['Iminus','K'],['SIGIminus','M']]:
            offset += 1
            item = CColumnGroupItem(parent = self, value = {'columnName' : name, 'columnType' : [coltype],
                                                            'partnerTo': 'Iplus', 'partnerOffset' : offset })
            self.__dict__['_qualifiers']['columnGroup'].append(item)
        CProgramColumnGroup.build(self)


class CHLColumnGroup(CProgramColumnGroup):
    QUALIFIERS = { 'guiLabel' : 'Hendrickson-Lattmann coefficients' }

    def build(self, qualifiers = {}):
        self.__dict__['_qualifiers']['columnGroup'] = []
        item = CColumnGroupItem( parent=self, value ={'columnName' : 'HLA', 'columnType' : ['A']})
        self.__dict__['_qualifiers']['columnGroup'].append(item)
        offset = 0
        for name in ['HLB','HLC','HLD']:
            offset += 1
            item =  CColumnGroupItem(parent = self, value = {'columnName' : name, 'columnType' : ['A'],
                                                             'partnerTo': 'HLA', 'partnerOffset' : offset})
            self.__dict__['_qualifiers']['columnGroup'].append(item)
        CProgramColumnGroup.build(self)


class CPhiFomColumnGroup(CProgramColumnGroup):
    QUALIFIERS = { 'guiLabel' : 'Phase and figure of merit' }

    def build(self,qualifiers={}):
        self.__dict__['_qualifiers']['columnGroup'] = []
        item = CColumnGroupItem( parent = self, value = {'columnName' : 'PHI', 'columnType' : ['P']})
        self.__dict__['_qualifiers']['columnGroup'].append(item)
        item = CColumnGroupItem( parent = self, value = {'columnName' : 'FOM', 'columnType' : ['W'],
                                                         'partnerTo': 'PHI', 'partnerOffset' : 1 })
        self.__dict__['_qualifiers']['columnGroup'].append(item)
        CProgramColumnGroup.build(self)


class CMapColumnGroup(CProgramColumnGroup):
    QUALIFIERS = {'guiLabel' : 'Structure factor and phase to define a map'}

    def build(self, qualifiers={}):
        # print("BUILD Map", )
        self.__dict__['_qualifiers']['columnGroup'] = []
        item = CColumnGroupItem(parent=self, value={'columnName' : 'F', 'columnType' : ['F']})
        self.__dict__['_qualifiers']['columnGroup'].append(item)
        item = CColumnGroupItem(parent=self, value={'columnName' : 'PHI', 'columnType' : ['P'],
                                                    'partnerTo': 'F', 'partnerOffset' : 1 })
        self.__dict__['_qualifiers']['columnGroup'].append(item)
        CProgramColumnGroup.build(self)


class CFreeRColumnGroup(CProgramColumnGroup):
    QUALIFIERS = {'guiLabel' : 'Set of FreeR flags'}

    def build(self,qualifiers={}):
        self.__dict__['_qualifiers']['columnGroup'] = []
        item = CColumnGroupItem( parent = self, value = {'columnName' : 'FREER', 'columnType' : ['I']})
        self.__dict__['_qualifiers']['columnGroup'].append(item)
        CProgramColumnGroup.build(self)


class CAnomalousIntensityColumnGroup(CProgramColumnGroup):
    '''
    Selection of I and AnomI columns from MTZ.
    Expected to be part of ab initio phasing dataset ( CDataset)
    '''
    QUALIFIERS = {'toolTipList' : ['The real part of the experimental intensity',
                                    'The anomalous part of the experimental intensity'],
                                        'guiLabel' : 'Intensity and anomalous intensity'}
    def build(self,qualifiers={}):
        self.__dict__['_qualifiers']['columnGroup'] = []
        item = CColumnGroupItem( parent=self, value ={'columnName' : 'I', 'columnType' : ['J']})
        self.__dict__['_qualifiers']['columnGroup'].append(item)
        item =  CColumnGroupItem( parent=self, value ={'columnName' : 'AnomI', 'columnType' : ['J']})
        self.__dict__['_qualifiers']['columnGroup'].append(item)
        CProgramColumnGroup.build(self)

    def columnGroupNames(self):
        return ['I','AnomI']

    # The CProgramColumnGroup base class expects a qualifier mtzFileKey which is the
    # name of an CMtzDatFile object - how does it get this if it is part of a CDataset?


class CAnomalousColumnGroup(CProgramColumnGroup):
    '''
    Selection of F/I and AnomF/I columns from MTZ.
    Expected to be part of ab initio phasing dataset ( CDataset)
    '''
    QUALIFIERS = {'toolTipList' : ['The real part of the experimental structure factors',
                                   'The anomalous part of the experimental structure factors']}

    def build(self,qualifiers={}):
        self.__dict__['_qualifiers']['columnGroup'] = []
        item = CColumnGroupItem(parent = self, value = {'columnName' : 'F/I', 'columnType' : ['F', 'J']})
        self.__dict__['_qualifiers']['columnGroup'].append(item)
        item = CColumnGroupItem( parent=self, value ={'columnName' : 'AnomF/I', 'columnType' : ['F', 'J']})
        self.__dict__['_qualifiers']['columnGroup'].append(item)
        CProgramColumnGroup.build(self)

    def columnGroupNames(self):
        return ['F/I', 'AnomF/I']

    # Use the generic getEtree/setEtree
    def getEtree(self):
        CCP4Data.CData.getEtree(self)

    def setEtree(self, element, checkValidity = True):
        CCP4Data.CData.setEtree(self, element)

class CExperimentalDataType(CCP4Data.CString):
    '''Experimental data type e.g. native or peak'''
    QUALIFIERS = {'onlyEnumerators' : True, 'enumerators' : ['native', 'derivative', 'SAD', 'peak', 'inflection',
                                                             'high_remote', 'low_remote', ''], 'default' : 'SAD' }

class CFormFactor(CCP4Data.CData):
    '''
    The for factor (Fp and Fpp) for a giving element and wavelength
    '''
    CONTENTS = {'Fp' : {'class' : CCP4Data.CFloat, 'qualifiers' : {'toolTip' : "Form factor F' for element at given wavelength"}},
                'Fpp' : {'class' : CCP4Data.CFloat, 'qualifiers' : {'toolTip' : "Form factor F'' for element at given wavelength"}}}
    CONTENTS_ORDER = ['Fp', 'Fpp']

    def validity(self,arg):
        err = CErrorReport()
        return err

    def guiValue(self):
        return self.__dict__['_value']['Fp'].__str__() + ','+self.__dict__['_value']['Fpp'].__str__()

    # Do we need lookup/calculate values if not provided?


class CAnomalousScatteringElement(CCP4ModelData.CElement):
    '''Definition of a anomalous scattering element'''
    QUALIFIERS = {'onlyEnumerators' : False, 'enumerators' : ['Br','Fe','Pt','Se'],
                  'charWidth' : 4, 'default' : 'Se' }


class CShelxLabel(CCP4Data.CString):
    QUALIFIERS = {'onlyEnumerators' : True, 'default' : 'UNDEFINED',
                  'enumerators' : [ 'UNDEFINED','HREM','LREM','PEAK','INFL','NAT','DERI'],
                  'menuText' : ['undefined','high remote','low remote','peak','inflection',
                                'native','derivative' ],
                  'toolTip' : 'Hint to Shelx for the use of the dataset' }


class CAsuComponent(CCP4Data.CData):
    '''A component of the asymmetric unit. This is for use in MR, defining
         what we are searching for. '''
    CONTENTS = {'moleculeType' : {'class' : CCP4Data.CString,
                'qualifiers' : {'onlyEnumerators' : True, 'enumerators' : ['PROTEIN', 'NUCLEIC'],
                                'menuText' : ['protein', 'nucleic acid'],
                                'default' : 'PROTEIN',
                                'toolTip' : 'Molecule type'}},
                   'seqFile' : {'class' : CCP4ModelData.CSeqDataFile,
                                'qualifiers' : {'jobCombo' : False,
                                                'mustExist' : True,
                                                'allowUndefined' : False } },
                   'numberOfCopies' : {'class' : CCP4Data.CInt,
                                       'qualifiers' : {'allowUndefined' : False,
                                                       'toolTip' : 'Number of copies of sequence',
                                                       'min' : 0, 'max' : 999, 'default': 1,
                                                       'enumerators' : [1,2,3,4,5,6,7,8,9,10,11,12]}}}

    def getTextItem(self):
        if self.__dict__['_value']['seqFile'].isSet():
            if self.__dict__['_value']['seqFile'].annotation.isSet():
                return str(self.__dict__['_value']['numberOfCopies']) + ' x '+str(self.__dict__['_value']['seqFile'].annotation)
            else:
                return str(self.__dict__['_value']['numberOfCopies']) + ' x '+str(self.__dict__['_value']['seqFile'])
        else:
            return '--'

    def getTableTextItem(self):
        return self.getTextItem()


class CAsuComponentList(CCP4Data.CList):
    SUBITEM = {'class' : CAsuComponent}
    QUALIFIERS = {'listMinLength' : 1, 'guiLabel' : 'Contents of asymmetric unit'}

    def validity(self,arg):
        try:
            mode = self.parent().parent().find('COMP_BY').__str__()
        except:
            mode = 'ASU'
        if mode != 'ASU':
            return CErrorReport()
        else:
            err = CCP4Data.CList.validity(self,arg)
            if len(arg)<1:
                err.append(self.__class__, 101, name = self.objectPath(), label = self.qualifiers('guiLabel'), stack = False)
        return err

    def saveToDb(self):
        saveList = []
        for obj in self.__dict__['_value']:
            if obj.seqFile.isSet():
                saveList.append(obj.seqFile)
        return saveList,None,{}

    def molecularWeight(self):
        wt = 0.0
        for ii in range(self.__len__()):
            seqAnalysis = self.__dict__['_value'][ii].seqFile.fileContent.getAnalysis()
            wt = wt + (seqAnalysis * float(self.__dict__['_value'][ii].numberOfCopies))
        return wt


class CMiniMtzDataFile(CMtzDataFile):

    SUBTYPE_MENU = None
    QUALIFIERS = {'mimeTypeName' : 'application/CCP4-mtz-mini', 'fileExtensions' : ['mtz','cif','ent'],
                  'fileContentClassName' : 'CMtzData', 'saveToDb' : True,
                  'correctColumns' : ['FQ','JQ','GLGL','KMKM','AAAA','PW','FP','I'],
                  'columnGroupClassList' : NotImplemented,
                  'toolTip' : 'Mini-MTZ file containing reflection,phases,FreeR set or map coefficients',
                  'helpFile' : 'data_files#MTZ'}
    QUALIFIERS_ORDER = ['fileExtensions', 'mimeTypeName', 'mimeTypeDescription', 'allowUndefined', 'mustExist',
                        'fromPreviousJob', 'jobCombo', 'fileContentClassName', 'isDirectory', 'saveToDb', 'requiredSubType',
                        'requiredContentFlag', 'correctColumns', 'columnGroupClassList', 'sameCrystalAs']
    QUALIFIERS_DEFINITION = {'correctColumns' : {'type' : list, 'listItemType' : str,
                                                 'description' : 'A list of coloumn data types expected in the file'}}

    ERROR_CODES = {201 : {'description' : 'Wrong number of columns'},
                   202 : {'description' : 'Wrong column types'},
                   203 : {'description' : 'No correct column types found in file'},
                   204 : {'description' : 'Duplicate or additional column types found in file'},
                   205 : {'description' : 'Columns in file have non-standard labels'},
                   206 : {'description' : 'File contains unmerged data'},
                   210 : {'description' : 'Failed creating mini-MTZ'},
                   211 : {'description' : 'Insufficient columns selected from imported MTZ'},
                   212 : {'description' : 'Data already imported as', 'severity' : SEVERITY_WARNING},
                   220 : {'description' : 'Can not convert file content, file does not exist'},
                   221 : {'description' : 'Can not convert file content, existing content insufficiently rich'},
                   222 : {'description' : 'Can not convert file content, bad input for target content'},
                   223 : {'description' : 'Can not recognise file content'},
                   224 : {'description' : 'Not possible to convert to required content - no mechanism implemented'},
                   225 : {'description' : 'Failed importing from an mmcif file - failed running cif2mtz'},
                   226 : {'description' : 'Failed importing from an mmcif file - no output from cif2mtz'}}

    def __init__(self,value={},qualifiers={},parent=None,name=None,fullPath=None,**kw):
        qualis = {}
        qualis.update(qualifiers)
        qualis.update(kw)
        CCP4File.CDataFile.__init__(self, value = value, qualifiers = qualis, parent = parent,
                                    name = name, fullPath = fullPath, keywords = kw)
        self.__dict__['sourceFileName'] = None

    def getSourceFileName(self):
        return self.__dict__['sourceFileName']

    def conversion(self,targetContent):
        return 'ok',targetContent

    def miniMtzType(self):
        '''
        Test for mini-MTZ type and contentFlag
        '''
        self.loadFile()
        signature = self.fileContent.columnSignature()
        for cls in [CObsDataFile, CPhsDataFile, CMapCoeffsDataFile, CFreeRDataFile]:
            contentFlag = 0
            for correctSig in cls.QUALIFIERS['correctColumns']:
                contentFlag += 1
                if signature == correctSig:
                    return cls,contentFlag
        return None, None

    def validColumns(self,correctColumns=None):
        rv = CErrorReport()
        self.loadFile()
        if not self.fileContent.merged:
            rv.append(self.__class__, 206, name=self.objectPath())
        if correctColumns is None:
            correctColumns = self.qualifiers('correctColumns')
        signature = self.getFileContent().columnSignature()
        nMap = 0
        nDup = 0
        badColumnNames = False
        sigListIndx = -1
        for items in correctColumns:
            sigListIndx += 1
            n = signature.count(items)
            if n > 1:
                nDup += 1
            if n >= 1 :
                # Do the columns have the correct column names?
                columnNames = self.fileContent.columnNames(signature.index(items), len(items))
                if columnNames != self.CONTENT_SIGNATURE_LIST[sigListIndx]:
                    badColumnNames = True
            nMap += n
            signature = re.sub(items, '', signature)
        # None of required columns found
        if nMap == 0:
            rv.append(self.__class__, 203, name=self.objectPath())
        # Ambiguity of required column or additional columns
        # Column selection will be displayed
        elif nDup > 1 or len(signature) > 0:
            rv.append(self.__class__, 204, name=self.objectPath())
        elif badColumnNames:
            """
            # Set the columnGroup data that woud normally be set by the column selection
            # dialog so we can just use the same code to call splitMtz
            # NO-- just handle this case the usual way for 'monster' mtz so it will go through the dialog
            columnGroup = self.columnGroup()
            for name in columnNames:
                columnGroup[sigListIndx].columnName.set(name)
            print 'CMiniMtzDataFile.validColumns setting columnGroup',columnGroup[sigListIndx]
            """
            rv.append(self.__class__, 205, name=self.objectPath())
        return rv

    def columnGroup(self):
        if 'columnGroupList' not in self.__dict__:
            self.__dict__['columnGroupList'] = []
            clsList = self.qualifiers('columnGroupClassList')
            if not isinstance(clsList,list): clsList = [clsList]
            for cls in clsList:
                self.__dict__['columnGroupList'].append(cls(parent=self.parent(),name=self.objectName()+'_COLUMNS',qualifiers={'mtzFileKey': self.objectName()}))
        return self.__dict__['columnGroupList']

    def defaultName(self,jobId=None):
        jobDirectory = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId = jobId)
        return  os.path.normpath(os.path.join(jobDirectory, self.objectName() + CCP4File.CDataFile.SEPARATOR + self.qualifiers('fileLabel') + '.mtz'))

    def set(self,value={},**kw):
        CCP4File.CDataFile.set(*[self,value], **kw)
        self.__dict__['sourceFileName'] = None

    def splitMtz(self,jobId=None,projectId=None,contentFlag=None,i2Labels=[],columnLabels=[]):
        errorReport = CErrorReport()
        # Set name for new split file and if it already exists remove previous refernce from db
        jobDirectory = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=jobId)
        filename = self.importFileName(jobDirectory=jobDirectory)
        #Set up the colin/colout
        colin = ','.join(columnLabels)
        if len(i2Labels) == 0:
            if contentFlag is not None:
                i2Labels = self.CONTENT_SIGNATURE_LIST[contentFlag-1]
            else:
                i2Labels = self.CONTENT_SIGNATURE_LIST[0]
        colout = ','.join(i2Labels)
        if isinstance(self,CFreeRDataFile):
            fileType = 10
        elif isinstance(self,CObsDataFile):
            fileType = 11
        elif isinstance(self,CPhsDataFile):
            fileType = 12
        elif isinstance(self,CMapCoeffsDataFile):
            fileType = 13
        else:
            fileType = None
        #Have we done the same import before?
        dbFileId, importId, checksum, dbAnnotation = CCP4Modules.PROJECTSMANAGER().alreadyImportedId(sourceFileName=self.__str__(),
                                                                                                     projectId=projectId,
                                                                                                     contentFlag=contentFlag,
                                                                                                     sourceFileReference=colin,
                                                                                                     fileType=fileType)
        print('CMiniMtzDataFile.splitMtz alreadyImportedId in', self.__str__(), projectId, fileType, contentFlag, colin)
        print('Testing if file previously imported previous db id', dbFileId)
        if dbFileId is not None:
            self.setDbFileId(dbFileId)
            errorReport.append(self.__class__, 212, name=self.objectPath(), details=dbAnnotation)
            return errorReport
        # Run cmtzsplit
        logFile =  os.path.normpath(os.path.join(jobDirectory, self.objectName() + '_mtzsplit.log'))
        cbin = os.path.normpath(os.path.join(CCP4Utils.getOSDir(), 'bin', 'cmtzsplit'))
        if not os.path.exists(cbin):
            cbin = os.path.normpath(os.path.join(CCP4Utils.getCCP4Dir(), 'bin', 'cmtzsplit'))
        arglist = ['-mtzin', self.__str__()]
        arglist.extend(['-mtzout', filename])
        arglist.extend(['-colin', colin,'-colout', colout])
        pid = CCP4Modules.PROCESSMANAGER().startProcess(cbin, arglist, logFile=logFile)
        status = CCP4Modules.PROCESSMANAGER().getJobData(pid)
        exitCode = CCP4Modules.PROCESSMANAGER().getJobData(pid, 'exitCode')
        if status != 0:
            errorReport.append(self.__class__, 210, 'Exit status:' + str(status), name=self.objectPath())
            return errorReport
        #Replace 'dataset' in the annotation with the correct dataset name
        annotation = self.annotation.__str__()
        if annotation.count('/dataset'):
            mtzColumn = self.fileContent.getColumn(columnLabels) # KJS - Put an s on the end here.
            dataset = ''
            if mtzColumn is not None and mtzColumn.dataset is not None and mtzColumn.dataset != 'HKL_base':
                dataset = '/' + mtzColumn.dataset.__str__()
            else:
                dataset = ''
            annotation = re.sub('/dataset', dataset, annotation)
        # Set sourceFileName and new file name for handling by PROJECTSMANAGER.importFiles()
        sourceFileName = self.__str__()
        self.setFullPath(filename)
        self.unsetFileContent()
        self.contentFlag.set(contentFlag)
        self.annotation.set(annotation)
        self.__dict__['sourceFileName'] = sourceFileName
        self.__dict__['sourceFileReference'] = colin[0:-1]
        return errorReport

    def contentSignature(self):
        return []

    def setContentFlag(self,reset=False):
        # Test file content to determine the contentFlag
        # is None => don't know or is not valid column label signature
        #print 'CMiniMtzDataFile.setContentFlag',self.objectName(),
        #print str(self.getFullPath()),self.isSet(),self.exists()
        #print 'CMiniMtzDataFile.setContentFlag',self.objectName(),self,fileInfo
        if (not reset) and self.__dict__['_value']['contentFlag'].isSet():
            return int(self.__dict__['_value']['contentFlag'])
        self.__dict__['_value']['contentFlag'].unSet()
        if (not self.isSet()) or (not self.exists()):
            return None
        # Try is info in Db
        if self.dbFileId.isSet() and not reset:
            flag = CCP4Modules.PROJECTSMANAGER().db().getFileInfo(fileId=str(self.dbFileId), mode='filecontent')
            if flag is not None:
                self.__dict__['_value']['contentFlag'].set(flag)
                return flag
        sigList = self.contentSignature()
        columnList = self.getFileContent().getListOfColumns()
        labelList=[]
        for item in columnList: labelList.append(item.columnLabel.__str__())
        flag = 1
        while flag <= len(sigList):
            if labelList == sigList[flag-1]:
                self.__dict__['_value']['contentFlag'].set(flag)
                return flag
            else:
                flag += 1
        return None

    def updateData(self):
        CMtzDataFile.updateData(self)

    def importFromCif(self, jobId=None):
        errorReport = CErrorReport()
        # Set name for new split file and if it already exists remove previous reference from db
        if jobId is not None:
            jobDirectory = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=jobId)
        else:
            jobDirectory = CCP4Utils.makeTmpFile(cdir=True)
        # Run cif2mtz
        logFile = os.path.normpath(os.path.join(jobDirectory,self.objectName()+'_cif2mtz.log'))
        cbin = os.path.normpath(os.path.join( CCP4Utils.getCCP4Dir(), 'bin', 'cif2mtz' ))
        hklout = os.path.normpath(os.path.join(jobDirectory,self.stripedName()+'.mtz'))
        arglist = ['hklin', self.__str__(), 'hklout', hklout]
        inputText = '''END\n'''
        pid = CCP4Modules.PROCESSMANAGER().startProcess(cbin, arglist, logFile=logFile, inputText=inputText)
        status = CCP4Modules.PROCESSMANAGER().getJobData(pid)
        exitCode = CCP4Modules.PROCESSMANAGER().getJobData(pid,'exitCode')
        if status != 0:
            errorReport.append(self.__class__, 225, 'Exit status:'+str(status), name=self.objectPath())
            return errorReport
        if os.path.exists(hklout):
            self.setFullPath(hklout)
        else:
            errorReport.append(self.__class__, 226, name=self.objectPath())
        return errorReport

    def saveAs(self, cformat='SHELX', fileName=None, mergeWith=[]):
        # Save as non-MTZ format probably by running CMtzDataFile.runMtz2various
        # mergeWith - list of other CMiniMtzDataFiles to be merged into the output file
        #Keep a list of the file labels (assumes we want mtz2various to ouput all)
        labin = self.columnNames()
        # Merge in any other required data objects using cmtzjoin
        if len(mergeWith) > 0:
            tmpFile = CCP4Utils.makeTmpFile(name=self.baseName.__str__(), extension='mtz')
            infiles = [[self.__str__(), self.columnNames(ifString=True)]]
            for mtzObj in mergeWith:
                labin.extend(mtzObj.columnNames())
                infiles.append ( [ mtzObj.__str__(), mtzObj.columnNames(ifString=True)] )
            hklin,err = self.runMtzjoin(tmpFile, infiles)
            if hklin is None or err.maxSeverity() > SEVERITY_WARNING:
                return None, err
        else:
            hklin = self.__str__()
        # Run mtz2various and return hklout,errorReport  (hklout is None if job failed)
        if cformat in ['SHELX']:
            return self.runMtz2various(hklin=hklin, labin=labin, output=cformat, hklout=fileName)
        # Simple test from pyi2.

    def datasetName(self):
        if not self.isSet():
            return ''
        for d in self.fileContent.datasets:
            if d != 'HKL_base':
                return str(d)
        return ''


class CMiniMtzDataFileList(CCP4Data.CList):
    SUBITEM={'class' : CMiniMtzDataFile}


class CObsDataFile(CMiniMtzDataFile):
    SUBTYPE_OBSERVED = 1
    SUBTYPE_DERIVED = 2
    SUBTYPE_REFERENCE = 3
    CONTENT_FLAG_IPAIR = 1
    CONTENT_FLAG_FPAIR = 2
    CONTENT_FLAG_IMEAN = 3
    CONTENT_FLAG_FMEAN = 4
    CONTENT_ANNOTATION = ['Anomalous Is', 'Anomalous SFs', 'Mean Is' ,'Mean SFs']
    # 4=not possible, 0 = same, 1=using mtzjoin, 2=using ctruncate, 3=ctruncate & mtzsplit
    #                               TO
    CONTENT_CONVERSION = [[0, 2, 2, 2],
                          [4, 0, 4, 2],
                          [4, 4, 0, 2],
                          [4, 4, 4, 0]]
    CONTENT_SIGNATURE_LIST = [['Iplus','SIGIplus','Iminus','SIGIminus'],
                              ['Fplus','SIGFplus','Fminus','SIGFminus'],
                              ['I','SIGI'], ['F','SIGF']]
    CONTENTS = {}
    CONTENTS.update(CMiniMtzDataFile.CONTENTS)
    CONTENTS['subType'] = {'class' : CCP4Data.CInt,
                           'qualifiers' : {'default' :SUBTYPE_OBSERVED,
                                           'enumerators' : [SUBTYPE_OBSERVED, SUBTYPE_DERIVED, SUBTYPE_REFERENCE],
                                           'onlyEnumerators' : True,
                                           'menuText' : ['observed data','derived data','reference data']}}
    QUALIFIERS = {'mimeTypeName' : 'application/CCP4-mtz-observed', 'mimeTypeDescription' : 'MTZ observed',
                  'fileExtensions' : ['mtz','cif','ent'],
                  'fileContentClassName' : 'CMtzData',
                  'fileLabel' : 'observed_data', 'guiLabel' : 'Reflections',
                  'toolTip' : "Observed structure factors or intensities",
                  'correctColumns' : ['KMKM','GLGL','JQ','FQ'],
                  'columnGroupClassList' : [CIPairColumnGroup, CFPairColumnGroup, CISigIColumnGroup, CFSigFColumnGroup],
                  'downloadModes' : ['ebiSFs'], 'helpFile' : 'data_files#Obs'}
    ERROR_CODES = {301 : { 'description' : 'Running ctruncate failed'},
                   302 : { 'description' : 'Running cmtzsplit to convert observed data type failed'},
                   303 : { 'description' : 'Running sftools failed'}}

    def conversion(self, targetContent):
        if not self.contentFlag.isSet():
            return 'ok',targetContent
        if isinstance(targetContent, list):
            import copy
            targetContentList = copy.deepcopy(targetContent)
        else:
            targetContentList = [targetContent]
        for targetContent in targetContentList:
            c = ['ok','mtzjoin','convert','','no'][self.CONTENT_CONVERSION[self.contentFlag-1][targetContent-1]]
            if c != 'no':
                return c,targetContent
        return c, targetContentList[0]

    def requiredContent(self):
        # Return the allowed contentFlag values for CDataFileView.getJobsWithOutputFiles to call to db getJobsWithOutputFiles
        contentList = self.qualifiers('requiredContentFlag')
        if contentList is None or contentList is NotImplemented or contentList == [self.CONTENT_FLAG_FMEAN]:
            return None
        retList = []
        for content in contentList:
            for n in range(1, 5):
                if self.CONTENT_CONVERSION[n-1][content-1] != 4:
                    if retList.count(n) == 0:
                        retList.append(n)
        if len(retList) == 1:
            return retList[0]
        else:
            return retList

    def contentSignature(self):
        return CObsDataFile.CONTENT_SIGNATURE_LIST

    def columnNames(self, ifString=False, content=None):
        if content is None:
            if not self.__dict__['_value']['contentFlag'].isSet():
                self.setContentFlag()
            if not self.__dict__['_value']['contentFlag'].isSet():
                return []
            else:
                content = int(self.contentFlag)
        if content in [CObsDataFile.CONTENT_FLAG_IPAIR,
                       CObsDataFile.CONTENT_FLAG_IMEAN,
                       CObsDataFile.CONTENT_FLAG_FPAIR,
                       CObsDataFile.CONTENT_FLAG_FMEAN]:
            if ifString:
                return ','.join(CObsDataFile.CONTENT_SIGNATURE_LIST[content-1])
            else:
                return CObsDataFile.CONTENT_SIGNATURE_LIST[content-1]
        else:
            if ifString:
                return ''
            else:
                return []

    def convert(self, targetContent=None, targetFile=None, parentPlugin=None):
        error = CErrorReport()
        if not self.isSet() or not self.exists():
            error.append(self.__class__, 220, name=self.objectPath())
            return None, error
        if targetContent is None or not isinstance(targetContent, int) or targetContent < 1 or targetContent > 4:
            error.append(self.__class__, 222, name=self.objectPath())
            return None, error
        if not self.__dict__['_value']['contentFlag'].isSet():
            flag = self.setContentFlag()
            if flag is None:
                error.append(self.__class__, 223, name=self.objectPath())
                return None, error
        if self.__dict__['_value']['contentFlag'] == targetContent:
            return self.fullPath.__str__(),error
        elif self.__dict__['_value']['contentFlag'] > targetContent:
            error.append(self.__class__, 221, name=self.objectPath())
            return None, error
        if targetFile is None:
            fsplit = os.path.splitext(self.fullPath.__str__())
            mode = ('asIPAIR', 'asFPAIR', 'asIMEAN', 'asFMEAN')[targetContent - 1]
            targetFile = fsplit[0] + '_'  + mode + fsplit[1]
            if os.path.exists(targetFile):
                return targetFile, error
        if self.__dict__['_value']['contentFlag'] == CObsDataFile.CONTENT_FLAG_IPAIR:
            if targetContent in (CObsDataFile.CONTENT_FLAG_FPAIR, CObsDataFile.CONTENT_FLAG_FMEAN, CObsDataFile.CONTENT_FLAG_IMEAN):
                return self.runTruncate(targetContent=targetContent, targetFile=targetFile, parentPlugin=parentPlugin)
        if self.__dict__['_value']['contentFlag'] == CObsDataFile.CONTENT_FLAG_IMEAN:
            if targetContent == CObsDataFile.CONTENT_FLAG_FMEAN:
                return self.runTruncate(targetContent=targetContent, targetFile=targetFile, parentPlugin=parentPlugin)
        if self.__dict__['_value']['contentFlag'] == CObsDataFile.CONTENT_FLAG_FPAIR:
            if targetContent == CObsDataFile.CONTENT_FLAG_FMEAN:
                return self.runSftools(targetContent=targetContent, targetFile=targetFile, parentPlugin=parentPlugin, mode='Amplitudes')
        # Its broke - no conversion possible
        error.append(self.__class__, 224, name=self.objectPath())
        return None, error

    def runSftools(self, targetContent=None, targetFile=None, parentPlugin=None, mode='Amplitudes'):
        error = CErrorReport()
        # Use SFTOOLS to calculate Fmean, SIGFMean
        cbin =  os.path.normpath(os.path.join( CCP4Utils.getCCP4Dir(), 'bin', 'sftools' ))
        if parentPlugin is not None:
            myDir =  os.path.normpath(os.path.join(parentPlugin.workDirectory, 'sftools'))
            os.mkdir(myDir)
        logFile =  os.path.normpath(os.path.join(myDir, 'sftools.log'))
        argList = []
        colPlus, colSigPlus, colMinus, colSigMinus = self.columnNames()
        if mode == 'Amplitudes':
            colOut, colTypeOut, colSigOut, colTypeSigOut = ('F', 'F', 'SIGF', 'Q')
        elif mode == 'Intensities':
            colTypeOut, colSigOut, colTypeSigOut = ('I', 'J', 'SIGI', 'Q')
        inputText = "READ " + self.fullPath.__str__() + "\n"
        #Evaluate weighted mean only if Fplus and Fminus both present
        inputText += "SELECT ONLY COL " + colPlus + " PRESENT\n"
        inputText += "SELECT MINUS COL " + colMinus + " ABSENT\n"
        inputText += "CALC COL SIGSQm = COL " + colSigMinus + " COL " + colSigMinus + " *\n" #Col 5
        inputText += "CALC COL FpXSIGSQm = COL " + colPlus + " COL SIGSQm *\n" #Col 6
        inputText += "CALC COL SIGSQp = COL " + colSigPlus + " COL "+colSigPlus+" *\n" #Col 7
        inputText += "CALC COL FmXSIGSQp = COL " + colMinus + " COL SIGSQp *\n" #Col 8
        inputText += "CALC COL SIGSQpPSIGSQm = COL SIGSQp COL SIGSQm +\n" #Col 9
        inputText += "CALC COL SIGSQpXSIGSQm = COL SIGSQp COL SIGSQm *\n" #Col 10
        inputText += "CALC COL FpXSIGSQmPFmXSIGSQp = COL FpXSIGSQm COL FmXSIGSQp +\n" #Col 11
        inputText += "CALC COL " + colOut + " = COL FpXSIGSQmPFmXSIGSQp COL SIGSQpPSIGSQm /\n"#Col 12
        inputText += "CALC COL SIGSQF = COL SIGSQpXSIGSQm COL SIGSQpPSIGSQm /\n"#Col 13
        inputText += "CALC COL " + colSigOut + " = COL SIGSQF LN 2 / EXP\n"#Col 14
        #Now deal with missing Fplus or Fminus
        inputText += "SELECT ONLY COL " + colMinus + " ABSENT\n"
        inputText += "CALC COL " + colOut + " = COL " + colPlus + " 0 +\n"
        inputText += "CALC COL " + colSigOut + " = COL " + colSigPlus + " 0 +\n"
        inputText += "SELECT ONLY COL " + colPlus + " ABSENT\n"
        inputText += "CALC COL " + colOut + " = COL " + colMinus + " 0 +\n"
        inputText += "CALC COL " + colSigOut + " = COL " + colSigMinus + " 0 +\n"
        inputText += "SELECT ALL\n"
        inputText += "SET TYPE COL " + colOut + "\n"
        inputText += colTypeOut + "\n"
        inputText += "SET TYPE COL " + colSigOut + "\n"
        inputText += colTypeSigOut + "\n"
        inputText += "WRITE " + targetFile + " COL " + colOut + " " + colSigOut + "\n"
        inputText += "STOP\n"
        pid = CCP4Modules.PROCESSMANAGER().startProcess(cbin, argList, inputText=inputText, logFile=logFile, cwd=myDir)
        status = CCP4Modules.PROCESSMANAGER().getJobData(pid)
        exitCode = CCP4Modules.PROCESSMANAGER().getJobData(pid,'exitCode')
        if status != 0 or not os.path.exists(targetFile):
            error.append(self.__class__, 303, str(targetFile), name=self.objectName())
            return None, error
        else:
            print("CObsDataFile.runSftools targetFile error", targetFile, error)
            return targetFile, error

    def runTruncate(self, targetContent=None, targetFile=None, parentPlugin=None):
        from wrappers.ctruncate.script import ctruncate
        from core import CCP4PluginScript
        error = CErrorReport()
        wrapper = ctruncate.ctruncate(self)
        if parentPlugin is not None:
            myDir = os.path.normpath(os.path.join(parentPlugin.workDirectory, 'ctruncate'))
            try:
                os.mkdir(myDir)
            except OSError as e:
                import errno
                if e.errno == errno.EEXIST:
                    import shutil
                    os.rename(myDir, os.path.join(parentPlugin.workDirectory, 'ctruncate_previous'))
                    shutil.rmtree(os.path.join(parentPlugin.workDirectory, 'ctruncate_previous'))
                    os.mkdir(myDir)
            wrapper.workDirectory = myDir
        inp = wrapper.container.inputData
        inp.HKLIN.setFullPath(self.fullPath.__str__())
        wrapper.container.controlParameters.OUTPUTMINIMTZ.set(True)
        wrapper.container.controlParameters.OUTPUTMINIMTZCONTENTFLAG.set(targetContent)
        if self.contentFlag == CObsDataFile.CONTENT_FLAG_IPAIR:
            inp.ISIGIanom.Ip = self.columnNames()[0]
            inp.ISIGIanom.SIGIp = self.columnNames()[1]
            inp.ISIGIanom.Im = self.columnNames()[2]
            inp.ISIGIanom.SIGIm = self.columnNames()[3]
        elif self.contentFlag == CObsDataFile.CONTENT_FLAG_IMEAN:
            inp.ISIGI.I = self.columnNames()[0]
            inp.ISIGI.SIGI = self.columnNames()[1]
        elif self.contentFlag == CObsDataFile.CONTENT_FLAG_FPAIR:
            wrapper.container.controlParameters.AMPLITUDES.set(True)
            inp.FSIGFanom.Fp = self.columnNames()[0]
            inp.FSIGFanom.SIGFp = self.columnNames()[1]
            inp.FSIGFanom.Fm = self.columnNames()[2]
            inp.FSIGFanom.SIGFm = self.columnNames()[3]
        wrapper.container.outputData.OBSOUT.setFullPath(targetFile)
        #MN Trying to get Ipair->Imean to work: have to tell ctruncate to output intensities
        if targetContent in [1, 3]:
            wrapper.container.controlParameters.OUTPUT_INTENSITIES.set(True)
        status = wrapper.process()
        if status != CCP4PluginScript.CPluginScript.SUCCEEDED:
            error.append(self.__class__, 301, name=self.objectName())
            return None, error
        else:
            return targetFile, error

    def convertObsMtz( self, infile=None, outfile=[], parentPlugin=None ):
        # Expect to use this to convert Fpairs to Fmean
        error = CErrorReport()
        cbin = os.path.normpath(os.path.join(CCP4Utils.getOSDir(), 'bin', 'cmtzsplit'))
        if not os.path.exists(cbin):
            cbin =  os.path.normpath(os.path.join(CCP4Utils.getCCP4Dir(), 'bin', 'cmtzsplit'))
        arglist = ['-mtzin', infile]
        if len(outfile) == 2:
            name, colin = outfile
            colout = ''
        else:
            name, colin, colout = outfile
        if parentPlugin is not None:
            logFile = os.path.normpath(os.path.join(self.parentPlugin.workDirectory, self.objectName() + '_splitMtz.log'))
        else:
            logFile = os.path.normpath(os.path.join(os.path.split(name)[0], self.objectName() + '_splitMtz.log'))
        arglist.append('-mtzout')
        arglist.append(name)
        arglist.append('-colin')
        arglist.append(colin)
        arglist.append('-colout')
        if len(colout) > 0:
            arglist.append(colout)
        else:
            arglist.append(colin)
        pid = CCP4Modules.PROCESSMANAGER().startProcess(cbin, arglist, logFile=logFile)
        status = CCP4Modules.PROCESSMANAGER().getJobData(pid)
        exitCode = CCP4Modules.PROCESSMANAGER().getJobData(pid, 'exitCode')
        if status == 0 and os.path.exists(outfile[0]):
            return outfile[0],error
        else:
            error.append(self, 302, self.__str__())
            return None, error


class CPhsDataFile(CMiniMtzDataFile):
    SUBTYPE_UNBIASED = 1
    SUBTYPE_BIASED = 2
    CONTENT_FLAG_HL = 1
    CONTENT_FLAG_PHIFOM = 2
    CONTENT_SIGNATURE_LIST = [['HLA','HLB','HLC','HLD'], ['PHI','FOM']]
    CONTENT_ANNOTATION = ['Hendrickson-Lattmann coeffs', 'Phi,FOM']
    CONTENTS = {}
    CONTENTS.update(CMiniMtzDataFile.CONTENTS)
    CONTENTS['subType'] = {'class' : CCP4Data.CInt, 'qualifiers' : {'default' : SUBTYPE_UNBIASED,
                                                                    'enumerators' : [SUBTYPE_UNBIASED,SUBTYPE_BIASED],
                                                                    'onlyEnumerators':True,
                                                                    'menuText' : ['unbiased data', 'biased data']}}
    QUALIFIERS = {'mimeTypeName' : 'application/CCP4-mtz-phases', 'mimeTypeDescription' : 'MTZ phases',
                  'fileExtensions' : ['mtz','cif','ent'], 'fileContentClassName' : 'CMtzData',
                  'guiLabel' : 'Phases', 'fileLabel' : 'phases',
                  'toolTip' : "Phases in Hendrickson-Lattmann or Phi/FOM form",
                  'correctColumns' : ['AAAA','PW'],
                  'columnGroupClassList' : [CHLColumnGroup, CPhiFomColumnGroup],
                  'helpFile' : 'data_files#Phs'}

    def contentSignature(self):
        return CPhsDataFile.CONTENT_SIGNATURE_LIST

    def columnNames(self,ifString=False,content=None):
        self.setContentFlag()
        if content is None:
            content = self.__dict__['_value']['contentFlag']
        if content == CPhsDataFile.CONTENT_FLAG_HL:
            if ifString:
                return 'HLA,HLB,HLC,HLD'
            else:
                return ['HLA', 'HLB', 'HLC', 'HLD']
        elif content == CPhsDataFile.CONTENT_FLAG_PHIFOM:
            if ifString:
                return 'PHI,FOM'
            else:
                return ['PHI', 'FOM']

    def conversion(self,targetContent):
        if isinstance(targetContent,list):
            import copy
            targetContentList = copy.deepcopy(targetContent)
        else:
            targetContentList = [targetContent]
        for targetContent in targetContentList:
            if targetContent == self.contentFlag:
                return 'ok', targetContent
            else:
                return 'convert', targetContent

    def convert(self, targetContent=None, targetFile=None, **kw):
        error = CErrorReport()
        if not self.isSet() or not self.exists():
            error.append(self.__class__, 220, name=self.objectPath())
            return None,error
        if targetContent is None or not isinstance(targetContent,int) or targetContent<1 or targetContent>4:
            error.append(self.__class__, 222, name=self.objectPath())
            return None,error
        if self.__dict__['_value']['contentFlag'] == targetContent:
            return self.fullPath.__str__(),error
        #MN CHECKME CHECK ME
        #I am going to disable this test...in fact Phi FOM can be "up-converted" to H L Coeffs
        #elif self.__dict__['_value']['contentFlag'] > targetContent:
        #  error.append(self.__class__,221,name=self.objectPath())
        #  return None,error
        #print 'CObsDataFile OK',self.__dict__['_value'],type(self.__dict__['_value'])
        if targetFile is None:
            fsplit = os.path.splitext(self.fullPath.__str__())
            mode = 'PHIFOM'
            targetFile = fsplit[0] + '_' + mode + fsplit[1]
        return self.runChltofom(targetContent=targetContent, targetFile=targetFile)

    def runChltofom(self, targetContent=None, targetFile=None):
        error = CErrorReport()
        from wrappers.chltofom.script import chltofom
        from core import CCP4PluginScript
        wrapper = chltofom.chltofom(self)
        wrapper.container.inputData.HKLIN.setFullPath(self.fullPath.__str__())
        wrapper.container.inputData.HKLIN.setContentFlag()
        wrapper.container.outputData.HKLOUT.setFullPath(targetFile)
        wrapper.container.controlParameters.OUTPUTMINIMTZ.set(True)
        if self.contentFlag == 2 and targetContent == 1:
            wrapper.container.controlParameters.DIRECTION = 'FOMTOHL'
        status = wrapper.process()
        if status != CCP4PluginScript.CPluginScript.SUCCEEDED:
            error.append(self.__class__, 301, name=self.objectName())
            return None, error
        else:
            return targetFile, error


class CMapCoeffsDataFile(CMiniMtzDataFile):
    SUBTYPE_NORMAL = 1
    SUBTYPE_DIFFERENCE = 2
    SUBTYPE_ANOM_DIFFERENCE = 3
    CONTENTS = {}
    CONTENTS.update(CMiniMtzDataFile.CONTENTS)
    CONTENT_FLAG_FPHI = 1
    CONTENT_SIGNATURE_LIST = [['F','PHI']]
    CONTENT_ANNOTATION = ['FPhi']
    CONTENTS['subType'] = {'class' : CCP4Data.CInt,
                           'qualifiers' : {'default': SUBTYPE_NORMAL, 'enumerators' : [SUBTYPE_NORMAL,SUBTYPE_DIFFERENCE, SUBTYPE_ANOM_DIFFERENCE],
                                           'onlyEnumerators' : True, 'menuText' : ['normal map', 'difference map', 'anomalous difference map']}}
    QUALIFIERS = {'mimeTypeName' : 'application/CCP4-mtz-map', 'mimeTypeDescription' : 'MTZ F-phi',
                  'fileExtensions' : ['mtz','cif','ent'], 'fileContentClassName' : 'CMtzData',
                  'fileLabel' : 'map_coefficients', 'guiLabel' : 'Map coefficients',
                  'toolTip' : "Electron density map coefficients: F,Phi",
                  'correctColumns' : ['FP','FQP'], 'columnGroupClassList' : [CMapColumnGroup],
                  'downloadModes' : ['Uppsala-EDS'], 'helpFile' : 'data_files#MapCoeffs'}

    def contentSignature(self):
        return CMapCoeffsDataFile.CONTENT_SIGNATURE_LIST

    def columnNames(self,ifString=False,content=None):
        if ifString:
            return 'F,PHI'
        else:
            return ['F','PHI']


class CFreeRDataFile(CMiniMtzDataFile):
    CONTENTS = {}
    CONTENTS.update(CMiniMtzDataFile.CONTENTS)
    CONTENT_SIGNATURE_LIST = [['FREER']]
    CONTENT_ANNOTATION = ['FreeR']
    CONTENTS['subType'] = {'class' : CCP4Data.CInt, 'qualifiers' : {'enumerators' : [], 'onlyEnumerators':True}}
    QUALIFIERS = {'mimeTypeName' : 'application/CCP4-mtz-freerflag', 'mimeTypeDescription' : 'FreeR flag',
                  'fileExtensions' : ['mtz','cif','ent'], 'fileContentClassName' : 'CMtzData',
                  'fileLabel' : 'freeRflag', 'guiLabel' : 'Free R set',
                  'toolTip' : "Set of reflections used for FreeR calculation", 'correctColumns' : ['I'],
                  'columnGroupClassList' : [CFreeRColumnGroup], 'helpFile' : 'data_files#FreeR'}

    def contentSignature(self):
        return CFreeRDataFile.CONTENT_SIGNATURE_LIST

    def columnNames(self,ifString=False,content=None):
        if ifString:
            return 'FREER'
        else:
            return ['FREER']

    def sameCrystal(self,other=None,testLevel=None):
        if self.fileContent is None:
            self.loadFile()
        if testLevel is None:
            testLevel = 1
        return self.getFileContent().sameCrystal(other.getFileContent(),testLevel)


class CShelxFADataFile(CCP4File.CDataFile):
    QUALIFIERS = {'mimeTypeName' : 'application/CCP4-shelx-FA',
                  'mimeTypeDescription' : 'Shelx FA',
                  'fileExtensions' : [ 'hkl' ],
                  'fileContentClassName' : None,
                  'fileLabel' : 'shelx_FA',
                  'guiLabel' : 'Shelx FA',
                  'toolTip' : "Data used by Shelx programs",
                  'helpFile' : 'data_files#shelxfa'}


class CPhaserSolDataFile(CCP4File.CDataFile):
    QUALIFIERS = {'mimeTypeName' : 'application/phaser-sol',
                  'mimeTypeDescription' : 'Phaser solution file',
                  'fileExtensions' : [ 'phaser_sol.pkl' ],
                  'fileContentClassName' : None,
                  'fileLabel' : 'phaser_sol',
                  'guiLabel' : 'Phaser solutions',
                  'toolTip' : "Possible solutions passed between runs of the Phaser program",
                  'helpFile' : 'data_files#phasersol'}

    def assertSame(self,other,diagnostic=False,**kw):
        try:
            import pickle
            selfStr = pickle.load(open(self.__str__())).unparse()
            print('CPhaserSolDataFile.assertSame selfStr',selfStr)
            otherStr = pickle.load(open(other.__str__())).unparse()
            print('CPhaserSolDataFile.assertSame otherStr',otherStr)
        except:
            return CErrorReport(self.__class__,301,name=self.objectPath(),details=str(self)+' : '+str(other))
        # Unsophisicated diff
        if sys.version_info > (3,0):
            from googlecode import diff_match_patch_py3
            dmp =  diff_match_patch_py3.diff_match_patch()
        else:
            from googlecode import diff_match_patch
            dmp =  diff_match_patch.diff_match_patch()
        diffs = dmp.diff_main(selfStr,otherStr)
        print('CPhaserSolDataFile.assertSame diffs', diffs)
        if len(diffs) > 1:
            return CErrorReport(self.__class__,313,name=self.objectPath(),details=str(self)+' : '+str(other))
        else:
            return CErrorReport(self.__class__,300,name=self.objectPath())

class CPhaserRFileDataFile(CCP4File.CDataFile):
    QUALIFIERS = {'mimeTypeName' : 'application/phaser-rfile',
                  'mimeTypeDescription' : 'Phaser rotation solution file',
                  'fileExtensions' : ['phaser_rlist.pkl'],
                  'fileContentClassName' : None,
                  'fileLabel' : 'phaser_rfile',
                  'guiLabel' : 'Phaser rotation solution',
                  'toolTip' : "Phaser rfile solutions for rotation search"}

class CRefmacKeywordFile(CCP4File.CDataFile):
    QUALIFIERS = {'mimeTypeName' : 'application/refmac-keywords',
                  'mimeTypeDescription' : 'Refmac keyword file',
                  'fileExtensions' : ['txt'],
                  'fileContentClassName' : None,
                  'fileLabel' : 'refmac_keywords',
                  'guiLabel' : 'Refmac keyword file',
                  'toolTip' : "A file containing keywords as they are meant to be read by refmac5"}

class CDialsJsonFile(CCP4File.CDataFile):
    QUALIFIERS = {'mimeTypeName' : 'application/dials-jfile',
                  'mimeTypeDescription' : 'Dials json data file',
                  'fileExtensions' : ['json','expt','jsn'],
                  'fileContentClassName' : None,
                  'fileLabel' : 'dials_jdata',
                  'guiLabel' : 'json data',
                  'toolTip' : "json data files"}

class CDialsPickleFile(CCP4File.CDataFile):
    QUALIFIERS = {'mimeTypeName' : 'application/dials-pfile',
                  'mimeTypeDescription' : 'Dials pickle data file',
                  'fileExtensions' : ['pickle','refl'],
                  'fileContentClassName' : None,
                  'fileLabel' : 'dials_pdata',
                  'guiLabel' : 'Xia2/Dials pickle data',
                  'toolTip' : "Xia2/Dials pickle data files"}

class CImosflmXmlDataFile(CCP4File.CDataFile):
    '''An iMosflm data file'''
    QUALIFIERS = {'fileLabel' : 'imosflm',
                  'mimeTypeName' : 'application/iMosflm-xml',
                  'mimeTypeDescription' : 'iMosflm data',
                  'guiLabel' : 'iMosflm data',
                  'fileExtensions' : ['imosflm.xml'],
                  'fileContentClassName' : None }

class CMergeMiniMtz(CCP4Data.CData):
    CONTENTS = {'fileName' : {'class' : CMiniMtzDataFile, 'qualifiers' : { 'fromPreviousJob' : False}},
                'columnTag' : {'class' : CCP4Data.CString}, 'columnNames' : {'class' : CCP4Data.CString}}
    CONTENTS_ORDER = ['fileName', 'columnTag', 'columnNames']
    ERROR_CODES = {201 : { 'description' : "Selected file is not a suitable 'mini' MTZ containing experimental data object"},
                   202 : { 'description' : 'Output column name list does not have correct number of names'}}

    def qualifiers(self,name=None,default=True,custom=True,contentQualifiers=True):
        if name == 'mimeTypeName':
            return "application/CCP4-mtz-mini"
        else:
            return  CCP4Data.CData.qualifiers(self,name=name,default=default,custom=custom,contentQualifiers=contentQualifiers)

    def validity(self,arg):
        v = CErrorReport()
        if self.__dict__['_value']['fileName'].isSet() and self.__dict__['_value']['fileName'].exists():
            cls,contentFlag = self.__dict__['_value']['fileName'].miniMtzType()
            if cls is None:
                v.append(self.__class__, 201, name = self.objectPath(), label = self.qualifiers('guiLabel'), stack = False)
            elif self.__dict__['_value']['columnNames'].isSet() :
                stdColumnNames = cls().columnNames(True,contentFlag)
                if  self.__dict__['_value']['columnNames'].__str__().count(',') != stdColumnNames.count(','):
                    v.append(self.__class__, 202, name = self.objectPath(),label = self.qualifiers('guiLabel'), stack = False)
        return v

    def getTableTextItems(self):
        return [self.__dict__['_value']['fileName'].guiLabel(useAnnotation=False,useObjectName=False),self.__dict__['_value']['columnNames'].__str__() ]

    def setColumnTag(self,overwrite=False):
        if not self.__dict__['_value']['fileName'].isSet():
            self.__dict__['_value']['columnTag'].unSet()
        if not overwrite and self.__dict__['_value']['columnTag'].isSet():
            return
        if not self.__dict__['_value']['fileName'].dbFileId.isSet():
            self.__dict__['_value']['columnTag'].unSet()
            return
        else:
            fileInfo = CCP4Modules.PROJECTSMANAGER().db().getFileInfo(fileId = self.__dict__['_value']['fileName'].dbFileId.__str__(),
                                                                      mode=['jobnumber','jobparamname','taskname'])
            if fileInfo['taskname'] is not None:
                self.__dict__['_value']['columnTag'].set(fileInfo['jobnumber'] + '_' + fileInfo['taskname'][0:20])
            else:
                self.__dict__['_value']['columnTag'].set(fileInfo['jobnumber'][0:20])

    def setColumnNames(self,mode='fromFile',overwrite=False):
        if not self.__dict__['_value']['fileName'].isSet():
            self.__dict__['_value']['columnNames'].unSet()
            return
        if not overwrite and self.__dict__['_value']['columnNames'].isSet():
            return
        fileColumnNames= self.__dict__['_value']['fileName'].fileContent.columnNames()
        if mode == 'applyTag' and self.__dict__['_value']['columnTag'].isSet():
            cls,contents = self.__dict__['_value']['fileName'].miniMtzType()
            if cls is not None:
                columnNames = cls.CONTENT_SIGNATURE_LIST[contents-1]
            else:
                columnNames = fileColumnNames
            if cls == CMapCoeffsDataFile and self.__dict__['_value']['fileName'].subType.isSet() and self.__dict__['_value']['fileName'].subType==2:
                columnNames = [ 'DIFF_'+columnNames[0],'DIFF_'+columnNames[1]]
            text = ''
            for name in columnNames:
                text = text + name+'_'+self.__dict__['_value']['columnTag'].__str__()+','
            self.__dict__['_value']['columnNames'].set(text[0:-1])
        else:
            text = str(fileColumnNames[0])
            for name in fileColumnNames[1:]:
                text = text + ',' + str(name)
            self.__dict__['_value']['columnNames'].set(text)


class CMergeMiniMtzList(CCP4Data.CList):
    SUBITEM = {'class' : CMergeMiniMtz}
    QUALIFIERS = {'listMinLength' : 2, 'saveToDb' : True }

    def saveToDb(self):
        fileObjList = []
        for item in self.__dict__['_value']:
            fileObjList.append(item.fileName)
        return fileObjList, None, {}


class CRunBatchRange(CCP4Data.CData):
    CONTENTS = {'runNumber' : {'class' : CCP4Data.CInt, 'qualifiers' : {'allowUndefined' : True, 'min' : 1} },
                'batchRange0' : {'class' : CCP4Data.CInt, 'qualifiers' : {'allowUndefined' : True, 'min' : 1}},
                'batchRange1' : {'class' : CCP4Data.CInt, 'qualifiers' : {'allowUndefined' : True, 'min' : 1}},
                'resolution' : {'class' : CCP4Data.CFloat , 'qualifiers' : {'min' : 0.0, 'allowUndefined' : True}},
                'fileNumber' : {'class' : CCP4Data.CInt, 'qualifiers' : {'allowUndefined' : True, 'min' : 1}}}
    QUALIFIERS = {'toolTip' : 'Specify range of reflections to treat as one run'}
    ERROR_CODES = {101 : {'description' : 'End of batch range less than start'},
                   102 : {'description' : 'All items must be set'}}

    def validity(self,arg):
        # Check validity of the individual items
        v = self.itemValidity(arg)
        if v.maxSeverity() > SEVERITY_WARNING:
            return v
        if not isinstance(arg,dict):
            arg = arg.get()
        #All or none must be set
        nSet = 0
        for item in ['runNumber','batchRange0','batchRange1']:
            if arg[item] is not None:
                nSet+=1
        if nSet not in [0,3]:
            v.append(self.__class__, 102, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        if v.maxSeverity() > 0:
            return v
        # Batch range
        if arg.get('batchRange1').__cmp__(arg.get('batchRange0')) < 0:
            v.append(self.__class__, 101, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        return v

    def getTableTextItems(self):
        ret = []
        for key in ['runNumber','batchRange0','batchRange1','resolution']:
            if not self.__dict__['_value'][key].isSet():
                ret.append('')
            else:
                ret.append( str(self.__dict__['_value'][key]))
        return ret

class CRunBatchRangeList(CCP4Data.CList):
    SUBITEM = {'class' : CRunBatchRange}
    QUALIFIERS = {'listMinLength' : 1}


class CDataset(CCP4Data.CData):
    '''
    The experimental data model for ab initio phasing
    '''
    CONTENTS = {'selected' : {'class' : CCP4Data.CBoolean }, 'obsDataFile' : {'class' : CObsDataFile },
                'crystalName' : {'class' : CCrystalName }, 'datasetName' : {'class' : CDatasetName },
                'formFactors' : {'class' : CFormFactor },
                'formFactorSource' : {'class' : CCP4Data.CString,
                                      'qualifiers' : {'onlyEnumerators' : True,
                                                      'enumerators' : [ 'no' , 'composition', 'xia2'],
                                                      'menuText' : [ 'user input' , 'atomic composition', 'from XIA2'],
                                                      'default' : 'no'}}}

    def getTextItem(self):
        return self.obsDataFile.baseName.__str__()

    def getTableTextItems(self):
        formText = ''
        if self.formFactors.Fp.isSet:
            formText+=str(self.formFactors.Fp)
        formText+=','
        if self.formFactors.Fpp.isSet:
            formText+=str(self.formFactors.Fpp)
        return [self.obsDataFile.baseName.__str__(),self.crystalName.__str__(),self.datasetName.__str__(),formText]


class CDatasetList(CCP4Data.CList):
    SUBITEM = {'class' : CDataset}


#===========================================================================================================
import unittest
def TESTSUITE():
    '''
    suite = unittest.TestLoader().loadTestsFromTestCase(testAssorted)
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testCell))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testMtz))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testSpaceGroup))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testComponents))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testComposition))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testCObsDataFile))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testCPhsDataFile))
    '''
    suite= unittest.TestLoader().loadTestsFromTestCase(testCPhsDataFile)
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)

class testAssorted(unittest.TestCase):

    def setUp(self):
        if QT():
            from PySide2 import QtCore
            self.app = CCP4Modules.QTAPPLICATION()
        from core import CCP4Container
        self.mummy = CCP4Container.CContainer()

    def testMtzColumn(self):
        t = CMtzColumn(columnLabel='foo',columnType='A',parent=self.mummy)
        self.assertEqual(t.columnLabel,'foo','Error instantiating CMtzColumn')
        try:
            t = CMtzColumn(columnLabel='foo',columnType='Z',parent=self.mummy)
        except CException as e:
            self.assertEqual(e[0]['code'],103,'Wrong error when instantiating CMtzColumn with bad column type')
        except:
            self.fail('Unexpected exception when instantiating CMtzColumn with bad column type')

    def testMtzData(self):
        filename =  os.path.normpath(os.path.join(CCP4Utils.getCCP4I2Dir(),'test','data','gere_nat.mtz'))
        t = CMtzData(parent=self.mummy, name='foo')
        t.loadFile(filename)
        self.assertEqual(19, t.getNColumns(), 'CMtzData loaded MTZ reports wrong number of columns')

class testMtz(unittest.TestCase):
    def setUp(self):
        self.testDataDir = os.path.normpath(os.path.join(CCP4Utils.getCCP4I2Dir(),'test','data'))
        # make all background jobs wait for completion
        CCP4Modules.PROCESSMANAGER().setWaitForFinished(10000)
        if QT():
            from PySide2 import QtCore
            self.app = CCP4Modules.QTAPPLICATION()
            self.mummy = QtCore.QObject(self.app)
        else:
            self.mummy = None

    def tearDown(self):
        CCP4Modules.PROCESSMANAGER().setWaitForFinished(-1)

    def test_1(self):
        self.mtz = CMtzDataFile( os.path.normpath(os.path.join(self.testDataDir,'gere_nat.mtz'),parent=self.mummy))
        print("test_1 getNColumns",self.mtz.getFileContent().getNColumns())
        self.assertEqual( 19, self.mtz.getFileContent().getNColumns(),'CMtzDataFile loaded MTZ reports wrong number of columns')

    def test_2(self):
        from core import CCP4Container
        self.dataContainer = CCP4Container.CContainer()
        self.dataContainer.loadContentsFromXml( os.path.normpath(os.path.join(self.testDataDir,'test_mtz_2.def.xml')))
        columns = self.dataContainer.testCProgramColumnGroup.F_SIGF.qualifiers('columnGroup')
        self.assertEqual(str(columns[0].columnName),'F','CProgramColumnGroup failed to load from test_mtz_2.def.xml')
        self.assertEqual(str(columns[1].columnName),'SIGF','CProgramColumnGroup failed to load from test_mtz_2.def.xml')
        self.dataContainer.loadDataFromXml( os.path.normpath(os.path.join(self.testDataDir,'test_mtz_2.params.xml')))
        #print 'test_2',self.dataContainer
        dataF = self.dataContainer.inputData.F_SIGF.F
        self.assertEqual(dataF,'F_nat','CProgramColumnGroup failed to load from test_mtz_2.params.xml')

    def test_3(self):
        from core import CCP4Container
        self.dataContainer = CCP4Container.CContainer()
        self.dataContainer.loadContentsFromXml( os.path.normpath(os.path.join(self.testDataDir,'test_mtz_2.def.xml')))
        # def file does not have SIGF defined
        self.dataContainer.loadDataFromXml( os.path.normpath(os.path.join(self.testDataDir,'test_mtz_3.params.xml')))
        fixed = self.dataContainer.inputData.F_SIGF.fix(self.dataContainer.testCProgramColumnGroup.F_SIGF.get())
        #print 'test_3 fixed',fixed,self.dataContainer.testCProgramColumnGroup.F_SIGF
        self.assertEqual(str(self.dataContainer.testCProgramColumnGroup.F_SIGF.SIGF),'SIGF_nat','CProgramColumnGroup.Partner failed to find unset column')


class testCell(unittest.TestCase):
    def testLength1(self):
        l = CCellLength(56.8)
        m = CCellLength()
        m.nm = 5.69
        self.assertEqual(l > m,False, 'Comparison of CCellLength failed')
        self.assertEqual(l + 0.2 > m,True, 'Addition and comparison of CCellLength failed')

    def testAngle1(self):
        t = CCellAngle(90.0)
        if t.rad<math.pi/2.0-0.001 or t.rad>math.pi/2.0+0.001:
            self.fail('Error return CCellAngle as radians')

    def testAngle2(self):
        try:
            t = CCellAngle(-56.0)
        except CException as e:
            self.assertEqual(len(e),1,'Unexpected exception length in setting  CCellAngle')
            self.assertEqual(e[0]['code'],101,'Unexpected exception in setting  CCellAngle')
        except:
            self.fail('Unexpected exception in setting  CCellAngle')
        else:
            self.fail('No exception in setting CCellAngle')

    def testCell1(self):
        c = CCell(a=78.0,b=56.0,c=13.5,alpha=89.0,beta=93.8,gamma=103.4)
        if c.a.nm > 7.81 or c.a.nm < 7.79:
            self.fail('Error returning cell length as nm')

class testSpaceGroup(unittest.TestCase):

    def setUp(self):
        self.symMan = SYMMETRYMANAGER()

    def test1(self):
        # test that the hard-coded chiral space groups match to xHM names in syminfo.lib
        for xSys in self.symMan.crystalSystems:
            for sgp in self.symMan.chiralSpaceGroups[xSys]:
                status,newSgp = self.symMan.spaceGroupValidity(sgp)
                if status == 5:
                    newSgpChiral = []
                    for item in newSgp:
                        ii = self.symMan.hmSpaceGroupList.index(item)
                        if self.symMan.pointGroupList[ii].count('-') == 0:
                            newSgpChiral.append(item)
                    print(sgp,'*',status,'*',newSgpChiral)
                elif status != 0:
                    self.fail('SYMMETRYMANAGER chrial space group name not found in syminfo.lib:'+sgp)

    def test2(self):
        s = CSpaceGroup()
        for sgp,expectedErr,expectedFix in [['P 21 21 21' , None, 'P 21 21 21'],
                                            ['P -1', 102,'P -1' ],
                                            ['P4/n b m', 103, 'P 4/n b m :1'],
                                            ['P21 1 1', 105 ,'P 21 1 1']]:
            rv = s.validity( sgp )
            if len(rv) == 0:
                if expectedErr is not None:
                    self.fail('No validity fail for CSpaceGroup:'+sgp)
            elif len(rv) > 1:
                self.fail('CErrorReport for CSpaceGroup longer than 1:'+sgp)
            elif rv[0]['code'] != expectedErr:
                self.fail('CErrorReport for CSpaceGroup does not give expected error:'+sgp)
            fix = s.fix(sgp)
            #print 'test2',sgp,fix
            if fix != expectedFix:
                self.fail('Incorrect CSpaceGroup.fix() for:'+sgp)

    def test3(self):
        s = CSpaceGroupCell()
        for sgp,cell,expectedErr in [
          ['P 21 21 21', {'a': 64.900, 'b':78.320, 'c':38.790, 'alpha':90.00, 'beta':90.00, 'gamma': 90.00}, None],
          ['P 21 21 21', {'a': 64.900, 'b':78.320, 'c':38.790, 'alpha':90.00, 'beta':91.00, 'gamma': 90.00}, 103],
          ['P 21 21 21', {'a': 64.900, 'b':78.320, 'c':64.900, 'alpha':90.00, 'beta':90.00, 'gamma': 90.00}, 101],
          ['P 21 1 1', {'a': 64.900, 'b':78.320, 'c':38.790, 'alpha':90.00, 'beta':90.00, 'gamma': 90.00}, 104]]:
            rv = s.validity({'spaceGroup' : sgp, 'cell' : cell})
            #print rv.report()
            if len(rv) == 0:
                if expectedErr is not None:
                    self.fail('No validity fail for CSpaceGroupCell:'+sgp)
            elif len(rv)>1:
                self.fail('CErrorReport for CSpaceGroupCell longer than 1:'+sgp)
            elif rv[0]['code'] != expectedErr:
                self.fail('CErrorReport for CSpaceGroupCell does not give expected error:'+sgp)


class testCObsDataFile(unittest.TestCase):
    def setUp(self):
        self.testDataDir =  os.path.normpath(os.path.join(CCP4Utils.getCCP4I2Dir(),'test','data'))
        self.app = CCP4Modules.QTAPPLICATION()
        from PySide2 import QtCore
        self.mummy = QtCore.QObject(self.app)

    def test1(self):
        self.obs = CObsDataFile(parent=self.mummy,fullPath= os.path.normpath(os.path.join(self.testDataDir,'rnase_obs_fpair.mtz')))
        self.obs.setContentFlag()
        print('testCObsDataFile.test1 contentFlag',self.obs.contentFlag)
        outfile = os.path.normpath(os.path.join(CCP4Utils.getTestTmpDir(),'testCObsDataFile.mtz'))
        if os.path.exists(outfile):
            os.remove(outfile)
        self.obs.convert(targetContent=CObsDataFile.CONTENT_FLAG_FMEAN,targetFile=outfile)
        self.assertTrue(os.path.exists(outfile),'CObsDataFile.convert failed')

class testCPhsDataFile(unittest.TestCase):
    def setUp(self):
        self.testDataDir =  os.path.normpath(os.path.join(CCP4Utils.getCCP4I2Dir(),'test','data'))
        self.app = CCP4Modules.QTAPPLICATION()
        from PySide2 import QtCore
        self.mummy = QtCore.QObject(self.app)

    def test1(self):
        self.phs = CPhsDataFile(parent=self.mummy,fullPath= os.path.normpath(os.path.join(self.testDataDir,'rnase25_mini_HL.mtz')))
        self.phs.setContentFlag()
        print('testCPhsDataFile.test1 contentFlag',self.phs.contentFlag)
        outfile =  os.path.normpath(os.path.join(CCP4Utils.getTestTmpDir(),'testCPhsDataFile.mtz'))
        if os.path.exists(outfile):
            os.remove(outfile)
        self.phs.convert(targetContent=CPhsDataFile.CONTENT_FLAG_PHIFOM,targetFile=outfile)
        self.assertTrue(os.path.exists(outfile),'CPhsDataFile.convert failed')
        # beware self.phs suddenly is something else..
        self.phs.setFullPath(outfile)
        self.phs.loadFile()
        columns = self.phs.fileContent.getListOfColumns()
        print('testCPhsDataFile.test1',columns)
        self.assertEqual(len(columns),2,'Output from testCPhsDataFile has wrong number of columns')
        self.assertTrue(columns[0].columnLabel.__str__() in ['PHI','FOM'] and columns[1].columnLabel.__str__() in ['PHI','FOM'],'Output from testCPhsDataFile has wrong column labels')
