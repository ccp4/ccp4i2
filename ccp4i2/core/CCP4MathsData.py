"""
Copyright (C) 2011 University of York
Liz Potterton Feb 2011 - CData subclasses for basic maths/geometry
"""

import math

from . import CCP4Data
from .CCP4ErrorHandling import CErrorReport, CException


class CXyz(CCP4Data.CData):

    CONTENTS = {'x' : {'class' : CCP4Data.CFloat}, 'y' : {'class' : CCP4Data.CFloat}, 'z' : {'class' : CCP4Data.CFloat}}
    ERROR_CODES = {201 : {'description' :'Attempting arithmetic with inappropriate data type'},
                   202 : {'description' :'Attempting arithmetic in unset data object'},
                   203 : {'description' :'Attempting arithmetic with unset data object as argument'},}

    def add(self,arg):
        if not isinstance(arg,CXyz):
            raise CException(self.__class__, 201, type(arg), name=self.objectPath())
        if not self.isSet():
            raise CException(self.__class__, 202, name=self.objectPath())
        if not arg.isSet():
            raise CException(self.__class__, 203, name=self.objectPath())
        self.x = self.x + arg.x
        self.y = self.y + arg.y
        self.z = self.z + arg.z
    
    def sub(self, arg):
        if not isinstance(arg,CXyz):
            raise CException(self.__class__, 201, type(arg))
        if not self.isSet():
            raise CException(self.__class__, 202, name=self.objectPath())
        if not arg.isSet():
            raise CException(self.__class__, 203, name=self.objectPath())   
        self.x = self.x - arg.x
        self.y = self.y - arg.y
        self.z = self.z - arg.z

    def div(self, arg):
        try:
            a = float(arg)
        except:
            raise CException(self.__class__, 201, type(arg))
        if not self.isSet(): raise CException(self.__class__, 202, name=self.objectPath())
        self.x = self.x / a
        self.y = self.y / a
        self.z = self.z / a

    def mul(self,arg):
        try:
            a = float(arg)
        except:
            raise CException(self.__class__,201,type(arg),name=self.objectPath())
        if not self.isSet(): raise CException(self.__class__,202,name=self.objectPath())
        self.x = self.x * a
        self.y = self.y * a
        self.z = self.z * a

    def __add__(self,arg):
        if not isinstance(arg,CXyz):
            raise CException(self.__class__,201,type(arg),name=self.objectPath())
        if not self.isSet(): raise CException(self.__class__,202,name=self.objectPath())
        if not arg.isSet(): raise CException(self.__class__,203,name=self.objectPath())     
        rv = CXyz(x=self.x+arg.x, y=self.y+arg.y, z=self.z+arg.z)
        return rv

    def __sub__(self,arg):
        if not isinstance(arg,CXyz):
            raise CException(self.__class__,201,type(arg),name=self.objectPath())
        if not self.isSet(): raise CException(self.__class__,202,name=self.objectPath())
        if not arg.isSet(): raise CException(self.__class__,203,name=self.objectPath())     
        rv = CXyz(x=self.x-arg.x, y=self.y-arg.y, z=self.z-arg.z)
        return rv

    def __div__(self,arg):
        try:
            a = float(arg)
        except:
            raise CException(self.__class__,201,type(arg),name=self.objectPath())
        if not self.isSet(): raise CException(self.__class__,202,name=self.objectPath())
        rv = CXyz(x=self.x/a, y=self.y/a, z=self.z/a)
        return rv

    def __mul__(self,arg):
        try:
            a = float(arg)
        except:
            raise CException(self.__class__,201,type(arg),name=self.objectPath())
        if not self.isSet(): raise CException(self.__class__,202,name=self.objectPath())
        rv = CXyz(x=self.x*a, y=self.y*a, z=self.z*a)
        return rv

    def length(self):
        if not self.isSet(): raise CException(self.__class__,202,name=self.objectPath())
        l2 = (self.x**2) +  (self.x**2) +  (self.x**2)
        return l2**0.5

class CXyzBox(CCP4Data.CData):

    CONTENTS = {'xMin' : {'class' : CCP4Data.CFloat},
                'yMin' : {'class' : CCP4Data.CFloat},
                'zMin' : {'class' : CCP4Data.CFloat},
                'xMax' : {'class' : CCP4Data.CFloat},
                'yMax' : {'class' : CCP4Data.CFloat},
                'zMax' : {'class' : CCP4Data.CFloat}}
    ERROR_CODES = {201 : {'description' :'Maximum x,y or z value less than minimum'}}

    def validity(self,arg):
        v = CErrorReport()
        if self.isSet():
            if (arg['xMin'] > arg['xMax']) or (arg['yMin'] > arg['yMax']) or (arg['zMin'] > arg['zMax']):
                v.append(self.__class__,201,name=self.objectPath(),label=self.qualifiers('guiLabel'),stack=False)
        return v

    def fix(self,arg):
        ret = {}
        for x in ['x','y','z']:
            if arg.get(x+'Min') > arg.get(x+'Max'):
                ret[x+'Min'] = arg[x+'Max']
                ret[x+'Max'] = arg[x+'Min']
            else:
                ret[x+'Min'] = arg[x+'Min']
                ret[x+'Max'] = arg[x+'Max']
        return ret

class CAngle(CCP4Data.CFloat):
    '''An angle'''

    def getRadians(self):
        if self._value is None:
            return 0.0
        else:
            return self._value * (math.pi/180.0)

    def setRadians(self,value):
        if value is None:
            self.set(value)
        else:
            self.set(value * (180.0/math.pi))

    PROPERTIES = {'rad' : {'fget' : getRadians , 'fset' : setRadians}}


class CEulerRotation(CCP4Data.CData):

    CONTENTS = {'alpha' : {'class' : CAngle},
                'beta' : {'class' : CAngle},
                'gamma' : {'class' : CAngle}}
    CONTENTS_ORDER = ['alpha', 'beta', 'gamma']


class CTransformation(CCP4Data.CData):
    CONTENTS = {'translation' : {'class' : CXyz},
                'rotation' : {'class' : CEulerRotation}}
    CONTENTS_ORDER = ['translation', 'rotation']

    def __getattr__(self, name):
        if name in ['x','y','z']:
            return self.__dict__['_value']['translation'].__getattr__(name)
        if name in ['alpha','beta','gamma']:
            return self.__dict__['_value']['rotation'].__getattr__(name)
        return CCP4Data.CData.__getattr__(self,name)

    def __setattr__(self,name,value):
        if name in ['x', 'y', 'z']:
            return self.__dict__['_value']['translation'].__setattr__(name,value)
        if name in ['alpha', 'beta', 'gamma']:
            return self.__dict__['_value']['rotation'].__setattr__(name,value)
        return CCP4Data.CData.__setattr__(self,name,value)

class CMatrix33(CCP4Data.CData):
    pass
