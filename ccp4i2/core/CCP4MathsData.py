"""
     CCP4MathsData.py: CCP4 GUI Project
     Copyright (C) 2011 University of York

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
   Liz Potterton Feb 2011 - CData subclasses for basic maths/geometry
"""


from core import CCP4Data
from core.CCP4ErrorHandling import *
from core.CCP4Config import QT, XMLPARSER
if QT():
    from core.CCP4QtObject import CObject
else:
    from core.CCP4Object import CObject
if XMLPARSER() == 'lxml':
    from lxml import etree

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
        import math
        if self._value is None:
            return 0.0
        else:
            return self._value * (math.pi/180.0)

    def setRadians(self,value):
        import math
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


#===========================================================================================================
import unittest
def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testXyz)
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testXyzBox))
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)

class testXyz(unittest.TestCase):

    def test1(self):
        a = CXyz(x=1.1,y=-2.2,z=3.3)
        a.mul(2)
        self.assertAlmostEqual(float(a.x),2.2,3,'Failed CXyz.mul x value')
        self.assertAlmostEqual(float(a.y),-4.4,3,'Failed CXyz.mul y value')
        self.assertAlmostEqual(float(a.z),6.6,3,'Failed CXyz.mul z value')

    def test2(self):
        a = CXyz(x=1.1,y=-2.2,z=3.3)
        b = a*2
        self.assertAlmostEqual(float(b.x),2.2,3,'Failed CXyz.__mul__ x value')
        self.assertAlmostEqual(float(b.y),-4.4,3,'Failed CXyz.__mul__ y value')
        self.assertAlmostEqual(float(b.z),6.6,3,'Failed CXyz.__mul__ z value')

    def test3(self):
        a = CXyz(x=1.1,y=-2.2,z=3.3)
        c = CXyz(x=10.0,y=-20.0,z=30.0)
        b = a + c
        self.assertAlmostEqual(float(b.x),11.1,3,'Failed CXyz.__add__ x value')
        self.assertAlmostEqual(float(b.y),-22.2,3,'Failed CXyz.__add__ y value')
        self.assertAlmostEqual(float(b.z),33.3,3,'Failed CXyz.__add__ z value')

    def test4(self):
        a = CXyz(x=1.1,y=-2.2,z=3.3)
        c = CXyz(x=10.0,y=-20.0,z=30.0)
        b = a + c
        self.assertAlmostEqual(float(b.x),11.1,3,'Failed CXyz.__add__ x value')
        self.assertAlmostEqual(float(b.y),-22.2,3,'Failed CXyz.__add__ y value')
        self.assertAlmostEqual(float(b.z),33.3,3,'Failed CXyz.__add__ z value')

    def test5(self):
        a = CXyz()
        c = CXyz(x=10.0,y=-20.0,z=30.0)
        e = None
        try:
            b = a + c
        except CException as e:
            pass
        if e is None: self.fail('No error code when expecting 202')
        self.assertEqual(e[0]['code'],202,'Wrong error code when expecting 202')

    def test6(self):
        a = CXyz()
        c = CXyz(x=10.0,y=-20.0,z=30.0)
        e = None
        try:
            b = c + a
        except CException as e:
            pass
        if e is None: self.fail('No error code when expecting 203')
        self.assertEqual(e[0]['code'],203,'Wrong error code when expecting 203')

    def test7(self):
        a = CXyz(x=10.0,y=-20.0,z=30.0)
        e = None
        try:
            b = a / 'foo'
        except CException as e:
            pass
        if e is None: self.fail('No error code when expecting 203')
        self.assertEqual(e[0]['code'],201,'Wrong error code when expecting 201')


class testXyzBox(unittest.TestCase):

    def test1(self):
        d = { 'xMin' : 12, 'xMax' : 6, 'yMin' : -3, 'yMax' : 27, 'zMin' : 0, 'zMax': 20 }
        try:
            a = CXyzBox(d)
        except CException as e:
            self.assertEqual(e[0]['code'],201,'Wrong validity return expecting 201')
        except:
            self.fail('Unexpected exception in setting CXyzBox')
        else:
            self.fail('No exception in setting CXyzBox')
        a = CXyzBox()
        a.set(a.fix(d))
        self.assertEqual(a.xMin,6,'Error in CXyzBox.fix()')

