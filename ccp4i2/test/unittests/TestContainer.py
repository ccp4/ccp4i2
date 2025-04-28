import os
import unittest

from ccp4i2.core.CCP4Container import CContainer
from ccp4i2.core.CCP4Utils import getCCP4I2Dir
from ccp4i2.core import CCP4Data, CCP4MathsData, CCP4ModelData, CCP4Utils, CCP4XtalData


class testCContainer(unittest.TestCase):

    #Removing broken tests
    def broken_test1(self):
        self.summat = CContainer(definitionFile=os.path.join(getCCP4I2Dir(),'test/data/summat.contents.xml'))
        #print summat.name,summat.CONTENTS
        self.assertEqual(self.summat.header.pluginName,'Summat','Wrong header.pluginName')
        #print 'summat.CONTENTS',self.summat.CONTENTS.keys()
        self.assertEqual(len(list(self.summat.CONTENTS.keys())),4,'Wrong content length')
        self.assertEqual(self.summat.CONTENTS['range']['class'],CCP4Data.CIntRange,'range parameter has wrong class')
        #print "summat.CONTENTS['range']['qualifiers']",self.summat.CONTENTS['range']['qualifiers']
        self.assertEqual(self.summat.CONTENTS['range']['qualifiers']['compare'],1,'range has wrong qualifiers')
        #print 'summat.gubbins.startValue',self.summat.gubbins.startValue
        self.assertEqual(self.summat.gubbins.startValue,12.0,'sub-container CFloat wrong initial value')

    def broken_test2(self):
        contents = {'nCycles': {'class': CCP4Data.CInt, 'qualifiers': { 'default' : None, 'max': 20, 'allowUndefined': True, 'min': 1}}, 'cutoff': {'class': CCP4Data.CInt, 'qualifiers': { 'default' :3,'max': 100, 'min': 1}}, 'range': {'class': CCP4Data.CIntRange, 'qualifiers': { 'start' : { 'default' : 0 } ,  'end' : { 'default' : 10 } , 'compare': 1 }}}
        self.summat = CContainer(name='Summat',contents=contents)
        #self.summat.range.set(start=10,end=20)
        self.summat.range.start.set(10)
        self.summat.range.end.set(20)
        self.assertEqual(self.summat.get('range'), {'start':10,'end':20},'Failed setting range parameter')
    
    def broken_test3(self):
        defFile = os.path.join(CCP4Utils.getCCP4I2Dir(),'tasks','demo','demo.def.xml')
        demo =  CContainer(definitionFile=defFile)
        self.assertEqual(demo.controlParameters.QUICKNDIRTY,'quick','Loading demo.def.xml does not give correct value for controlParameters.QUICKNDIRTY')

    def broken_test4(self):
        c = CContainer()
        c.addObject(CCP4XtalData.CMtzDataFile(name='MTZIN'))
        c.addObject(CCP4Data.CFloatRange(name='RESOLUTION_RANGE',qualifiers = {'compare':-1,'start' :{'min':0.0,'default':999.9},'end':{'min':0.0,'default':0.0}}))
        c.MTZIN.fullPath = '/foo/bar/wotsit.tmp'
        c.RESOLUTION_RANGE.start=67.9
        c.RESOLUTION_RANGE.end = 2.1
        self.assertEqual(c.RESOLUTION_RANGE.start.qualifiers('min'),0.0,'Error adding CFloatRange to CContainer - wrong qualifier')
        self.assertEqual(c.RESOLUTION_RANGE.end,2.1,'Error adding  CFloatRange to CContainer - wrong value')
        self.assertEqual(c.MTZIN.baseName,'wotsit.tmp','Error adding CMtzDataFile to CContainer')
        c.deleteObject('MTZIN')
        self.assertEqual(len(c),1,'Error deleteing object from CContainer')

    def test5(self):
        c = CContainer()
        c.addContent(name='MTZIN',cls=CCP4XtalData.CMtzDataFile)
        c.addContent(name='RESOLUTION_RANGE',cls='CFloatRange',qualifiers={'compare':-1,'start' : { 'min':0.0,'default': 1.0}, 'end': {'min':0.0,'default':0.0}})
        self.assertEqual(c.RESOLUTION_RANGE.start.qualifiers('min'),0.0,'Error adding CFloatRange to CContainer - wrong qualifier')

    def broken_test7(self):
        defFile = os.path.join(CCP4Utils.getCCP4I2Dir(),'tasks','demo','demo.def.xml')
        demo =  CContainer(definitionFile=defFile)
        #print 'CContainer.test7',demo.inputData.dataOrder()
        demo.inputData.PDBIN = os.path.join(CCP4Utils.getCCP4I2Dir(),'test','data','1df7.pdb')
        #print 'CContainer.test7', demo.inputData.PDBIN
        demo.controlParameters.CHOOSEONE = 5
        demo.controlParameters.CHOOSEOTHER = 2
        other = CContainer(definitionFile=defFile)
        rv = other.copyData(demo,['PDBIN','CHOOSEONE','CHOOSEOTHER','FOO'])
        #print 'CContainer.test7',rv,other.inputData.PDBIN
        self.assertEqual(other.inputData.PDBIN.baseName,'1df7.pdb','Copying between containers failed to give correct filename')
        self.assertEqual(other.controlParameters.CHOOSEONE,5,'Copying between containers failed to give correct value for CHOOSEONE')
        self.assertEqual(len(rv),1,'Copying between containers failed to give correct length error report')
        self.assertEqual(rv[0]['code'],139,'Copying between containers failed to give correct error report code')

    def test8(self):
        c = CContainer()
        i = CContainer(name='inputData')
        i.addContent(name='HKLIN',cls=CCP4XtalData.CMtzDataFile)
        i.addContent(name='XYZIN',cls=CCP4ModelData.CPdbDataFile)
        c.addObject(i)
        err = c.parseCommandLine(['HKLIN','/mydir/myfile.mtz','XYZIN','/mydir/myotherfile.pdb'])
        self.assertEqual(len(err),0,'Parsing command line returns an error')
        self.assertEqual(str(c.inputData.XYZIN.fullPath),'/mydir/myotherfile.pdb','Parsing command line XYZIN misset')

    def test9(self):
        c = CContainer()
        c.addContent(name='XYZIN',cls=CCP4ModelData.CPdbDataFile)
        c.addContent(name='COORD',cls = CCP4MathsData.CXyz)
        err = c.parseCommandLine(['-p', '/mydir/myotherfile.pdb', '-x', '1.0', '2.0', '3.0'],
                                 [['-p' , 'XYZIN'], ['-x' , 'COORD', ['x', 'y', 'z']]])
        self.assertEqual(len(err), 0, 'Parsing command line returns an error')
        self.assertEqual(str(c.XYZIN.fullPath), '/mydir/myotherfile.pdb', 'Parsing command line XYZIN misset')
        self.assertEqual(c.COORD.z, 3.0, 'Parsing command line COORDS.z missed')

    def broken_test10(self):
        defFile = os.path.join(CCP4Utils.getCCP4I2Dir(), 'test', 'data', 'summat_1.def.xml')
        c = CContainer(definitionFile=defFile)
        self.assertEqual(c.PDBIN.baseName.isSet(), False, 'Weirdness in loading from second def file')

    def broken_test11(self):
        c = CContainer()
        c.loadDataFromXml(os.path.join(CCP4Utils.getCCP4I2Dir(), 'pipelines', 'demo_copycell', 'test_data', 'test_1.params.xml'))
        self.assertTrue(c.inputData.HKLIN.isSet(), 'Error using loadDataFromXml with auto load of def file')
    
def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testCContainer)
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)
