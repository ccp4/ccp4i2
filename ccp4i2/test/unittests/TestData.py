import sys
import unittest
import xml.etree.ElementTree as ET

from ccp4i2.core import QtCore
from ccp4i2.core.CCP4Data import CList, CInt, CException, CIntRange, CFloat, CDict, CDataFile


def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testCListAppend)
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testCListAssorted))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testQObject))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testDict))
    # suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testTable))  # KJS : This is broken
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)

class testCListAppend(unittest.TestCase):
    def setUp(self):
        self.l = CList([2,3,4],listMinLength=3,listCompare=1, subItem= {'class' : CInt, 'qualifiers' : {'default':0,'min':0}})

    def testAppend0(self):
        self.l.append(5)
        self.assertEqual(self.l.get(),[2,3,4,5])

    def testAppend1(self):
        self.l.append(CInt(5))
        self.assertEqual(self.l.get(),[2,3,4,5])

    def testAppend2(self):
        # Switch off the validity checks to allow an 'uninitialised' CInt to append
        self.l.setQualifiers(min=NotImplemented,listCompare= NotImplemented)
        print('testAppend2 qualifiers',self.l.qualifiers(),self.l.subItemQualifiers())
        self.l.append(CInt(default=0))
        self.assertEqual(self.l.get(),[2,3,4,0])

    def testAppend3(self):
        # Should not append value < min
        # Expect CInt 101
        #self.failUnlessRaises(CException,self.l.append,-1)
        try:
            self.l.append(-1)
        except CException as e:
            self.assertEqual(len(e),1,'Unexpected exception length in CList.append that should fail item quailfier test')
            self.assertEqual(e[0]['code'],101,'Unexpected exception in CList.append that should fail item qualifier test')
        else:
            self.fail('No exception in CList.append that should fail item quailfier test')

    def testAppend4(self):
        # Should not append value < last item in list
        #self.failUnlessRaises(CException,self.l.append,3)
        try:
            self.l.append(3)
        except CException as e:
            self.assertEqual(len(e),1,'Unexpected exception length in CList.append that should fail comparison test')
            self.assertEqual(e[0]['code'],103,'Unexpected exception in CList.append that should fail comparison test')
        else:
            self.fail('No exception in CList.append that should fail comparison test')

class testCListAssorted(unittest.TestCase):
    def setUp(self):
        self.l = CList(subItemClass=CIntRange,listMaxLength=4)

    def testList1(self):
        self.l.set({'start' : 2, 'end':6})
        self.l.append({'start' : 12, 'end':16})
        self.assertEqual(self.l.get(),[{'start' : 2, 'end':6},{'start' : 12, 'end':16}])
        self.assertEqual(self.l[1],{'start' : 12, 'end':16})

    def testList2(self):
        self.l.set({'start' : 2, 'end':6})
        n =  CList({'start' : 12, 'end':16},subItemClass=CIntRange,listMaxLength=4)
        m = self.l + n
        self.assertEqual(m.get(),[{'start' : 2, 'end':6},{'start' : 12, 'end':16}])

    def testList3(self):
        self.l.set({'start' : 2, 'end':6})
        self.l.insert(0,{'start' : 12, 'end':16})
        self.l.reverse()
        self.assertEqual(self.l.get(),[{'start' : 2, 'end':6},{'start' : 12, 'end':16}])

    def testList4(self):
        self.l.set({'start' : 2, 'end':6})
        m =  [  CIntRange(start=4,end=14),CIntRange(start=5,end=15)]
        n = [CIntRange(start=6,end=16),CIntRange(start=7,end=17)]
        self.l.extend(m)
        self.assertEqual(self.l.get(),[{'start' : 2, 'end':6},{'start' : 4, 'end':14},{'start' : 5, 'end':15}])
        try:
            self.l.extend(m)
        except CException as e:
            self.assertEqual(len(e),1,'Unexpected exception length in CList.extend')
            self.assertEqual(e[0]['code'],108,'Unexpected exception in CList.extend')
        else:
            self.fail('No exception in CList.extend that should fail')
        self.l.remove({'start' : 4, 'end':14})
        self.assertEqual(self.l[1].end,15,'Fail after CList.remove')

    def testList5(self):
        testXML = '''<CList>
  <CIntRange>
    <start>2</start>
    <end>6</end>
  </CIntRange>
  <CIntRange>
    <start>4</start>
    <end>14</end>
  </CIntRange>
  <CIntRange>
    <start>5</start>
    <end>15</end>
  </CIntRange>
</CList>
'''
        self.l.set([{'start' : 2, 'end':6},{'start' : 4, 'end':14},{'start' : 5, 'end':15}])
        element = self.l.getEtree()
        ET.indent(element)
        text = ET.tostring(element)
        #print text
        self.assertEqual(text,testXML,'Failed writing XML comparison')
        m = CList(subItemClass=CIntRange,listMaxLength=4)
        m.setEtree(element)
        self.assertEqual(self.l.get(),m.get(),'Failed write/read XML etree')

    def testList6(self):
        self.l.set({'start' : 2, 'end':6})
        j = CList(subItemClass=CIntRange,listMaxLength=4)
        j.append({'start' : 4, 'end':8})
        self.assertEqual(len(j),1,'With two lists - second list wrong length')

class testQObject(unittest.TestCase):
#class testQObject():
    def setUp(self):
        self.bleeped = False
        self.app = QtCore.QCoreApplication(sys.argv)
        self.master = QtCore.QObject(self.app)

    @QtCore.Slot()
    def bleep(self):
        print('BLEEP!!!')
        self.bleeped = True

    def test1(self):
        f = CFloat(parent=self.master)
        f.dataChanged.connect(self.bleep)
        f.set(12.0)
        self.assertTrue(self.bleeped,'dataChanged signal not connected')

class testDict(unittest.TestCase):
    def test1(self):
        d = CDict(subItem={'class' : CIntRange, 'qualifiers' : {'compare' : -1, 'start' : {'min':0}}})
        e = d.set({'foo' : {'start' : 10, 'end' : 5}})
        #print 'd.foo',d.foo
        #print 'error',e
        self.assertEqual(d.foo.start, 10, 'Failed to set CDict')

    def test2(self):
        d = CDict(subItem = {'class' : CIntRange, 'qualifiers' : {'compare' : -1, 'start' : {'min':0}}})
        e = d.set({'foo' : {'start' : -10, 'end' : -20}})
        print('testDict.test2',d,e)
        if len(e) > 0:
            self.assertEqual(e[0]['code'],101,'Setting incorrect dict item does not give correct error code 101')
        else:
            self.fail('Setting incorrect dict item does not give error')

    def test3(self):
        d = CDict(subItem = { 'class' : CIntRange, 'qualifiers' : { 'compare' : -1, 'start' : {'min':0}}})
        e = d.set( { 'foo' : { 'start' : 10, 'end' : 20}} )
        self.assertEqual(e[0]['code'],102,'Setting incorrect dict item does not give correct error code 102')

    def test4(self):
        d = CDict(subItemClass=CDataFile)
        d.PDBIN = { 'project' : 'FOO' , 'baseName' : 'bar.pdb' }
        d.PDBINX = '/foo/bar_x.pdb'
        self.assertEqual(d.dataOrder(),['PDBIN','PDBINX'],'Failed creating Dict using __setattr__')
        self.assertEqual(str(d.PDBINX),'/foo/bar_x.pdb','Failed creating Dict using __setattr__ - 2')
