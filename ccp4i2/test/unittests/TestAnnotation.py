import time
import unittest

from ccp4i2.core.CCP4Annotation import CAnnotation, CAnnotationList, CTime
from ccp4i2.core.CCP4Utils import getUserId


def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testAnnotation)
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)

class testAnnotation(unittest.TestCase):

    def test1(self):
        a = CAnnotation('Try this text')
        self.assertEqual(a.text,'Try this text','Failed setting CAnnotation text')
        self.assertEqual(a.author,getUserId(),'Failed setting CAnnotation author')
        if a.time < int(time.time())-10 or a.time > int(time.time()):
            self.fail('Failed setting CAnnotation time')

    def test2(self):
        a = CAnnotation(text='Try this text')
        t = CTime()
        t.setCurrentTime()
        d = a.time-t
        self.assertEqual(d.__class__, CTime, 'Time difference not a CTime object')
        self.assertEqual(0 <= d < 10, True, 'Time difference notin expected range')

    def test3(self):
        annoList = CAnnotationList()
        annoList.append('Test this string')
        annoList.append(CAnnotation('Test another string'))
        self.assertEqual(len(annoList),2,'Annotation list wrong length')
        self.assertEqual(annoList[1].text,'Test another string','Annotation wrong text')
