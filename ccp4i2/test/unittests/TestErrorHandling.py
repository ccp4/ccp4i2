import unittest

from ccp4i2.core.CCP4Data import CData
from ccp4i2.core.CCP4ErrorHandling import CException


def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testError)
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)

class testError(unittest.TestCase):

    def test1(self):
        e = CException(CData, 1, 'foo')
        tree = e.getEtree()
        f = CException()
        f.setEtree(tree)
        self.assertEqual(f[0]['code'],1,'Error save/restore CException to etree = wrong code')
