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
