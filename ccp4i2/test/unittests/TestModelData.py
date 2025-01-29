def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testPdbData)
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testRange))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testSeqData))
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)


class testPdbData(unittest.TestCase):

    def test1(self):
        p = CPdbData()
        p.loadFile(CPdbDataFile(project='CCP4I2_TOP', relPath='test/data', baseName='1df7.pdb'))
        self.assertEqual(p.composition.chains,['A'], 'Error loading CPDBData - wrong chains')
        self.assertEqual(p.composition.moleculeType, ['PROTEIN', 'MONOMER'], 'Error loading CPDBData - wrong moleculeType')

    def test2(self):
        p = CPdbData()
        resNum, altLoc = p.interpretResidueInput('123.1')
        self.assertEqual(resNum, 123, 'Error in CPdbData.interpretResidueInput - wrong resNum')
        self.assertEqual(altLoc, '1', 'Error in CPdbData.interpretResidueInput - wrong altLoc')

    def test3(self):
        p = CPdbDataFile(project='CCP4I2_TOP', relPath='test/data', baseName='1df7.pdb')
        p.selection = 'A/10-20'
        self.assertEqual(p.isSelectionSet(), True, 'Error - CPdbDataFile.isSelectionSet() wrong')
        p.getSelectedAtomsPdbFile(fileName=os.path.join(CCP4Utils.getTestTmpDir(), '1df7_selection.pdb'))


class testRange(unittest.TestCase):

    def test1(self):
        c = CCP4Container.CContainer()
        c.addContent(name='XYZIN', cls=CPdbDataFile)
        c.addContent(name='DOMAINLIST', cls=CCP4Data.CList, subItem={'class':CResidueRangeList } )
        c.DOMAINLIST.append([])
        c.DOMAINLIST[0].append({'chainId':'A', 'firstRes':'1', 'lastRes':'20'})
        c.DOMAINLIST[0].append({'chainId':'B', 'firstRes':'1', 'lastRes':'40'})
        #print 'testRange.test1',c.DOMAINLIST[0].subItemClass()
        #print 'testRange.test1',c.DOMAINLIST[0],type(c.DOMAINLIST[0])
        #print 'testRange.test1',c.DOMAINLIST[0][1],type(c.DOMAINLIST[0][1])
        self.assertEqual(c.DOMAINLIST[0][1].chainId, 'B', 'Failed to set CResidueRangeList data')


class testSeqData(unittest.TestCase):

    def test1(self):
        seq = CSequence()
        seq.loadFile(CSeqDataFile(project='CCP4I2_TOP', relPath='test/data', baseName='1mzr.fasta'))
        print(seq.__dict__['_value'])
        seq.loadFile(CSeqDataFile(project='CCP4I2_TOP', relPath='test/data', baseName='1mzr.pir'))
        print(seq.__dict__['_value'])

class testCAsuContent(unittest.TestCase):

    def testMolecularWeight(self):
        asuFile = CAsuDataFile(project='CCP4I2_TOP', relPath='demo_data/gamma', baseName='gamma.asu.xml')
        asuFile.loadFile()

    def testCreateCAsuContent(self):
        sequence = CAsuContentSeq({'sequence':'QWERTY', 'nCopies':1, 'polymerType':'PROTEIN', 'name':'DUMMY', 'desccription':'', 'source':''})
        asuContent = CAsuContent()
        asuContent.seqList.addItem(sequence)
        print(asuContent)
        self.assertAlmostEqual(asuContent.molecularWeight(), second=881.93, places=2)
