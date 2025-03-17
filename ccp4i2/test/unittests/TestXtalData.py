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
        self.app = QTAPPLICATION()
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
        PROCESSMANAGER().setWaitForFinished(10000)
        self.app = QTAPPLICATION()
        self.mummy = QtCore.QObject(self.app)

    def tearDown(self):
        PROCESSMANAGER().setWaitForFinished(-1)

    def test_1(self):
        self.mtz = CMtzDataFile( os.path.normpath(os.path.join(self.testDataDir,'gere_nat.mtz'),parent=self.mummy))
        print("test_1 getNColumns",self.mtz.getFileContent().getNColumns())
        self.assertEqual( 19, self.mtz.getFileContent().getNColumns(),'CMtzDataFile loaded MTZ reports wrong number of columns')

    def test_2(self):
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
        self.app = QTAPPLICATION()
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
        self.app = QTAPPLICATION()
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
