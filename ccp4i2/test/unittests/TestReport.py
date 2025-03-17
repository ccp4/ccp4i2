def TESTSUITE():
  suite = unittest.defaultTestLoader.loadTestsFromTestCase(testReport)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)

class testReport(unittest.TestCase):

  def test1(self):
    testFile = CCP4File.CDataFile(project='CCP4I2_TOP',relPath='test/data',baseName='test_report.html')
    r = CReport()
    r.loadFromXmlFile(str(testFile))
    x = r.getDataObject(id='final_results')
    self.assertEqual(x.cell.a,123.4,'Failed to load container and cell data')
    x = r.getDataObject(id='data_table_1')
    self.assertEqual(x.nColumns,9,'Failed to load CReportTableModel number of columns')
