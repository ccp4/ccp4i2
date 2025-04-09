import os
import unittest

from ccp4i2.core.CCP4Utils import getTestTmpDir, getCCP4I2Dir
from ccp4i2.core.CCP4Annotation import CUserId


def TESTSUITE():
  suite = unittest.defaultTestLoader.loadTestsFromTestCase(testDbConnect)
  #suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testFilePath))
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)

class testDbConnect(unittest.TestCase):

  def test1(self):    
    fileName = os.path.join(getTestTmpDir(),'sqlite.db')
    if os.path.exists(fileName): os.remove(fileName)
    
    userId = CUserId()
    userId.setCurrentUser()
    userName = str(userId)
    print('testDbUser',fileName,userName)
    objUIdb = UIdb(filename =fileName ,hostname='localhost',user=userName)
    cur=objUIdb.establishConnection()
    self.assertFalse(cur is None,'Failed to open Db connection')

    schemaFile = os.path.join(getCCP4I2Dir(),'dbapi','database_schema.sql')
    print('Attempting to read schema',schemaFile)
    objUIdb.readFile(schemaFile)

    #cur.execute("SELECT UserName FROM Users WHERE UserName = '" + username + "'")
    cur.execute("SELECT UserName FROM Users")
    for eachUser in cur.fetchall():
      print('Initial users:'+eachUser)

    cur.createUser('lizp','ysbl')

    cur.execute("SELECT UserName FROM Users")
    for eachUser in cur.fetchall():
      print('Added users:'+eachUser)

    cur.close()
