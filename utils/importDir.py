from __future__ import print_function

"""
An even more simplified version of the stuff in CCP4TestDb.py
This reads a directory's DATABASE.db.xml and attempts to build a db entry for it.

TODO:
Does not seem quite perfect yet. Imports 198 out of 259 for one of my preojects. Certainly ignores all "Pending" jobs. Hmm. Possibly failed too.
So probably need os.walk through CCP4_JOBS instead of looking at files in XML file?
"""
import sys
import os
if __name__ == "__main__":
    sys.path.append(os.path.join(os.path.dirname(__file__),".."))
import shutil

from PySide6 import QtCore
from lxml import etree

from dbapi import CCP4DbApi,CCP4DbUtils
from core import CCP4File
from core.CCP4ErrorHandling import *
from core.CCP4Modules import PROJECTSMANAGER,TASKMANAGER,QTAPPLICATION,JOBCONTROLLER
from core.CCP4QtObject import CObject
from qtcore import CCP4Export   

def parse_from_unicode(unicode_str):
    utf8_parser = etree.XMLParser(encoding='utf-8')
    s = unicode_str.encode('utf-8')
    return etree.fromstring(s, parser=utf8_parser)

def importFilesFromDirXML(dbFileName,dirName,importProjectComments=False):
    """
    dbFileName: the database file
    dirName: the directory we are scanning.
    """
    pm = PROJECTSMANAGER()
    db = CCP4DbApi.CDbApi(parent=pm, fileName=dbFileName, mode='sqlite', createDb=True)
    pm.setDatabase(db)
    pm.startCheckForFinishedJobs()
    
    #xmlFile = PROJECTSMANAGER().extractDatabaseXml(compressedFile)
    xmlFile = os.path.join(dirName,"DATABASE.db.xml")
    dbImport = CCP4DbApi.CDbXml(db=db,xmlFile=xmlFile,diagnostic=False)
    importProjectInfo = dbImport.loadProjectInfo()
    print('importProject importProjectInfo',importProjectInfo)
    
    dbImport.projectDirectory = dirName
    dbImport.createProject()
    dbImport.createTempTables()
    dbImport.loadTempTable()
    dbImport.setExclImportedFiles()
    
    with open(xmlFile) as f:
        t = f.read()
    
        if sys.version_info < (3,0):
            xmlparser = etree.XMLParser()
            tree = etree.fromstring(t, xmlparser)
        else:
            tree = parse_from_unicode(t)
        
        doneJobs = []
        if len(tree.xpath("ccp4i2_body"))>0:
            #Import jobs with input/output files and their files.
            if len(tree.xpath("ccp4i2_body")[0].xpath("fileTable"))>0:
                ft = tree.xpath("ccp4i2_body")[0].xpath("fileTable")[0]
                for f in ft:
                    if len(tree.xpath("//job[@jobid='"+f.attrib["jobid"]+"']"))>0:
                        job = tree.xpath("//job[@jobid='"+f.attrib["jobid"]+"']")[0]
                        jobNo = job.attrib["jobnumber"]
                        if not jobNo in doneJobs:
                            ifImport,newJobNumber = dbImport.importThisFile( jobNumber = jobNo, fileName = None )
                            doneJobs.append(jobNo)
                        ifImport,newJobNumber = dbImport.importThisFile( jobNumber = jobNo, fileName = f.attrib["filename"] )
 
            #Import jobs with no input/output files.
            if len(tree.xpath("ccp4i2_body")[0].xpath("jobTable"))>0:
                jt = tree.xpath("ccp4i2_body")[0].xpath("jobTable")[0]
                for j in jt:
                    if len(tree.xpath("//job[@jobid='"+j.attrib["jobid"]+"']"))>0:
                        job = tree.xpath("//job[@jobid='"+j.attrib["jobid"]+"']")[0]
                        jobNo = job.attrib["jobnumber"]
                        if not jobNo in doneJobs:
                            ifImport,newJobNumber = dbImport.importThisFile( jobNumber = jobNo, fileName = None )
                            doneJobs.append(jobNo)

    dbImport.db.commit()
    pm.makeProjectDirectory(directory=dirName)
        
    dbImport.cleanupTempTables()
    stats = dbImport.importStats()
    for key,value in stats.items():
        print('importProject stats', key,value)
    
    dbImport.importTempTables()
    if importProjectComments: dbImport.importProjectCommentsTempTables()
    dbImport.removeTempTables()
    dbImport.db.projectReset.emit({'projectId':dbImport.projectId})
    dbImport = None

if __name__ == "__main__":
    fileName = sys.argv[1]
    dirName = sys.argv[2]

    importFilesFromDirXML(fileName,dirName)

