"""
An even more simplified version of the stuff in CCP4TestDb.py
This reads a directory's DATABASE.db.xml and attempts to build a db entry for it.

TODO:
Does not seem quite perfect yet. Imports 198 out of 259 for one of my preojects. Certainly ignores all "Pending" jobs. Hmm. Possibly failed too.
So probably need os.walk through CCP4_JOBS instead of looking at files in XML file?
"""

import os
import sys
import xml.etree.ElementTree as ET

from ..core.CCP4Modules import PROJECTSMANAGER
from ..dbapi import CCP4DbApi


def parse_from_unicode(unicode_str):
    return ET.fromstring(unicode_str)


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
