"""
This program constructs a CCP4I2 database from a set of ZIP files in a folder.
It can be used in two ways:
  1) Construct a new DB "from scratch" for just the newly imported projects
  2) Append to an existing database (-a option). This will only work if the project name does not already exist in the database.

  Usage:
        python3 ImportAllProjects.py input_folder_containing_zip_files output_folder_where_project_files_will_go [-a|--append] [-d|--dbFile=outputdbname]
"""

import argparse
import glob
import os
import shutil
import sqlite3
import sys
import tempfile

from ..core import CCP4Utils
from ..core.CCP4ErrorHandling import Severity
from ..core.CCP4Modules import PROJECTSMANAGER
from .QApp import QTAPPLICATION


def ImportZipFile(compressedFile,destDirName):
    from ..dbapi import CCP4DbApi
    from ..qtcore import CCP4Export
    try:
      xmlFile = PROJECTSMANAGER().extractDatabaseXml(compressedFile)
    except CException as e:
      print('Import project: Failed extracting database XML file from compressed file',compressedFile)
      return

    try:
      dbImport = CCP4DbApi.CDbXml(db=PROJECTSMANAGER().db(),xmlFile=xmlFile)
      importProjectInfo = dbImport.loadProjectInfo()
    except:
      print('Import project: Error attempting to read database file in\n'+str(compressedFile))
      return

    try:
      projectInfo =  PROJECTSMANAGER().db().getProjectInfo(projectId=dbImport.projectId)
    except:
      projectInfo = None

    print('handleImportProject1 importProjectInfo',importProjectInfo)
    print('handleImportProject1 projectInfo',projectInfo)
    print('handleImportProject1 dbImport.projectName',dbImport.projectName)

    dirName = os.path.join(destDirName,dbImport.projectName)

    if projectInfo is None:
      ret = dbImport.createProject(projectDirectory=dirName)
      if ret.maxSeverity()>Severity.WARNING:
          print('Error creating project',dbImport.projectName,'in database')
          return
      dbImport.createTempTables()
      dbImport.loadTempTable()
      dbImport.setExclImportedFiles()
      if dbImport.errReport.maxSeverity()>Severity.WARNING:
          print('Error report from the import process..')
          print(dbImport.errReport.report())
      if not os.path.exists(dirName):
        try:
          os.mkdir(dirName)
        except:
          print('Import project - Failed to create directory:',dirName)
          return

      print('ImportAllProjects creating import thread')
      importThread = CCP4Export.ImportProjectThread(None,projectDir=dirName,compressedFile=compressedFile,
                                                         dbImport=dbImport,diagnostic=True)
      errReport = importThread.run()

      print('ImportAllProjects.importProject import thread running',errReport.report())

      dbImport.cleanupTempTables()
      stats = dbImport.importStats()

      for key,value in list(stats.items()):
          if key == 'failedFiles':
            if len(value)>0: print('Failed to import files..(probably already present)')
            for item in value:
                print('Job_'+str(item[4]),item[2])
          else:
            print('ImportAllProjects stats', key,value)

      if errReport.maxSeverity()>Severity.WARNING:
        dbImport.removeTempTables()
        text = 'ERRORS UNPACKING DATA FILES\n'
        for err in errReport: text = text + err['details'] + '\n'
        print('Error unpacking project files to '+dirName)
        print(text)
        return

      print('Loading project data to database')
      dbImport.importTempTables()
      dbImport.removeTempTables()
      print('Finished importing project')

      if stats['jobsTotal']>stats['jobsImported'] or stats['filesTotal']>stats['filesImported']:
        text = 'Some of the jobs/files are already in the database in project '+str(dbImport.projectName)+'\n' + \
        'Imported '+str(stats['jobsImported'])+' new jobs from '+str(stats['jobsTotal'])+' in file \n' + \
       'and '+str(stats['filesImported'])+' new data files from '+str(stats['filesTotal'])+' in file\n'
      else:
        text = 'Successfully imported '+str(stats['jobsImported'])+' jobs and '+str(stats['filesImported'])+' data files'

      if stats.get('incrJobNumber',0) > 0:
          text = text +'\nImporting jobs '+str(stats['importMin'])+' to '+str(stats['importMax'])+' have been renumbered\n'+str(int(stats['importMin'])+int(stats['incrJobNumber']))+' to '+str(int(stats['importMax'])+int(stats['incrJobNumber'])) +' to avoid clash with existing jobs'

      print('Import complete',text)

      #PROJECTSMANAGER().backupDBXML()

    else:
      print("Project '"+dbImport.projectName+"' already exists")

def ImportAll(zipDir,dbFile,destDirName,appendDB):
    from . import startup
    try:
        os.mkdir(destDirName)
    except FileExistsError as e:
        pass
    except OSError as e:
        print(error)
        return

    QTAPPLICATION(graphical=False)

    zipFiles = glob.glob(os.path.join(zipDir,"*.zip"))
    
    if not appendDB:
        tfile = tempfile.NamedTemporaryFile(delete=False)
        tfn = tfile.name+".sqlite"
        tfile2 = tempfile.NamedTemporaryFile(delete=False)
        tfn2 = tfile2.name
        tfile.close()
        tfile2.close()
        pm = startup.startProjectsManager(dbFileName=tfn)
    else:
        pm = startup.startProjectsManager(dbFileName=dbFile)

    for z in zipFiles:
        try:
            ImportZipFile(z,destDirName)
        except:
            print("----------------------------------------")
            print("----------------------------------------")
            print("Importing",z,"failed")
            print("----------------------------------------")
            print("----------------------------------------")

    pm.db().close()

    if not appendDB:
        #Now we use sqlite api to make fresh, "clean" hopefully copy of new db...
        con = sqlite3.connect(tfn)
        sql = "".join([s+"\n" for s in con.iterdump()])

        f2trunc = open(tfn2,"w")
        f2trunc.close()
        print("Writing to",tfn2)

        conbak = sqlite3.connect(tfn2)
        cur = conbak.cursor()

        try:
                cur.executescript(sql)
        except:
                print("Fail",com)
                conbak.close()
                raise

        conbak.commit()
        #... and then copy to desited location.
        shutil.copy(tfn2,dbFile)

        print()
        print("############################################################")
        print("Saved database",dbFile)
        print()
        print("If you wish to use this as a 'real' CCP4I2 database, then you")
        print("Need to copy it to (e.g. Mac and Linux) $HOME/.CCP4I2/db :")
        print()
        print("cp "+dbFile+" $HOME/.CCP4I2/db/database.sqlite")
        print()
        print("And remove any ealier -wal and -shm files")
        print()
        print("rm $HOME/.CCP4I2/db/database.sqlite-wal")
        print("rm $HOME/.CCP4I2/db/database.sqlite-shm")
        print()
        print("You *must* close any instance of CCP4I2 before doing the above.")
        print("############################################################")
        conbak.close()
    else:
        print()
        print("############################################################")
        print("Added projects to database",dbFile)
        print()
        print("If you wish to use this as a 'real' CCP4I2 database, then you")
        print("Need to copy it *and any corresponding -wal and -shm files")
        print("to (e.g. Mac and Linux) $HOME/.CCP4I2/db :")
        print()
        print("cp "+dbFile+" $HOME/.CCP4I2/db/database.sqlite")
        print("cp "+dbFile+"-wal $HOME/.CCP4I2/db/database.sqlite-wal")
        print("cp "+dbFile+"-shm $HOME/.CCP4I2/db/database.sqlite-shm")
        print()
        print("You *must* close any instance of CCP4I2 before doing the above.")
        print("############################################################")

    return

if __name__ == "__main__":

    parser = argparse.ArgumentParser( prog='ImportAllProjects')
    parser.add_argument('inputfolder',help="Folder containing input zip files")
    parser.add_argument('outputfolder',help="Output folder where project files will go. e.g. $HOME/CCP4I2_PROJECTS")
    parser.add_argument('-a', '--append', action='store_true',help="This flag will attempt to add to existing database rather than creating new")
    parser.add_argument('-d', '--dbFile',help="Output database file name, defaults to 'database.sqlite'")
    args = parser.parse_args()
    print(args.inputfolder)
    print(args.outputfolder)
    print(args.append)
    print(args.dbFile)

    destDirName = CCP4Utils.getProjectDirectory()

    if len(sys.argv)>2:
        destDirName = args.outputfolder

    if args.dbFile:
        dbName = args.dbFile
    else:
        dbName = "database.sqlite"

    ImportAll(args.inputfolder,dbName,destDirName,args.append)
