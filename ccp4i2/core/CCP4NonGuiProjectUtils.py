from __future__ import print_function


"""
     CCP4NonGuiProjectUtils.py: CCP4 GUI Project
     Copyright (C) 2015 University of Newcastle

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the sstatusW
     license to address the requirements of UK law.
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php

     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
"""

"""
     Martin Noble Jul 2015
     
     Create a class to allow for programatic (non gui) import of a .ccp4_project.zip 
     There is a distressing degree of code in common with CCP4ProjectManagerGui...would
     be better if one were calling the other more.
     
     Initially offers one call
"""
from ccp4i2.core.CCP4Modules import PROJECTSMANAGER
from ccp4i2.core.CCP4ErrorHandling import CException, SEVERITY_WARNING
from PySide2 import QtCore
from ccp4i2.core import CCP4Utils

DIAGNOSTIC=True

class CCP4NonGuiProjectUtils(QtCore.QObject):
  def __init__(self, compressedArchive=None):
    super(CCP4NonGuiProjectUtils,self).__init__()
    if compressedArchive is not None:
        self.importCompressedArchive(compressedArchive)
    
  def importCompressedArchive(self, compressedArchive=None):
    import os
    if not os.path.isfile(compressedArchive): return

    try:
      xmlFile = PROJECTSMANAGER().extractDatabaseXml(compressedArchive)
    except CException as e:
      print('Import project - Failed extracting database XML file from compressed file')
      return
    except:
      print('Import project - Error extracting database xml file from ')
      return

    try:
      from dbapi import CCP4DbApi
      self.dbImport = CCP4DbApi.CDbXml(db=PROJECTSMANAGER().db(),xmlFile=xmlFile)
      #self.dbImport.setDiagnostic(True)
      importProjectInfo = self.dbImport.loadProjectInfo()  
    except:
      print('Import project - Error attempting to read database file in\n'+str(compressedArchive))
      return
    try:
      projectInfo =  PROJECTSMANAGER().db().getProjectInfo(projectId=self.dbImport.projectId)
    except:
      projectInfo = None
   
    #print 'handleImportProject1 importProjectInfo',importProjectInfo

    if projectInfo is None:
      print('Is a new project',self.dbImport.projectId)
      import os
      dirNameRoot = os.path.normpath(os.path.join(CCP4Utils.getProjectDirectory(),self.dbImport.projectName))
      dirName = dirNameRoot
      incrementCounter = 0
      while os.path.exists(dirName):
          dirName = dirNameRoot + str(incrementCounter)
          incrementCounter += 1
      print(dirName)
      self.importProject(compressedArchive, dirName)
    else:
      print('Is an existing project', self.dbImport.projectId,' In', projectInfo['projectdirectory'])
      self.importProject(compressedArchive, projectInfo['projectdirectory'], existingProject=self.dbImport.projectId)

  def importProject(self,compressedArchive,dirName,existingProject=None):
    import os
    if DIAGNOSTIC: print('CProjectManagerDialog.importProject',compressedArchive,dirName,existingProject)
    # Load the database.xml into temporary tables in db
    self.dbImport.projectDirectory = dirName
    if existingProject is None:
      ret = self.dbImport.createProject()
      if ret.maxSeverity()>SEVERITY_WARNING:
        print(ret.report())
        print('Error creating project in database')
        return
      print('Succeeded in creating project')
    
    self.dbImport.createTempTables()
    self.dbImport.loadTempTable()
    # If loading jobs to an existing project flag up jobs in temp tables that
    # are already in db
    if existingProject is not None:
      self.dbImport.setExclInTempTables()
    # Flag imported files to be imported (there is no checking yet that they exist)
    self.dbImport.setExclImportedFiles()

    if DIAGNOSTIC: print('CProjectManagerDialog.importProject setting Temp Tables',self.dbImport.errReport.report())
    
    if self.dbImport.errReport.maxSeverity()>SEVERITY_WARNING:
      if DIAGNOSTIC: print('Error report from the import process..')
      if DIAGNOSTIC: print(self.dbImport.errReport.report())
      print('Error loading data from project export file')
    
    # Make project directory if necessary
    if not os.path.exists(dirName):
      try:
        os.mkdir(dirName)
      except:
        print('Import project - Failed to create directory:', dirName)
        return

    print('Unpacking project files to '+dirName)
    from qtcore import CCP4Export
    # Unpack project files from the tar file (possibly in separate thread) 
    # Pass import thread dbImport to enable query database and flagging loaded jobs/files
    if DIAGNOSTIC: print('CProjectManagerDialog.importProject creating import thread')
    self.importThread = CCP4Export.ImportProjectThread(self,projectDir=dirName,compressedFile=compressedArchive,
                                                       dbImport=self.dbImport,diagnostic=DIAGNOSTIC)
    errReport = self.importThread.run()
    if DIAGNOSTIC: print('CProjectManagerDialog.importProject import thread running',errReport.report())
    
    self.dbImport.cleanupTempTables()
    self.dbImport.listTempJobs('TempJobs after cleanup')
    self.dbImport.listTempFiles('TempFiles after cleanup')
    self.dbImport.listTempFileUses('TempFileUses after cleanup')
    stats = self.dbImport.importStats()
    if DIAGNOSTIC:
      for key,value in list(stats.items()):
        if key == 'failedFiles':
          if len(value)>0: print('Failed to import files..(probably already present)')
          for item in value:
              print('Job_'+str(item[4]),item[2])
        else:
          print('CProjectManagerDialog.importProject stats', key,value)
    if errReport.maxSeverity()>SEVERITY_WARNING:
      self.dbImport.removeTempTables()
      text = 'ERRORS UNPACKING DATA FILES\n'
      for err in errReport: text = text + err['details'] + '\n'
      print('Import failed',text)
      print('Error unpacking project files to '+dirName)
      return

    print('Loading project data to database')
    self.dbImport.importTempTables()
    self.dbImport.removeTempTables()
    print('Finished importing project')

    self.dbImport.db.projectReset.emit({'projectId':self.dbImport.projectId})
    
    if stats['jobsTotal']>stats['jobsImported'] or stats['filesTotal']>stats['filesImported']:
      text = 'Some of the jobs/files are already in the database in project '+str(self.dbImport.projectName)+'\n' + \
      'Imported '+str(stats['jobsImported'])+' new jobs from '+str(stats['jobsTotal'])+' in file \n' + \
     'and '+str(stats['filesImported'])+' new data files from '+str(stats['filesTotal'])+' in file\n'
    else:
      text = 'Successfully imported '+str(stats['jobsImported'])+' jobs and '+str(stats['filesImported'])+' data files'

    if stats.get('incrJobNumber',0) > 0:
        text = text +'\nImporting jobs '+str(stats['importMin'])+' to '+str(stats['importMax'])+' have been renumbered\n'+str(int(stats['importMin'])+int(stats['incrJobNumber']))+' to '+str(int(stats['importMax'])+int(stats['incrJobNumber'])) +' to avoid clash with existing jobs'
    if len(text)>0:  print('Import complete',text)
        
  def createImportProgress(self):
    pass

  def handleImportProgress(self,ret):
    #print 'handleImportProgress',ret
    jobNo,done = ret

  def doneImportProgress(self):
    pass
