from __future__ import print_function

"""
     CCP4DbMackup.py: CCP4 GUI Project
     Copyright (C) 2012 STFC

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
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
   Liz Potterton June 2012 - create and benchmark mock database
"""

from dbapi import CCP4DbApi,CCP4ModelData,CCP4XtalData

class BuildDb:

  TASKS = { 
  'SCA' : { 'in' : ['UNM'],
          'out' : ['OBS'] },
  'MR_1' : { 'in' : ['OBS', 'PDB', 'SEQ'],
           'out' : ['PDB','PHZ'] },
  'MR_2' : { 'in' : ['OBS', 'PHZ', 'SEQ'],
           'out' : ['PDB_6','PDB_7','PDB_8','PDB_9','PDB_10'],
           'sub' : [  'MR_2_A' , 'MR_2_B','MR_2_C', 'MR_2_D', 'MR_2_E','MR_2_F' ]  },
  'BUCREF' : { 'in' : ['OBS', 'PHZ', 'PDB'],
             'out' : ['PDB'],
             'sub' : [ 'BUC','REF','BUC','REF','BUC','REF',] },
  'REF_1' : { 'in' : ['OBS', 'PHZ', 'PDB'],
            'out' : ['PHZ','MAP','PDB' ] },
  'BLD_1' :  { 'in' : ['OBS', 'PHZ', 'PDB'],
             'out' : ['PDB'] },
  'MR_2_A' : { 'in' : ['OBS', 'PHZ', 'SEQ'],
              'out' : ['PDB_1','PDB_2','PDB_3','PDB_4','PDB_5'] },
  'MR_2_B' : { 'in' : ['OBS', 'PHZ', 'PDB_1'], 'out' : ['PDB_6']  },
  'MR_2_C' : { 'in' : ['OBS', 'PHZ', 'PDB_2'],  'out' : ['PDB_7'] },
  'MR_2_D' : { 'in' : ['OBS', 'PHZ', 'PDB_3'],  'out' : ['PDB_8'] },
  'MR_2_E' : { 'in' : ['OBS', 'PHZ', 'PDB_4'],  'out' : ['PDB_9'] },
  'MR_2_F' : { 'in' : ['OBS', 'PHZ', 'PDB_5'],  'out' : ['PDB_10'] },
  'BUC' : { 'in' : ['OBS', 'PHZ', 'PDB'], 'out' : ['PDB'] },
  'REF' : { 'in' : ['OBS', 'PHZ', 'PDB'], 'out' : ['PDB','MAP'] }
  }

  PATTERNS =  { 1 : [ 'SCA','SCA','MR_1','MR_2','BUCREF','REF_1','BLD_1','REF_1','BLD_1','REF_1','REF_1','BLD_1','REF_1' ] }

  FILECLASS = { 'PDB' : CCP4ModelData.CPdbDataFile,
                'SEQ' : CCP4ModelData.CSeqDataFile,
                'UNM' : CCP4XtalData.CUnmergedMtzDataFile,
                'OBS' : CCP4XtalData.CObsDataFile,
                'PHZ' : CCP4XtalData.CPhsDataFile,
                'MAP' : CCP4XtalData.CMapCoeffsDataFile }
  FILETYPE = { 'PDB' : 2,
                'SEQ' : 1,
                'UNM' : 5,
                'OBS' : 11,
                'PHZ' : 12,
                'MAP' : 13 }
  
  def __init__(self,fileName=None,mode='sqlite'):
    self.mode=mode
    self.db = CCP4DbApi.CDbApi(mode=mode,fileName=fileName,userName='me')
    self.db.setPreference('testFileExists',False)
    self.projectNumber = 0


  def makeProject(self,pattern=1):
    self.projectNumber = self.projectNumber + 1
    topJobList = []
    projectId = self.db.createProject(projectName='test_project_'+str(self.projectNumber),projectDirectory='/Users/lizp/Desktop/test_projects')
    for taskName in self.PATTERNS.get(pattern):
      if len(taskName)>0:
        if len(topJobList)>0:
          previousJobId = topJobList[-1]
        else:
          previousJobId = None
        jobId = self.makeJob(projectId=projectId,taskName=taskName,previousJobId=previousJobId)
        topJobList.append(jobId)
      else:
        topJobList = []
    return projectId


  def makeJob(self,projectId,taskName=None,parentJobId=None,previousJobId=None,details=None):
    if details is None: details = self.TASKS.get(taskName)
    jobId = self.db.createJob(projectId,taskName,jobTitle='whatever',status=CCP4DbApi.JOB_STATUS_FINISHED,parentJobId=parentJobId)
    self.makeFilesIn(projectId,jobId,details.get('in',[]),previousJobId=previousJobId)
    self.makeFilesOut(projectId,jobId,details.get('out',[]))
    if details.get('sub',None) is not None:
      subJobId = None
      for subName in details['sub']:
        subJobId = self.makeJob(projectId,taskName=subName,parentJobId=jobId,previousJobId=subJobId)
    return jobId

  def makeFilesIn(self,projectId,jobId,fileInList,previousJobId=None):
    for fileIn in fileInList:
      if fileIn.count('_'):
        ftype,num = fileIn.split('_')
      else:
        ftype = fileIn
        num = 0
      if previousJobId is not None:
        fileIdList = self.db.getFileByJobContext(contextJobId=previousJobId,fileType=self.FILETYPE[ftype],projectId=projectId)
      else:
        fileIdList = []
      if len(fileIdList)>0:
        self.db.createFileUse(jobId=jobId,fileId=fileIdList[0],role=CCP4DbApi.FILE_ROLE_IN)
      else:
        fileObject = self.FILECLASS[ftype]()
        fileObject.baseName=fileIn+'.'+fileObject.qualifiers('fileExtensions')[0]
        self.db.createFile(jobId=jobId,projectId=projectId,fileObject=fileObject,sourceFileName='/foo/bar/'+fileIn)

  def makeFilesOut(self,projectId,jobId,fileOutList):
    for fileOut in fileOutList:
      if fileOut.count('_'):
        ftype,num = fileOut.split('_')
      else:
        ftype = fileOut
        num = 0
      fileObject = self.FILECLASS[ftype]()
      fileObject.baseName=fileOut+'.'+fileObject.qualifiers('fileExtensions')[0]
      
      self.db.createFile(jobId=jobId,projectId=projectId,fileObject=fileObject)
      

def makeDb():
  import os
  try:
    os.remove('/Users/lizp/Desktop/mockDb.sqlite')
    os.remove('/Users/lizp/Desktop/mockDb.xml')
  except:
    pass
  mDb = BuildDb('/Users/lizp/Desktop/mockDb.sqlite')
  projectId = mDb.makeProject()
  err = mDb.db.exportProjectXml(projectId=projectId,fileName='/Users/lizp/Desktop/mockDb.xml')
  print(err.report())
  
