from __future__ import print_function

"""
     CCP4databaseDef.py: CCP4 GUI Project
     Copyright (C) 2010 University of York

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
   Liz Potterton May 2010: created
"""

##@package CCP4DatabaseDef  Python representation of old ccp4i database.def files

"""
Example use:
  t = CDatabaseDef('..../database.def')
  for iJob in range(1,t.nJobs):
    if t.jobs.has_key(iJob):
      print 'taskname',t.jobs[iJob]['TASKNAME']

"""

class CDatabaseDef:
  def __init__(self,filename):
    self.headerArray = {}
    self.dataArray = {}
    self.project_name =''
    self.project_dir = ''
    self.jobs = {}
    self.headerArray,typeArray,self.dataArray = self.read(filename)
    #print 'dataArray',self.dataArray

    # Get project name and dir from header
    if 'PROJECT' in self.headerArray:
      self.project_name = self.headerArray['PROJECT'][0]
      if len(self.headerArray['PROJECT'])>1:
        self.project_dir = self.headerArray['PROJECT'][1]
    print('Read project',self.project_name,self.project_dir)

    self.nJobs = int(self.dataArray.get('NJOBS',0))

    self.jobs = self.extractJobs(self.dataArray)

  def read (self,filename='',headerOnly=0):
      '''
      Read the def file into dicts for header, types and data
      A database.def should not have any type information
      This is left in the code in case want to reuse elsewhere
      '''
      headerArray = {}
      typeArray = {}
      dataArray = {}
      try:
        f = open(filename)
      except:
        print('ERROR opening',filename)
        return headerArray,typeArray,dataArray
      if not f:
        print('ERROR opening',filename)
        return headerArray,typeArray,dataArray
      try:
        content = f.readlines()
        f.close()
      except:
        print('ERROR reading ',filename)
        return headerArray,typeArray,dataArray

      
      for line in content:
        words = self.splitDefLine(line)
        #print 'CCP4Data.parseDefFile',words
        if len(words)>=2:
          if words[0] == '#CCP4I' and len(words)>=3:
            #print 'header',words
            headerArray[words[1]] = words[2:]
          elif words[0][0] == '#':
            pass
          elif not headerOnly:
           if len(words)>2:
             if  words[1][0]=='_':
               typeArray[words[0]] = words[1][1:]
             else:
               typeArray[words[0]] = words[1]
             idata = 2
           else:
             idata = 1
           if words[idata].strip('"')=='':
             dataArray[words[0]] = ''
           else:
             dataArray[words[0]]=words[idata]
     
      return  headerArray,typeArray,dataArray

  def splitDefLine(self,line=''):
    import re
    m = re.search(r'(.*)\"(.*)\"(.*)',line)
    if not m:
      return line.split()
    else:
      a,b,c = m.groups()
      a = a.strip()
      c = c.strip()
      rv = []
      if a:rv.extend(a.split())
      rv.append(b)
      if c: rv.append(c.split())
      return rv
    
  def splitDefList(self,line=''):
    '''
    Parse the INPUT_FILE and OUTPUT_FILE from database.def
    This is a space-separated list
    If the file name includes spaces then the name is enclosed in curly brace
    '''
    import re
    rv = []
    while len(line)>0:
      m = re.search(r'(.*?)\{(.*?)\}(.*)',line)
      if not m:
        rv.extend(line.split())
        line = ''       
      else:
        a,b,c = m.groups()
        #print 'splitDefList',a,'*',b,'*',c
        a = a.strip()
        c = c.strip()
        if a:rv.extend(a.split())
        rv.append(b)
        line = c.strip()
        #print 'new line',line
    return rv


  def extractJobs(self,dataArray):

    jobs = {}

    for key,value in list(dataArray.items()):
      try:
        if key.count(',') == 1:
          dataType,jobId = key.split(',')
          try:
            jobId = int(jobId)
          except:
            jobId = -1
          if jobId>0:
            # Initialise the data for one job
            if jobId not in jobs:
              jobs[jobId] = { 'STATUS' : '',
                              'DATE' : 0,
                              'LOGFILE' : '',                             
                              'TASKNAME' : '', 
                              'TITLE' : '',
                              'INPUT_FILES' : [],
                              'INPUT_FILES_DIR' : [],
                              'OUTPUT_FILES' : [],
                              'OUTPUT_FILES_DIR' : []   }
            # Check that the dataType is one of the recognised properties for a job
            if dataType in jobs[jobId]:
              if ['INPUT_FILES','OUTPUT_FILES'].count(dataType):
                jobs[jobId][dataType] = self.splitDefList(value)
              elif ['INPUT_FILES_DIR','OUTPUT_FILES_DIR'].count(dataType):
                jobs[jobId][dataType] = value.split(' ')
              else:
                jobs[jobId][dataType] = value
      except:
         print('ERROR interpreting ',key)

    return jobs 


