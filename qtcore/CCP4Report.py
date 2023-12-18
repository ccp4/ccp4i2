from __future__ import print_function

"""
     CCP4Report.py: CCP4 GUI Project
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

##@package CCP4Report (QtCore) Data types for tables,graphs etc found in reports

from PySide2 import QtCore
import os
from lxml import etree
from xml.etree import ElementTree as ET
from core.CCP4ErrorHandling import *
from core.CCP4DataManager import DATAMANAGER
from core.CCP4Data import CData
from report.CCP4ReportParser import CCP4NS

NSMAP = { 'ccp4' : CCP4NS }


class CReportTable(QtCore.QAbstractTableModel):

  def __init__(self,parent=None,name=None,qualifiers={},**kw):
    QtCore.QAbstractTableModel.__init__(self,parent)
    if name is not None: self.setObjectName(name)
    self.nColumns = 0
    self.nRows = 0
    self.title = None
    self.columnLabels = []
    self.columns = []
    self.columnTypes = []

      
  def setEtree(self,element=None,checkValidity=True):
    self.setObjectName(str(element.get('id')))
    e = element.find('Ncolumns')
    if e is not None: 
      self.nColumns = int(e.text)
      # 'CReportTableModel.setEtree Ncolumns',self.nColumns
    e = element.find('title')
    if e is not None: self.title = e.text
    e = element.find('labels')
    if e is not None:
      labels_text = e.text
      if labels_text: self.columnLabels = labels_text.split()
      if len(self.columnLabels) != self.nColumns:
         print('CReportTableModel.setEtree number of labels',len(self.columnLabels),'does not equal Ncolumns',self.nColumns)
    while len(self.columnLabels)<self.nColumns:
        self.columnLabels.append(None)
    #for ic in range(0,len(self.columnLabels)):
    #  self.setHeaderData(ic, QtCore.Qt.Horizontal,self.columnLabels[ic] )    
        
    e = element.find('data')
    if e is not None:
      try:
        eleList = e.text.split()
      except:
        eleList = []

    # Read first row and determine data type (int or float)
    # Add approariately typed array to self.columns list
    if len(eleList)> self.nColumns:
      import array
      self.columns = []
      self.columnTypes = []
      for ic in range(0,self.nColumns):
        try:
          ivalue = int(eleList[ic])
          self.columns.append(array.array('i'))
          self.columnTypes.append('int')
        except:
           try:
             fvalue = float(eleList[ic])
             self.columns.append(array.array('f'))
             self.columnTypes.append('float')
           except:             
             self.columnTypes.append('text')
             self.columns.append(None)
      #print 'CReportTableModel.setEtree columnTypes',self.columnTypes

      self.nRows = len(eleList)/self.nColumns

      ie = -1
      for ir in range(0,self.nRows):
        for ic in range(0,self.nColumns):
          ie = ie + 1
          try:
            if self.columnTypes[ic] == 'int':
              self.columns[ic].append(int(eleList[ie]))
            elif self.columnTypes[ic] == 'float':
              self.columns[ic].append(float(eleList[ie]))
          
          except:
            if self.columnTypes[ic] == 'int':
              self.columns[ic].append(0)
            elif self.columnTypes[ic] == 'float':        
              self.columns[ic].append(0.0)
    
    return 0
      
  def rowCount(self,idx):
    return self.nRows

  def columnCount(self,idx):
    return self.nColumns

  def data(self,idx=None,role=QtCore.Qt.DisplayRole):
    #print 'CReportTableModel.data index',idx.column(),idx.row(),self.columns[idx.column()]
    if self.columnTypes[idx.column()] != 'text':
      if role == QtCore.Qt.DisplayRole:
        return self.columns[idx.column()][idx.row()]

#FIXME PYQT - or maybe None? This used to return QVariant.
    return None

  def columnDataType(self,column=-1):
    if column>=0 and column<len(self.columnTypes):
      return self.columnTypes(ic)
      
  def headerData(self,section,orientation=QtCore.Qt.Horizontal,role=QtCore.Qt.DisplayRole):
    if role==QtCore.Qt.DisplayRole:
      if orientation == QtCore.Qt.Horizontal:
        if section>=0 and section<len(self.columnLabels):
          #print 'headerData',section,self.columnLabels[section]
          v = self.columnLabels[section]
          return v
#FIXME PYQT - or maybe None? This used to return QVariant.
      return None

  def qualifiers(self):
    return {}


class CPlotDirectives(CData):

  pass  
    
    
class CReport(QtCore.QObject):

  ERROR_CODES =  { 101 : {'description' : 'Report XML file does not found' },
                   102 : { 'description' : 'Error parsing report file' },
                   103 : { 'description' : 'Error converting input to Python string' },
                   104 : {'description' :'Error converting to text to save report' },
                   105 : { 'description' :'Error saving report to file' },
                   106 : { 'description' :'Can not find data class corresponding to type in file' }
                   }
  def __init__(self,parent=None):
      QtCore.QObject.__init__(self,parent)
      self.ccp4_data = {}
      self.linkIds = {}
      self.dataArray = {}
      self._jobInfo = {}
      self.baseRef= None
      self.resetBaseHref = None
      self.topFolds = []
      
  def objectPath(self):
    return str(self.objectName())

  def loadFromXmlFile(self,filename):
    self.filename = os.path.normpath(filename)
    #print 'CReport.loadFromXmlFile',self.filename
    if not os.path.exists(self.filename):
      raise CException(self.__class__,101,self.filename)
    try:
      from core import CCP4Utils
      root = CCP4Utils.openFileToEtree(self.filename)
      #root = tree.getroot()
      #print 'CReport.loadFromXmlFile',root
    except ET.ParseError as e:
      #print 'CReport.loadFromXmlFile',e.filename,e.lineno,e.msg,e.offset,e.print_file_and_line,e.text
      raise CException(self.__class__,102,e.msg+' in file '+filename)
    except:
      raise CException(self.__class__,102,filename)

    self.loadFromEtree(root)
    #print 'loadFromXmlFile loaded tree' , self.toString(root)
    
  def loadFromEtree(self,root):
    #try:
      self.extractBaseHref(root)
      self.extractCCP4Data(root)
      self.extractLinkIds(root)
      self.extractTopFolds(root)
    #except:
    #  raise CException('CReport.loadFromEtree','Error extracting ccp4_data elements')
    

  def containsCCP4Data(self):
    if self.ccp4_data:
      return True
    else:
      return False

  def toString(self,root):
    text = ''
    try:
      text = ET.tostring(root,pretty_print=True)
    except:
      raise CException(self.__class__,103)
    return text
  
  def saveToXmlFile(self,root,filename):
    try:
      text = ET.tostring(root,pretty_print=True)
    except:
      raise CException(self.__class__,104,filename)
    #print 'CReport.saveXMLFile',filename,text
    try:
      from core import CCP4Utils
      CCP4Utils.saveFile(self,fileName=filename,text=text)
    except:
      raise CException(self.__class__,105,filename)
    

  def tostring(self,root):
    if root is not None: return ET.tostring(root)
    return ''

  def extractBaseHref(self,root):
    import sys
    self.baseHref = None
    self.resetBaseHref= None
    if sys.platform[0:3] == 'win': return
    ifSame = False
    eleList = root.findall('.//html/head/base')
    if len(eleList)>0 and eleList[0].get('href') is not None:
      self.baseHref = QtCore.QUrl(eleList[0].get('href')).toLocalFile().__str__()
      if self.baseHref is None: return
      from core import CCP4Utils
      localBasePath = os.path.join(CCP4Utils.getCCP4I2Dir(),'docs')
      try:
        ifSame = CCP4Utils.samefile(self.baseHref,localBasePath)
      except:
        ifSame = False
      #print 'extractBaseHref samepath', localBasePath,ifSame
    if not ifSame and self.baseHref is not None and self.baseHref.count('ccp4i2')>0:
      eleList[0].set('href','file://'+localBasePath+'/')
      import tempfile
      path,base0 = os.path.split(self.filename)
      path,dir0 = os.path.split(path)
      self.resetBaseHref = os.path.join(tempfile.gettempdir(),dir0+'_'+base0)
      CCP4Utils.saveEtreeToFile(root,self.resetBaseHref)
      
      #print 'extractBaseHref',self.resetBaseHref

  def extractCCP4Data(self,root):
    import copy
    #for element in root.iterdescendants(): print element.tag
    body = root.find('body')
    if body is not None:
      #for element in body.findall('ccp4:ccp4_data',namespaces=NSMAP):
      for tag in ('ccp4_data','{http://www.ccp4.ac.uk/ccp4ns}ccp4_data'):
        for element in body.findall('.//'+tag):
          #print 'CReport.extractCCP4Data found ccp4_data',etree.tostring(element,pretty_print=True)
          data_type = element.get('type','')
          id = element.get('id','')
          #print 'CReport.extractCCP4Data',data_type,id
          if data_type not in self.ccp4_data: self.ccp4_data[data_type] = {}
          self.ccp4_data[data_type][id] = copy.deepcopy(element)
      

  def parseCCP4Data(self,element):
     data_type = str(element.get('type'))
     name = str(element.get('id'))
     dataClass = DATAMANAGER().getCss(data_type)
     if dataClass is None:
       raise CException(self.__class__,106,data_type)
     #print 'CReport.parseCCP4Data',data_type,name,dataClass
     dataobj = dataClass(parent=self,name=name)
     #print 'CReport.parseCCP4Data',name,data_type,etree.tostring(element,pretty_print=True)
     if data_type == 'CContainer':
       e = dataobj.loadFromEtree(element)
     else:
       #print 'CReport.parseCCP4Data',name,data_type,etree.tostring(element,pretty_print=True)
       e = dataobj.setEtree(element)      
     self.dataArray[name] = dataobj
     return e
     #print 'CReport.parseCCP4Data dataArray ',self.dataArray


  def extractTopFolds(self,root):
    topFolds = []
    foldElements = root.findall(".//span[@class='folder']")
    #print 'into CCP4Report.listTopFolds',len(foldElements)
    for ele in foldElements:
      isTop = True
      parent = ele.findall('..')[0]
      while parent.tag != 'body' and isTop:
        if parent.tag == 'div' and parent.get('class','').count('hidesection'):
          isTop = False
        else:
          parent = parent.findall('..')[0]
      if isTop:
        #print ET.tostring(ele,pretty_print=True)
        if ele.text.count( u"\u25BC") or  ele.text.count( u"\u25B6"):
          #Confused - seems to work if chop first two chrs rather than first seven
          txt = str(ele.text[2:])
          topFolds.append([txt,ele.get('title',None)])
        else:
          topFolds.append([str(ele.text),ele.get('title',None)])     
    #print 'CCP4Report.listTopFolds',topFolds
    self.topFolds = topFolds

  def extractLinkIds(self,root):
    # In reports a link to a potentially non-existant sub-job report has an id that is the jobId
    # of the sub-job. This will be helpful for calling CReportGenerator to create the sub-job report
    # path = 'body/div/a'
    print("############################################################")
    print("extractLinkIds")
    print(root)
    print("############################################################")
    path = 'body/div[@class="sub-job-list"]/span/[@class="folder_link"]/a'
    for ele in root.findall(path):
      print('CReport.extractLinkIds',ele.get('href'),ele.get('id'))
      if ele.get('href') is not None and ele.get('id') is not None and ele.get('id')[0:5] == 'jobId':
        self.linkIds[str(ele.get('href'))] = int(ele.get('id')[5:])
        

  def getLinkId(self,path):
    print('CReport.getLinkId',path,type(path))
    print(self.linkIds)
    return self.linkIds.get(path,None)
     
  def getEtree( self,dataType=None,id=None):
    if dataType is not None:
      if dataType in self.ccp4_data and id in self.ccp4_data[dataType]:
        rv = self.ccp4_data[dataType][id]
        #del self.ccp4_data[dataType][id]
        return rv
    else:
      for dataType in list(self.ccp4_data.keys()):
        if id in self.ccp4_data[dataType]:
          rv = self.ccp4_data[dataType][id]
          #del self.ccp4_data[dataType][id]
          return rv

  def setJobInfo(self,jobInfo={}):
    # A kludge to inject info from core program into web page
    self._jobInfo.update(jobInfo)
    #print 'CReport.setJobInfo',self._jobInfo

  def getJobInfo(self):
    return self._jobInfo
    

  def getDataObject(self,dataType=None,id=None,subElement=None):
    #print 'CReport.getDataObject',data_type,id,'*',sub_element,'*',self.ccp4_data.keys()
    # Check for sane input
    if id is None:
      return None
    elif dataType is not None:
      # Ensure the ccp4_data element has been parsed
      if dataType in self.ccp4_data and id in self.ccp4_data[dataType]:
        self.parseCCP4Data(self.ccp4_data[dataType][id])
        del self.ccp4_data[dataType][id]
    else:
      for dataType in list(self.ccp4_data.keys()):
        if id in self.ccp4_data[dataType]:
          self.parseCCP4Data(self.ccp4_data[dataType][id])
          del self.ccp4_data[dataType][id]
          break

    #print 'CReport.getDataObject dataArray',self.dataArray
    # Return the appropriate CData object
    if id in self.dataArray:
      if dataType == 'CContainer' and subElement is not None:
        return self.dataArray[id].__getattr__(subElement)
      else:
        return self.dataArray[id]
    else:
      return None

#========================================================================================   
import unittest
def TESTSUITE():
  suite = unittest.defaultTestLoader.loadTestsFromTestCase(testReport)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)

class testReport(unittest.TestCase):

  def test1(self):
    from core import CCP4File
    testFile = CCP4File.CDataFile(project='CCP4I2_TOP',relPath='test/data',baseName='test_report.html')
    r = CReport()
    r.loadFromXmlFile(str(testFile))
    x = r.getDataObject(id='final_results')
    self.assertEqual(x.cell.a,123.4,'Failed to load container and cell data')
    x = r.getDataObject(id='data_table_1')
    self.assertEqual(x.nColumns,9,'Failed to load CReportTableModel number of columns')
    
