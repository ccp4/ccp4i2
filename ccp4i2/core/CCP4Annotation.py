from __future__ import print_function

"""
     CCP4Annotation.py: CCP4 GUI Project
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
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See theShostna
     GNU Lesser General Public License for more details.
"""

"""
   Liz Potterton Aug 2010 - Generic annotation object
"""
import os
import re
import time
import types
from core import CCP4Data
from core.CCP4Data import CBaseData
from core import CCP4Utils
from core.CCP4ErrorHandling import *

class CTime(CCP4Data.CInt):
    '''The time. Uses Python time module'''
    QUALIFIERS = {'min' : 0, 'label' : 'Time', 'toolTip' : 'Time and date as hh:mm dd/mm/yyyy', 'format' : '%H:%M %d/%b/%y'}
    QUALIFIERS_ORDER = ['format']
    QUALIFIERS_DEFINITION = {'format' : {'type' : str,
                                         'description' : 'Argument to Python time.strftime to display time in human readable format'}}

    def setCurrentTime(self):
        self.set(int(time.time()))
    
    def setEtree(self,element,checkValidity=True):
        #print 'CTime.setEtree',element.text
        try:
            localtime = time.strptime(str(element.text),self.qualifiers('format'))
        except:
            try:
                localtime = time.strptime(str(element.text), '%H:%M %d/%b/%Y')
            except:
                localtime = time.gmtime()
        epochTime = time.mktime(localtime)
        self.__dict__['_value'] = int(epochTime)

    def __str__(self):
        if self._value is None:
            return ''
        else:
            return time.strftime(self.qualifiers('format'),time.localtime(self._value))

    def __int__(self):
        if self._value is None:
            return 0
        else:
            return self.__dict__['_value']

    def date(self):
        return time.strftime("%d/%b/%y", time.localtime(self._value))

    def set(self,value):
        # Try setting assuming value is an int
        try:
            CBaseData.set(self,value)
        except CException as e:
            try:
                localtime = time.strptime(value,self.qualifiers('format'))
                epochTime = time.mktime(localtime)
                self.__dict__['_value'] = epochTime
            except:
                raise e

    def getQDateTime(self):
        from PySide2 import QtCore
        q = QtCore.QDateTime()
        if self.__dict__['_value'] is not None:
            q.setTime_t(int(self.__dict__['_value']))
            q.addMSecs(int(1000*(self.__dict__['_value']-int(self.__dict__['_value']))))
        return q

    def setQDateTime(self,QDateTime=None):
        t = QDateTime.toTime_t()
        print('CTime.setQDateTime',t)
        self.__dict__['_value'] = t

    PROPERTIES = {'QDateTime' : {'fget' : getQDateTime, 'fset' : setQDateTime}}

class CUserId(CCP4Data.CString):
    """A user ID"""
    QUALIFIERS = {}
    QUALIFIERS_DEFINITION = {}
    QUALIFIERS = { 'label' : 'User id',
                     'toolTip' : 'User id as me@myplace.ac.uk' }

    def setCurrentUser(self):
        self.set(CCP4Utils.getUserId())


class CHostName(CCP4Data.CString):
    '''Computer name'''
    QUALIFIERS = {'label' : 'Machine name',
                  'toolTip' : 'Hostname as mycomputer.myplace.ac.uk'}

    def setCurrentHost(self):
        self.set( CCP4Utils.getHostName() )


class CHostname(CHostName):
    pass


class CUserAddress(CCP4Data.CData):
    """User id and platform node"""
    CONTENTS = {'platformNode' : {'class' : CCP4Data.CString},
                'userId' : {'class' : CUserId}}
    CONTENTS_ORDER = ['platformNode','userId']
    QUALIFIERS = {'label' : 'User id and current machine',
                  'toolTip' : 'User id as me@myplace.ac.uk and machine name'}

    def setCurrent(self):
        import platform
        self._value['userId'].setCurrentUser()
        self._value['platformNode'].set(platform.node())

    def __str__(self):
        return str(self.userId)+'@'+str(self.platformNode)


def currentUserAddress():
    userAddress = CUserAddress()
    userAddress.setCurrent()
    return str(userAddress)


class CAnnotation(CCP4Data.CData):
    """Annotation text with user id and time"""
    CONTENTS = {'text' : {'class' : CCP4Data.CString, 'qualifiers' : {'allowUndefined' : True, 'charWidth' : -1}},
                'time' : {'class' : CTime , 'qualifiers' : {'allowUndefined' : True, 'default' : None}},
                'author' : {'class' : CUserId, 'qualifiers' : {'allowUndefined' : True, 'default' : None}}}
    QUALIFIERS = {'label' : 'Annotation', 'toolTip' : 'Enter your comments'}

    def __init__(self,value={}, qualifiers= {}, parent=None, name=None,**kw):
        # Enable initialisiation with just a text string
        if isinstance(value,(str,CCP4Data.CString)):
            valu = {'text' : value}
        else:
            valu = {}
            valu.update(value)
        qualis = {}
        qualis.update(qualifiers)
        if len(kw) > 0:
            vs, qs, us = CCP4Data.sortArguments(self.__class__, kw)
            if len(us) > 0:
                for key in list(us.keys()):
                    print(self.className(), 'input argument unrecognised:', key)
            qualis.update(qs)
            valu.update(vs)
        CCP4Data.CData.__init__(self,value=valu,qualifiers=qualis,parent=parent,name=name)

    def set(self,value={},**kw):
        # Reimplement to allow input of a Python string
        # and to automatically add time and author
        #print 'CAnnotation.set value',value,kw
        if isinstance(value,str):
            valu = {'text' : value}
        elif isinstance(value, CAnnotation):
            valu = value
        elif value is None:
            valu = {}
            valu.update(kw)
        else:
            valu = {}
            valu.update(value)
            valu.update(kw)
        #print 'CAnnotation.set valu',valu
        CCP4Data.CData.set(self,valu)
        if not self.time.isSet():
            self.time.setCurrentTime()
        if not self.author.isSet():
            self.author.setCurrentUser()


class CAnnotationList(CCP4Data.CList):
    """A list of annotation"""
    SUBITEM = { 'class' : CAnnotation } 

    def append(self,arg):
        if isinstance(arg,(str,CCP4Data.CString)):
            anno = CAnnotation(text=arg)
            #print 'CAnnotationList.append',anno
            CCP4Data.CList.append(self,anno)
        else:
            CCP4Data.CList.append(self,arg)

'''
class CUsersList(CCP4Data.CList):
    SUBITEM = { 'class' : CUserId }
    QUALIFIERS = { 'default' : [CUserId('liz')] }
    pass
'''

class CFont(CCP4Data.CData):
    '''Simplified Qt font options'''
    #QtGui.QFont styles
    StyleNormal = 0
    StyleItalic = 1
    StyleOblique = 2

    CONTENTS = {'family' : {'class' : CCP4Data.CString , 'qualifiers' : {'default' : 'Helvetica'}},
                'style'  : {'class' : CCP4Data.CInt , 'qualifiers' : {'onlyEnumerators' : True,
                                                                      'default' : StyleNormal,
                                                                      'enumerators' : [StyleNormal, StyleItalic, StyleOblique],
                                                                      'menuText' : ['normal','italic','oblique']}},
               'pointSize' : {'class' : CCP4Data.CInt, 'qualifiers' : {'min' : 1, 'default' : 12}},
               'weight' :{'class' : CCP4Data.CInt,
                          'qualifiers' : {'min' : 0, 'max' : 99, 'default' : 50, 'allowUndefined' : False,
                                          'enumerators' : [25,50,63,75,87],
                                          'menuText' : ['light', 'normal', 'demi-bold', 'bold', 'black']}}}

    def getQFont(self):
        from PySide2 import QtGui, QtWidgets
        italic = (self.style == QtGui.QFont.StyleItalic)
        f = QtGui.QFont(str(self.family),int(self.pointSize),int(self.weight),italic)
        #print 'CFont.getQFont', f, str(f.family()),f.pointSize()
        return f

    def setQFont(self,font=None):
        if font is None:
            return
        self.family = str(font.family())
        self.style = font.style()
        self.pointSize = font.pointSize()
        self.weight = font.weight()
        self.dataChanged.emit()


class CAuthor(CCP4Data.CString):
    '''Placeholder for bibliographic author'''
    pass


class CBibReference(CCP4Data.CData):
    '''Bibliographic reference'''
    CONTENTS_ORDER = ['pmid', 'title', 'authorList', 'source', 'url' , 'selected']
    CONTENTS = {'pmid' : {'class' : CCP4Data.CInt}, 'title' : {'class' : CCP4Data.CString},
                'authorList' : {'class' : CCP4Data.CList, 'subItem' :{'class' : CAuthor}},
                'source' : {'class' : CCP4Data.CString}, 'url' : {'class' : CCP4Data.CString},
                'selected' : {'class' : CCP4Data.CBoolean}}
    ERROR_CODES = {101 : {'description' : 'Failed to load Medline data'}}
    '''
    def loadFromMedline(self,text):
        #print 'loadFromMedline',text
        import re
        m = re.search(r'TI  -(.*)\n[A-Z]',text,re.DOTALL)
        if m is not None:
          self.title = m.groups()[0].strip()
        else:
          return CErrorReport(self.__class__,101)
        m = re.search(r'SO  -(.*)',text)
        if m is not None:
          self.source = m.groups()[0].strip()
        else:
          return CErrorReport(self.__class__,101)
        m = re.findall(r'AU  -(.*)',text)
        for item in m:
          self.authorList.addItem()
          self.authorList[-1].set(item.strip())
        m = re.search(r'URL -(.*)',text)
        if m is not None: self.url = m.groups()[0].strip()
        return CErrorReport()
        #print 'CBibReferenceGroup.loadFromMedline',self.get()
    '''

    def loadFromMedline(self,text):
        lineList = text.split('\n')
        idx = 0
        while idx < len(lineList):
            if lineList[idx][0:5] == 'TI  -':
                self.title.set(lineList[idx][5:])
                while not '-' in lineList[idx+1][0:5]:
                    idx += 1
                    self.title.set(str(self.title) + lineList[idx][5:])
            elif lineList[idx][0:5] == 'SO  -':
                self.source.set(lineList[idx][5:])
            elif lineList[idx][0:5] == 'AU  -':
                print('author', lineList[idx])
                self.authorList.addItem()
                self.authorList[-1].set(lineList[idx][5:])
            elif lineList[idx][0:5] == 'URL -':
                self.url.set(lineList[idx][5:])
            idx += 1
        if not self.source.isSet() and not self.url.isSet():
            return CErrorReport(self.__class__, 101)
        return CErrorReport()

    def get0(self):
        authorText = ''
        if len(self.authorList)>0:
            for authorObj in self.authorList:
                authorText += authorObj.get0()+', '
        rv = {'title' : self.__dict__['_value']['title'], 'authorList' : authorText[0:2],
              'source' : self.__dict__['_value']['source'],
              'selected' : self.__dict__['_value']['selected'] }
        return rv


class CBibReferenceGroup(CCP4Data.CData):
    '''Set of bibliographic references for a task'''
    CONTENTS_ORDER = ['taskName', 'version', 'title', 'references']
    CONTENTS = {'taskName' : {'class' : CCP4Data.CString},
                'version' : {'class' : CCP4Data.CString},
                'title' : {'class' : CCP4Data.CString},
                'references' : {'class' : CCP4Data.CList , 'subItem' : { 'class' : CBibReference}}}
    ERROR_CODES = {100 : { 'description' : 'Failed attempting to load MedLine file - file not found'},
                   101 : { 'description' : 'Failed attempting to find references file'},
                   102 : { 'description' : 'Error copying file'}}

    def loadFromMedline(self, fileNameList=[], taskName=None, version=None):
        from core import CCP4TaskManager
        if len(fileNameList) == 0 and taskName is not None:
            self.__dict__['fileNameList'] = CCP4TaskManager.TASKMANAGER().searchReferenceFile(taskName,version=version)
            if len(fileNameList) == 0:
                self.__dict__['fileNameList'] = CCP4TaskManager.TASKMANAGER().searchReferenceFile(taskName,version=version,drillDown=True)
        else:
            self.__dict__['fileNameList'] = fileNameList
        #print 'CBibReferenceGroup.loadFromMedline',fileNameList
        self.__dict__['_value']['taskName'].set(taskName)
        self.__dict__['_value']['version'].set(version)
        self.__dict__['_value']['title'].set('References for '+CCP4TaskManager.TASKMANAGER().getTitle(taskName,version=version))
        for fileName in self.__dict__['fileNameList']:
            text = CCP4Utils.readFile(fileName=fileName)
            textList = text.split('\nPMID- ')
            for text in textList:
                self.__dict__['_value']['references'].addItem()
                rv = self.__dict__['_value']['references'][-1].loadFromMedline(text)
                if len(rv)>0:
                    del self.__dict__['_value']['references'][-1]
        if len(self.__dict__['_value']['references'])>0:
            self.__dict__['_value']['references'][0].selected = True

    def export(self, cformat, fileName):
        import shutil
        err = CErrorReport()
        #print 'CBibReferenceGroup.export',cformat,fileName
        if os.path.splitext(fileName)[1] == '':
            fileName = fileName +'.' + cformat + '.txt'
        #fileNameList = CCP4TaskManager.TASKMANAGER().searchReferenceFile(name=self.taskName.get(),version=self.version.get(),cformat=format)
        #if len(fileNameList)==0:
        #  err.append(self.__class__,101,details='For taskname: '+str(self.taskName)+' version: '+str(self.version)+' cformat: '+str(format))
        #else:
        if cformat == 'bibtex':
            source = re.sub(r'\.medline','.bibtex',self.__dict__['fileNameList'][0])
        else:
            source = self.__dict__['fileNameList'][0]       
        try:
            shutil.copyfile(source,fileName)
        except:
            err.append(self.__class__,102,details='From: '+str(source)+' to: '+str(fileName))
        return err

class CMetaDataTag(CCP4Data.CData):
    '''This class will extend list of enumerators if new value for string is entered'''
    CONTENTS_ORDER = ['tag']
    CONTENTS = {'tag' : {'class' : CCP4Data.CString}}
    QUALIFIERS = {'enumeratorsFunction' : None,
                  'addEnumeratorFunction' : None}
    QUALIFIERS_ORDER = ['enumeratorsFunction', 'addEnumeratorFunction']
    QUALIFIERS_DEFINITION = {'enumeratorsFunction' : {'type' : types.MethodType ,
                                                      'definition' : 'Function returning list of enumerators'},
                             'addEnumeratorFunction' : {'type' : types.MethodType ,
                                                        'definition' : 'Function to add to list of enumerators'}}

    def __init__(self, **kw):
        CCP4Data.CData.__init__(self, **kw)
        '''
        print 'CMetaDataTag.__init__',repr(self.parent()),self.parent
        try:
          print self.parent().index(self)
        except:
          print 'no index'
        import traceback
        traceback.print_stack(limit=6)
        '''

    def validity(self, arg):
        if isinstance(arg,CMetaDataTag):
            arg=arg.get()
        return self.tag.validity(arg['tag'])

    def getTextItem(self):
        if self.tag.isSet() and len(self.tag) > 0:
            return self.tag.__str__()
        else:
            return '--'

    def getEnumerators(self):
        func = self.qualifiers('enumeratorsFunction')
        if func is None:
            return []
        else:
            return func()

    def addEnumerator(self, tag):
        if tag is None or len(tag)==0: return -1
        if tag in self.qualifiers('enumeratorsFunction')(): return -1
        return self.qualifiers('addEnumeratorFunction')(tag)


class CMetaDataTagList(CCP4Data.CList):
    SUBITEM = {'class' : CMetaDataTag}
    QUALIFIERS = {'listMinLength' : 1}

    def __init__(self,**kw):
        CCP4Data.CList.__init__(self,**kw)
        #print 'CMetaDataTagList.__init__',len(self)

    def set(self, arg, **kw):
        if isinstance(arg,list) and len(arg) > 0 and (arg[0] is None or isinstance(arg[0], str)):
            arg0 = []
            for item in arg:
                arg0.append({'tag' : item})
            #print 'CMetaDataTagList.set',arg0
            return CCP4Data.CList.set(self, arg0, **kw)
        else:
            return CCP4Data.CList.set(self, arg, **kw)


class CServerGroup(CCP4Data.CData):
    from core.CCP4File import CDataFile
    '''One or more compute servers used in "remote" running'''
    CONTENTS = {'name' : {'class' : CCP4Data.CString},
                'mechanism' : {'class' : CCP4Data.CString , 'qualifiers' :
                              {'enumerators' : ['ssh', 'ssh_shared', 'qsub_local', 'qsub_remote','slurm_remote', 'custom'],
                               'menuText' : ['ssh', 'ssh with shared filesystem', 'qsub queue', 'qsub on another machine','Slurm on another machine', 'custom'],
                               'onlyEnumerators' : True, 'default' : 'ssh'}},
                'serverList' : {'class' : CCP4Data.CList, 'qualifiers': {'minLength' : 1},
                                'subItem' : {'class' : CHostname , 'qualifiers' : { 'default' : None , 'allowUndefind' : False}}},
                'userExtensible' : {'class' : CCP4Data.CBoolean, 'qualifiers' : { 'default' : False}},
                'customCodeFile' : {'class' : CDataFile, 'qualifiers' : {'allowUndefind' : True, 'fileExtensions' : ['py']}},
                'queueOptionsFile' : {'class' : CDataFile, 'qualifiers' : {'allowUndefind' : True}},
                'ccp4Dir' : {'class' : CCP4Data.CString, 'qualifiers' : {'allowUndefind' : True}},
                'tempDir' : {'class' : CCP4Data.CString, 'qualifiers' : {'allowUndefind' : True}},
                'sge_root' : {'class' : CCP4Data.CString, 'qualifiers' : {'allowUndefind' : True}},
                'keyFilename' : {'class' : CCP4Data.CString, 'qualifiers' : {'allowUndefind' : True}},
                'validate' : {'class' : CCP4Data.CString, 
                              'qualifiers' : {'onlyEnumerators' : True, 'default' : 'password',
                                              'enumerators' : ['password', 'key_filename', 'pass_key_filename']}},
                'timeout' : {'class' : CCP4Data.CFloat},
                'maxTries' : {'class' : CCP4Data.CInt, 'qualifiers' : {'default' : 2}}}
    CONTENTS_ORDER = ['name', 'mechanism', 'serverList', 'userExtensible', 'ccp4Dir', 'tempDir',
                      'sge_root', 'keyFilename', 'validate', 'customCodeFile', 'queueOptionsFile']

    #'enumerators' : ['ssh', 'ssh_shared', 'qsub_local', 'qsub_remote', 'qsub_shared', 'custom'],
    #                          'menuText' : ['ssh', 'ssh with shared filesystem', 'qsub on this machine', 'qsub on another machine', 'qsub with shared filesystem', 'custom' ],


class CDateRange(CCP4Data.CData):
    '''A date range - may be on a scale of years,months or days'''
    CONTENTS = { 'year' : {'class' : CCP4Data.CInt , 'qualifiers' : {'enumerators' : []}},
                 'month' : {'class' : CCP4Data.CString ,  'qualifiers' : {'onlyEnumerators' : True,
                            'enumerators' : ['January', 'February' ,'March', 'April', 'May', 'June',
                                             'July', 'August' ,'September', 'October', 'November', 'December'],
                            'default' : 'January'}},
                 'day' : { 'class' : CCP4Data.CInt, 'qualifiers' : {'default' : 1, 'min' : 1, 'max' : 31}},
                 'yearRange' :  {'class' : CCP4Data.CInt, 'qualifiers' : {'default' : 0, 'min' : 0, 'max' : 100}},
                 'monthRange' : {'class' : CCP4Data.CInt, 'qualifiers' : {'default' : 0, 'min' : 0, 'max' : 12}},
                 'dayRange' : {'class' : CCP4Data.CInt, 'qualifiers' : {'default' : 0, 'min' : 0, 'max' : 30}}}
    CONTENTS_ORDER = ['year', 'month', 'day', 'yearRange', 'monthRange', 'dayRange']
    MONTHLENGTH = [ 0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]  # KJS - Err.. Gregorian calendar & all that ?
    MONTHLIST = ['January', 'February' ,'March', 'April', 'May', 'June',
                 'July', 'August' ,'September', 'October', 'November', 'December']
    
    def __init__(self, value={}, qualifiers= {}, parent=None, name=None, **kw):
        CCP4Data.CData.__init__(self,value=value,qualifiers= qualifiers,parent=parent,name=name,**kw)
        # Setting defaults for year enumerator etc on assumption that we a looking to the recent past
        year, month, day = time.strftime("%Y %B %d",time.localtime(time.time())).split()
        print('CDateRange.__init__',day,month,year)
        yearList = []
        for y in range(2010,int(year)+1):
            yearList.insert(0,y)
        self.__dict__['_value']['year'].setQualifier('enumerators', yearList)
        self.__dict__['_value']['year'].set(int(year))
        self.__dict__['_value']['month'].set(month)
        self.__dict__['_value']['day'].set(day)
        dayList = list(range(1,self.MONTHLENGTH[self.MONTHLIST.index(month)+1]+1))
        self.__dict__['_value']['day'].setQualifier('enumerators', dayList)

    def dateString(self, addYear=0, addMonth=0, addDay=0):
        iYear = int(self.year)+addYear
        if not self.month.isSet():
            iMon = addMonth
        else:
            iMon = (self.MONTHLIST.index(str(self.month)) + addMonth)
        if self.day.isSet():
            iDay = int(self.day)+addDay
        else:
            iDay = addDay
        iMon = iMon + iDay/30
        iDay = max(1,iDay%30)
        iYear = iYear + iMon/12
        iMon = iMon%12
        date = str(iYear) + ' ' + self.MONTHLIST[iMon] + ' ' + ('0' + str(iDay))[-2:]
        #print 'dateString',date
        return date

    def epochTime(self,dateString=None):
        localtime = time.strptime(dateString,"%Y %B %d")
        #print 'CDateRange.epochTime localtime',localtime
        epochTime = time.mktime(localtime)
        #print 'CDateRange.epochTime',epochTime
        return epochTime

    def epochRange(self):
        start = self.epochTime(self.dateString())
        nY =nM = nD= 0
        if self.yearRange.isSet():
            nY = int(self.yearRange)
        if self.monthRange.isSet():
            nM = int(self.monthRange)
        if self.dayRange.isSet():
            nD = int(self.dayRange)
        if not self.month.isSet():
            # We are dealing in whole years - so the end is start of next year
            end = self.epochTime(self.dateString(addYear=1+nY,addMonth=nM,addDay=nD))
        elif not self.day.isSet():
            # We are dealing in whole month so end is start of next month
            end = self.epochTime(self.dateString(addYear=nY,addMonth=1+nM,addDay=nD))
        else:
            end = self.epochTime(self.dateString(addYear=nY,addMonth=nM,addDay=1+nD))
        if nY > 0 or nM > 0 or nD > 0:
            start = self.epochTime(self.dateString(addYear=-nY,addMonth=-nM,addDay=-nD))
        #print 'epochRange',start,end,time.strftime("%Y %B %d",time.localtime(end))
        return start,end


#===========================================================================================================
import unittest

def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testAnnotation)
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)

class testAnnotation(unittest.TestCase):

    def test1(self):
        a = CAnnotation('Try this text')
        self.assertEqual(a.text,'Try this text','Failed setting CAnnotation text')
        self.assertEqual(a.author,CCP4Utils.getUserId(),'Failed setting CAnnotation author')
        if a.time < int(time.time())-10 or a.time > int(time.time()):
            self.fail('Failed setting CAnnotation time')

    def test2(self):
        a = CAnnotation(text='Try this text')
        t = CTime()
        t.setCurrentTime()
        d = a.time-t
        self.assertEqual(d.__class__, CTime, 'Time difference not a CTime object')
        self.assertEqual(0 <= d < 10, True, 'Time difference notin expected range')

    def test3(self):
        annoList = CAnnotationList()
        annoList.append('Test this string')
        annoList.append(CAnnotation('Test another string'))
        self.assertEqual(len(annoList),2,'Annotation list wrong length')
        self.assertEqual(annoList[1].text,'Test another string','Annotation wrong text')
