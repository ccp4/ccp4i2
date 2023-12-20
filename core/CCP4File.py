from __future__ import print_function


"""
     CCP4File.py: CCP4 GUI Project
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
   Liz Potterton Aug 2010 - File handling classes
"""

import os
import re
import sys
import types
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from PySide2 import QtCore

from core import CCP4Data
from core import CCP4Config
from core.CCP4Modules import MIMETYPESHANDLER, PROJECTSMANAGER
from core.CCP4ErrorHandling import *
from report.CCP4ReportParser import CCP4NS

from xml.etree import ElementTree as ET

# Version number of form n.m or n.m.i
class CVersion(CCP4Data.CString):
    '''A (string) version number of the form n.m.i'''

    PYTHONTYPE = str
    QUALIFIERS = {'allowUndefined' : True, 'default' : None, 'charWidth' : 10}
    QUALIFIERS_ORDER = ['allowUndefined', 'default', 'charWidth']
    QUALIFIERS_DEFINITION = {'allowUndefined' :{'type' : bool,
                                                'description' : 'Flag if allow an unset value at run time'},
                             'default' : {'description' : 'A default value'},
                             'charWidth' : {'type' : int,
                                            'description' : 'Number of characters allowed for widget in GUI'}}
    ERROR_CODES = {101 : {'description' : 'Version is not of form n.m or n.m.i'}}

    def validity(self, arg):
        validityObj = CErrorReport()
        if arg is None:
            if self.qualifiers('allowUndefined'):
                validityObj.append(self.__class__, 1, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            else:
                validityObj.append(self.__class__, 2, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            return validityObj
        if not isinstance(arg, str):
            validityObj.append(self.__class__, 5, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            return validityObj
        if arg.count('alpha_rev'):
            pass
        else:
            items = arg.split('.')
            if len(items) < 2 or len(items) > 3:
                validityObj.append(self.__class__, 101, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            else:
                for i in items:
                    try:
                        ii = int(i)
                    except:
                        validityObj.append(self.__class__, 101, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
                        break
        return validityObj

    def setCurrent(self):
        try:
            #MN Change...This uses CCP4Update module to get_revno....importing CCP4Update
            #is not safe, and in some circumstances does not work, so I have put it within the try except
            #clause
            from core import CCP4Update
            rev = CCP4Update.get_revno()
        except:
            rev = 0
        if rev > 0:
            self.__dict__['_value'] = 'alpha_rev_' + str(rev)
        else:
            try:
                self.__dict__['_value'] = CCP4Config.VERSION()
            except:
                self.__dict__['_value'] = 'Unknown'


class CProjectName(CCP4Data.CString):
    '''The name of a CCP4i project or directory alias'''

    PYTHONTYPE = str
    QUALIFIERS = {'allowUndefined' : True, 'allowAlias' : True, 'allowUnfound' : True, 'default' : None}
    QUALIFIERS_ORDER = ['allowUndefined','allowAlias','allowUnfound','default']
    QUALIFIERS_DEFINITION = {'allowUndefined' : {'type' : bool,
                                                 'description' : 'Flag if allow undefined value at run time' },
                             'allowAlias' : {'type' : bool,
                                             'description' : 'Flag if allow project to be directory alias at run time'},
                             'allowUnfound' : {'type' : bool,
                                               'description' : 'Flag if allow unfound project at run time'},
                             'default' : {'type' :str}}
    ERROR_CODES = {101 : {'description' : 'Invalid project name'},
                   102 : {'description' : 'Project does not have directory set'},
                   103 : {'description' : 'Project directory does not exist'},
                   104 : {'severity' : SEVERITY_WARNING,
                          'description' : 'Warning - Project name is a directory alias'},
                   105 : {'severity' : SEVERITY_WARNING,
                          'description' : 'Warning - Project does not have directory set'},
                   106 : {'severity' : SEVERITY_WARNING,
                          'description' : 'Warning - Project directory does not exist'}}

    def validity(self, arg):
        v = CErrorReport()
        #print 'CProjectName.validity',arg,type(arg),self.parent()
        #print 'CProjectName.validity returning OK',self.parent().objectName()
        # NOT FULLY IMPLEMENTED
        # Need to test that it is valid project name
        if arg is None:
            if self.qualifiers('allowUndefined'):
                v.append(self.__class__, 1, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            else:
                v.append(self.__class__, 2, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
                #print 'arg None', len(v)
            return v
        elif arg == 'FULLPATH':
            return v

        d = PROJECTSMANAGER().getProjectDirectory(projectName=arg, testAlias=True)
        #print 'CProjectName.validity',arg,d
        if d is None:
            if not self.qualifiers('allowAlias'):
                v.append(self.__class__, 101, arg, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            return v
        if not os.path.exists(d):
            if self.qualifiers('allowUnfound'):
                v.append(self.__class__, 106, arg, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            else:
                v.append(self.__class__, 103, arg, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        return v

    def directory(self):
        if self.__dict__['_value'] == 'workDirectory':
            from core import CCP4PluginScript
            obj = self
            while isinstance(obj,CCP4Data.CData):
                obj = obj.parent()
                if isinstance(obj,CCP4PluginScript.CPluginScript):
                    return obj.getWorkDirectory()
            return None
        d = PROJECTSMANAGER().getProjectDirectory(projectName=str(self._value), testAlias=True)
        if d is None:
            return CFilePath()
        else:
            return CFilePath(d)

class CProjectId(CCP4Data.CUUID):
    '''The CCP4i2 database project id - a global unique id'''
    PYTHONTYPE = str
    QUALIFIERS = {'allowUndefined' : True, 'allowUnfound' : True, 'default' : None}
    QUALIFIERS_ORDER = ['allowUndefined','allowUnfound','default']
    QUALIFIERS_DEFINITION = {'allowUndefined' : {'type' : bool,
                                                 'description' : 'Flag if allow undefined value at run time'},
                             'allowUnfound' : {'type' : bool ,
                                               'description' : 'Flag if allow unfound project at run time'},
                             'default' : {'type' :str}}
    ERROR_CODES = {201 : {'description' : 'Unrecognised projectId'},
                   202 : {'description' : 'Project does not have directory set'},
                   203 : {'description' : 'Project directory does not exist'},
                   205 : {'severity' : SEVERITY_WARNING,
                          'description' : 'Warning - Project does not have directory set'},
                   206 : {'severity' : SEVERITY_WARNING,
                          'description' : 'Warning - Project directory does not exist'}}

    def validity(self, arg):
        v = CErrorReport()
        # NOT FULLY IMPLEMENTED
        # Need to test that it is valid project name
        if arg is None:
            if self.qualifiers('allowUndefined'):
                v.append(self.__class__, 1, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            else:
                v.append(self.__class__, 2, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            return v
        arg = str(arg)
        d = PROJECTSMANAGER().getProjectDirectory(projectId=arg, testAlias=True)
        if d is None:
            if not self.qualifiers('allowUnfound'):
                v.append(self.__class__, 201, arg, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
                return v
            else:
                return v
        if not os.path.exists(d):
            if self.qualifiers('allowUnfound'):
                v.append(self.__class__, 206, arg, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            else:
                v.append(self.__class__, 203, arg, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        return v

    def directory(self):
        #print 'CProjectId.directory', self._value,PROJECTSMANAGER().getProjectDirectory(projectId=str(self._value))
        if self.__dict__['_value'] == 'workDirectory':
            print("CProjectId.directory IS workDirectory")
            from core import CCP4PluginScript
            obj = self
            while isinstance(obj, CCP4Data.CData):
                obj = obj.parent()
                if isinstance(obj, CCP4PluginScript.CPluginScript):
                    return obj.getWorkDirectory()
            return None
        d = PROJECTSMANAGER().getProjectDirectory(projectId=str(self._value))
        if d is None:
            return CFilePath()
        else:
            return CFilePath(d)


class CFilePath(CCP4Data.CString):
    '''A file path'''
    ALLOWED_CHARACTERS_IGNORE = 0
    ALLOWED_CHARACTERS_WARN = 1
    ALLOWED_CHARACTERS_FAIL = 2
    ALLOWED_CHARACTERS_FIX = 3
    QUALIFIERS = {'allowUndefined' : True, 'allowedCharacters' : '',
                  'allowedCharactersMode': ALLOWED_CHARACTERS_WARN, 'default' : None }
    QUALIFIERS_ORDER = ['allowUndefined', 'allowedCharacters', 'allowedCharactersMode', 'default']
    QUALIFIERS_DEFINITION = {'allowUndefined' : {'type' : bool,
                                                 'description' : 'Flag if allow undefined value at run time'},
                             'allowedCharacters' : {'type' : str,
                                                    'description' : 'Set of characters allowed in file name'},
                             'allowedCharactersMode': {'type' : int,
                                                       'description' : 'Handling of violation of allowed characters'},
                             'default' : {'type' : str, 'description' : 'Default file path'}}
    ERROR_CODES = {101 : {'description' : 'Invalid characters in file name'},
                   102 : {'severity' : SEVERITY_WARNING, 'description' : 'Invalid characters in file name'}}

    def coerce(self, arg):
        # Convert \ to / as Python and Qt work with that internally
        # and deal with it on Windows (we shall see!)
        arg = CCP4Data.CString.coerce(self, arg)
        if arg is None:
            return arg
        else:
            arg = re.sub(r'\\', '/', arg)
        return arg

    def validity(self, arg):
        # Dr. Google says checking if file path is valid is very difficult
        # and dependent on the underlying filesystem.
        # Insist on no weird characters - ???What is part of 'standard' Windows paths?
        v = CErrorReport()
        if arg is None:
            if self.qualifiers('allowUndefined'):
                v.append(self.__class__, 1, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            else:
                v.append(self.__class__, 2, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            return v
        charMode = self.qualifiers('allowedCharactersMode')
        if charMode > 0:
            invalidChars = re.findall('[^a-zA-Z0-9/_.: \-+' + self.qualifiers('allowedCharacters') + ']', arg)
            if len(invalidChars) > 0:
                if [CFilePath.ALLOWED_CHARACTERS_WARN,CFilePath.ALLOWED_CHARACTERS_FIX].count(charMode):
                    v.append(self.__class__, 102, invalidChars.__repr__().lstrip('[').rstrip(']'), stack=False, name=self.objectPath(), label=self.qualifiers('guiLabel'))
                else:
                    v.append(self.__class__, 101, invalidChars.__repr__().lstrip('[').rstrip(']'), stack=False, name=self.objectPath(), label=self.qualifiers('guiLabel'))
        return v

    def fix(self, arg):
        # NOT FULLY IMPLEMENTED
        # ??? Is the str.translate() function useful here
        # Replace an unallowed char with underscore
        return re.sub('[^a-zA-Z0-9/_.: \-' + self.qualifiers('allowedCharacters') + ']', '_', arg)

    def abspath(self):
        if self._value is None:
            return CFilePath()
        else:
            return CFilePath(os.path.abspath(self._value))

    def dirname(self):
        if self._value is None:
            return CFilePath()
        else:
            return CFilePath(os.path.dirname(self._value))

    def basename(self):
        if self._value is None:
            return CFilePath()
        else:
            return CFilePath(os.path.basename(self._value))

    def exists(self):
        if self._value is None:
            return False
        else:
            return os.path.exists(os.path.normpath(self._value))

    def isfile(self):
        if self._value is None:
            return False
        else:
            return os.path.isfile(self._value)

    def isdir(self):
        if self._value is None:
            return False
        else:
            return os.path.isdir(self._value)

    def samefile(self, arg):
        other = self.getValue(arg)
        if other is None or self._value is None:
            return False
        else:
            from core import CCP4Utils
            return CCP4Utils.samefile(self._value, other)

    def pathsplit(self):
        if self._value is None:
            return []
        else:
            return os.path.split(self._value)

    def splitext(self):
        return os.path.splitext(self._value)

    def __add__(self, arg):
        other = self.getValue(arg)
        #print 'in __add__', other
        if other is None:
            return CFilePath(self._value)
        elif self._value is None:
            return CFilePath(other)
        else:
            return CFilePath(os.path.normpath(os.path.join(self._value, other)))

    def __radd__(self, arg):
        #print 'in __radd__'
        other = self.getValue(arg)
        if other is None:
            return CFilePath(self._value)
        elif self._value is None:
            return CFilePath(other)
        else:
            return CFilePath(os.path.normpath(os.path.join(other, self._value)))

    """
    def __cmp__(self,arg):
        other =  self.getValue(arg)
        if self._value is None or other is None:
          if self._value is not None:
            return -1
          elif other is not None:
            return 1
          else:
            return 0
        else:
          #MN 'str' and 'unicode' types do not have __cmp__, but equality can be tested using ==
          if isinstance(self._value,basestring) and isinstance(other,basestring): return self._value == other
          return self._value.__cmp__(other)
    """

    def __str__(self):
        if self._value is None:
            return ''
        else:
            return os.path.normpath(self._value.__str__())


class CDataFile(CCP4Data.CData):
    dataLoaded = QtCore.Signal()
    '''A data file - expected  to have associated class for file contents'''
    SEPARATOR = '-'
    CONTENTS = {'project' : {'class' : CProjectId }, 'baseName' : {'class' : CFilePath},
                'relPath' : {'class' : CFilePath }, 'annotation' : {'class' : CCP4Data.CString},
                'dbFileId' : {'class' : CCP4Data.CUUID},
                'subType' : {'class' : CCP4Data.CInt, 'qualifiers' : {'default' : None}},
                'contentFlag' : {'class' : CCP4Data.CInt, 'qualifiers': {'min' : 0, 'default' : None}}}
    QUALIFIERS = {'allowUndefined' : True, 'mustExist' : False, 'fromPreviousJob' : False,
                  'jobCombo' : True, 'mimeTypeName' : '', 'mimeTypeDescription' : '',
                  'fileLabel' : None, 'fileExtensions' : [], 'fileContentClassName' : None,
                  'isDirectory' : False, 'saveToDb' : True, 'requiredSubType' : None, 'requiredContentFlag' : None }

    QUALIFIERS_ORDER = ['fileExtensions', 'mimeTypeName', 'mimeTypeDescription', 'fileLabel', 'allowUndefined',
                        'mustExist', 'fromPreviousJob', 'jobCombo', 'fileContentClassName', 'isDirectory',
                        'saveToDb', 'requiredSubType', 'requiredContentFlag']
    QUALIFIERS_DEFINITION = {'allowUndefined' : {'type' : bool,
                                                 'description' : 'Flag if data file can be undefined at run time'},
                             'mustExist' : {'type' : bool,
                                            'description' :'Flag if data file must exist at run time'},
                             'fromPreviousJob' : {'type' : bool,
                                                  'description' :'Flag if input data file can be inferred from preceeding jobs'  },
                             'jobCombo' : {'type' : bool,
                                           'description' :'Flag if data widget should be a combo box ' },
                             'mimeTypeName' : {'type' : str, 'description' : ''},
                             'mimeTypeDescription' : {'type' : str, 'description' : ''},
                             'fileLabel' : {'type' : str, 'description' : 'Label for file'},
                             'fileExtensions' : {'type' : list, 'listItemType' : str,
                                                 'description' : 'A list of strings containing allowed file extensions (no dot)'},
                             'fileContentClassName' : {'type' : str , 'editable' : False,
                                                       'description' : 'A string containing the name of a class which will hold the file contents' },
                             'isDirectory' : {'type' : bool,
                                              'description' : 'Flag if the data is a directory' },
                             'ifInfo' : {'type' : bool,
                                         'description' :'Flag if gui widget should have info icon'},
                             'saveToDb' : {'type' : bool,
                                           'description' :'Save the name of this file in the database' },
                             'requiredSubType' : {'type' : list, 'listItemType' : int,
                                                   'description' : 'A list of allowed sub types'},
                             'requiredContentFlag' : {'type' : list, 'listItemType' : int,
                                                      'description' : 'A list of allowed content flags'}}

    ERROR_CODES = {101 : {'description' : 'File does not exist'},
                   102 : {'description' : 'No mime type for data file'},
                   103 : {'description' : 'Attempting to set file content with inappropriate data'},
                   104 : {'description' : 'There is no file content class specified for this type of file'},
                   105 : {'description' : 'The file content class specified for this type of file can not be found'},
                   300 : {'description' : 'Passed' , 'severity' : SEVERITY_OK},
                   305 : {'description' : 'Neither original nor test file exists', 'severity' : SEVERITY_OK},
                   306 : {'description' : 'Original file does not exists'},
                   307 : {'description' : 'Test file does not exist '},
                   308 : {'description' : 'Files failed checksum comparison'},
                   309 : {'description' : 'Files failed size comparison'},
                   310 : {'description' : 'No comparison testing implemented for this file type', 'severity' : SEVERITY_WARNING},
                   311 : {'description' : 'Failed loading original file for comparison'},
                   312 : {'description' : 'Failed loading test file for comparison'},
                   313 : {'description' : 'Files failed simple text diff comparison'},
                   320 : {'description' : 'Unrecognised error attempting to load file'}}

    BLOCK_LOAD_FILES = False    # Might be set true for testing (eg CCP4DbApi)

    def __init__(self, value={}, qualifiers={}, parent=None, name=None, fullPath=None, checkDb=True, keywords={}, **kw):
        #print 'CDataFile.__init__',value,fullPath,kw
        qualis = {}
        qualis.update(qualifiers)
        if value is not None and fullPath is None and isinstance(value, str):
            fPath = value
            value = {}
        else:
            fPath = fullPath
        valu = {}
        valu .update(value)
        kwords = {}
        kwords.update(keywords)
        kwords.update(kw)
        if len(kwords) > 0:
            vs,qs,us = CCP4Data.sortArguments(self.__class__, kwords, aliases=['projectName'])
            if len(us) > 0:
                for key in list(us.keys()):
                    print(self.className(), 'input argument unrecognised:', key)
            qualis.update(qs)
            valu.update(vs)
        self.__dict__['_fileContent'] = None
        #self.__dict__['_dbFileId'] = None
        #print 'CDataFile.__init__',valu,fPath
        CCP4Data.CData.__init__(self, value=valu, qualifiers=qualis, parent=parent, name=name)
        for item in ['dbFileId','annotation']:
            if item not in self.__dict__['_value']:
                self.__dict__['_value'][item] = CCP4Data.CString()
        for item in ['contentFlag','subType']:
            if item not in self.__dict__['_value']:
                self.__dict__['_value'][item] = CCP4Data.CInt()
        #if not self.__dict__['_value']['annotation'].isSet(): self.__dict__['_value']['annotation'].set(self.qualifiers('mimeTypeDescription'))
        if fPath is not None:
            self.setFullPath(fPath, checkDb=checkDb)

    def buildItem(self, key=None, cls=None, qualifiers={}, subItem={}):
        # Reimplement to prevent the propagation of dataChanged signal
        #print 'CData.buildItem',key,qualifiers,subItem
        if len(subItem) > 0:
            self.__dict__['_value'][key] = cls(qualifiers=qualifiers, parent=self, name=key, subItem=subItem)
        else:
            self.__dict__['_value'][key] = cls(qualifiers=qualifiers, parent=self, name=key)

    def getExt(self):
        # Return the actual extension of file
        if not self.__dict__['_value']['baseName'].isSet():
            return None
        base, ext = os.path.splitext(self.__dict__['_value']['baseName'].get())
        return ext

    def getFullPath(self):
        path = self.makeFullPath(projectId=self.__dict__['_value']['project'].get(),
                                 relPath=self.__dict__['_value']['relPath'].get(),
                                 baseName=self.__dict__['_value']['baseName'].get())
        if len(path) > 0:
            return CFilePath(path, qualifiers=self.__dict__['_value']['relPath'].qualifiers(custom=True))
        else:
            return CFilePath()

    def makeFullPath(self, projectId=None, relPath=None, baseName=None, projectName=None):
        if projectId is None and projectName is not None:
            projectId = self.getProjectId(projectName=projectName)
        if projectId is None:
            path = ''
        elif projectId == 'workDirectory':
            #print 'CDataFile.makeFullPath project IS workDirectory'
            from core import CCP4PluginScript
            obj = self
            while isinstance(obj,CCP4Data.CData):
                obj = obj.parent()
                if isinstance(obj,CCP4PluginScript.CPluginScript):
                    path = obj.getWorkDirectory()
                    obj = None
        else:
            path = PROJECTSMANAGER().getProjectDirectory(projectId=projectId, testAlias=True)
            #print 'CDataFile.makeFullPath project path',projectId,type(projectId),path
            if path is None:
                path = ''
        if relPath is not None:
            path = os.path.join(path, relPath)
        if baseName is not None:
            path = os.path.join(path, baseName)
        if len(path) == 0:
            return ''
        else:
            return os.path.normpath(path)

    def setAnnotation(self, value=''):
        self.annotation.set(value)

    def getAnnotation(self):
        return self.annotation.get()

    def setFullPath(self, value='', checkDb=True):
        #print 'CDataFile.setFullPath', value
        if value is None or value == '':
            self.setDefault()
            self.updateData()
            return
        if isinstance(value, CFilePath):
            value = value.abspath()
            path = value.dirname()
            baseName = value.basename()
        else:
            value = os.path.abspath(value)
            path = os.path.dirname(value)
            baseName = os.path.basename(value)
        #print 'CDataFile.setFullPath',value,path,baseName,'project:',projectId,relPath,os.path.exists(str(value))
        if self.qualifiers('mustExist') and not os.path.exists(str(value)):
            raise CException(self.__class__, 101, 'Filename: ' + str(value), name=self.objectPath())
        if checkDb:
            projectName, relPath, projectId = PROJECTSMANAGER().interpretDirectory(path)
        else:
            projectId = None
        if projectId is not None:
            self.__dict__['_value']['project'].set(projectId)
            self.__dict__['_value']['relPath'].set(relPath)
        else:
            self.__dict__['_value']['project'].setDefault()
            self.__dict__['_value']['relPath'].set(path)
        self.__dict__['_value']['baseName'].set(baseName)
        self.__dict__['_value']['dbFileId'].unSet()
        self.__dict__['_value']['annotation'].unSet()
        self.updateData()

    def setDbFileId(self, fileId):
        self.unSet()
        try:
            modelist = ['filename', 'relpath', 'annotation', 'projectid', 'filesubtype', 'filecontent']
            fileInfo = PROJECTSMANAGER().db().getFileInfo(fileId=fileId, mode=modelist)
        except Exception as e:
            print('Error in CDataFile.setDbFileId', self.objectName())
            print(e)
            return
        self.__dict__['_value']['project'].set(fileInfo['projectid'])
        self.__dict__['_value']['relPath'].set(fileInfo['relpath'])
        self.__dict__['_value']['baseName'].set(fileInfo['filename'])
        self.__dict__['_value']['annotation'].set(fileInfo['annotation'])
        self.__dict__['_value']['subType'].set(fileInfo['filesubtype'])
        self.__dict__['_value']['contentFlag'].set(fileInfo['filecontent'])
        self.__dict__['_value']['dbFileId'].set(fileId)
        #print 'CDataFile.setDbFileId',self.objectName(),self.__dict__['_value']
        self.updateData()

    def stripedName(self):
        if self.baseName.get() is None:
            return None
        else:
            return os.path.splitext(self.baseName.__str__())[0]

    def isSet(self, allowUndefined=False, allowDefault=True, allSet=True):
        try:
            return self.__dict__['_value']['baseName'].isSet()
        except:
            #print 'Error in CDataFile.isSet()',self.__dict__['_value'].keys()
            return False

    def guiLabel(self, useAnnotation=True, useObjectName=True):
        #print 'CDataFile.guiLabel', self.objectName(), self, self.annotation, self.dbFileId
        from core import CCP4TaskManager
        if self.isSet():
            if useAnnotation and self.__dict__['_value']['annotation'].isSet():
                return self.__dict__['_value']['annotation'].__str__()
            if self.__dict__['_value']['dbFileId'].isSet():
                fileInfo = PROJECTSMANAGER().db().getFileInfo(fileId=str(self.__dict__['_value']['dbFileId']), mode=['jobnumber', 'taskname'])
                #print 'CDataFile.guiLabel fileInfo',fileInfo
                gLabel = self.qualifiers('mimeTypeDescription') + ' from ' + str(fileInfo['jobnumber']) + ' ' \
                          + CCP4TaskManager.TASKMANAGER().getTitle(fileInfo['taskname'])
                return gLabel
            else:
                return os.path.splitext(self.baseName.__str__())[0]
        elif useObjectName:
            return str(self.objectName())
        else:
            return ''

    def getProjectId(self, projectName):
        try:
            return PROJECTSMANAGER().db().getProjectId(projectName=projectName)
        except:
            return None

    def set(self, value={}, **kw):
        #print 'CDataFile.set',value,kw
        if isinstance(value, str):
            #print 'CDataFile.set from string',value
            self.setFullPath(value)
        elif isinstance(value, CDataFile):
            #print 'CDataFile.set from CDataFile',value.get()
            CCP4Data.CData.set(self, value.get())
        elif isinstance(value, CFilePath):
            self.setFullPath(str(value))
        elif isinstance(value, dict):
            if 'fullPath' in value:
                self.setFullPath(value['fullPath'])
            else:
                if value.get('project', None) is None and kw.get('project', None) is None:
                    if kw.get('projectName', None) is not None:
                        kw['project'] = self.getProjectId(kw['projectName'])
                        del kw['projectName']
                    elif value.get('projectName', None) is not None:
                        value['project'] = self.getProjectId(value['projectName'])
                        del value['projectName']
                #print 'CDataFile.set updated', value, kw
                valu = {'dbFileId' : None, 'project' : None, 'baseName' :None, 'relPath' :None, 'annotation': None}
                valu.update(value)
                valu.update(kw)
            CCP4Data.CData.set(self, valu)
        self.unsetFileContent()

    def unSet(self):
        CCP4Data.CData.unSet(self)
        for key in ['sourceFileName', 'sourceFileAnnotation']:
            if key in self.__dict__:
                del self.__dict__[key]
        self.unsetFileContent()

    def validity(self, arg={}):
        #v = self.itemValidity(arg)
        #print 'CDataFile.validity',self.objectPath(),arg,type(arg)
        v = CErrorReport()
        if arg is None:
            if self.qualifiers('allowUndefined'):
                v.append(self.__class__, 1, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            else:
                v.append(self.__class__, 2, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            return v
        else:
            if isinstance(arg, (CDataFile, CFilePath)):
                fullPath = arg.__str__()
            elif isinstance(arg, str):
                fullPath = arg
            elif isinstance(arg, dict):
                if isinstance(arg.get('baseName'), CFilePath):
                    fullPath = self.makeFullPath(projectId=arg['project'].get(), relPath=arg['relPath'].get(), baseName=arg['baseName'].get()).__str__()
                else:
                    fullPath = self.makeFullPath(projectId=arg.get('project',None), relPath=arg.get('relPath',None), baseName=arg.get('baseName',None)).__str__()
            #print 'CDataFile.validity fullPath',self.objectPath(),fullPath
            if fullPath == '':
                if self.qualifiers('allowUndefined'):
                    v.append(self.__class__, 1, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
                else:
                    v.append(self.__class__, 2, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
                #print 'CDataFile.validity fullPath',self.objectName(),arg,type(arg),fullPath,type(fullPath),self.qualifiers('allowUndefined'),v.maxSeverity()
                return v
            if self.qualifiers('mustExist') and not os.path.exists(str(fullPath)):
                v.append(self.__class__, 101, name=self.objectPath(),details='Filename: '+str(fullPath), label=self.qualifiers('guiLabel'), stack=False)
        #print 'CDataFile.validity',len(v),v.maxSeverity()
        return v

    def abspath(self):
        p = self.fullPath
        return p.abspath()

    def exists(self):
        p = self.fullPath
        #print 'CDataFile.exists',self.__dict__['_value']['project'],p
        return p.exists()

    def dirname(self):
        return self.project.directory() + self.relPath

    def samefile(self,arg):
        return self.fullPath.samefile(self.getValue(arg))

    def getValue(self,arg):
        if isinstance(arg, CDataFile):
            return arg.fullPath
        else:
            return arg

    def fileExtensions(self):
        return self.qualifiers('fileExtensions')

    def __str__(self):
        return self.getFullPath().__str__()

    def getTextItem(self):
        #print 'CDataFile.getTextItem',self.annotation.__str__(),self.fullPath.__str__()
        if self.annotation.isSet():
            return self.annotation.__str__()
        else:
            return self.fullPath.__str__()

    def loadFile(self, initialise=False):
        rv = CErrorReport()
        #print 'CDataFile.loadFile',str(self),self.BLOCK_LOAD_FILES,self.fileContentClass()
        if self.BLOCK_LOAD_FILES:
            return
        if self.qualifiers('fileContentClassName') is None:
            return
        if self.__dict__['_fileContent'] is None or initialise:
            self.__dict__['_fileContent'] = self.fileContentClass()(parent=self, name='fileContent')
        #print 'CDataFile.loadFile', self.fullPath, self.fullPath.exists(), self.fullPath.isfile()
        if self.fullPath.exists() and self.fullPath.isfile():
            try:
                #print 'CDataFile.loadFile _fileContent', type(self.__dict__['_fileContent'])
                self.__dict__['_fileContent'].loadFile(self.fullPath)
            except CException as e:
                rv.extend(e)
            except Exception as e:
                rv.append(self.__class__, 320, self.__str__() + '\n\n' + str(e), exc_info=sys.exc_info())
        else:
            self.__dict__['_fileContent'].unSet()
            self.__dict__['_lastLoadedFile'] = None
        self.dataLoaded.emit()
        return rv

    def fileContentClass(self):
        className = self.qualifiers('fileContentClassName')
        if className is None:
            raise CException(self.__class__, 104, name=self.objectPath())
        from core import CCP4DataManager
        cls = CCP4DataManager.DATAMANAGER().getClass(className)
        if cls is None:
            raise CException(self.__class__, 105, 'Contents class name:' + str(className), name=self.objectPath())
        return cls

    def getFileContent(self):
        #print 'CDataFile.getFileContent',self.objectName()
        if self.__dict__['_fileContent'] is None:
            self.loadFile()
        return self.__dict__['_fileContent']

    def setFileContent(self, fileContent=None):
        if fileContent is None or not isinstance(fileContent, self.fileContentClass()):
            raise CException(self.__class__, 103, name=self.objectPath())
        else:
            self.__dict__['_fileContent'] = fileContent

    PROPERTIES = {'fullPath' : {'fget' : getFullPath, 'fset' : setFullPath},
                  'fileContent' : {'fget' : getFileContent, 'fset' : setFileContent},
                  'annotation' : {'fget' : getAnnotation, 'fset' : setAnnotation}
                  }

    def unsetFileContent(self):
        self.__dict__['_fileContent'] = None

    def requiredContent(self):
        return None

    def setContentFlag(self, reset=False):
        # Set the contentFlag - mostly relevant and reimplemented CMiniMtzDataFile
        return None

    def getTableTextItems(self):
        return [self.__str__()]

    def saveToDb(self):
        ifSave = self.qualifiers('saveToDb')
        if ifSave is not NotImplemented and ifSave:
            return [self], None, {}
        else:
            return [], None, {}

    def setOutputPath(self, jobName='', projectId=None, relPath=None):
        name = self.objectName().__str__()
        if len(name) == 0:
            try:
                name = self.parent().objectName().__str__()
            except:
                pass
        baseName = jobName + name
        fileLabel = self.qualifiers('label')
        if fileLabel is not None and fileLabel is not NotImplemented:
            baseName += CDataFile.SEPARATOR + fileLabel
        baseName += '.' + self.fileExtensions()[0]
        self.set({'project' : projectId, 'relPath' : relPath, 'baseName' : baseName } )

    def checksum(self, blockSize=256*128, ifHex=True, filePath = None):
        '''From http://stackoverflow.com/questions/1131220/get-md5-hash-of-big-files-in-python
        Block size directly depends on the block size of your filesystem
        to avoid performances issues
        Here I have blocks of 4096 octets (Default NTFS)'''
        if filePath is None:
            try:
                filePath = self.__str__()
                #print 'CDataFile.checksum path', path, os.path.exists(path)
            except:
                return None
        if not os.path.exists(filePath):
            return None
        import hashlib
        md5 = hashlib.md5()
        with open(filePath,'rb') as f:
            for chunk in iter(lambda: f.read(blockSize), b''):
                md5.update(chunk)
        if ifHex:
            return md5.hexdigest()
        else:
            return md5.digest()

    def assertSame(self, arg, testPath=False, testChecksum=True, testSize=False, testDiff=False, diagnostic=False, fileName=None):
        from core import CCP4Utils
        if isinstance(arg, self.__class__):
            other = arg.__str__()
        else:
            other = arg
        # Expect to use file with name given by self.__str__() but can input an alternate fileName for example
        # if have create temp file with changed end-of-line
        if fileName is None:
            fileName=self.__str__()
        name = self.objectPath(False)
        if not self.isSet():
            if len(other) == 0:
                return CErrorReport(self.__class__, 315, name=name)
            else:
                return CErrorReport(self.__class__, 303, name=name, details=str(self) + ' : ' + str(other))
        elif len(other) ==0:
            return CErrorReport(self.__class__, 302, name=name, details=str(self) + ' : ' + str(other))
        report = CErrorReport()
        if testPath:
            if other != self.__str__():
                report.append(self.__class__, 304, name=name, details=str(self) + ' : ' + str(other))
        if testChecksum:
            # Do a checksum comparison
            otherSum = self.checksum(filePath=other)
            selfSum = self.checksum(filePath=fileName)
            if selfSum is None:
                if otherSum is None:
                    report.append(self.__class__, 305, name=name, details=str(self) + ' : ' + str(other))
                else:
                    report.append(self.__class__, 306, name=name, details=str(self) + ' : ' + str(other))
            elif otherSum is None:
                report.append(self.__class__, 307, name=name, details=str(self) + ' : ' + str(other))
            if otherSum != selfSum:
                report.append(self.__class__, 308, name=name, details=str(self) + ' : ' + str(other))
        if testSize:
            try:
                selfSize = os.path.getsize(fileName)
            except:
                selfSize = None
            try:
                otherSize = os.path.getsize(str(other))
            except:
                otherSize = None
            if selfSize is None:
                if otherSize is None:
                    report.append(self.__class__, 305, name=name, details=str(self) + ' : ' + str(other))
                else:
                    report.append(self.__class__, 306, name=name, details=str(self) + ' : ' + str(other))
            elif otherSize is None:
                report.append(self.__class__, 307, name=self.objectName(), details=str(self) + ' : ' + str(other))
            diff = (abs(float(selfSize - otherSize)))/min(selfSize, otherSize)
            #print 'CDataFile.assertSame testSize diff',selfSize,otherSize,diff
            if diff > 0.01:
                report.append(self.__class__, 309, name=name, details='by ' + str(int(100*diff)) + '% : ' + str(self) + ' : ' + str(other))
        if testDiff:
            # Unsophisicated diff
            if sys.version_info > (3,0):
                from googlecode import diff_match_patch_py3
                dmp =  diff_match_patch_py3.diff_match_patch()
            else:
                from googlecode import diff_match_patch
                dmp =  diff_match_patch.diff_match_patch()
            diffs = dmp.diff_main(CCP4Utils.readFile(fileName), CCP4Utils.readFile(other))
            #print 'CDataFile.assertSame diffs', diffs
            if len(diffs) > 1:
                report.append(self.__class__, 313, name=self.objectPath(), details=str(self) + ' : ' + str(other))
        if len(report) == 0:
            report.append(self.__class__, 300, name=name)
        return report

    def importFileName(self, jobId=None, jobDirectory=None, ext=None):
        if ext is None:
            ext = '.' + self.qualifiers('fileExtensions')[0]
        #print 'CDataFile.importFileName jobId', jobId
        if jobDirectory is None:
            jobDirectory = PROJECTSMANAGER().jobDirectory(jobId=jobId)
        objName = self.objectName()
        if len(objName) == 0:
            try:
                objName = self.parent().objectName()
                if isinstance(self.parent(), CCP4Data.CList):
                    try:
                        idx = self.parent().index(self)
                        if idx >= 0:
                            objName = objName + '_' + str(idx)
                    except:
                        pass
            except:
                pass
        # MN Here we try to avoid overwriting if a single task imports multiple instances of a particular file type
        filename = os.path.normpath(os.path.join(jobDirectory, objName + '-' + self.qualifiers('fileLabel') + ext))
        iExtra = 0
        while os.path.isfile(filename):
            filename =  os.path.normpath(os.path.join(jobDirectory, objName + '-' + self.qualifiers('fileLabel') + '_' + str(iExtra) + ext))
            iExtra += 1
        #print 'CMiniMtzDataFile.importFileName', filename
        return filename

    def importFile(self, jobId=None, sourceFileName=None, ext=None, annotation=None, validatedFile=None, jobNumber=None):
        import shutil
        if sourceFileName is None:
            sourceFileName = self.__str__()
        if ext is None:
            ext=os.path.splitext(self.__str__())[1]
        filename = self.importFileName(jobId=jobId, ext=ext)
        #print 'CDataFile.importFile copy', sourceFileName, filename
        if validatedFile is not None and os.path.exists(validatedFile):
            shutil.copyfile(validatedFile, filename)
            self.setFullPath(filename)
        else:
            shutil.copyfile(sourceFileName, filename)
            self.setFullPath(filename)
        if not self.annotation.isSet():
            if annotation is None:
                annotation = self.qualifiers('guiLabel') + ' imported from ' + os.path.split(sourceFileName)[1]
                if jobNumber is not None and not annotation.count(' by job '):
                    annotation = annotation + ' by job ' + str(jobNumber)
            self.annotation.set(annotation)
        self.setContentFlag(reset=True)
        #print 'CDataFile.importFile setContentFlag',self.contentFlag
        # Set the hidden sourceFileName to be picked up by PROJECTMANAGER.importFiles() at run time
        # - the original file is copied to CCP4_IMPORTED_FILES when user has confirmed the selection by clicking 'Run'
        self.__dict__['sourceFileName'] = sourceFileName

    def isDosFile(self):
        from core import CCP4Utils
        text = CCP4Utils.readFile(self.__str__())
        isDFile = (text.find('\r\n') >= 0)
        return isDFile

    def resetLineEnd(self, outputFilename=None, toDos=False):
        from core import CCP4Utils
        text = CCP4Utils.readFile(self.__str__())
        isDos = text.find('\r\n') >= 0
        if isDos == toDos:
            if outputFilename is not None:
                CCP4Utils.saveFile(text=text, fileName=outputFilename)
                return outputFilename
            else:
                return self.__str__()
        if outputFilename is None:
            import tempfile
            outputFilename = tempfile.mktemp(suffix='.txt')
        if toDos:
            textOut = re.sub(r'\n', '\r\n', text)
        else:
            textOut = re.sub(r'\r\n', '\n', text)
        CCP4Utils.saveFile(fileName=outputFilename, text=textOut)
        return outputFilename


# The header info in a file
class CFileFunction(CCP4Data.CString):
    '''List of recognised XML file functions'''

    QUALIFIERS = {'enumerators' : ['DEF', 'PARAMS', 'LOG', 'PROJECTDIRECTORIES', 'COM', 'REFMAC', 'OUTPUT',
                                  'STATUS', 'PROJECTDATABASE', 'MGSCENE', 'JOBSERVERSTATUS', 'WORKFLOW',
                                  'COMFILEPATCH', 'CUSTOMTASK', 'IMPORTEDJOB', 'I1SUPPLEMENT',
                                  'ASUCONTENT','UNKNOWN'], 'onlyEnumerators' : True}
    QUALIFIERS_DEFINITION = {'enumerators' : {'type' : list},
                            'onlyEnumerators' : {'type' : bool, 'editable' : False}}


class CI2XmlHeader(CCP4Data.CData):
    '''Container for header info from XML file'''
    from core import CCP4Annotation

    CONTENTS = {'function' : {'class' : CFileFunction}, 'userId' : {'class' : CCP4Annotation.CUserId},
                'hostName' : {'class' : CCP4Annotation.CHostName}, 'creationTime' : {'class' : CCP4Annotation.CTime},
                'ccp4iVersion' : {'class' : CVersion}, 'pluginName' : {'class' : CCP4Data.CString},
                'pluginVersion' : {'class' : CVersion}, 'pluginTitle' : {'class' : CCP4Data.CString},
                'projectName' : {'class' : CProjectName},
                'projectId' : {'class' : CProjectId, 'qualifiers': {'allowUnfound' : True}},
                'jobId' : {'class' : CCP4Data.CUUID}, 'jobNumber' : {'class' : CCP4Data.CString},
                'comment' : {'class' : CCP4Data.CString}, 'OS' : {'class' : CCP4Data.CString}}

    ERROR_CODES = {101 : {'description' : 'Attempting to read header from non-existant Xml file'},
                   102 : {'description' : 'Error loading file to read header'},
                   103 : {'description' : 'Error finding <ccp4i2_header> in file'},
                   104 : {'description' : 'Error interpreting header from file'},
                   105 : {'description' : 'File does not have <ccp4i2> root node'}}

    def setCurrent(self):
        self.userId.setCurrentUser()
        self.creationTime.setCurrentTime()
        self.hostName.setCurrentHost()
        self.ccp4iVersion.setCurrent()
        if sys.platform == 'darwin':
            self.OS = 'MacOSX'
        elif sys.platform[0:5] == 'linux':
            self.OS = 'Linux'
        elif sys.platform[0:3] == 'win':
            self.OS = 'Windows'

    def loadFromXml(self, fileName=None, checkValidity=True):
        from core import CCP4Utils
        if fileName is None:
            raise CException(self.__class__, 101, name = self.objectPath())
        if not os.path.exists(fileName):
            raise CException(self.__class__, 101, 'Filename: ' + fileName, name = self.objectPath())
        root = None
        try:
            root = CCP4Utils.openFileToEtree(fileName)
        except:
            pass
        if root is None:
            raise CException(self.__class__, 102, 'Filename: ' + fileName, name=self.objectPath())
        if not root.tag in [CI2XmlDataFile.ROOT_TAG, '{'+CCP4NS+'}' + CI2XmlDataFile.ROOT_TAG]:
            raise CException(self.__class__, 105, 'Filename: ' + fileName, name=self.objectPath())
        try:
            headerEtree = root.find('ccp4i2_header')
        except:
            headerEtree = None
        if headerEtree is None:
            raise CException(self.__class__, 103, 'Filename: ' + fileName, name=self.objectPath())
        #try:
        rv = self.setEtree(headerEtree, checkValidity = checkValidity)
        #except:
        #  raise CException(self.__class__,104,'Filename: '+fileName)
        #print 'CI2XmlHeader.loadFromXml',rv.report(),str(self)
        return rv


class CDataFileContent(CCP4Data.CData):
    '''Base class for classes holding file contents'''
    pass


class CTextDataFile(CDataFile):
    '''A text data file'''

    QUALIFIERS = {'mimeTypeName' : '"text/plain"', 'mimeTypeDescription' : 'Standard plain text',
                  'fileLabel' : None,'fileExtensions' : ['txt','log']}

class CMmcifData(CDataFileContent):
    '''Generic mmCIF data.
     This is intended to be a base class for other classes
     specific to coordinates, reflections or geometry data.'''
    pass


class CMmcifDataFile(CDataFile):
    '''A generic mmCIF format file.
     This is intended to be a base class for other classes
     specific to coordinates, reflections or geometry data.'''

    QUALIFIERS = { 'fileExtensions' : ['cif','ent'] }


class CXmgrDataFile(CDataFile):
    '''An xmgr format file. This is the input format for xmgrace, as output by scala or aimless'''
    QUALIFIERS = {'mimeTypeName' : 'application/grace', 'fileExtensions' : ['xmgr'] }


class CPostscriptDataFile(CDataFile):
    '''A postscript format file'''
    QUALIFIERS = {'mimeTypeName' : 'application/postscript','fileExtensions' : ['ps'], 'guiLabel' : 'Postscript file' }

class CPDFDataFile(CDataFile):
    '''An PDF format file'''
    QUALIFIERS = {'mimeTypeName' : 'application/x-pdf','fileExtensions' : ['pdf'], 'guiLabel' : 'PDF file' }

class CSceneDataFile(CDataFile):
    '''An xml format file for defining scene in CCP4mg.'''
    QUALIFIERS = {'fileLabel' : 'scene', 'mimeTypeName' : 'application/CCP4-scene',
                  'mimeTypeDescription' : 'CCP4mg scene file', 'guiLabel' : 'CCP4mg scene',
                  'fileExtensions' : ['scene.xml'], 'fileContentClassName' : NotImplemented }


class CXmlDataFile(CDataFile):
    '''A reference to an XML file'''
    QUALIFIERS = {'fileExtensions' : ['xml'], 'saveToDb' : False, 'mimeTypeName' : 'application/xml' }
    ERROR_CODES = {1001 : {'description' : 'Unknown error reading XML file'},
                   1002 : {'description' : 'Error trying to find root node in XML'},
                   1006 : {'description' : 'Attempting to save XML file with incorrect body'},
                   1007 : {'description' : 'Error creating XML text'},
                   1008 : {'description' : 'Error saving XML text to file'},
                   1009 : {'description' : 'Error reading XML file'},
                   1010 : {'description' : 'XML file does not exist'},
                   1011 : {'description' : 'No file name given for making I2XMlDataFile'},
                   1012 : {'description' : 'Error creating I2XMlDataFile object'},
                   1013 : {'description' : 'Error creating I2XMlDataFile file'}}

    def loadFile(self, printout=False):
        if CCP4Config.XMLPARSER() == 'lxml':
            from lxml import etree
        if not self.fullPath.exists():
            raise CException(self.__class__, 1010, 'Filename: ' + str(self.fullPath), name=self.objectPath())
        else:
            fileName = self.fullPath.get()
        root = self.getEtreeRoot(fileName)
        if printout:
            print(ET.tostring(root)) # KJS - problem here. etree not defined. Fix in.
        return root
        
        
    def assertSame(self, arg, testPath=False, testChecksum=True, testSize=False, testDiff=False, diagnostic=False, fileName=None):
        #MN Now here we have an issue that assertSame on an XML file is fraught with difficulties, and an identical checkSum is probably far too stringent.  I'm gonna suggest that we should remove testChecksum for now, with a view to putting in a more intelligent comparison later
        report = CDataFile.assertSame(self, arg, testPath, False, testSize, testDiff, diagnostic, fileName)
        return report


    def getEtreeRoot(self,fileName=None,useLXML=True):
        from core import CCP4Utils
        if fileName is None:
            fileName = self.fullPath.get()
        if CCP4Config.XMLPARSER() == 'lxml':
            from lxml import etree
        try:
            tree = CCP4Utils.openFileToEtree(fileName,useLXML=useLXML)
        except ET.LxmlError as e:
            raise CException(self.__class__, 1009, fileName + ' : ' + str(e), name=self.objectPath())
        except Exception as e:
            raise CException(self.__class__, 1001, fileName + ' : ' + str(e), name=self.objectPath())
        return tree

    '''
    Think this is now redundant as openFileToEtree does it
    try:
      root = tree.getroot()
    except Exception as e:
      raise CException(self.__class__,1002,fileName+' : '+str(e))
    return root
    '''

    def saveFile(self, bodyEtree=None):
        from core import CCP4Utils
        if CCP4Config.XMLPARSER() == 'lxml':
            from lxml import etree
        fileName = self.fullPath.get()
        print('bodyEtree is', bodyEtree)
        if bodyEtree is None:
            pass
        elif not isinstance(bodyEtree, ET.Element):
            import traceback
            traceback.print_stack()
            raise CException(self.__class__, 1006, fileName, name=self.objectPath())
        if CCP4Config.DEVELOPER():
            ET.indent(bodyEtree)
            text = ET.tostring(bodyEtree, xml_declaration=True)
            CCP4Utils.saveFile(fileName=fileName, text=text, overwrite=1)
        else:
            try:
                text = ET.tostring(bodyEtree, xml_declaration=True) # KJS - Code is broken here. Fix in.
            except:
                raise CException(self.__module__, 1007, fileName)
            try:
                CCP4Utils.saveFile(fileName=fileName, text=text, overwrite=1)
            except:
                raise CException(self.__module__, 1008, fileName)

    def makeI2XmlDataFile(self, fileName=None, overWrite=False, header=None, **kw):
        if not self.fullPath.exists():
            raise CException(self.__class__, 1010, 'Filename: ' + str(self.fullPath), name=self.objectPath())
        else:
            myFileName = self.fullPath.get()
        if fileName is None and not overWrite:
            raise CException(self.__class__, 1011, name=self.objectPath())
        if fileName is None:
            fileName = str(self.getFullPath())
        if isinstance(fileName,(CDataFile, CFilePath)):
            fileName = str(fileName)
        try:
            c = CI2XmlDataFile(fileName)
        except:
            raise CException(self.__class__, 1012, fileName, name=self.objectPath())
        if header is not None:
            c.header.set(header)
        # Expect kw to be header data
        if len(kw) > 0:
            c.header.set(kw)
        if not c.header.userId.isSet():
            c.header.userId.setCurrentUser()
        if not c.header.creationTime.isSet():
            c.header.creationTime.setCurrentTime()
        if CCP4Config.XMLPARSER() == 'lxml':
            from lxml import etree
        body = ET.Element(CI2XmlDataFile.BODY_TAG)
        body.append(self.getEtreeRoot(myFileName))
        c.saveFile(bodyEtree=body)

class CEBIValidationXMLDataFile(CXmlDataFile):
    '''An XLM file returned from the EBI validation server '''
    QUALIFIERS = {'mimeTypeName' : 'application/EBI-validation-xml','fileExtensions' : ['xml'], 'guiLabel' : 'EBI Validation XML' }

class CI2XmlDataFile(CXmlDataFile):
    '''A reference to an XML file with CCP4i2 Header'''
    ROOT_TAG = 'ccp4i2'
    HEADER_TAG = 'ccp4i2_header'
    BODY_TAG = 'ccp4i2_body'
    CONTENTS = {'project' : {'class' : CProjectId}, 'baseName' : {'class' : CFilePath},
                'relPath' : {'class' : CFilePath}, 'header' : {'class' : CI2XmlHeader}}
    QUALIFIERS = {'fileExtensions' : ['xml'], 'autoLoadHeader' : True}
    QUALIFIERS_ORDER = [ 'autoLoadHeader']
    QUALIFIERS_DEFINITION = {'autoLoadHeader' : { 'type' : bool}}
    ERROR_CODES = {1003 : {'description' : 'XML does not have <ccp4i2> root node'},
                   1004 : {'severity' : SEVERITY_WARNING, 'description' : 'XML does not have <ccp4i2_header> section'},
                   1005 : {'description' : 'XML does not have <ccp4i2_body> section'}}

    def __init__(self, value={}, qualifiers={}, parent=None, name=None, fullPath=None, keywords={}, **kw):
        kw['checkDb'] = False
        CXmlDataFile.__init__(self, value=value, qualifiers=qualifiers, parent=parent, name=name, fullPath=fullPath, keywords=keywords, **kw)

    def updateData(self):
        if self.qualifiers('autoLoadHeader'):
            if self.getFullPath().exists():
                self.loadFile()

    def loadFile(self, printout=False):
        if not self.fullPath.exists():
            raise CException(self.__class__, 1010, 'Filename: ' + str(self.fullPath), name=self.objectPath())
        else:
            fileName = self.fullPath.get()
        root = self.getEtreeRoot(fileName)
        if printout:
            if CCP4Config.XMLPARSER() == 'lxml':
                ET.indent(root)
                print(ET.tostring(root))
        return self.loadHeader(root)

    def getEtreeRoot(self, fileName=None,useLXML=True):
        from core import CCP4Utils
        if CCP4Config.XMLPARSER() == 'lxml':
            from lxml import etree
        if fileName is None:
            fileName = self.fullPath.get()
        #print 'getEtreeRoot fileName', fileName, type(fileName)
        try:
            root = CCP4Utils.openFileToEtree(fileName,useLXML=useLXML)
        except ET.LxmlError as e:
            raise CException(self.__class__, 1009, fileName + ' : ' + str(e), name=self.objectPath())
        except Exception as e:
            raise CException(self.__class__, 1001, fileName + ' : ' + str(e), name=self.objectPath())
        '''
        try:
          root = tree.getroot()
        except Exception as e:
          raise CException(self.__class__, 1002, fileName+' : ' + str(e))
        '''
        if not root.tag in [CI2XmlDataFile.ROOT_TAG, '{' + CCP4NS + '}' + CI2XmlDataFile.ROOT_TAG]:
            raise CException(self.__class__, 1003, fileName, name=self.objectPath())
        return root

    def loadHeader(self, root=None):
        err = CErrorReport()
        root = self.getEtreeRoot()
        header_etree = root.find(CI2XmlDataFile.HEADER_TAG)
        if header_etree is None:
            header_etree = root.find('{' + CCP4NS + '}' + CI2XmlDataFile.HEADER_TAG)
        if header_etree is None:
            err.append(self.__class__, 1004, self.fullPath.get(), name = self.objectPath())
            self.__dict__['_value']['header'].setDefault()
        else:
            #try:
            err.extend(self.__dict__['_value']['header'].setEtree(header_etree))
            #except:
            #  print 'Error loading header'
            #  from lxml import etree
            #  print ET.tostring(header_etree)
        return err

    def getBodyEtree(self):
        if not self.fullPath.exists():
            raise CException(self.__class__, 101, self.fullPath.get(), name=self.objectPath())
        else:
            fileName = self.fullPath.get()
        root = self.getEtreeRoot(fileName)
        body_etree = root.find(CI2XmlDataFile.BODY_TAG)
        if body_etree is None:
            raise CException(self.__class__, 1005, fileName, name=self.objectPath())
        return body_etree

    def saveFile(self, bodyEtree=None,useLXML=False):
        from core import CCP4Utils
        if CCP4Config.XMLPARSER() == 'lxml' and useLXML:
            from lxml import etree
            testType = etree._Element
        else:
            from xml.etree import ElementTree as ET
            testType = ET.Element
        fileName = self.getFullPath()
        if bodyEtree is None:
            pass
        elif not isinstance(bodyEtree,testType):
            raise CException(self.__class__, 1006, fileName, name=self.objectPath())
        else:
            if bodyEtree.tag != CI2XmlDataFile.BODY_TAG:
                bodyEtree.tag = CI2XmlDataFile.BODY_TAG
        doc = ET.parse(StringIO( '<ccp4:ccp4i2 xmlns:ccp4="' + CCP4NS + '"></ccp4:ccp4i2>'))
        root = doc.getroot()
        headerEtree = self.__dict__['_value']['header'].getEtree(useLXML=useLXML)
        headerEtree.tag = CI2XmlDataFile.HEADER_TAG
        root.append(headerEtree)
        if bodyEtree is not None:
            if len(root.findall('./' + CI2XmlDataFile.BODY_TAG)) > 0:
                root.remove(root.findall('./' + CI2XmlDataFile.BODY_TAG))
            root.append(bodyEtree)
        if CCP4Config.DEVELOPER():
            if useLXML:
                text = etree.tostring(doc, xml_declaration=True)
                CCP4Utils.saveFile(fileName=fileName, text=text,overwrite=1)
            else:
                doc.write(os.path.normpath(str(fileName)), encoding = "ASCII", xml_declaration = True)  
        else:
            if useLXML:
                try:
                    ET.indent(doc)
                    text = ET.tostring(doc, xml_declaration=True)
                except:
                    raise CException(self.__module__, 1007, fileName)
                try:
                    CCP4Utils.saveFile(fileName=fileName, text=text, overwrite=1)
                except:
                    raise CException(self.__module__, 1008, fileName)
            else:
                try:
                    doc.write(os.path.normpath(str(fileName)), encoding = "ASCII", xml_declaration = True)  
                except:
                    raise CException(self.__module__, 1008, fileName)

def compareXmlFiles(xmlFile1, xmlFile2):
    from core import CCP4Utils
    tree1 = CCP4Utils.openFileToEtree(xmlFile1)
    tree2 = CCP4Utils.openFileToEtree(xmlFile2)
    return compareEtreeNodes(tree1, tree2)

def compareEtreeNodes(node1, node2):
    errList = []
    if node1.tag != node2.tag:
        errList.append([getNodePath(node1), 'tag', node1.tag, node2.tag])
    if node1.text != node2.text:
        errList.append([getNodePath(node1), 'text', node1.text, node2.text])
    for key, value1 in list(node1.attrib.items()):
        if key not in node2.attrib:
            errList.append([getNodePath(node2), 'missingAttrib2', key, None])   # KJS - tree1 not defined, change to node1
        if node2.attrib.get(key) != value1:
            errList.append([getNodePath(node1, key), 'diffAttrib', value1, node2.attrib.get(key)])
    for key in list(node2.attrib.keys()):
        if key not in node1.attrib:
            errList.append([getNodePath(node1), 'missingAttrib1', None, key])
    for child1 in node1.getchildren():
        hits = getEquivalentNode(child1, node2)
        if len(hits) == 1:
            errList.extend(compareEtreeNodes(child1, hits[0]))
        else:
            errList.append([getNodePath(child1), 'missingChild2', str(child1.tag), None] )
    for child2 in node2.getchildren():
        hits = getEquivalentNode(child2, node1)
        if len(hits) != 1:
            errList.append([getNodePath(child2), 'missingChild1', None, str(child2.tag)])
    return errList

DBIDSFORUNIQUENODE = {'job' : ['jobid'], 'file' : ['fileid'], 'importfile' : ['importid'],
                      'fileuse' : ['fileid','jobid'], 'exportfile' : ['exportid'],
                      'comment' : ['commentid' ], 'xdata' : ['xdataid'],
                      'projectexport' : ['projectexportid'], 'projectimport' : ['projectimportid']}

def getEquivalentNode(node, otherParent):
    idItems = DBIDSFORUNIQUENODE.get(node.tag, [])
    path = './' + str(node.tag)
    if len(idItems) > 0:
        path = path + "[@" + idItems[0]+"='" + node.get(idItems[0]) + "'"
        if len(idItems)>1:
            path = path + " and @"+idItems[1]+"='" + node.get(idItems[1]) + "'"
        path = path + ']'
    #print path
    return otherParent.findall(path)

def getNodePath(node, key=None):
    name = getUniqueNodeLabel(node)
    if key is not None:
        name = name + '.' + key
    return name
    '''
      path.append(getUniqueNodeLabel(node))
      for p in node.iterancestors():
        path.insert(0,getUniqueNodeLabel(p))
      return path
    '''

def getUniqueNodeLabel(node):
    name = re.sub('{http://www.ccp4.ac.uk/ccp4ns}', '', str(node.tag))
    idItems = DBIDSFORUNIQUENODE.get(node.tag, [])
    for item in idItems:
        name = name + ':' + node.get(item)
    return name

def printCompareEtreeNodesErrors(errList, xmlFile1=None, xmlFile2=None):
    for path, errType, value1, value2 in errList:
        print('{:15}{}'.format(errType, path))

def xmlFileHeader(fileName):
    h = CI2XmlHeader()
    try:
        h.loadFromXml(fileName)
    except CException as e:
        raise e
    except Exception as e:
        raise e # CException(self.__class__, 115, fileName, name=self.objectPath()) # KJS - A global function does not have a class instance ....
    return h

def cloneI2XmlFile(sourceFile, targetFile, header={}, current=True, taskFrame=None, taskName=None, suggestedParams=None):
    import getpass
    from core import CCP4ModelData, CCP4TaskManager, CCP4Utils
    xFile = CI2XmlDataFile(sourceFile)
    if current:
        xFile.header.setCurrent()
    xFile.header.set(header)
    #This stops the header being overwritten with that from targetFile
    xFile.setQualifiers({'autoLoadHeader' : False})
    body = xFile.getBodyEtree()
    #FIXME - And then there are the lists. Not considered yet.
    #This sees if we have a sequence in old task which has become ASU in new version of task and so creates and runs a ProvideAsuContents.
    if len(body.findall("inputData"))>0:
        if len(body.findall("inputData/SEQIN")) > 0 or len(body.findall("inputData/SEQUENCE_LIST")) > 0 or len(body.findall("inputData/AWA_SEQIN"))>0:
            if len(body.findall("inputData/SEQIN")) > 0:
                oldTag = "SEQIN"
            elif len(body.findall("inputData/SEQUENCE_LIST")) > 0:
                oldTag = "SEQUENCE_LIST"
            elif len(body.findall("inputData/AWA_SEQIN")) > 0:
                oldTag = "AWA_SEQIN"
            #OK, old task had a SEQIN
            print("--------------------------------------------------")
            print("OK, old task had a ",oldTag)
            print("Now I need to know the def.xml for",taskName)
            defFile = CCP4TaskManager.TASKMANAGER().lookupDefFile(taskName)
            if defFile is None:
                defFile = CCP4TaskManager.TASKMANAGER().searchDefFile(taskName)
                if defFile is None:
                    print("I did not find a def file")
            print("Found def file",defFile)
            parser = ET.XMLParser()
            f = open(defFile)
            s = f.read()
            f.close()
            tree = ET.fromstring(s, parser)
            inputs = []
            input_els = []
            for inputCont in tree.findall("ccp4i2_body/container[@id='inputData']"):
                for inputEl in inputCont.findall("content"):
                    try:
                        inputs.append(inputEl.attrib["id"])
                        input_els.append(inputEl)
                    except:
                        print("Failed with",inputEl)
                        exc_type, exc_value,exc_tb = sys.exc_info()[:3]
                        sys.stderr.write(str(exc_type)+'\n')
                        sys.stderr.write(str(exc_value)+'\n')

            for additionalFile in tree.findall("ccp4i2_body/file"):
                try:
                    add_xml = additionalFile.findall("CI2XmlDataFile")[0]
                    add_relPath = add_xml.findall("relPath")[0].text
                    add_baseName = add_xml.findall("baseName")[0].text
                    add_defFile = os.path.join(CCP4Utils.getCCP4I2Dir(),add_relPath,add_baseName)
                    parser2 = ET.XMLParser()
                    f2 = open(add_defFile)
                    s2 = f2.read()
                    f2.close()
                    tree2 = ET.fromstring(s2, parser2)
                except:
                    print("Failed to read additional file")
                else:
                    print("Read additional file",add_defFile,tree2)
                for inputCont in tree2.findall("ccp4i2_body/container[@id='inputData']"):
                    for inputEl in inputCont.findall("content"):
                        try:
                            inputs.append(inputEl.attrib["id"])
                            input_els.append(inputEl)
                        except:
                            print("Failed with",inputEl)
                            exc_type, exc_value,exc_tb = sys.exc_info()[:3]
                            sys.stderr.write(str(exc_type)+'\n')
                            sys.stderr.write(str(exc_value)+'\n')
            print("--------------------------------------------------")
            doStuff = False
            if len(body.findall("inputData/"+oldTag+"/dbFileId")) > 0 and len(body.findall("inputData/"+oldTag+"/baseName")) > 0 and len(body.findall("inputData/"+oldTag+"/project")) > 0 and len(body.findall("inputData/"+oldTag+"/annotation")) > 0 and len(body.findall("inputData/"+oldTag+"/relPath")) > 0 and len(body.findall("inputData/"+oldTag+"/selection")) == 0:
                print("OLD JOB "+oldTag+" IS THE OLD TYPE")

                if "ASUIN" in inputs and not oldTag in inputs:
                    print("Looks like a simple "+oldTag+" -> ASUIN translation")
                    doStuff = True
                    newTag = "ASUIN"
                if "SEQIN" in inputs and not "ASUIN" in inputs:
                    print("SEQIN to SEQIN, change may or may not be required....")
                    theElement = input_els[inputs.index("SEQIN")]
                    if len(theElement.findall("className")) > 0 and theElement.findall("className")[0].text == "CSeqDataFile":
                        print("New task is a sequence too, so do nothing")
                    if len(theElement.findall("className")) > 0 and theElement.findall("className")[0].text == "CAsuDataFile":
                        print("New task is a CASU, so perhaps ought to do something...")
                        doStuff = True
                        newTag = "SEQIN"
                if "AWA_SEQIN" in inputs and not "ASUIN" in inputs:
                    theElement = input_els[inputs.index("AWA_SEQIN")]
                    if len(theElement.findall("className")) > 0 and theElement.findall("className")[0].text == "CSeqDataFile":
                        print("New task is a sequence too, so do nothing")
                    if len(theElement.findall("className")) > 0 and theElement.findall("className")[0].text == "CAsuDataFile":
                        print("New task is a CASU, so perhaps ought to do something...")
                        doStuff = True
                        newTag = "AWA_SEQIN"
            asu_inps = []
            if len(body.findall("inputData/"+oldTag+"/CSeqDataFile")) > 0:
                #Again, we need to check how this compares with current task spec.
                print("We have a list of sequences ...")
                for df in body.findall("inputData/"+oldTag+"/CSeqDataFile"):
                    seq_dbFileId = df.findall("dbFileId")[0].text
                    seq_baseName = df.findall("baseName")[0].text
                    seq_project = df.findall("project")[0].text
                    seq_annotation = df.findall("annotation")[0].text.split()[-1]
                    seq_relPath = df.findall("relPath")[0].text
                    print(seq_dbFileId,seq_baseName,seq_project,seq_annotation,seq_relPath)
                    projectName=PROJECTSMANAGER().db().getProjectInfo(projectId=seq_project,mode='projectdirectory')
                    seqFileName = os.path.join(projectName,seq_relPath,seq_baseName)
                    seqFile = CCP4ModelData.CSeqDataFile(parent=None)
                    seqFile.setFullPath(seqFileName)
                    seqFile.__dict__['format'] = 'internal'
                    seqFile.fileContent.loadInternalFile(str(seqFile))
                    seqContent = seqFile.fileContent.sequence
                    seq1 = re.sub(r'[\r\n ]', '', str(seqContent))
                    seq2 = re.sub(r'[^A-Z]', '', seq1).upper()
                    seq3 = re.sub(r'[BJOXZ]', '', seq2)
                    asu_inps.append({'name':seq_annotation,'sequence':seq3})
                if "ASUIN" in inputs and not oldTag in inputs:
                    print("Looks like a list SEQIN -> ASUIN translation")
                    newTag = "ASUIN"
                    inpD = body.findall("inputData")[0]
                    SEQIN_el = inpD.findall(oldTag)[0]
                    inpD.remove(SEQIN_el)
                    print("REMOVED OLD")
                else:
                    pass #Do nothing at present, may be more complicated in the future.
                    asu_inps = []
            print("Do stuff is", doStuff)
            if doStuff:
                seq_dbFileId = body.findall("inputData/" + oldTag + "/dbFileId")[0].text
                seq_baseName = body.findall("inputData/" + oldTag + "/baseName")[0].text
                seq_project = body.findall("inputData/" + oldTag + "/project")[0].text
                seq_annotation = body.findall("inputData/" + oldTag + "/annotation")[0].text.split()[-1] #Remove possible leading job number and space.
                seq_relPath = body.findall("inputData/" + oldTag + "/relPath")[0].text
                inpD = body.findall("inputData")[0]
                SEQIN_el = inpD.findall(oldTag)[0]
                inpD.remove(SEQIN_el)
                print("REMOVED OLD")
                projectName=PROJECTSMANAGER().db().getProjectInfo(projectId=seq_project, mode='projectdirectory')
                seqFileName = os.path.join(projectName, seq_relPath, seq_baseName)
                seqFile = CCP4ModelData.CSeqDataFile(parent=None)
                seqFile.setFullPath(seqFileName)
                seqFile.__dict__['format'] = 'internal'
                seqFile.fileContent.loadInternalFile(str(seqFile))
                seqContent = seqFile.fileContent.sequence
                seq1 = re.sub(r'[\r\n ]', '', str(seqContent))
                seq2 = re.sub(r'[^A-Z]', '', seq1).upper()
                seq3 = re.sub(r'[BJOXZ]', '', seq2)
                asu_inps = {'name' : seq_annotation, 'sequence' : seq3}
                print(asu_inps)
            if len(asu_inps) > 0:
                print(asu_inps)
                if taskFrame is not None:
                    from core import CCP4Container
                    asuJob = taskFrame.openTask(taskName='ProvideAsuContents')
                    container = CCP4Container.CContainer(parent=taskFrame, definitionFile=CCP4TaskManager.TASKMANAGER().lookupDefFile('ProvideAsuContents'), guiAdmin=True)
                    #This needs to contain list of stuff when list input!
                    container.inputData.ASU_CONTENT.set(asu_inps)
                    container.outputData.ASUCONTENTFILE.set({'project' : seq_project, 'baseName' : "ASUCONTENTFILE.asu.xml", "relPath":"CCP4_JOBS/job_" + str(asuJob.info["jobnumber"])})
                    fName = PROJECTSMANAGER().makeFileName(jobId=asuJob.jobId, mode='JOB_INPUT')
                    if os.path.isfile(fName):
                        os.remove(fName)
                    projectName = PROJECTSMANAGER().db().getProjectInfo(str(seq_project),'projectName')
                    f = CI2XmlDataFile(fullPath=fName)
                    cHeader = container.getHeader()
                    if cHeader is not None:
                        f.header.set(cHeader)
                    f.header.function.set('PARAMS')
                    f.header.jobId.set(asuJob.jobId)
                    f.header.projectId.set(seq_project)
                    f.header.projectName.set(projectName)
                    f.header.userId.set(getpass.getuser())
                    f.header.pluginName = 'ProvideAsuContents'
                    f.header.ccp4iVersion.set(CCP4Config.VERSION())
                    if sys.platform == 'darwin':
                        f.header.OS.set('MacOSX')
                    elif sys.platform[0:5] == 'linux':
                        f.header.OS.set('Linux')
                    elif sys.platform[0:3] == 'win':
                        f.header.OS.set('Windows')
                    f.header.jobNumber.set(str(asuJob.info["jobnumber"]))
                    print("SETTING JOB NUMBER IN HEADER TO",asuJob.info["jobnumber"])
                    f.header.hostName.set(CCP4Utils.getHostName())
                    f.saveFile(bodyEtree=container.getEtree())
                    PROJECTSMANAGER().runInternalTask(jobId = asuJob.jobId, projectId=seq_project, taskName = 'ProvideAsuContents')
                print("--------------------",taskName,"--------------------")
                if taskName != "clamshelx":
                    jobFiles = PROJECTSMANAGER().db().getJobFilesInfo(jobId = asuJob.jobId)
                    ASUIN = ET.SubElement(body.findall("inputData")[0],newTag)
                    ASUIN_dbFileId = ET.SubElement(ASUIN,"dbFileId")
                    ASUIN_dbFileId.text = "UNKNOWN"
                    if len(jobFiles)>0:
                        ASUIN_dbFileId.text = jobFiles[0]['fileId']
                    ASUIN_baseName = ET.SubElement(ASUIN,"baseName")
                    ASUIN_baseName.text = "ASUCONTENTFILE.asu.xml"
                    ASUIN_project = ET.SubElement(ASUIN,"project")
                    ASUIN_project.text = seq_project
                    ASUIN_annotation = ET.SubElement(ASUIN,"annotation")
                    ASUIN_annotation.text = str(asuJob.info["jobnumber"]) + " Asu content file from Define AU contents"
                    print("SETTING JOB NUMBER IN ASUIN_annotation TO",asuJob.info["jobnumber"])
                    ASUIN_relPath = ET.SubElement(ASUIN,"relPath")
                    ASUIN_relPath.text = "CCP4_JOBS/job_"+str(asuJob.info["jobnumber"])
                    print("SETTING JOB NUMBER IN ASUIN_relPath TO",asuJob.info["jobnumber"])
                    ASUIN_selection = ET.SubElement(ASUIN,"selection")
                    ASUIN_selection_item = ET.SubElement(ASUIN_selection,"item")
                    ASUIN_selection_item_key = ET.SubElement(ASUIN_selection_item,"key")
                    ASUIN_selection_item_key.text = seq_annotation
                    ASUIN_selection_item_value = ET.SubElement(ASUIN_selection_item,"value")
                    ASUIN_selection_item_value.text = "True"
                    #print "A NEW HOPE"
                    #print ET.tostring(ASUIN)
                    #print ASUIN
                    #print ET.tostring(body.findall("inputData/"+newTag)[0])

    if suggestedParams is not None:
        print("cloneI2XmlFile verdict suggestions")
        ctlParams =  body.findall("controlParameters")
        if len(ctlParams)>0:
            try:
                for suggestion in suggestedParams:
                    print(suggestion.tag)
                    poss = body.findall("controlParameters/"+suggestion.tag)
                    if len(poss)>0:
                        poss[0].text = suggestion.text
                    else: #Let's hope it's a real parameter....
                        newEl = ET.SubElement(ctlParams[0],suggestion.tag)
                        newEl.text = suggestion.text
            except:
                import traceback
                print("Some problem with cloning with verdict suggestions...."); sys.stdout.flush()
                exc_type, exc_value, exc_tb = sys.exc_info()[:3]
                sys.stderr.write(str(exc_type) + '\n')
                sys.stderr.write(str(exc_value) + '\n')
                traceback.print_tb(exc_tb)
                raise

    #Exclude interruptStatus container
    interruptContainer = body.findall("./interruptStatus")
    if len(interruptContainer) > 0:
        body.remove(interruptContainer[0])
    xFile.setFullPath(targetFile)
    xFile.saveFile(bodyEtree=body)
    """
    print "**************************************************"
    parser = ET.XMLParser()
    f = open(targetFile)
    s = f.read()
    f.close()
    tree = ET.fromstring(s, parser)
    print "The new file"
    print ET.tostring(tree)
    print "**************************************************"
    """

class CSearchPath(CCP4Data.CData):
    CONTENTS = {'name' : {'class' : CCP4Data.CString}, 'path' : {'class' : CDataFile},}
    #           'searchPath' : { 'class' : CCP4Data.CList, 'subItem' : {'class' : CDataFile}}

    def getTextItem(self):
        return self.name.__str__()


class CSearchPathList(CCP4Data.CList):
    SUBITEM = {'class' : CSearchPath}


class CExportedFile(CCP4Data.CData):
    CONTENTS = {'exportId' : {'class' : CCP4Data.CUUID}}

    def __init__(self,value=[], qualifiers={}, parent=None, name=None):
        CCP4Data.CData. __init__(self, value, parent=parent, qualifiers=qualifiers, name=name)
        #import traceback
        #print 'CExportedFile.__init__'
        #traceback.print_stack()

    def getTextItem(self):
        # This needs to tally with CExportedFileCombo.load
        from core import CCP4TaskManager
        from dbapi import CCP4DbApi
        if self.exportId.isSet():
            exFile = PROJECTSMANAGER().db().getExportFileInfo(exportId=CCP4DbApi.UUIDTYPE(self.__dict__['_value']['exportId']),mode='label')
            if len(exFile) > 0:
                title = CCP4TaskManager.TASKMANAGER().getTitle(exFile['taskname'])
                if len(exFile['exportfilename']) > 22:
                    return '..' + exFile['exportfilename'][-20:] + ' '+ exFile['jobnumber'] + ' ' + title
                else:
                    return exFile['exportfilename']+' ' + exFile['jobnumber'] + ' ' + title
        return '--'


class CExportedFileList(CCP4Data.CList):
    SUBITEM = {'class' : CExportedFile}


class CExePath(CCP4Data.CData):
    CONTENTS = {'exeName' : {'class' : CCP4Data.CString},
                'exePath' : {'class' : CDataFile, 'qualifiers' : {'mustExist' : True, 'allowUndefined' : False}}}
    CONTENTS_ORDER = ['exeName','exePath']
    QUALIFIERS = {}

    def __init__(self, value=[], qualifiers={}, parent=None, name=None):
        CCP4Data.CData. __init__(self, value, parent=parent, qualifiers=qualifiers, name=name)

    def getTextItem(self):
        if self.__dict__['_value']['exeName'].isSet():
            ret0 = str(self.__dict__['_value']['exeName'])
        else:
            ret0 = '-'
        if self.__dict__['_value']['exePath'].isSet():
            ret1 = str(self.__dict__['_value']['exePath'])
        else:
            ret1 = '-'
        return '{0} {1}'.format(ret0, ret1)


class CExePathList(CCP4Data.CList):
    SUBITEM = {'class' : CExePath}
    QUALIFIERS = {'listMinLength' : 1}

    def __init__(self, value=[], qualifiers={}, parent=None, name=None, build=True, **kw):
        CCP4Data.CList.__init__(self, value=value, qualifiers=qualifiers, parent=parent, name=name, build=build, **kw)
        self.__dict__['exeLookup'] = None

    def setupExeLookup(self):
        self.__dict__['exeLookup'] = {}
        for item in self.__dict__['_value']:
            if os.path.exists(str(item.exePath)):
                self.__dict__['exeLookup'][item.exeName] = str(item.exePath)
        #print 'CExePathList.setupExeLookup',self.__dict__['exeLookup']

    def getExecutable(self, name):
        if self.__dict__['exeLookup'] is None:
            self.setupExeLookup()
        return self.__dict__['exeLookup'].get(name, name)


#===========================================================================================================
import unittest
def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testProject)
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testFilePath))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testDataFile))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testI2XmlDataFile))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testXmlDataFile))
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)

class testProject(unittest.TestCase):

    def testCProjectName1(self):
        p = CProjectName()
        p.set('CCP4I2_TEST')
        self.assertEqual(p,'CCP4I2_TEST','Failed to set good project name')
    """
    def testCProjectName2(self):
        p = CProjectName()
        try:
          p.set('rubbish')
        except CException as e:
          self.assertEqual(len(e),1,'Unexpected exception length in setting bad CProjectName')
          self.assertEqual(e[0]['code'],104,'Unexpected exception in setting bad CProjectName')
        except:
          self.fail('Unexpected exception in setting bad CProjectName')
        else:
          self.fail('No exception in setting bad CProjectName')

    def testCProjectName3(self):
        import shutil
        from core.CCP4Utils import getHOME
        PROJECTSMANAGER().createProject('DUMMY_CCP4I2_TEST',os.path.join(getHOME(),'DUMMY_CCP4I2_TEST'))
        shutil.rmtree(os.path.join(getHOME(),'DUMMY_CCP4I2_TEST'))

        p = CProjectName()
        try:
          p.set('DUMMY_CCP4I2_TEST')
        except CException as e:
          self.assertEqual(len(e),1,'Unexpected exception length in setting project with nonexistant directory')
          self.assertEqual(e[0]['code'],106,'Unexpected exception in setting project with nonexistant directory')
        except:
          self.fail('Unexpected exception in setting bad project with nonexistant directory')
        else:
          self.fail('No exception in setting project with nonexistant directory')
        PROJECTSMANAGER().deleteProject('DUMMY_CCP4I2_TEST')
    """

class testFilePath(unittest.TestCase):

    def testFilePath1(self):
        f = CFilePath('foo/bar/fo_ehue.tmp')
        self.assertEqual(f, 'foo/bar/fo_ehue.tmp', 'Error setting CFilePath')

    def testFilePath2(self):
        try:
            f = CFilePath('foo/bar/fo_e*ue.tmp', allowedCharactersMode=CFilePath.ALLOWED_CHARACTERS_FAIL)
        except CException as e:
            self.assertEqual(len(e), 1, 'Unexpected exception length in setting bad CFilePath')
            self.assertEqual(e[0]['code'], 101, 'Unexpected exception in setting project bad CFilePath')
        except:
            self.fail('Unexpected exception in setting bad CFilePath')
        else:
            self.fail('No exception in setting bad CFilePath')

    def testFilePath3(self):
        try:
            f = CFilePath()
            f.set(f.fix('foo/bar/fo_e*ue.tmp'))
        except:
            self.fail('Unexpected exception in fixing bad CFilePath')
        self.assertEqual(f, 'foo/bar/fo_e_ue.tmp', 'Error setting CFilePath')

class testDataFile(unittest.TestCase):

    def setUp(self):
        PROJECTSMANAGER().makeDefaultProject('CCP4I2_TEST')

    def testDataFile0(self):
        f = CDataFile()
        self.assertEqual(f.baseName.get(), None,'Initialising CDataFile does not have baseName as None')
        self.assertEqual(f.fullPath.get(), None,'Initialising CDataFile does not have fullPath as None')

    def testDataFile1(self):
        f = CDataFile(projectName='CCP4I2_TEST', relPath='foo/bar', baseName='myfile.txt')
        #print 'testDataFile1',f.fullPath,f.relPath,f.project.directory()
        #print 'testDataFile1',PROJECTSMANAGER().getProjectDirectory(projectName='CCP4I2_TEST')
        #print 'testDataFile1',PROJECTSMANAGER().db().listProjects(toTerm=True)
        self.assertEqual(f.fullPath,os.path.join(PROJECTSMANAGER().getProjectDirectory(projectName='CCP4I2_TEST'),'foo/bar/myfile.txt'))

    def testDataFile2(self):
        f = CDataFile(mustExist=True)
        projectId=PROJECTSMANAGER().db().getProjectId(projectName = 'CCP4I2_TEST')
        try:
            f.set(projectId=projectId,relPath = 'foo/bar',baseName = 'myfile.txt')
        except CException as e:
            self.assertEqual(len(e),1,'Unexpected exception length in setting non-existant CDataFile')
            self.assertEqual(e[0]['code'],101,'Unexpected exception in setting non-existant CDataFile')
        except:
            self.fail('Unexpected exception in setting non-existant CDataFile')
        else:
            self.fail('No exception in setting non-existant CDataFile')

    def testDataFile3(self):
        f = CDataFile()
        f.fullPath = os.path.join(PROJECTSMANAGER().getProjectDirectory(projectName='CCP4I2_TEST'),'foo/bar/myfile.txt')
        projectName = PROJECTSMANAGER().db().getProjectInfo(str(f.project),'projectName')
        self.assertEqual(projectName,'CCP4I2_TEST','CDataFile setting full path gives wrong project')
        self.assertEqual(f.relPath,'foo/bar','CDataFile setting full path gives wrong relPath')
        self.assertEqual(f.baseName,'myfile.txt','CDataFile setting full path gives wrong baseName')

    def testDataFile4(self):
        f = CDataFile(fullPath = os.path.join(PROJECTSMANAGER().getProjectDirectory(projectName='CCP4I2_TEST'),'foo/bar/myfile.txt'))
        projectName = PROJECTSMANAGER().db().getProjectInfo(str(f.project),'projectName')
        self.assertEqual(projectName,'CCP4I2_TEST','CDataFile setting full path gives wrong project')
        self.assertEqual(f.relPath,'foo/bar','CDataFile setting full path gives wrong relPath')
        self.assertEqual(f.baseName,'myfile.txt','CDataFile setting full path gives wrong baseName')

    def testDataFile5(self):
        bindir = os.path.join(os.environ['CCP4I2_TOP'],'bin')
        f = CDataFile(relPath=bindir,baseName='browser')
        #print 'testDataFile5',bindir,f.fullPath
        h = CFilePath(os.path.join(bindir,'browser'))
        self.assertTrue(f.samefile(os.path.join(bindir,'browser')),'Failed CDataFile.samefile(Python string)')
        self.assertTrue(f.samefile(h),'Failed CDataFile.samefile(CFilePath)')

    """
      def testDataFile5(self):
        try:
          f = CDataFile(project='rubbish',baseName='whatever')
        except CException as e:
          self.assertEqual(len(e),1,'Unexpected exception length in setting CDataFile with bad project')
          self.assertEqual(e[0]['code'],104,'Unexpected exception in setting CDataFile with bad project')
        except:
          self.fail('Unexpected exception in setting CDataFile with bad project')
        else:
          self.fail('No exception in setting CDataFile with bad project')
    """


class testI2XmlDataFile(unittest.TestCase):

    def test1(self):
        c = CI2XmlDataFile( projectName='CCP4I2_TOP',relPath='test/data',baseName='pdbset.def.xml')
        self.assertEqual(c.header.pluginTitle,'PDBSet','Error loading header of I2XmlDataFile')

    def test2(self):
        c = CI2XmlDataFile(projectName='CCP4I2_TOP',relPath='test/data',baseName='pdbset.def.xml')
        e = c.getBodyEtree()
        self.assertEqual(str(e.find('container').get('id')),'inputData','Error in getBodyEtree')

    def test3(self):
        c = CI2XmlDataFile(projectName='CCP4I2_TEST',baseName='testXmlDataFile.xml')
        if c.fullPath.exists() : os.remove(c.fullPath.get())
        from core import CCP4Config
        if CCP4Config.XMLPARSER() == 'lxml': from lxml import etree
        ele = ET.Element(CI2XmlDataFile.BODY_TAG)
        c.saveFile(bodyEtree=ele)
        self.assertTrue(os.path.exists(c.fullPath.get()),'No file written by saveFile')


class testXmlDataFile(unittest.TestCase):

    def test1(self):
        c = CXmlDataFile(projectName='CCP4I2_TOP',relPath='test/data',baseName='refmac_as2m1_n.xml')
        e = c.loadFile()
        self.assertEqual(e.tag,'REFMAC','Failed to load refmac xml to CXmlDataFile')

        f = CDataFile(projectName='CCP4I2_TEST',baseName='refmac_as2m1_n.refmac.xml')
        if f.fullPath.exists():  os.remove(str(f.fullPath))
        d = c.makeI2XmlDataFile(fileName=f,function='REFMAC',pluginVersion='5.1.1')
        self.assertTrue( f.fullPath.exists(),'Failed to create I2 version of refmac xml file')


def testTiming():
    f = CDataFile()
    start = time.perf_counter()
    for n in range(0,1000):
        v = f.validity('/foo')
        print('testTiming',time.perf_counter()-start)
