from __future__ import print_function

"""
     CCP4Config.py: CCP4 GUI Project
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
   Liz Potterton Sept 2010 - Separate CCP4Config out from core.CCP4Modules
"""

import os
import sys
import re
import glob
from lxml import etree
from xml.etree import ElementTree as ET

def DEFCONFIG():
    if  CConfig.insts is None:
        from core import CCP4Utils
        CConfig(os.path.join(CCP4Utils.getDotDirectory(),'configs','ccp4i2_config.params.xml'))
    return CConfig.insts

def CONFIG(fileName=None, mode='ccp4i2'):
    if CConfig.insts is None:
        CConfig(fileName, mode=mode)
    return CConfig.insts

utf8_parser = etree.XMLParser(encoding='utf-8')

def parse_from_unicode(unicode_str):
    s = unicode_str.encode('utf-8')
    return ET.fromstring(etree.tostring(etree.fromstring(s, parser=utf8_parser)))

class CConfig:
    '''This class should be instantiated when process (browser,jobcontroller) is started'''
    insts = None

    def __init__(self, filename=None, mode='ccp4i2', **kw):
        from core import CCP4Utils
        if CConfig.insts is None:
            CConfig.insts = self
        self._xmlMode = 'lxml'
        self.developer = True
        self.graphical = False
        self.qt = True
        self.jobControllerMode = 'server'
        self.dbFile = None
        self.dbMode ='sqlite'
        self.dbUser = None
        self.maxRunningProcesses = 4
        # search paths for external programs - convention: program name is lower case
        self.searchPath = {}
        self.loadVersion()
        # Load local installation config file
        localFile = os.path.join(CCP4Utils.getCCP4I2Dir(), 'local_setup','ccp4i2_config.params.xml')
        if os.path.exists(localFile):
            self.loadDataFromXml(localFile)
        # Load users config file
        if filename is None and mode is not None:
            filename = os.path.join(CCP4Utils.getDotDirectory(), 'configs', mode + '_config.params.xml')
        if filename is not None:
            if os.path.exists(filename):
                self.loadDataFromXml(filename)
            else:
                self.saveDataToXml(filename)
        if len(kw) > 0:
            for key, value in list(kw.items()):
                if hasattr(self, key):
                    setattr(self, key, value)
        #print 'CConfig.init', self.searchPath

    def loadDataFromXml(self, fileName):
        errList = []
        from core import CCP4Utils
        text = CCP4Utils.readFile(fileName)
        #print 'CConfig.load text', fileName, text
        root = parse_from_unicode(text)
        body = root.find('ccp4i2_body')
        for element in body.iter():
            tag = element.tag
            value = element.text
            #print 'CConfig.loadDataFromXml', tag, value, type(tag), type(value)
            if tag in ['developer', 'graphical', 'qt']:
                if value.lower() == 'true':
                    setattr(self, tag, True)
                elif value.lower() == 'false':
                    setattr(self, tag, False)
                else:
                    errList.append(tag)
            elif tag == 'dbMode':
                if value in ['sqlite', 'qt_sqlite']:
                    self.dbMode = value
                else:
                    errList.append(tag)
            elif tag in ['dbFile', 'dbUser']:
                if value.lower() == 'none':
                    setattr(self, tag, None)
                else:
                    setattr(self, tag, value)
            elif tag in ['maxRunningProcesses']:
                if value.lower() == 'none':
                    setattr(self, tag, 1)
                else:
                    setattr(self, tag, int(value))
            # Assume that there could be more than one file providing this config info
            # and aim to give the contents of the last (user-specific) file priority in the search path
            elif tag in ['searchPath']:
                for exe in element:
                    if exe.tag not in self.searchPath:
                        self.searchPath[exe.tag] = {}
                    for platform in exe:
                        if platform.tag not in self.searchPath[exe.tag]:
                            self.searchPath[exe.tag][platform.tag] = []
                        itemList = []
                        for item in platform:
                            itemList.append(item.text)
                        itemList.reverse()
                        for item in itemList:
                            self.searchPath[exe.tag][platform.tag].insert(0, item)
        if len(errList) > 0:
            print('ERROR loading from config file:', fileName)
            print('ERROR with data:', end=' ')
            for item in errList:
                print(item, end=' ')
            print(' ')

    def saveDataToXml(self, fileName):
        from core import CCP4File
        root = ET.Element('configs')
        for tag in ['developer', 'graphical', 'qt', 'dbMode', 'dbFile', 'dbUser', 'maxRunningProcesses']:
            ele = ET.Element(tag)
            ele.text = str(getattr(self, tag))
            root.append(ele)
        ele = ET.Element('searchPath')
        root.append(ele)
        #print 'CConfig.saveDataToXml',self.searchPath.keys()
        for exe in list(self.searchPath.keys()):
            exeEle = ET.Element(exe)
            ele.append(exeEle)
            for platform in list(self.searchPath[exe].keys()):
                platformEle = ET.Element(platform)
                exeEle.append(platformEle)
                for item in self.searchPath[exe][platform]:
                    pathEle = ET.Element('path')
                    pathEle.text = item
                    platformEle.append(pathEle)
        # Initialise CI2XmlDataFile thisway to avoid calling PROJECTSMANAGER().db to
        # attempt to interpret a fullPath.  This CConfig may be instantiated before
        # the PROJECTSMANAGER has been
        relPath,baseName = os.path.split(fileName)
        f = CCP4File.CI2XmlDataFile({'relPath' :relPath, 'baseName' : baseName })
        f.header.setCurrent()
        f.header.function.set('PARAMS')
        f.header.pluginName = 'ccp4i2_config'
        f.header.pluginTitle = 'CCP4i2 Configuration'
        ET.indent(root)
        f.saveFile(root)

    def loadVersion(self):
        # Do not use CCP4File here - config not set - broken without Qt 
        from core import CCP4Utils
        text = CCP4Utils.readFile(os.path.join(CCP4Utils.getCCP4I2Dir(), 'core', 'version.params.xml'))
        _tree = ET.fromstring(text)
        tree=ET.Element('ccp4i2root')
        tree.append(_tree)
        self.ccp4iVersion = tree.findall('./ccp4i2/ccp4i2_header/ccp4iVersion')[0].text
        print('ccp4i2 version', self.ccp4iVersion)
        revision = CCP4Utils.getProgramVersion('ccp4i2',mode='revision')
        print('ccp4i2 source revision', revision)

    def set(self, key='', value=''):
        if ['qt', 'developer', 'graphical'].count(key) and [True, False].count(value):
            setattr(self, key, value)
        elif key == 'dbMode' and value in ['sqlite']:
            setattr(self, key, value)
        elif key in ['dbFile', 'dbUser']:
            setattr(self, key, value)

    def setXmlMode(self):
        try:
            self._xmlMode = 'lxml'
        except:
            try:
                self._xmlMode = 'elementtree'
            except:
                self._xmlMode = ''

    def xmlMode(self):
        if self._xmlMode is None:
            self.setXmlMode()
        return self._xmlMode

#==============================================FUNCTIONS

def PATH(exe, firstOnly = True):
    if not CConfig.insts:
        CConfig()
    exe = exe.lower()
    if sys.platform in ['win32']:
        platform = 'windows'
    elif  sys.platform in ['darwin']:
        platform = 'macosx'
    else:
        platform = 'linux'
    #print 'CConfig.PATH',exe,platform,CConfig.insts.searchPath
    if not hasattr(CConfig.insts,'searchPath') or exe not in CConfig.insts.searchPath or \
                         platform not in CConfig.insts.searchPath[exe]:
        if firstOnly:
            return None
        else:
            return []
    for path in CConfig.insts.searchPath[exe][platform]:
        # Test for envvar that glob does not seem to handle
        m = re.match('\$([^/]*)(.*)',path)
        if m is not None:
            env,relpath = m.groups()
            if env in os.environ:
                exeList = glob.glob(os.path.join(os.environ[env], relpath))
        else:
            exeList = glob.glob(path)
        #print 'CConfig.PATH',path,m,exeList
        if len(exeList) == 1:
            return exeList[0]
        elif len(exeList) > 1:
            # Sort to reverse data order (ie most recent first)
            exeList.sort(cmpFileData)   # : KJS This doesn't look it will work.... (no cmpFileData !)
            if firstOnly:               # : KJS Looks like this file is never used, barring some commented out lines
                return exeList[0]
            else:
                return exeList
    if firstOnly:
        return None
    else:
        return []

def cmpFileDate(file1, file2):
    if os.path.getmtime(file1) > os.path.getmtime(file1):
        return 1
    else:
        return -1
  
def DEVELOPER():
    if not CConfig.insts:
        CConfig()
    return CConfig.insts.developer

def VERSION():
    if not CConfig.insts:
        CConfig()
    return CConfig.insts.ccp4iVersion

def GRAPHICAL():
    if not CConfig.insts:
        CConfig()
    return CConfig.insts.graphical

def QT():
    # Beware this loads config params file with dependencies
    # on CCP4Data which may not yet be properly loaded
    if not CConfig.insts:
        CConfig()
    return CConfig.insts.qt

def XMLPARSER():
    # Ditto comments in QT()
    if  not CConfig.insts:
        CConfig()
    return CConfig.insts.xmlMode()

def DBMODE():
    if  not CConfig.insts:
        CConfig()
    return CConfig.insts.dbMode

def DBFILE():
    if  not CConfig.insts:
        CConfig()
    return CConfig.insts.dbFile

def DBUSER():
    if  not CConfig.insts:
        CConfig()
    return CConfig.insts.dbUser 
