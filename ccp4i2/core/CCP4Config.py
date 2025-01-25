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

   Liz Potterton Sept 2010 - Separate CCP4Config out from core.CCP4Modules
"""

import os

from lxml import etree

from . import CCP4File, CCP4Utils
from .. import __version__


def CONFIG(fileName=None):
    if CConfig.insts is None:
        CConfig(fileName)
    return CConfig.insts


utf8_parser = etree.XMLParser(encoding='utf-8')


def parse_from_unicode(unicode_str):
    s = unicode_str.encode('utf-8')
    return etree.fromstring(s, parser=utf8_parser)


class CConfig:
    '''This class should be instantiated when process (browser,jobcontroller) is started'''
    insts = None

    def __init__(self, filename=None, **kw):
        if CConfig.insts is None:
            CConfig.insts = self
        self.developer = True
        self.graphical = False
        self.dbFile = None
        self.dbUser = None
        self.maxRunningProcesses = 4
        # search paths for external programs - convention: program name is lower case
        self.searchPath = {}
        print('ccp4i2 version', __version__)
        # Load local installation config file
        localFile = os.path.join(CCP4Utils.getCCP4I2Dir(), 'local_setup','ccp4i2_config.params.xml')
        if os.path.exists(localFile):
            self.loadDataFromXml(localFile)
        # Load users config file
        if filename is None:
            filename = os.path.join(CCP4Utils.getDotDirectory(), 'configs', 'ccp4i2_config.params.xml')
        if os.path.exists(filename):
            self.loadDataFromXml(str(filename))
        else:
            self.saveDataToXml(str(filename))
        for key, value in list(kw.items()):
            if hasattr(self, key):
                setattr(self, key, value)

    def loadDataFromXml(self, fileName):
        errList = []
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
        root = etree.Element('configs')
        for tag in ['developer', 'graphical', 'qt', 'dbFile', 'dbUser', 'maxRunningProcesses']:
            ele = etree.Element(tag)
            ele.text = str(getattr(self, tag))
            root.append(ele)
        ele = etree.Element('searchPath')
        root.append(ele)
        #print 'CConfig.saveDataToXml',self.searchPath.keys()
        for exe in list(self.searchPath.keys()):
            exeEle = etree.Element(exe)
            ele.append(exeEle)
            for platform in list(self.searchPath[exe].keys()):
                platformEle = etree.Element(platform)
                exeEle.append(platformEle)
                for item in self.searchPath[exe][platform]:
                    pathEle = etree.Element('path')
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
        f.saveFile(root)

    def set(self, key='', value=''):
        if ['qt', 'developer', 'graphical'].count(key) and [True, False].count(value):
            setattr(self, key, value)
        elif key in ['dbFile', 'dbUser']:
            setattr(self, key, value)

#==============================================FUNCTIONS

def DEVELOPER():
    if not CConfig.insts:
        CConfig()
    return CConfig.insts.developer

def GRAPHICAL():
    if not CConfig.insts:
        CConfig()
    return CConfig.insts.graphical

def DBFILE():
    if  not CConfig.insts:
        CConfig()
    return CConfig.insts.dbFile

def DBUSER():
    if  not CConfig.insts:
        CConfig()
    return CConfig.insts.dbUser
