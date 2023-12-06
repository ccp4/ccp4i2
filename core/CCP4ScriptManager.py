from __future__ import print_function

"""
     CCP4ScriptManager.py: CCP4 GUI Project
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
   Liz Potterton Oct 2010 - Class to keep track of all CCP4Scripts
"""
import os
import re
import glob
from core.CCP4ErrorHandling import *

def SCRIPTMANAGER():
    if CScriptManager.insts is None:
        CScriptManager.insts = CScriptManager()
    return CScriptManager.insts


class CScriptManager:

    ERROR_CODES = {101 : {'description' : 'No definintion (def.xml) file for script'},
                   102 : {'description' : 'Target directory for new plugin/pipeline does not exist'},
                   103 : {'description' : 'Directory for wrapper plugin already exists'},
                   104 : {'description' : 'No name provided for new wrapper plugin'},
                   105 : {'description' : 'Error creating plugin/pipeline directory'},
                   106 : {'description' : 'Error creating plugin/pipeline sub-directories'},
                   107 : {'description' : 'Error copying file to wrapper plugin sub-directory'},
                   108 : {'description' : 'Error editing wrapper plugin template file'},
                   109 : {'description' : 'Wrapper plugin template file not found or not readable'},
                   110 : {'description' : 'License template file not found or not readable'},
                   111 : {'description' : 'No name provided for new pipeline'},
                   112 : {'description' : 'Directory for the pipeline already exists'},
                   120 : {'description' : 'Error importing Python script'}}
    insts = None

    def __init__(self):
        from core.CCP4Utils import getCCP4I2Dir
        self._searchPath = [os.path.join(getCCP4I2Dir(),'wrappers'),os.path.join(getCCP4I2Dir(),'pipelines')]
        plines = glob.glob(os.path.join(getCCP4I2Dir(),'pipelines','*','wrappers'))
        for item in plines:
            self._searchPath.append(item)
        self.clsLookup= {}

    def searchPath(self):
        return self._searchPath

    def appendSearchPath(self,path):
        # Append a path which could include wildcard * - must work
        # with glob
        path = os.path.abspath(path)
        if not self._searchPath.count(path):
            self._searchPath.append(path)

    def buildClsLoookup(self):
        import inspect
        from core.CCP4Utils import globSearchPath,importFileModule
        from qtgui.CCP4TaskWidget import CTaskWidget
        myErrorReport = CErrorReport()
        pyFileList = globSearchPath(self.searchPath(), '*/*.py')
        #print 'pyFileList',pyFileList
        for pyFile in pyFileList:
            module,err = importFileModule(pyFile)
            if err is not None:
                myErrorReport.append(self.__class__, 999, stack=False, details=str(err))
                break
            #print 'module',module
            clsList = inspect.getmembers(module, inspect.isclass)
            for className, cls in clsList:
                if issubclass(cls,CTaskWidget):
                    # The taskName used for internal reference is the name of directory
                    # containing CTask*.py and *.def.xml files
                    taskName =  os.path.split(os.path.split(pyFile)[0])[1]
                    taskModule = getattr(cls, 'TASKMODULE', 'test')
                    self.clsLookup[taskName] = cls
                    if taskModule not in self.moduleLookup: self.moduleLookup[taskModule] = []
                    self.moduleLookup[taskModule].append(taskName)

    def printLookup(self):
        classNameList = list(self.clsLookup.keys())
        classNameList.sort()
        for item in classNameList:
            print("{0:20}".format(item, self.clsLookup[item]))

    def getClass(self, className=''):
        if className in self.clsLookup:
            return self.clsLookup[className]
        else:
            return None


    def makeWrapperPlugin(self, targetDir=None, pipeline=None, name=None, template='simple_wrapper', license='No license'):
        from core.CCP4Utils import getCCP4I2Dir, saveFile
        print('SCRIPTMANAGER.makeWrapperPlugin', name, pipeline, template, license)
        report = CErrorReport()
        if name is None:
            raise CException(self.__class__, 104)
        if pipeline is not None:
            targetDir = os.path.join(getCCP4I2Dir(), 'pipelines', pipeline, 'wrappers')
        if targetDir is None:
            targetDir = os.path.join(getCCP4I2Dir(), 'wrappers')
        if not os.path.exists(targetDir):
            raise CException(self.__class__, 101, targetDir)
        newPluginDir = os.path.join(targetDir, name)
        if os.path.exists(newPluginDir):
            report.append(self.__class__, 103, newPluginDir)
        else:
            try:
                os.mkdir(newPluginDir)
            except:
                raise CException(self.__class__, 105, newPluginDir)
        for subDirName in ['script', 'test_data']:
            subDir = os.path.join(newPluginDir, subDirName)
            if not os.path.exists(subDir):
                try:
                    os.mkdir(subDir)
                except:
                    report.append(self.__class__, 106, subDir)
        try:
            wrapperText = self.makeWrapperScript(name=name, template=template, license=license)
            wrapperFile = os.path.join(newPluginDir, 'script', name + '.py')
            saveFile(fileName=wrapperFile, text=wrapperText)
        except CException as e:
            report.extend(e)
        except Exception as e:
            report.append(self.__class__, 108)
        return newPluginDir, report

    def makeWrapperScript(self, template='simple_wrapper', license=None, name=''):
        from core.CCP4Utils import getCCP4I2Dir, readFile
        if license is not None:
            licenseFile = os.path.join(getCCP4I2Dir(), 'data', 'plugin_templates', 'license_templates', license)
            licenseText = readFile(licenseFile)
            if len(licenseText) == 0:
                raise CException(self.__class__, 110, licenseFile)
        else:
            licenseText = ''
        templateFile = os.path.join(getCCP4I2Dir(), 'data', 'plugin_templates', template)
        if templateFile[-3:] != '.txt':
            templateFile = templateFile + '.txt'
        text = readFile(templateFile)
        if len(text) == 0:
            raise CException(self.__class__, 109)   # KJS : removed newPluginDir (not defined here)
        text = licenseText + text
        wrapperText = re.sub('\$NAME\$', name, text)
        return wrapperText


    def listOfPipelines(self):
        from core.CCP4Utils import getCCP4I2Dir
        dirList = glob.glob(os.path.join(getCCP4I2Dir(), 'pipelines', '*'))
        pipelineList = []
        for d in dirList:
            pipelineList.append(os.path.split(d)[1])
        return pipelineList

    def listOfPluginTemplates(self):
        from core.CCP4Utils import getCCP4I2Dir
        fileList = glob.glob(os.path.join(getCCP4I2Dir(), 'data', 'plugin_templates', '*'))
        templateList = []
        for f in fileList:
            templateList.append(os.path.split(f)[1])
        templateList.remove('license_templates')
        return templateList

    def listOfPluginLicenses(self):
        from core.CCP4Utils import getCCP4I2Dir
        fileList = glob.glob(os.path.join(getCCP4I2Dir(), 'data', 'plugin_templates', 'license_templates', '*'))
        licenseList = []
        for f in fileList:
            licenseList.append(os.path.split(f)[1])
        return licenseList

    def makePipeline(self, name=None, targetDir=None):
        from core.CCP4Utils import getCCP4I2Dir
        if name is None:
            raise CException(self.__class__, 111)
        if targetDir is None:
            targetDir = os.path.join(getCCP4I2Dir(), 'pipelines')
        if not os.path.exists(targetDir):
            raise CException(self.__class__, 102, targetDir)
        pipelineDir = os.path.join(targetDir, name)
        if os.path.exists(pipelineDir):
            raise CException(self.__class__, 112, pipelineDir)
        try:
            os.mkdir(pipelineDir)
        except:
            raise CException(self.__class__, 105, pipelineDir)
        for subDirName in ['script', 'test_data', 'wrappers']:
            subDir = os.path.join(pipelineDir, subDirName)
            try:
                os.mkdir(subDir)
            except:
                raise CException(self.__class__, 106, subDir)
        return pipelineDir
