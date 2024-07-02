from __future__ import print_function

"""
     CCP4DataManager.py: CCP4 GUI Project
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
   Liz Potterton Aug 2010 - Class to keep track of all CCP4Data and CCP4Widget classes
"""

import os
import re
import types
import json
from core.CCP4Config import GRAPHICAL,DEVELOPER
from core.CCP4ErrorHandling import *
from core import CCP4Utils
import importlib


TIMING=True
def DATAMANAGER():
    if CDataManager.insts is None:
        if TIMING:
            import time
            t1 = time.time()
        CDataManager.insts = CDataManager()
        if TIMING:
            t2 = time.time()
            print("Starting ",CDataManager.insts.__class__, t2-t1)
    return CDataManager.insts

class CDataManager:

    EXCLUDE_CLASSES = ['CData', 'CBaseData', 'CCollection', 'CList', 'CDict', 'CContainer']
    EXCLUDE_FILES = ['CCP4ComFilePatchManager.py', 'CCP4CustomTaskManager.py', 'CCP4ImportedJobManager.py', 'CCP4WorkflowManager.py', 'CCP4Preferences.py']
    ERROR_CODES = {101: { 'description' : 'No model given for widget()'},
                   102: { 'description' : 'The model is not a CCP4 CData'},
                   103: { 'description' : 'No parent widget given for widget()'},
                   104: { 'description' : 'Parent widget is not a Qt QWidget'},
                   105: { 'description' : 'No suitable widget class found for model'},
                   106: { 'description' : 'Error handling qualifiers for widget'},
                   107: { 'description' : 'Loading DATAMANAGER, failed attempting to import module'},
                   108: { 'description' : 'Undetermined error creating widget'}}

    insts = None

    def __init__(self):
        self.clsLookup= {}
        self.widgetLookup= {}
        self.viewLookup= {}
        self.toUpper = {}
        self._searchPath = [os.path.join(CCP4Utils.getCCP4I2Dir(), 'core')]
        if GRAPHICAL():
            self._searchPath.extend([os.path.join(CCP4Utils.getCCP4I2Dir(), 'qtgui')])
        rv = self.buildClassLookup()
        if len(rv) > 0:
            print(rv.report())

    def searchPath(self):
        return self._searchPath

    def reportClassSearchPath(self):
        return [os.path.join(CCP4Utils.getCCP4I2Dir(), 'report')]

    def appendSearchPath(self,path):
        # Append a path which could include wildcard * - must work
        # with glob
        path = os.path.abspath(path)
        if not self._searchPath.count(path):
            self._searchPath.append(path)

    def customClasses(self):
        if GRAPHICAL():
            return []
        else:
            return []

    def buildClassLookup(self):
        try:
            myErrorReport = CErrorReport()
            import os
            cacheDir = os.path.split(__file__)[0]
            cacheFilePath = os.path.join(cacheDir,"CDataManagerCache.json")
            with open(cacheFilePath,"r") as cacheFile:
                jsonText = cacheFile.read()
                compositeDict = json.loads(jsonText)
                for key in compositeDict:
                    #print ('Processing lookup', key)
                    lookup = compositeDict.get(key)
                    if key != "toUpper":
                        for lookupKey in lookup:
                            lookup[lookupKey]['class'] = None
                    setattr(self, key, lookup)
                    # Check for issues (duplicates in clsLookup)
                    chkset = set()
                    #chkset.add("CMtzColumn") # add deliberate error
                    if key == "clsLookup":
                        for key in lookup:
                            lookup_cont = lookup[key]
                            clsNam = lookup_cont.get("clsName")
                            if clsNam not in chkset :
                                chkset.add(clsNam)
                            else:
                                print ("DUPLICATION PROBLEM flagged in CDataManagerCache JSON : Class ", clsNam)
                                tmpM = "Potentially serious duplication problem flagged.\n" + \
                                       "    " +clsNam + "\n"
                                tmpM += "Recommend using grep to look for duplication in code"
                                if GRAPHICAL():
                                    QtGui.QMessageBox.warning(None, "Problem Loading CachedLookups JSON", tmpM,
                                                                    "Developer Action Recommended")
                            #print ("KEY:", key, "  CONT:", lookup[key])
                    #print (lookup)
            return myErrorReport
        except:
            print('Falling back to building CDataManager lookups from scratch')
            return self.buildClassLookupFromScratch()

    def buildClassLookupFromScratch(self):
        import inspect
        import core
        import qtgui
        from core import CCP4Data
        from qtgui import CCP4Widgets
        
        myErrorReport = CErrorReport()
        
        #print dir(core)
        #print dir(qtgui)
        modules = inspect.getmembers(core, inspect.ismodule)
        if GRAPHICAL():
            graphModules = inspect.getmembers(qtgui, inspect.ismodule)
            modules = (modules+graphModules)
        #print "CDataManager.buildClassLookupFromScratch modules:"
        #for module in modules:
        #    print module
        '''
        for moduleName, module in modules:
            clsList = inspect.getmembers(module, inspect.isclass)
            '''
        import os, glob
        coreDirSearch = os.path.join(CCP4Utils.getCCP4I2Dir(), 'core', '*.py')
        #print "#CDataManager.buildClassLookup",coreDirSearch
        coreFiles = glob.glob(coreDirSearch)
        #print "#CDataManager.buildClassLookup",coreFiles
        coreModules = []
        for coreFile in coreFiles:
            filename, fileextension = os.path.splitext(coreFile)
            fileroot = os.path.basename(filename)
            coreModules.append("core."+fileroot)
        #print coreModules
        import importlib
        
        qtguiDirSearch = os.path.join(CCP4Utils.getCCP4I2Dir(), 'qtgui', '*.py')
        #print "#CDataManager.buildClassLookup",qtguiDirSearch
        qtguiFiles = glob.glob(qtguiDirSearch)
        #print "#CDataManager.buildClassLookup",qtguiFiles
        qtguiModules = []
        for qtguiFile in qtguiFiles:
            filename, fileextension = os.path.splitext(qtguiFile)
            fileroot = os.path.basename(filename)
            qtguiModules.append("qtgui."+fileroot)
        #print qtguiModules
        
        allModules = coreModules + qtguiModules
        
        for moduleName in allModules:
            module = importlib.import_module(moduleName)
            clsList = inspect.getmembers(module, inspect.isclass)
            
            for name, cls in clsList:
                if issubclass(cls, CCP4Data.CData):
                    self.clsLookup[name] = {'class':cls,'clsModule':cls.__module__,'clsName':cls.__name__}
                    self.toUpper[name.lower()] = name
            if GRAPHICAL():
                for name,cls in clsList:
                    if issubclass(cls, CCP4Widgets.CViewWidget):
                        self.widgetLookup[name] = {'class':cls,'clsModule':cls.__module__,'clsName':cls.__name__}
                        model = getattr(cls, 'MODEL_CLASS', None)
                        if model is not None:
                            if isinstance(model, tuple):
                                for item in model:
                                    self.viewLookup[item.__name__] = {'class':cls,'clsName':cls.__name__,'clsModule':cls.__module__}
                            else:
                                self.viewLookup[model.__name__] = {'class':cls,'clsName':cls.__name__,'clsModule':cls.__module__}
        for name, cls, widgetCls in self.customClasses():
            self.clsLookup[name] = {'class':cls, 'clsModule':cls.__module__,'clsName':cls.__name__}
            self.toUpper[name.lower()] = name
            if widgetCls is not None:
                self.viewLookup[cls.__name__] = {'class':widgetCls,'clsName':widgetCls.__name__,'clsModule':widgetCls.__module__}
        
        compositeDict = {'viewLookup':self.viewLookup, 'widgetLookup':self.widgetLookup, 'clsLookup':self.clsLookup,'toUpper':self.toUpper}
        import os
        cacheDir = os.path.split(__file__)[0]
        #print "CacheDir is", cacheDir
        cacheFilePath = os.path.join(cacheDir,"CDataManagerCache.json")
        #print "cacheFilePath is",cacheFilePath
        import json
        class MyEncoder(json.JSONEncoder):
            def default(self,obj):
                #print "#CTaskManager.MyEncoder Handling case of",obj
                try:
                    return obj.__module__+obj.__name__
                except:
                    return "UnSerializable Object"
        with open(cacheFilePath,"w") as cacheFile:
            #print "file open"
            try:
                jsonText = json.dumps(compositeDict, indent=4, sort_keys=True, separators=(',', ': '), cls=MyEncoder)
                cacheFile.write(jsonText)
            except Exception as e:
                print(e)
        return myErrorReport

    def printLookup(self):
        classNameList = list(self.clsLookup.keys())
        classNameList.sort()
        for item in classNameList:
            print("{0:20} {1:30} {2:30}".format(item, self.clsLookup[item], self.viewLookup.get(self.clsLookup[item], 'No widget')))

    def getClass(self, className=''):
        import importlib
        #print "#CDataManager.getClass", className
        if className in self.clsLookup:
            if self.clsLookup[className].get('class',None) is None:
                module = importlib.import_module(self.clsLookup[className]['clsModule'])
                clsName = self.clsLookup[className]['clsName']
                self.clsLookup[className]['class'] = getattr(module, clsName)
            #print "#CDataManager.getClass", className, self.clsLookup[className]['class']
            return self.clsLookup[className]['class']
        elif className in self.toUpper:
            if self.clsLookup[self.toUpper[className]].get('class',None) is None:
                module = importlib.import_module(self.clsLookup[self.toUpper[className]]['clsModule'])
                clsName = self.clsLookup[self.toUpper[className]]['clsName']
                self.clsLookup[self.toUpper[className]]['class'] = getattr(module, clsName)
            #print "#CDataManager.getClass", className, self.clsLookup[self.toUpper[className]]['class']
            return self.clsLookup[self.toUpper[className]]['class']
        elif className in self.widgetLookup:
            if self.widgetLookup[className].get('class',None) is None:
                module = importlib.import_module(self.widgetLookup[className]['clsModule'])
                clsName = self.widgetLookup[className]['clsName']
                self.widgetLookup[className]['class'] = getattr(module, clsName)
            #print "#CDataManager.getClass", className, self.widgetLookup[className]['class']
            return self.widgetLookup[className]['class']
        else:
            return None

    def getWidgetClass(self, model=None, modelClass=None, name=''):
        from core import CCP4Data
        if modelClass is None:
            try:
                modelClass = model.__class__
            except:
                pass
        modelClassName = modelClass.__name__
        #raise CException(self.__class__, 102)
        widgetClassDict = self.viewLookup.get(modelClassName, None)
        #print 'getWidgetClass',name,modelClass,widgetClass
        if widgetClassDict is not None:
            if widgetClassDict.get('class',None) is None:
                import importlib
                module = importlib.import_module(widgetClassDict['clsModule'])
                widgetClassDict['class'] = getattr(module,widgetClassDict['clsName'])
            return widgetClassDict['class']

        # Try looking for a CView for a parent class
        cls = modelClass
        widgetClass = None
        while cls != CCP4Data.CData and widgetClass is None:
            try:
                cls = cls.__bases__[0]
                clsName = cls.__name__
                widgetClassDict = self.viewLookup.get(clsName, None)
                if widgetClassDict is not None:
                    if widgetClassDict.get('class',None) is None:
                        import importlib
                        module = importlib.import_module(widgetClassDict['clsModule'])
                        widgetClassDict['class'] = getattr(module,widgetClassDict['clsName'])
                        widgetClass =widgetClassDict['class']
                    return widgetClassDict['class']
            except:
                print("Fell back to CCP4Data.CData", cls)
                cls = CCP4Data.CData
            print('getWidgetClass cls, widgetClass', cls, widgetClass)
        if widgetClass is None:
            raise CException(self.__class__, 105, 'For data: ' + name)
        return widgetClass

    def widget(self, model=None, parentWidget=None, qualifiers={}, name='', modelClass=None):
        if not GRAPHICAL():
            return None
        from PySide6 import QtGui, QtWidgets
        widgetClass = self.getWidgetClass(model=model, modelClass=modelClass, name=name)
        if parentWidget is None:
            raise CException(self.__class__, 103)
        if not isinstance(parentWidget, QtWidgets.QWidget):
            raise CException(self.__class__, 104)
        try:
            if model is not None:
                qualis = {'dragType' : model.__class__.__name__[1:]}
                qualis.update(model.qualifiers())
            elif modelClass is not None:
                qualis = {'dragType' : modelClass.__name__[1:]}
            else:
                qualis = {}
            qualis.update(qualifiers)
        except:
            raise CException(self.__class__, 106)
        if DEVELOPER():
            widget = widgetClass(parent=parentWidget, model=model, qualifiers=qualis)
        else:
            try:
                widget = widgetClass(parent=parentWidget, model=model, qualifiers=qualis)
            except CException as e:
                raise e
            except:
                raise CException(self.__class__, 108, stack=False)
        return widget

    def buildQStandardItemModel(self, parent=None, mode=None):
        from PySide6 import QtGui, QtWidgets
        import inspect
        from core import CCP4Data
        myErrorReport = CErrorReport()
        model = QtGui.QStandardItemModel(parent)
        root = model.invisibleRootItem()
        if mode == 'repgen':
            pyFileList = CCP4Utils.globSearchPath(self.reportClassSearchPath(), '*.py')
        else:
            pyFileList = CCP4Utils.globSearchPath(self.searchPath(), '*.py')
        for pyFile in pyFileList:
            if os.path.split(pyFile)[1] not in self.EXCLUDE_FILES:
                module, err = CCP4Utils.importFileModule(pyFile)
                moduleItem = None
                if module is None:
                    myErrorReport.append(self.__class__, 107, pyFile, stack=False, details=str(err))
                clsList = inspect.getmembers(module, inspect.isclass)
                for name,cls in clsList:
                    if issubclass(cls, CCP4Data.CData) and not ['CData', 'CBaseData', 'CCollection', 'CList', 'CDict', 'CContainer'].count(name):
                        if moduleItem is None:
                            moduleItem = QtGui.QStandardItem(module.__name__)
                            root.appendRow(moduleItem)
                        clsItem = QtGui.QStandardItem(cls.__name__)
                        if cls.__doc__ is not None:
                            clsItem.setToolTip(cls.__doc__)
                        moduleItem.appendRow(clsItem)
        return model

    def makeHtmlClassListing(self):
        import inspect
        from qtgui import CCP4DefEd
        from core import CCP4Data
        from core import CCP4Modules
        from core import CCP4Annotation
        text = {}
        indexText = ''
        myErrorReport = CErrorReport()   # KJS : seems to have been missing
        pyFileList = CCP4Utils.globSearchPath(self.searchPath(), '*.py')
        for pyFile in pyFileList:
            if os.path.split(pyFile)[1] not in self.EXCLUDE_FILES:
                module, err = CCP4Utils.importFileModule(pyFile)
                if module is None:
                    myErrorReport.append(self.__class__, 107, pyFile, stack=False, details=str(err))
                moduleName = os.path.splitext(os.path.basename(pyFile))[0]
                text = '<html>\n'
                clsList = inspect.getmembers(module, inspect.isclass)
                nClass = 0
                for name, cls in clsList:
                    if issubclass(cls,CCP4Data.CData) and not name in self.EXCLUDE_CLASSES:
                        nClass = nClass + 1
                        info = CCP4DefEd.CClassInfo(name, CCP4Modules.QTAPPLICATION())
                        text = text + info.getInfo(html=False)
                text = text + '</html>\n'
                if nClass > 0:
                    fileName = os.path.join(CCP4Utils.getCCP4I2Dir(), 'docs', 'developers', 'modules', moduleName + '.html')
                    CCP4Utils.saveFile(fileName=fileName, text=text)
                    indexText = indexText + '<h4>' + moduleName + '</h4>\n'
                    for name,cls in clsList:
                        if issubclass(cls, CCP4Data.CData) and not name in self.EXCLUDE_CLASSES:
                            doc = cls.__doc__
                            if doc is None:
                                doc = ''
                            indexText = indexText + '<a href="./' + moduleName + '.html#' + name + '">' + name + '</a> ' + doc+ '<br/>\n'
        index0Text = CCP4Utils.readFile(fileName= os.path.join(CCP4Utils.getCCP4I2Dir(),'docs','developers','modules', 'index0.html'))
        allText = re.sub('INSERT_CLASS_LINKS_HERE', indexText, index0Text)
        t = CCP4Annotation.CTime()
        t.setCurrentTime()
        allText = re.sub('INSERT_DATE', str(t), allText)
        CCP4Utils.saveFile(fileName=os.path.join(CCP4Utils.getCCP4I2Dir(), 'docs', 'developers', 'modules', 'index.html'), text=allText)
        return

#=======================================================================================================
import unittest

# Weird - this is failing with Qt error message - can not create widget when no GUI
# But works OK when do same thing manually 

'''
def TESTSUITE():
  suite = unittest.defaultTestLoader.loadTestsFromTestCase(testDataManager)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)

'''


class testDataManager(unittest.TestCase):

    def setUp(self):
        '''
        if GRAPHICAL():
          from core.CCP4Modules import QTAPPLICATION
          self.app = QTAPPLICATION()
          print 'testDataManager.setUp',GRAPHICAL(),type(self.app)
          from PySide6 import QtGui, QtWidgets
          self.dialog = QtWidgets.QDialog()
        '''
        from core.CCP4Data import CFloat
        self.testData = CFloat()
        self.manager = DATAMANAGER()

    def test1(self):
        from qtgui import CCP4Widgets
        if GRAPHICAL():
            from core.CCP4Modules import QTAPPLICATION
            app = QTAPPLICATION()
            widgetClass = self.manager.getWidgetClass(self.testData)
            self.assertEqual(widgetClass, CCP4Widgets.CFloatView, 'Failed to return correct widget class')
        else:
            self.fail('Can not test CCP4DataMananger in non-graphical mode')

    def test2(self):
        if GRAPHICAL():
            from core.CCP4Modules import QTAPPLICATION
            app = QTAPPLICATION()
            from qtgui import CCP4Widgets
            from PySide6 import QtGui, QtWidgets
            dialog = QtWidgets.QDialog()
            widget = self.manager.widget(model=self.testData,parentWidget=dialog)
            self.assertTrue(isinstance(widget,CCP4Widgets.CFloatView),'Failed to create CFloatView widget')
        else:
            self.fail('Can not test CCP4DataMananger in non-graphical mode')

