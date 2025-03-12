"""
Liz Potterton Aug 2010 - Class to keep track of all CCP4Data and CCP4Widget classes
"""

import glob
import importlib
import inspect
import json
import os

from PySide2 import QtGui, QtWidgets

from . import CCP4Utils
from .CCP4Config import CONFIG
from .CCP4ErrorHandling import CException


def DATAMANAGER():
    if CDataManager.insts is None:
        CDataManager.insts = CDataManager()
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
        if CONFIG().graphical:
            self._searchPath.extend([os.path.join(CCP4Utils.getCCP4I2Dir(), 'qtgui')])
        self.buildClassLookup()

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
        return []

    def buildClassLookup(self):
        try:
            cacheDir = os.path.split(__file__)[0]
            cacheFilePath = os.path.join(cacheDir,"CDataManagerCache.json")
            with open(cacheFilePath,"r") as cacheFile:
                jsonText = cacheFile.read()
                compositeDict = json.loads(jsonText)
                for key in compositeDict:
                    lookup = compositeDict.get(key)
                    if key != "toUpper":
                        for lookupKey in lookup:
                            lookup[lookupKey]['class'] = None
                    setattr(self, key, lookup)
                    # Check for issues (duplicates in clsLookup)
                    chkset = set()
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
                                if CONFIG().graphical:
                                    QtGui.QMessageBox.warning(None, "Problem Loading CachedLookups JSON", tmpM,
                                                                    "Developer Action Recommended")
        except:
            print('Falling back to building CDataManager lookups from scratch')
            self.buildClassLookupFromScratch()

    def buildClassLookupFromScratch(self):
        from . import CCP4Data
        from .. import core
        from .. import qtgui
        from ..qtgui import CCP4Widgets

        modules = inspect.getmembers(core, inspect.ismodule)
        if CONFIG().graphical:
            graphModules = inspect.getmembers(qtgui, inspect.ismodule)
            modules = (modules+graphModules)
        coreDirSearch = os.path.join(CCP4Utils.getCCP4I2Dir(), 'core', '*.py')
        coreFiles = glob.glob(coreDirSearch)
        coreModules = []
        for coreFile in coreFiles:
            filename, fileextension = os.path.splitext(coreFile)
            fileroot = os.path.basename(filename)
            coreModules.append("core."+fileroot)
        qtguiDirSearch = os.path.join(CCP4Utils.getCCP4I2Dir(), 'qtgui', '*.py')
        qtguiFiles = glob.glob(qtguiDirSearch)
        qtguiModules = []
        for qtguiFile in qtguiFiles:
            filename, fileextension = os.path.splitext(qtguiFile)
            fileroot = os.path.basename(filename)
            qtguiModules.append("qtgui."+fileroot)
        
        allModules = coreModules + qtguiModules
        
        for moduleName in allModules:
            module = importlib.import_module(moduleName)
            clsList = inspect.getmembers(module, inspect.isclass)
            
            for name, cls in clsList:
                if issubclass(cls, CCP4Data.CData):
                    self.clsLookup[name] = {'class':cls,'clsModule':cls.__module__,'clsName':cls.__name__}
                    self.toUpper[name.lower()] = name
            if CONFIG().graphical:
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
        cacheDir = os.path.split(__file__)[0]
        cacheFilePath = os.path.join(cacheDir,"CDataManagerCache.json")

        class MyEncoder(json.JSONEncoder):
            def default(self,obj):
                #print "#CTaskManager.MyEncoder Handling case of",obj
                try:
                    return obj.__module__+obj.__name__
                except:
                    return "UnSerializable Object"
        with open(cacheFilePath,"w") as cacheFile:
            try:
                jsonText = json.dumps(compositeDict, indent=4, sort_keys=True, separators=(',', ': '), cls=MyEncoder)
                cacheFile.write(jsonText)
            except Exception as e:
                print(e)

    def printLookup(self):
        classNameList = list(self.clsLookup.keys())
        classNameList.sort()
        for item in classNameList:
            print("{0:20} {1:30} {2:30}".format(item, self.clsLookup[item], self.viewLookup.get(self.clsLookup[item], 'No widget')))

    def getClass(self, className=''):
        if className in self.clsLookup:
            if self.clsLookup[className].get('class',None) is None:
                module = importlib.import_module(self.clsLookup[className]['clsModule'])
                clsName = self.clsLookup[className]['clsName']
                self.clsLookup[className]['class'] = getattr(module, clsName)
            return self.clsLookup[className]['class']
        elif className in self.toUpper:
            if self.clsLookup[self.toUpper[className]].get('class',None) is None:
                module = importlib.import_module(self.clsLookup[self.toUpper[className]]['clsModule'])
                clsName = self.clsLookup[self.toUpper[className]]['clsName']
                self.clsLookup[self.toUpper[className]]['class'] = getattr(module, clsName)
            return self.clsLookup[self.toUpper[className]]['class']
        elif className in self.widgetLookup:
            if self.widgetLookup[className].get('class',None) is None:
                module = importlib.import_module(self.widgetLookup[className]['clsModule'])
                clsName = self.widgetLookup[className]['clsName']
                self.widgetLookup[className]['class'] = getattr(module, clsName)
            return self.widgetLookup[className]['class']
        else:
            return None

    def getWidgetClass(self, model=None, modelClass=None, name=''):
        from . import CCP4Data
        if modelClass is None:
            try:
                modelClass = model.__class__
            except:
                pass
        modelClassName = modelClass.__name__
        widgetClassDict = self.viewLookup.get(modelClassName, None)
        if widgetClassDict is not None:
            if widgetClassDict.get('class',None) is None:
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
        if not CONFIG().graphical:
            return None
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
        try:
            widget = widgetClass(parent=parentWidget, model=model, qualifiers=qualis)
        except CException as e:
            raise e
        except:
            raise CException(self.__class__, 108, stack=False)
        return widget

    def buildQStandardItemModel(self, parent=None, mode=None):
        from . import CCP4Data
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
