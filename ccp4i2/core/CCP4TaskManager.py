"""
Copyright (C) 2010 University of York
Liz Potterton Sept 2010 - Class to keep track of all CCP4Tasks
"""

import glob
import importlib
import inspect
import json
import os
import pkgutil
import re
import sys
import time
import traceback
import xml.etree.ElementTree as ET

from . import CCP4Utils
from ..dbapi import CCP4DbApi
from .CCP4Config import CONFIG
from .CCP4Modules import CUSTOMTASKMANAGER
from .CCP4Modules import PREFERENCES
from .CCP4Modules import PROJECTSMANAGER
from .CCP4Modules import WORKFLOWMANAGER
from .CCP4ErrorHandling import CErrorReport, Severity


MODULE_ORDER = ['data_entry', 'data_processing', 'data_reduction', 'bigpipes', 'alpha_fold', 'expt_phasing', 'bioinformatics',
                'molecular_replacement', 'density_modification', 'model_building', 'refinement', 'ligands',
                'validation', 'export', 'expt_data_utility', 'model_data_utility', 'developer_tools',
                'preferences', 'wrappers', 'test']
MODULE_TITLES = {'data_entry'            : 'Import merged data, AU contents, alignments or coordinates',
                 'bigpipes'              : 'Data to complete structure solution including ligand fitting',
                 'data_processing'       : 'Integrate X-ray images',
                 'data_reduction'        : 'X-ray data reduction and analysis',
                 'alpha_fold'            : 'AlphaFold and RoseTTAFold Utilities',
                 'expt_phasing'          : 'Experimental phasing',
                 'bioinformatics'        : 'Bioinformatics including model preparation for Molecular Replacement',
                 'molecular_replacement' : 'Molecular Replacement',
                 'density_modification'  : 'Density modification',
                 'model_building'        : 'Model building and Graphics',
                 'refinement'            : 'Refinement',
                 'ligands'               : 'Ligands',
                 'validation'            : 'Validation and analysis',
                 'export'                : 'Export and Deposition',
                 'expt_data_utility'     : 'Reflection data tools',
                 'model_data_utility'    : 'Coordinate data tools',
                 'developer_tools'       : 'Developer tools',
                 'preferences'           : 'Preferences',
                 'wrappers'              : 'Wrappers',
                 'demo'                  : 'Demo code for developers only',
                 'test'                  : 'Test code for developers only'}
                #'refln_data_analysis' : 'Reflection data analysis',
MODULE_DEFAULTS = {'molecular_replacement' : ['mrbump_basic', 'phaser_simple', 'phaser_pipeline', 'molrep_pipe', 'molrep_den','csymmatch', 'parrot', 'phaser_rnp_pipeline' ],
                   'data_reduction'        : ['aimless_pipe', 'freerflag', 'matthews', 'molrep_selfrot'],
                   'data_entry'            : ['import_merged', 'ProvideAsuContents', 'ProvideSequence', 'ProvideAlignment'],
                   'data_processing'       : ['xia2_dials', 'xia2_xds', 'imosflm'] ,
                   'expt_phasing'          : ['crank2', 'shelx', 'phaser_EP_AUTO'],
                   'alpha_fold'            : ['ccp4mg_edit_model','mrparse','editbfac'],
                   'refinement'            : ['servalcat_pipe', 'prosmart_refmac'],
                   'bioinformatics'        : ['ccp4mg_edit_model' , 'ccp4mg_edit_nomrbump', 'chainsaw', 'sculptor','phaser_ensembler','clustalw'],
                   'bigpipes'              : ['SubstituteLigand','dr_mr_modelbuild_pipeline'],
                   'model_building'        : ['modelcraft', 'coot_rebuild','coot_script_lines', 'coot_find_waters','arp_warp_classic','shelxeMR','dr_mr_modelbuild_pipeline'],
                   'validation'            : ['validate_protein', 'edstats', 'privateer', 'qtpisa'],
                   'export'                : ['PrepareDeposit'],
                   'expt_data_utility'     : ['pointless_reindexToMatch', 'phaser_EP_LLG', 'cmapcoeff', 'chltofom', 'cphasematch', 'ctruncate', 'splitMtz', 'scaleit'],
                   'model_data_utility'    : ['csymmatch', 'gesamt', 'coordinate_selector', 'qtpisa'],
                   'developer_tools'       : [],
                   'test'                  : ['demo_copycell']}

# The taskname (the key) maps to different directory and def.xml filename
TASKNAME2DIR = {'refmac' : 'refmac_i2'}

def TASKMANAGER():
    if CTaskManager.insts is None:
        CTaskManager.insts = CTaskManager()
    return CTaskManager.insts

def LISTTASKS():
    tm = TASKMANAGER()
    text = '\n\nThe current menu\n\n'
    tree = tm.taskTree()
    for moduleName,moduleTitle,taskList in tree:
        text += "{0:30} {1:50}\n".format(moduleName, moduleTitle)
        for taskName in taskList:
            if tm.getTaskAttribute(taskName, 'TASKTITLE') is not None:
                text += "                  {0:30} {1:50}\n".format(taskName, tm.getTaskAttribute(taskName, 'TASKTITLE'))
            else:
                text += "                  {0:30} {1:50}\n".format(taskName, "None")
    text += tm.printLookup()
    return text

class CTaskManager:

    ERROR_CODES = {101 : {'description' : 'No definintion (def.xml) file for task'},
                   102 : {'description' : 'Undefined error opening def file in task viewer'},
                   103 : {'description' : 'Undefined error loading task viewer in browser'},
                   104 : {'description' : 'Attempting to open task viewer in non-graphical mode'},
                   105 : {'description' : 'taskName not recognised in whatNext'},
                   106 : {'description' : 'whatNext error importing module'},
                   107 : {'description' : 'whatNext error calling module whatNext'},
                   108 : {'description' : 'Error importing Python module'},
                   120 : {'description' : 'Error loading local task customisation file'}}

    insts = None
    def __init__(self):
        self.timeFindingIcons = 0
        self.timeFindingDefXML = 0
        self._guiSearchPath = [os.path.join(CCP4Utils.getCCP4I2Dir(), 'tasks', '*')]
        self._searchPath = [os.path.join(CCP4Utils.getCCP4I2Dir(), 'wrappers', '*', 'script'), os.path.join(CCP4Utils.getCCP4I2Dir(), 'wrappers2', '*', 'script'), os.path.join(CCP4Utils.getCCP4I2Dir(), 'pipelines', '*', 'script'), os.path.join(CCP4Utils.getCCP4I2Dir(), 'pipelines', '*', 'wrappers', '*', 'script'), os.path.join(CCP4Utils.getCCP4I2Dir(), 'deprecated')]
        self.taskLookup= {}
        self.moduleLookup= {}
        self.scriptLookup = {}
        self.reportLookup = {}
        self.performanceClassLookup = {}
        self.taskPerformanceClassLookup = {}
        self.i1TaskLookup = None
        #Local customisations
        self.blockLocal = []
        errReport = self.buildLookup()
        # Do i2 task title lookup separately on demand
        if len(errReport) > 0:
            print('ERROR LOADING TASK MODULES')
            print(errReport.report())
        errReport = self.loadCustomisation()
        if len(errReport) > 0:
            print('ERROR LOADING LOCAL TASK CUSTOMISATIONS')
            print(errReport.report())

    def searchPath(self):
        return self._guiSearchPath + self._searchPath

    def appendSearchPath(self, path):
        # Append a path which could include wildcard * - must work with glob
        path = os.path.abspath(path)
        if not self._searchPath.count(path):
            self._searchPath.append(path)

    def buildLookup(self):
        try:
            return self.loadCachedClassLookups()
        except:
            print("Falling back to redetermining lookup")
            return self.buildLookupFromScratch()

    def loadCachedClassLookups(self):
        thisDir = os.path.split(__file__)[0]
        jsonPath = os.path.join(thisDir, "CachedLookups.json")
        try:
            with open(jsonPath,"r") as jsonFile:
                jsonText = jsonFile.read()
                compositeDict = json.loads(jsonText)
                self.taskLookup = compositeDict["tasks"]
                self.scriptLookup = compositeDict["scripts"]
                #Issue here is that this is not a full dictionary
                #so only get the module name squished into the
                #place wheres the class should be stored
                self.taskPerformanceClassLookup = compositeDict["taskPerformanceClasses"]
                self.reportLookup = compositeDict["reports"]
                self.moduleLookup = compositeDict["modules"]
                self.taskIconLookup = compositeDict["taskIconLookup"]
                #Do the imports for task dictionaties to produce the class mappings
                chkset = set()
                storeinfo = []
                for iTable, lookupTable in enumerate([self.taskLookup, self.scriptLookup, self.reportLookup]):
                    for entityName, versions in lookupTable.items():
                        for version, versionDict in versions.items():
                            clsNam = versionDict.get("clsName")
                            tskNam = versionDict.get("taskName")
                            modNam = versionDict.get("clsModule")
                            storeinfo.append([clsNam, tskNam, modNam])
                            # print(version, versionDict)
                            if clsNam not in chkset : # and clsNam != "crank2_phas_report" : # second part to fake failure
                                chkset.add(clsNam)
                            else:
                                print ("---------------------")
                                print("SERIOUS ISSUE: Trying to add duplicate ", clsNam)
                                print ("Duplicate Wrapper Defined in JSON for", clsNam)
                                tmpM = "Potentially serious duplication problem flagged.\n" + \
                                       "    " +clsNam + "\n"
                                tmpM += "Recommend checking the following wrappers : \n"
                                for x in storeinfo :
                                    if clsNam == x[0]:
                                        print (x[0], x[1], x[2])
                                        if x[1] is None:
                                            tmpM += "    " + str(x[2]) + "\n"
                                        else:
                                            tmpM += "    " + str(x[1]) + "\n"
                                #print (clsNam, chkset)
                                print ("---------------------")
                            #Don't bother to load class, but set it to "None", obliging it to be loaded when needed
                            versionDict["class"] = None
                #Set the class to none on the performance classses, so as to force lazy load
                for taskName, performanceClassDict in self.taskPerformanceClassLookup.items():
                    performanceClassDict["class"] = None
            return CErrorReport()
        except:
            raise

    def explore_package(self, module_name, module_name_list):
        loader = pkgutil.get_loader(module_name)
        fname = os.path.dirname(loader.get_filename())
        for sub_module in pkgutil.walk_packages([fname], prefix=module_name+"."):
            _, sub_module_name, _ = sub_module
            module_name_list.append(sub_module_name)

    def buildLookupFromScratch(self):
        print("#CCP4TaskManager.buildLookupFromScratch")
        from . import CCP4PerformanceData
        from . import CCP4PluginScript
        from ..report import CCP4ReportParser
        graphical = CONFIG().graphical
        myErrorReport = CErrorReport()
        loadTarget = [[self._searchPath, True]]
        if graphical:
            from ..qtgui import CCP4TaskWidget
            loadTarget.insert(0, [self._guiSearchPath, False])
        # Get the performance data classes
        clsList = inspect.getmembers(CCP4PerformanceData, inspect.isclass)
        for className,cls in clsList:
            if issubclass(cls, CCP4PerformanceData.CPerformanceIndicator) and not cls == CCP4PerformanceData.CPerformanceIndicator:
                self.performanceClassLookup[className] = {"class":cls, "clsName":cls.__name__, "clsModule":cls.__module__}
        moduleLookup = {'wrappers' : []}
        guiedTasks = []
        #Toplevel
        module_name_list = []
        self.explore_package("wrappers", module_name_list)
        self.explore_package("wrappers2", module_name_list)
        self.explore_package("pipelines", module_name_list)
        module_name_list = [module_name for module_name in
                            module_name_list if "crank2.crank2" not in module_name and 'PdbView' not in module_name]

        for module_name in module_name_list:
            pyFile = None
            mymodule = None
            try:
                #print("Trying",module_name); sys.stdout.flush()
                mymodule = CCP4Utils.importModule(module_name)
                pyFile = mymodule.__file__
            except:
                print("Failed to import module", module_name)
                errtype, value, traceback = sys.exc_info()
                print(errtype)
                print(value)
            if pyFile is not None:
                pyFile = mymodule.__file__
                if True:
                    clsList = inspect.getmembers(mymodule, inspect.isclass)
                    for className, cls in clsList:
                        if graphical and issubclass(cls, CCP4TaskWidget.CTaskWidget) and not cls == CCP4TaskWidget.CTaskWidget:
                            # The taskName used for internal reference is the name of directory
                            # containing CTask*.py and *.def.xml files
                            taskName = getattr(cls, 'TASKNAME', None)
                            if taskName is None:
                                taskName = os.path.split(os.path.split(pyFile)[0])[1]
                            taskVersion = getattr(cls, 'TASKVERSION', '0.0.0')
                            guiName = getattr(cls, 'GUINAME', taskName)
                            guiVersion = getattr(cls, 'GUIVERSION', '0.0.0')
                            taskModuleList = getattr(cls, 'TASKMODULE', 'test')
                            maintainer = getattr(cls, 'MAINTAINER', 'Nobody')
                            shortTitle = getattr(cls, 'SHORTTASKTITLE', None)
                            taskTitle = getattr(cls, 'TASKTITLE', None)
                            description = getattr(cls, 'DESCRIPTION', None)
                            rank = getattr(cls, 'RANK', None)
                            if shortTitle is None:
                                shortTitle = getattr(cls, 'TASKTITLE', None)
                            if not isinstance(taskModuleList, list):
                                taskModuleList = [taskModuleList]
                            if 'guiName' not in self.taskLookup:
                                self.taskLookup[guiName] = {}
                            self.taskLookup[guiName][guiVersion] = {'class' : cls, 'taskName' : taskName,
                                                                    'taskVersion' : taskVersion,
                                                                    'MAINTAINER' : maintainer, 'TASKTITLE':taskTitle,
                                                                    'DESCRIPTION': description, 'shortTitle' : shortTitle,
                                                                    'RANK': rank,
                                                                    'clsName':cls.__name__,'clsModule':cls.__module__}
                            for taskModule in taskModuleList:
                                if taskModule not in moduleLookup:
                                    moduleLookup[taskModule] = []
                                if not guiName in moduleLookup[taskModule]:
                                    moduleLookup[taskModule].append(guiName)
                            guiedTasks.append(taskName)
                        elif issubclass(cls,CCP4PluginScript.CPluginScript) and not cls==CCP4PluginScript.CPluginScript:
                            taskModuleList = getattr(cls, 'TASKMODULE', None)
                            if taskModuleList != 'deprecated':
                                taskName = getattr(cls, 'TASKNAME', None)
                                if taskName is None:
                                    taskName = os.path.split(os.path.split(os.path.split(pyFile)[0])[0])[1]
                                taskVersion = getattr(cls, 'TASKVERSION', '0.0.0')
                                maintainer = getattr(cls, 'MAINTAINER', 'Nobody')
                                taskTitle = getattr(cls, 'TASKTITLE', None)
                                description = getattr(cls, 'DESCRIPTION', None)
                                interruptLabel = getattr(cls, 'INTERRUPTLABEL', None)
                                if not isinstance(taskModuleList, list):
                                    taskModuleList = [taskModuleList]
                                subTasks = getattr(cls, 'SUBTASKS', [])
                                #print 'buildLookup CPluginScript',taskName,taskModule
                                internal = issubclass(cls, CCP4PluginScript.CInternalPlugin)
                                if 'taskName' not in self.scriptLookup:
                                    self.scriptLookup[taskName] = {}
                                self.scriptLookup[taskName][taskVersion] = {'class' : cls, 'subTasks': subTasks,
                                                                            'internal' : internal,
                                                                            'TASKTITLE':taskTitle,'DESCRIPTION':description,
                                                                            'MAINTAINER' : maintainer, 'INTERRUPTLABEL' : interruptLabel,
                                                                        'clsName':cls.__name__,'clsModule':cls.__module__}
                            if cls.PERFORMANCECLASS is not None and cls.PERFORMANCECLASS in self.performanceClassLookup:
                                self.taskPerformanceClassLookup[taskName] = self.performanceClassLookup[cls.PERFORMANCECLASS]
                            for taskModule in taskModuleList:
                                if taskModule is not None and not taskName in guiedTasks:
                                    if taskModule not in moduleLookup:
                                        moduleLookup[taskModule] = []
                                    if not taskName in moduleLookup[taskModule]:
                                        moduleLookup[taskModule].append(taskName)
                            #else:
                            #  if not taskName in moduleLookup['wrappers'] and not taskName in ['customtask','workflow']:
                            #    moduleLookup['wrappers'].append(taskName)
                        elif issubclass(cls, CCP4ReportParser.Report) and not cls == CCP4ReportParser.Report:
                            taskName = getattr(cls, 'TASKNAME', None)
                            taskVersion = getattr(cls, 'TASKVERSION', '0.0.0')
                            maintainer = getattr(cls, 'MAINTAINER', 'Nobody')
                            if taskName is None:
                                taskName = re.sub('_report', '', cls.__name__)
                            if 'taskName' not in self.reportLookup:
                                self.reportLookup[taskName] = {}
                            self.reportLookup[taskName][taskVersion] = {'class' : cls,
                                                                        'MAINTAINER' : maintainer, 'modes' : [],
                                                                        'clsName':cls.__name__,'clsModule':cls.__module__}
                            try:
                                if cls.RUNNING:
                                    self.reportLookup[taskName][taskVersion]['modes'].append('Running')
                                if cls.FAILED:
                                    self.reportLookup[taskName][taskVersion]['modes'].append('Failed')
                            except:
                                pass
        # Sort the task order according to MODULE_DEFAULTS - allow a wrapper to be
        # included if it is listed in the MODULE_DEFAULTS
        for module,taskList in list(moduleLookup.items()):
            self.moduleLookup[module] = []
            if module in MODULE_DEFAULTS:
                for item in MODULE_DEFAULTS[module]:
                    if item in taskList:
                        self.moduleLookup[module].append(item)
                    elif item in moduleLookup['wrappers']:
                        self.moduleLookup[module].append(item)
                for item in taskList:
                    if not item in self.moduleLookup[module]:
                        self.moduleLookup[module].append(item)
            else:
                self.moduleLookup[module].extend(taskList)
        #---------------------------------------------------------------
        class MyEncoder(json.JSONEncoder):
            def default(self, obj):
                try:
                    return obj.__module__+obj.__name__
                except:
                    return "UnSerializable Object"
        #---------------------------------------------------------------
        if not hasattr(self, "taskIconLookup"):
            self.taskIconLookup = {}
        for task in self.taskLookup:
            iconPath = self.searchIconFileNotInLookup(task)
            ccp4i2Root = CCP4Utils.getCCP4I2Dir()
            if str(iconPath).startswith(ccp4i2Root):
                self.taskIconLookup[task] = iconPath[len(ccp4i2Root)+1:]
        for moduleName in MODULE_ORDER:
            iconPath = self.searchIconFileNotInLookup(moduleName)
            ccp4i2Root = CCP4Utils.getCCP4I2Dir()
            if str(iconPath).startswith(ccp4i2Root):
                self.taskIconLookup[moduleName] = iconPath[len(ccp4i2Root)+1:]

        compositeDict = {"scripts":self.scriptLookup,
                         "tasks":self.taskLookup,
                         "taskPerformanceClasses":self.taskPerformanceClassLookup,
                         "reports":self.reportLookup,
                         "modules":self.moduleLookup,
                         "performanceClasses":self.performanceClassLookup,
                         "taskIconLookup":self.taskIconLookup}

        jsonBlob = json.dumps(compositeDict, cls=MyEncoder, indent=4, sort_keys=True, separators=(',', ': '), ensure_ascii=False)
        thisDir = os.path.split(__file__)[0]
        jsonPath = os.path.join(thisDir, "CachedLookups.json")
        with open(jsonPath,"w") as jsonFile:
            jsonFile.write(jsonBlob)
        return myErrorReport

    def loadCustomisation(self):
        err = CErrorReport()
        customFile = os.path.join(CCP4Utils.getCCP4I2Dir(), 'custom_code', 'task_customisation.xml')
        if not os.path.exists(customFile):
            return err
        try:
            from . import CCP4File
            xmlFile = CCP4File.CI2XmlDataFile(customFile)
            body = xmlFile.getBodyEtree()
        except:
            err.append(self.__class__, 120, customFile)
            return err
        self.blockLocal = []
        for node in body:
            if node.tag in ['blocklocal', 'blockLocal']:
                for item in node:
                    if item.tag in ['taskname', 'taskName']:
                        self.blockLocal.append(str(item.text))
        return err

    def printLookup(self):
        i2dir = CCP4Utils.getCCP4I2Dir()
        text = '\n\nGui (CTaskWidget) list:\n\n'
        classNameList = list(self.taskLookup.keys())
        classNameList.sort()
        for task in classNameList:
            for version in list(self.taskLookup[task].keys()):
                module = CCP4Utils.importModule(self.taskLookup[task][version]['clsModule'])
                pyFile = module.__file__
                path = os.path.relpath(pyFile, i2dir)
                maintainer = self.taskLookup[task][version]['MAINTAINER']
                if maintainer is None:
                    maintainer = "None"
                taskModuleName = self.getTaskAttribute(task,'TASKMODULE')
                if taskModuleName is None:
                    taskModuleName = "None"
                else:
                    taskModuleName = str(taskModuleName)
                taskModuleTitle = self.getTaskAttribute(task,'TASKTITLE')
                if taskModuleTitle is None:
                    taskModuleTitle = "None"
                text += "{0:30} {4:20} {5:30} {2}\n".format(task, version, path, maintainer, taskModuleName, taskModuleTitle)
        text += '\n\nScript (CPLuginScript) list:\n\n'
        classNameList = list(self.scriptLookup.keys())
        classNameList.sort()
        for task in classNameList:
            for version in list(self.scriptLookup[task].keys()):
                module = CCP4Utils.importModule(self.scriptLookup[task][version]['clsModule'])
                pyFile = module.__file__
                path = os.path.relpath(pyFile, i2dir)
                text += "{0:30} {1:10} {3:30} {2}\n".format(task, version, path, self.scriptLookup[task][version]['MAINTAINER'])
        text += '\n\nReport list:\n\n'
        classNameList = list(self.reportLookup.keys())
        classNameList.sort()
        for task in classNameList:
            for version in list(self.reportLookup[task].keys()):
                module = CCP4Utils.importModule(self.reportLookup[task][version]['clsModule'])
                pyFile = module.__file__
                path = os.path.relpath(pyFile, i2dir)
                text += "{0:30} {1:10} {3:30} {2}\n".format(task, version, path, self.reportLookup[task][version]['MAINTAINER'], str(self.reportLookup[task][version]['modes']))
        return text

    def sortVersions(self, inList):
        return inList

    def getTaskData(self, name='', version=None):
        if name in self.taskLookup:
            versions = self.sortVersions(list(self.taskLookup[name].keys()))
            if version is not None and version != '0.0.0':
                if version in versions:
                    return self.taskLookup[name][version]
                else:
                    return self.taskLookup[name][versions[0]]
            else:
                return self.taskLookup[name][versions[0]]
        return {}

    def getScriptData(self, name='', version=None):
        if name in self.scriptLookup:
            versions = self.sortVersions(list(self.scriptLookup[name].keys()))
            if version is not None and version != '0.0.0':
                if version in versions:
                    return self.scriptLookup[name][version]
                else:
                    return self.scriptLookup[name][versions[0]]
            else:
                return self.scriptLookup[name][versions[0]]
        return {}

    def getReportData(self, name='', version=None):
        rv = {}
        if name in self.reportLookup:
            versions = self.sortVersions(list(self.reportLookup[name].keys()))
            if version is not None and version != '0.0.0':
                if version in versions:
                    rv = self.reportLookup[name][version]
                    rv['version'] = version
            else:
                rv = self.reportLookup[name][versions[0]]
                rv['version'] = versions[0]
        return rv

    def getTaskWidgetClass(self, name='', version=None):
        data = self.getTaskData(name, version)
        if self.lazyLoadClassForDict(data) is not None:
            return data['class']
        if name in WORKFLOWMANAGER().getList():
            data = self.getTaskData('workflow')
            return self.lazyLoadClassForDict(data)
        else:
            return None

    def getReportClass(self, name='', version=None, jobStatus=None, doReload=False):
        from ..report import CCP4ReportParser
        data = self.getReportData(name, version)
        if jobStatus == 'Running remotely':
            jobStatus = 'Running'
        if jobStatus is None or jobStatus in CCP4DbApi.FINISHED_JOB_STATUS or jobStatus in data.get('modes',[]):
            correspondingClass = self.lazyLoadClassForDict(data)
            if doReload:
                # Attempt to reload the module
                cls = self.lazyLoadClassForDict(data)
                if cls is None:
                    return None
                module = __import__(cls.__module__)
                module = importlib.reload(module)
                clsList = inspect.getmembers(module, inspect.isclass)
                for className, cls1 in clsList:
                    if issubclass(cls1, CCP4ReportParser.Report) and getattr(cls1, 'TASKNAME') == name and getattr(cls1, 'TASKVERSION') == data['version']:
                        self.reportLookup[name][data['version']]['class'] = cls1
                        print('Reloading task report module for', className)
                        return cls1
            # Return the currently imported cls even if reload has failed
            return self.lazyLoadClassForDict(data)
        else:
            return None

    def lazyLoadClassForDict(self, data):
        if data.get('class', None) is not None:
            #print '#CTaskManager lazy Load of already loaded', data['class'].__name__
            return data['class']
        moduleName = data.get('clsModule', None)
        if moduleName is None:
            return None
        className = data.get('clsName',None)
        if className is None:
            return None
        try:
            importedModule = CCP4Utils.importModule(moduleName)
        except:
            traceback.print_exc()
            importedModule = None
            
        if importedModule is None: return None
        correspondingClass = getattr(importedModule,className,None)
        data['class'] = correspondingClass
        #print '#CTaskManager lazy Load of class', data['class'].__name__
        #traceback.print_stack(limit=3)
        return correspondingClass

    def getPluginScriptClass(self, name='', version=None):
        #print 'getPluginScriptClass',self.scriptLookup.keys(),self.scriptLookup.has_key(name)
        data = self.getScriptData(name, version)
        #print 'getPluginScriptClass',
        cls = self.lazyLoadClassForDict(data)
        if cls is not None:
            return cls
        # Try if is a 'derived' gui
        data = self.getTaskData(name, version)
        pluginName = data.get('taskName', None)
        #print 'getPluginScriptClass ',name,data,pluginName
        if pluginName is not None and pluginName!= name:
            data = self.getScriptData(pluginName, data.get('taskVersion', None))
            cls = self.lazyLoadClassForDict(data)
            if cls is not None:
                return cls
        if name in WORKFLOWMANAGER().getList():
            data = self.getScriptData('workflow', version)
            cls = self.lazyLoadClassForDict(data)
            return cls
        elif name in CUSTOMTASKMANAGER().getList():
            data = self.getScriptData('customtask', version)
            cls = self.lazyLoadClassForDict(data)
            return cls
        else:
            return None

    def getClass(self, name='', version=None):
        data = self.getTaskData(name, version)
        if len(data) == 0:
            data = self.getScriptData(name, version)
        return self.lazyLoadClassForDict(data)

    def isInternalPlugin(self, name='', version=None):
        data = self.getScriptData(name, version)
        return data.get('internal', False)

    def getTitle(self, taskName='', version=None):
        #traceback.print_stack()
        #print '\n#CTaskManager.getTitle', taskName, version
        if taskName is None:
            return 'UNKNOWN'
        taskTitle = None    # KJS : Look at further
        data = self.getTaskData(taskName, version)
        taskTitle = data.get('TASKTITLE', None)
        if taskTitle is None:
            data = self.getScriptData(taskName,version)
            taskTitle = data.get('TASKTITLE', None)
        if taskTitle is not None:
            return taskTitle
        else:
            # Try customisations - will return same input taskName if not recognised
            title = WORKFLOWMANAGER().getTitle(taskName)
            if title == taskName:
                title = CUSTOMTASKMANAGER().getTitle(taskName)
            if title is None:
                title = taskName
            return title

    def getShortTitle(self, taskName='', version=None, substitute=True):
        if taskName is None:
            return 'UNKNOWN'
        shortTitle = self.getTaskData(taskName, version).get('shortTitle', None)
        if shortTitle is not None:
            return shortTitle
        else:
            return taskName

    def getAllShortTitles(self):
        ''' Return list if taskName,shortTitle'''
        ret = []
        for taskName, versionDict in list(self.taskLookup.items()):
            for version, infoDict in list(versionDict.items()):
                ret.append((taskName,infoDict.get('shortTitle', taskName)))
                break
        return ret

    def getTaskLabel(self, taskName='', version=None, substitute=True):
        print('CTaskManager.getTaskLabel', taskName, version)
        if taskName is None:
            return 'UNKNOWN'
        label = None
        data = self.getTaskData(taskName,version)
        cls = self.lazyLoadClassForDict(data)
        if cls is not None:
            label = getattr(cls, 'TASKLABEL', None)
        if label is not None:
            return label
        else:
            return taskName

    def getTaskAttribute(self, taskName, attribute, version=None, default=None, script=False):
        if attribute in ['INTERRUPTLABEL','MAINTAINER','internal','subTasks']:
            return self.getScriptData(name=taskName, version=version).get(attribute, None)
        if attribute in ['blockLocal']:
            return (taskName in self.blockLocal)
        if script:
            data = self.getScriptData(taskName, version)
        else:
            data = self.getTaskData(taskName, version)
        if data is not None:
            if attribute in data:
                return data.get(attribute, default)
            else:
                cls = self.lazyLoadClassForDict(data)
                if hasattr(cls, attribute):
                    return getattr(cls, attribute)
                else:
                    return default
        else:
            return default

    def getReportAttribute(self, taskName, attribute, version=None, default=None):
        cls = self.getReportClass(taskName, version)
        print('CTaskManager.getReportAttribute',taskName,attribute,cls)
        if cls is not None:
            return getattr(cls, attribute, default)
        else:
            return default

    def getPerformanceClass(self, taskName):
        performanceClassDict = self.taskPerformanceClassLookup.get(taskName, None)
        if performanceClassDict is None:
            return None
        cls = self.lazyLoadClassForDict(performanceClassDict)
        return self.lazyLoadClassForDict(performanceClassDict)

    def taskTree(self, shortTitles=False):
        # Return list of [module_name,taskList]
        tree = []
        for module in MODULE_ORDER:
            if module not in ['preferences', 'test', 'demo', 'wrappers', 'deprecated']:
                if module in self.moduleLookup:
                    taskList = []
                    if shortTitles:
                        for item in self.moduleLookup[module]:
                            taskList.append([item,self.getShortTitle(item)])
                    else:
                        taskList.extend(self.moduleLookup[module])
                    tree.append([module, MODULE_TITLES.get(module, module), taskList])
        if bool(PREFERENCES().SHOW_WRAPPERS):
            for module in ['wrappers', 'test', 'demo']:
                if module in self.moduleLookup:
                    taskList = []
                    if shortTitles:
                        for item in self.moduleLookup[module]:
                            taskList.append([item, self.getShortTitle(item)])
                    else:
                        taskList.extend(self.moduleLookup[module])
                    tree.append([module, MODULE_TITLES.get(module, module), taskList])
        return tree

    def searchDefFile(self, name=None, version=None):
        if name is None:
            return None
        defFileList = CCP4Utils.globSearchPath(self.searchPath(), os.path.join(name + '.def.xml'))
        if len(defFileList) > 0:
            return defFileList[0]
        pluginName = self.getTaskAttribute(taskName=name, attribute='USEPLUGIN')
        if pluginName is not None:
            defFileList = CCP4Utils.globSearchPath(self.searchPath(), os.path.join(pluginName + '.def.xml'))
            if len(defFileList) > 0:
                return defFileList[0]
        return None

    def searchXrtFile(self, name='', version=None, jobStatus=None):
        for path in self.searchPath():
            if jobStatus == 'Running':
                fileName = os.path.join(re.sub(r'\*', name, path), name + '.running.xrt')
                if os.path.exists(fileName):
                    return fileName
            else:
                fileName = os.path.join(re.sub(r'\*', name, path), name + '.xrt')
                if os.path.exists(fileName):
                    return fileName
        return None

    def searchReferenceFile(self, name=None, version=None, cformat='medline', drillDown=False):
        if name is None:
            return []
        defFileList = CCP4Utils.globSearchPath(self.searchPath(), name + '.' + cformat + '.txt')
        if drillDown:
            subTaskList = self.getScriptData(name=name, version=version).get('subTasks', [])
            #print 'CTaskManager.searchReferenceFile subTaskList',name,subTaskList
            for subTask in subTaskList:
                defFileList.extend(self.searchReferenceFile(name=subTask, cformat=cformat, drillDown=True))
        return defFileList


    def searchHelpFile(self, name=None, version=None):
        if name is None:
            return None
        sph_pth = os.path.join(CCP4Utils.getCCP4I2Dir(), 'docs', 'sphinx', 'build', 'html')
        #helpFileList = CCP4Utils.globSearchPath(self.searchPath(),os.path.join(name+'.html'))
        #if len(helpFileList)>0: return helpFileList[0]
        helpFile = os.path.join(sph_pth, 'tasks', name, 'index.html')
        if os.path.exists(helpFile):
            return helpFile
        helpFile = os.path.join(sph_pth, 'tasks', name, name + '.html')
        if os.path.exists(helpFile):
            return helpFile
        helpFile = os.path.join(CCP4Utils.getCCP4I2Dir(), 'pipelines', name, 'docs',
                                'sphinx', 'build', 'html', 'index.html')
        if os.path.exists(helpFile):
            return helpFile
        return None

    def searchIconFile(self, name=None):
        if not hasattr(self, "taskIconLookup"):
            self.taskIconLookup = {}
        if name in self.taskIconLookup:
            ccp4i2Root = CCP4Utils.getCCP4I2Dir()
            iconPath = os.path.join(ccp4i2Root, self.taskIconLookup[name])
        else:
            ccp4i2Root = CCP4Utils.getCCP4I2Dir()
            iconPath = self.searchIconFileNotInLookup(name)
            if iconPath.startswith(ccp4i2Root):
                self.taskIconLookup[name] = iconPath[len(ccp4i2Root)+1:]
        return iconPath
        #return self.taskIconLookup[name] # KJS: Remove bogus double return

    def searchIconFileNotInLookup(self, name=None):
        t1 = time.time()
        if name is not None:
            for ext in ['.svg', '.png', '.ico']:
                pixFile = os.path.join(CCP4Utils.getCCP4I2Dir(), 'qticons', name + ext)
                if os.path.exists(pixFile):
                    self.timeFindingIcons += (time.time() - t1)
                    return pixFile
        if name is not None:
            for ext in ['.svg', '.png', '.ico']:
                pixFileList = CCP4Utils.globSearchPath(self.searchPath(),os.path.join(name + ext))
                if len(pixFileList) > 0:
                    self.timeFindingIcons += (time.time() - t1)
                    return pixFileList[0]
        if CONFIG().developer: print("Icon not found", name)
        self.timeFindingIcons += (time.time() - t1)
        return os.path.join(CCP4Utils.getCCP4I2Dir(), 'qticons', 'ccp4.png')

    def hasReportDefinition(self, name='', jobStatus=None):
        if jobStatus == 'Running remotely':
            jobStatus = 'Running'
        if self.getReportClass(name=name, jobStatus=jobStatus) is not None:
            return True
        elif self.searchXrtFile(name=name, jobStatus=jobStatus) is not None:
            return True
        defFile = WORKFLOWMANAGER().getDefFile(name)
        return (defFile is not None)

    def lookupDefFile(self, name='', version=None):
        if version == '':
            version = None

        data = self.getTaskData(name, version)
        if data.get('clsModule', None) is not None:
            module = CCP4Utils.importModule(data.get('clsModule', None))
            if module is not None:
                defFile = os.path.join(os.path.split(module.__file__)[0], name + '.def.xml')
                if os.path.exists(defFile):
                    return defFile

        data = self.getScriptData(name, version)
        if data.get('clsModule', None) is not None:
            module = CCP4Utils.importModule(data.get('clsModule', None))
            if module is not None:
                defFile = os.path.join(os.path.split(module.__file__)[0], name + '.def.xml')
                if os.path.exists(defFile):
                    return defFile
        # No py files so try just simple directory names = task name
        # check it is not one of the mismatched tasname/dir
        name = TASKNAME2DIR.get(name, name)
        defFile = os.path.join(CCP4Utils.getCCP4I2Dir(), 'tasks', name, name + '.def.xml')
        if os.path.exists(defFile):
            return defFile
        for d in ['wrappers', 'pipelines', 'wrappers2']:
            defFile = os.path.join(CCP4Utils.getCCP4I2Dir(), d, name, 'script', name + '.def.xml')
            if os.path.exists(defFile):
                return defFile
        fileList = glob.glob(os.path.join(CCP4Utils.getCCP4I2Dir(), 'pipelines', '*', 'wrappers', name,'script', name + '.def.xml'))
        if len(fileList) > 0:
            return fileList[0]
        # Try if it is a workflow - returns None if not recognised
        defFile = WORKFLOWMANAGER().getDefFile(name)
        if defFile is not None:
            return defFile
        # Try derived name - used by crank2 to load the crank2 parent for subjobs crank2_phas, crank2_comb etc
        if '_' in name:
            der_name = name.split('_')[0]
            defFile = os.path.join(CCP4Utils.getCCP4I2Dir(), 'pipelines', der_name, 'script', der_name + '.def.xml')
            if os.path.exists(defFile):
                return defFile
            defFile = os.path.join(CCP4Utils.getCCP4I2Dir(), 'tasks', der_name, der_name + '.def.xml')
            if os.path.exists(defFile):
                return defFile
        defFile = CUSTOMTASKMANAGER().getDefFile(name)
        return defFile

    def exportJobFiles(self, taskName=None, jobId=None, mode='menu'):
        data = self.getScriptData(taskName)
        if len(data) == 0:
            print('CTaskManager.exportMtz no script data for', taskName)
            return []
        try:
            moduleName=data['clsModule']
            mod = CCP4Utils.importModule(moduleName)
            if mod is None:
                print('CTaskManager.exportMtz no module for', taskName, data['clsModule'])
                return []
        except:
            print('CTaskManager.exportMtz exception importing  module',data['clsModule'])
            return []
        if mode == 'menu':
            if not 'exportJobFileMenu' in dir(mod):
                #print 'CTaskManager.exportMtz no exportJobFileMenu for', taskName, dir(mod)
                return []
            else:
                try:
                    return mod.exportJobFileMenu(jobId=jobId)
                except:
                    print('CTaskManager.exportJobFiles menu failed for', taskName)
                    return []
        else:
            try:
                return mod.exportJobFile(jobId=jobId, mode=mode)
            except:
                print('CTaskManager.exportJobFiles failed for', taskName)
                return []

    def exportMtzColumnLabels(self, taskName=None, jobId=None, paramNameList=[], sourceInfoList=[]):
        data = self.getTaskData(taskName)
        if len(data) == 0:
            print('CTaskManager.exportMtzColumnLabels no module for', taskName)
            return []
        try:
            mod = CCP4Utils.importModule(data['clsModule'])
        except:
            print('CTaskManager.exportMtzColumnLabels no module for', taskName)
            return []
        try:
            return mod.exportMtzColumnLabels(jobId=jobId, paramNameList=paramNameList, sourceInfoList=sourceInfoList)
        except:
            print('CTaskManager.exportMtzColumnLabels failed for', taskName)
            return []

    def whatNext(self, taskName, jobId=None):
        #print 'CCP4TaskManager.whatNext',taskName,jobId,type(jobId),self.taskLookup.has_key(taskName)
        # Get the jobs parent - the parent pipeline should provide the info for what to do next
        childJobs = [jobId]
        pipelineTaskName = taskName
        childJobNumber = None
        if jobId is None:
            parentJobId = None
        else:
            parentJobId = PROJECTSMANAGER().db().getJobInfo(jobId=jobId, mode='parentJobId')
        while parentJobId is not None:
            childJobs.insert(0, parentJobId)
            parentJobId = PROJECTSMANAGER().db().getJobInfo(jobId=parentJobId, mode='parentJobId')
        if len(childJobs) > 1:
            pipelineTaskName = PROJECTSMANAGER().db().getJobInfo(jobId=childJobs[0], mode='taskName')
            childJobNumber = PROJECTSMANAGER().db().getJobInfo(jobId=childJobs[-1], mode='jobNumber')
        #print 'CTaskManager.whatNext',taskName,jobId,childJobs,pipelineTaskName
        projectName = PROJECTSMANAGER().db().getJobInfo(jobId=jobId, mode='projectname')
        if parentJobId is not None:
            jobId = parentJobId
        whatList = []
        data = self.getTaskData(pipelineTaskName)
        if len(data) == 0:
            self.getScriptData(pipelineTaskName)
        if len(data) > 0:
            try:
                mod = CCP4Utils.importModule(data.get('clsModule'))
            except:
                pass
            else:
                try:
                    whatList = mod.whatNext(jobId=jobId, childTaskName=taskName, childJobNumber=childJobNumber, projectName=projectName)
                except:
                    #print 'CDbApi.whatNext failed in running whatNext from module'
                    pass
        if len(whatList) == 0:
            cls = self.getTaskWidgetClass(pipelineTaskName)
            if cls is not None:
                whatList = getattr(cls, 'WHATNEXT', [])
            if len(whatList)==0:
                cls = self.getPluginScriptClass(pipelineTaskName)
                if cls is not None:
                    whatList = getattr(cls, 'WHATNEXT', [])
        rv = []
        for name in whatList:
            if isinstance(name, list):
                rv.append((name[0], self.getShortTitle(name[0]), name[1]))
            else:
                rv.append((name, self.getShortTitle(name), None))
        return rv

    def getWhatNextPage(self, taskName=None, jobId=None):
        if taskName is None and jobId is not None:
            taskName = PROJECTSMANAGER().db().getJobInfo(jobId=jobId, mode='taskname')
        if taskName is None:
            return None
        #page = os.path.join(CCP4Utils.getCCP4I2Dir(),'docs','whatnext',taskName+'.html')
        page = os.path.join(CCP4Utils.getCCP4I2Dir(), 'tasks', taskName, taskName + '.whatnext.html')
        #print 'getWhatNextPage',page
        if os.path.exists(page):
            #return 'whatnext/'+taskName+'.html'
            return page
        else:
            return None

    def getI1TaskTitle(self, taskName):
        if self.i1TaskLookup is None:
            self.loadI1TaskTitles()
        return self.i1TaskLookup.get(taskName, taskName)

    def loadI1TaskTitles(self, fileName=None):
        from ..qtgui import CCP4I1Projects
        if fileName is None:
            fileName = os.path.join(CCP4Utils.getCCP4Dir(), 'share', 'ccp4i', 'etc', 'modules.def')
        metaData, params, err = CCP4I1Projects.readI1DefFile(fileName)
        if err.maxSeverity() > Severity.WARNING or 'TASK_NAME' not in params:
            return CErrorReport(self.__class__, 103, details=fileName, stack=False)
        self.i1TaskLookup = {}
        for nT in range(len(params['TASK_NAME'][1])):
            self.i1TaskLookup[params['TASK_NAME'][1][nT]] = params['TASK_TITLE'][1][nT]
        return CErrorReport()

class CMakeDocsIndex():

    ERROR_CODES = {101: {'description' : 'Error saving to file - do you have write permission?'},
                   102: {'description' : 'Error searching for task doc files'},
                   103: {'description' : 'Error creating html'}}

    def __init__(self):
        pass

    def run(self):
        from ..report import CCP4ReportParser
        # List of documented tasks
        try:
            dList = glob.glob(os.path.join(CCP4Utils.getCCP4I2Dir(), 'docs', 'sphinx', 'build', 'html', 'tasks', '*'))
            docList = []
            for f in dList:
                docList.append(os.path.split(f)[1])
        except:
            return CErrorReport(self.__class__, 102)
        try:
            # html document as etree
            docTree = CCP4ReportParser.htmlDoc(htmlBase='.', title='CCP4i2 Tasks', cssFile='task.css')
            body = docTree.getroot().xpath('./body')[0]
            #print 'makeDocsIndex',docList
            t = ET.Element('div')
            t.text = 'CCP4i2 Task Documentation'
            t.set('class', 'title')
            body.append(t)
            for module, title, taskList in TASKMANAGER().taskTree():
                modEle = ET.Element('div')
                modEle.set('class', 'module')
                imgEle = ET.Element('img')
                pixFile = TASKMANAGER().searchIconFile(module)
                if pixFile is not None:
                    imgEle.set('src', '../../qticons/' + os.path.split(pixFile)[1])
                imgEle.set('alt', module)
                modEle.append(imgEle)
                e = ET.Element('p')
                e.text = title
                modEle.append(e)
                body.append(modEle)
                listEle = ET.Element('div')
                listEle.set('class', 'taskList')
                body.append(listEle)
                for taskName in taskList:
                    rank = TASKMANAGER().getTaskAttribute(taskName, 'RANK')
                    pixFile = TASKMANAGER().searchIconFile(taskName)
                    imgEle = ET.Element('img')
                    if pixFile is not None:
                        imgEle.set('src', '../../qticons/' + os.path.split(pixFile)[1])
                    imgEle.set('alt', taskName)
                    imgEle.set('class', 'taskIcon')
                    taskEle = ET.Element('p')
                    taskEle.text = TASKMANAGER().getTitle(taskName)
                    if taskName in docList:
                        hrefEle = ET.Element('a')
                        hrefEle.set('href', './' + taskName + '/index.html')
                        hrefEle.append(taskEle)
                    else:
                        hrefEle = None
                    descEle = ET.Element('p')
                    descEle.set('class', 'taskDesc')
                    descEle.text = TASKMANAGER().getTaskAttribute(taskName, 'DESCRIPTION')
                    div1 = ET.Element('div')
                    div1.set('class', 'taskTitle')
                    div1.append(imgEle)
                    if hrefEle is not None:
                        div1.append(hrefEle)
                    else:
                        div1.append(taskEle)
                    listEle.append(div1)
                    listEle.append(descEle)
        except Exception as e:
            return CErrorReport(self.__class__, 101, str(e))
        sph_pth = os.path.join(CCP4Utils.getCCP4I2Dir(), 'docs', 'sphinx', 'build', 'html')
        try:
            ET.indent(docTree)
            text = ET.tostring(docTree)
            CCP4Utils.saveFile(os.path.join(sph_pth, 'tasks', 'index.html'), text, overwrite=True)
        except:
            return CErrorReport(self.__class__, 101, details=os.path.join(sph_pth, 'tasks', 'index.html'))
        return CErrorReport()
