from __future__ import print_function

"""
     CCP4Container.py: CCP4 GUI Project
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
   Liz Potterton Aug 2010 - Generic container for CCP4Data objects
"""
import os
import re
import types
from core import CCP4Data
from core.CCP4ErrorHandling import *
from core.CCP4DataManager import DATAMANAGER
from core.CCP4Config import DEVELOPER,XMLPARSER


class CContainer(CCP4Data.CData):
    ERROR_CODES = {101 : {'description' : 'Error parsing XML'},
                   102 : {'description' : 'Missing information'},
                   103 : {'description' : 'Unknown data class'},
                   104 : {'description' : 'Error creating data object'},
                   105 : {'description' : 'Error setting data object qualifiers'},
                   106 : {'description' : 'Error loading container definition'},
                   107 : {'description' : 'XML file does not have correct function defined in the header'},
                   108 : {'description' : 'XML undefined error interpreting sub-container'},
                   109 : {'description' : 'Error attempting to access unknown attribute','severity' : SEVERITY_WARNING },
                   110 : {'description' : 'Error creating sub-container'},
                   111 : {'description' : 'XML file does not have expected pluginName defined in the header'},
                   113 : {'description' : 'Attempting to add object that is not a CData'},
                   114 : {'description' : 'Attempting to add object without valid name'},
                   115 : {'description' : 'Attempting to add object with name that is already in container'},
                   116 : {'description' : 'Error while attempting to add object'},
                   117 : {'description' : 'Attempting to delete object with unrecognised name'},
                   118 : {'description' : 'Error while attempting to delete object'},
                   119 : {'description' : 'Error while attempting to set this container as object parent'},
                   120 : {'description' : 'Attempting to add object of unrecognised class to container contents'},
                   121 : {'description' : 'Error while attempting to add to container contents'},
                   122 : {'description' : 'Error while attempting to make object from new content in container'},
                   123 : {'description' : 'Unknown error while reading container header'},
                   124 : {'description' : 'Definition of sub-content for data of class that does not require sub-content'},
                   125 : {'description' : 'Unknown error while reading container content'},
                   126 : {'description' : 'No id for sub-container in XML file'},
                   127 : {'description' : 'Attempting to load container data from file that does not exist'},
                   128 : {'description' : 'Unknown error creating XML for sub-container'},
                   129 : {'description' : 'Error retieving data object for XML'},
                   130 : {'description' : 'Error saving data object to XML'},
                   131 : {'description' : 'Unknown error writing container contents to XML file'},
                   132 : {'description' : 'Error changing object name - no name given'},
                   133 : {'description' : 'Error changing object name - object with new name already exists'},
                   134 : {'description' : 'Error changing object name - no object with old name'},
                   135 : {'description' : 'Unknown error changing object name'},
                   136 : {'description' : 'Error inserting object in container data order'},
                   137 : {'description' : 'Unknown error restoring data from database'},
                   138 : {'description' : 'Attempting to copy from otherContainer which is not a CContainer'},
                   139 : {'severity' : SEVERITY_WARNING, 'description' : 'Attempting to copy data which is not in this container'},
                   140 : {'severity' : SEVERITY_WARNING, 'description' : 'Attempting to copy data which is not in the other container'},
                   141 : {'severity' : SEVERITY_WARNING, 'description' : 'Unknown error copying data'},
                   142 : {'description' : 'Unrecognised class name in file'},
                   143 : {'description' : 'Item in file does not have an id'},
                   144 : {'description' : 'Item id in file is not unique'},
                   145 : {'description' : 'Failed setting command line argument'},
                   146 : {'description' : 'Insufficient arguments at end of command line'},
                   147 : {'description' : 'Error handling XmlDataFile for file element in def xml'},
                   148 : {'description' : 'XmlDataFile for file element in def xml: file not found'},
                   149 : {'description' : 'XmlDataFile for file element in def xml: can not read xml'},
                   150 : {'description' : 'loadDataFromXml could not find plugin def file'},
                   160 : {'description' : 'Error in adding guiAdmin to CContainer'},
                   161 : {'description' : 'Error adding object to guiAdmin'},
                   162 : {'description' : 'Error adding guiAdmin to CContainer'},
                   310 : {'description' : 'Different number of file objects to compare'},
                   311 : {'description' : 'Different number of XData objects to compare'},
                   312 : {'description' : 'Different number of key-value pairs to compare'},
                   313 : {'description' : 'Different values of key-value pair'},
                   314 : {'description' : 'Error running comparison of object'}}

    def __init__(self, name=None, etreeElement=None, contents={}, definitionFile='', parent=None, guiAdmin=False, **kw):
        self.__dict__['CONTENTS'] = {}
        CCP4Data.CData.__init__(self, parent=parent, name=name)
        self.__dict__['_dataOrder'] = []  # Retain order of data from XML to use in auto-generated GUI
        report = CException()
        if etreeElement is not None:
            report = self.loadContentsFromEtree(etreeElement)
            if guiAdmin:
                try:
                    report.extend(self.addGuiAdminContents())
                except:
                    report.append(self.__class__, 160, str(name), name=self.objectPath())
        elif len(contents) > 0:
            report.extend(self.setContents(contents))
        elif len(definitionFile) > 0:
            report.extend(self.loadContentsFromXml(definitionFile, guiAdmin=guiAdmin))
        # Throw exception if any of these data loading methods
        # return error
        if  report.maxSeverity() >= SEVERITY_ERROR:
            myException = CException()
            myException.extend(report)
            raise myException

    def __getattr__(self, name, default=NotImplemented):
        #print 'in CData.__getattr__',name
        if name in self.CONTENTS:
            return self._value[name]
        elif name == 'header':
            return getattr(self.__dict__, 'header', None)
        elif name == 'paramsHeader':
            return getattr(self.__dict__, 'paramsHeader', None)
        elif name == 'version':
            return self.getVersion()
        if default is NotImplemented:
            raise CException(self.__class__, 109, name, name=self.objectPath())
        else:
            return default

    def contents(self, name=None):
        if name is None:
            return self.CONTENTS
        elif name in self.CONTENTS:
            return self.CONTENTS[name]
        else:
            return None

    def setContents(self,contents={}):
        #myException = CException()
        myException = CErrorReport()
        for name,defn in list(contents.items()):
            obj = None
            cls = None
            if 'class' in defn:
                cls = defn['class']
            elif 'className' in defn:
                cls = DATAMANAGER().getClass(defn['className'])
                if cls is None:
                    myException.append(self.__class__,103,'Class: '+defn['className']+' specified for: '+name,name=self.objectPath())
            else:
                myException.append(self.__class__,102,'For: '+name,name=self.objectPath())
            if cls is not None:
                try:
                    obj = cls(parent=self)
                    obj.setObjectName(name)
                except CException as e:
                    obj = None
                    myException.extend(e)
                #except:
                #  obj = None
                #  myException.append(self.__class__,104,'For: '+name,name=self.objectPath())
            if obj is not None and 'qualifiers' in defn:
                # Set the qualifiers - note that if this fails the obj is deleted
                try:
                    obj.setQualifiers(qualifiers=defn['qualifiers'])
                except CException as e:
                    myException.extend(e)
                    obj = None
                except:
                    myException.append(self.__class__, 105, 'Tag: ' + name, name=self.objectPath())
                    obj = None
            if obj is not None:
                self.__dict__['_value'][name] = obj
                self.__dict__['CONTENTS'][name] = {'class' : cls, 'qualifiers' : obj.qualifiers()}
                self.__dict__['_dataOrder'].append(name)
        # Not sure if we treat this as exception or return
        #if len(myException)>0: raise myException
        return myException

    def update(self, container=None):
        '''Update from another container'''
        for key in container.dataOrder():
            try:
                obj = self.__getattr__(key)
            except:
                pass
            else:
                try:
                    obj.set(container.__getattr__(key))
                except:
                    print('ERROR in CContainer.update for key', key)

    def loadContentsFromXml(self, fileName=None, guiAdmin=False):
        from core import CCP4File
        f = CCP4File.CI2XmlDataFile(fullPath=fileName)
        f.loadFile()
        #print 'CContainer.loadContentsFromXml', f.header.function, fileName
        if f.header.function != 'DEF':
            raise CException(self.__class__, 107, 'Filename: ' + fileName, name=self.objectPath())
        self.__dict__['header'] = f.header
        self.__dict__['header'].setParent(self)
        #print 'loadContentsFromXml', f.header.pluginName, f.header.pluginTitle,guiAdmin
        if f.header.pluginTitle is not None:
            self.setObjectName(str(f.header.pluginTitle))
        if f.header.pluginName is not None:
            self.setObjectName(str(f.header.pluginName))
        errorReport = self.loadContentsFromEtree(element=f.getBodyEtree())
        from lxml import etree
        #print 'CContainer.loadContentsFromXml after loadContentsFromEtree',self.objectName()
        if guiAdmin:
            self.addGuiAdminContents()
        return errorReport

    def addParamsSubContainers(self):
        for name in ['inputData', 'controlParameters', 'outputData']:
            self.addObject(CContainer(parent=self, name=name), name=name)

    def addGuiAdminContents(self):
        err = CErrorReport()
        if self.__dict__['_dataOrder'].count('guiAdmin'):
            return err
        from core import CCP4DataManager
        container = CContainer(name='guiAdmin')
        for id, className,qualifiers in [['followFrom', 'CFollowFromJob', {'guiLabel' : 'Follow from'}],
                                         ['jobTitle', 'CJobTitle', {'default' : '', 'guiLabel' : 'Job title'}],
                                         ['jobStatus', 'CJobStatus', {}], ['patchSelection', 'CPatchSelection', {}]]:
            cls = CCP4DataManager.DATAMANAGER().getClass(className)
            if cls is not None:
                try:
                    obj = cls(parent=container, name=id, qualifiers=qualifiers)
                    #print 'CContainer.addGuiAdminContents', id, obj
                    container.addObject(obj)
                except:
                    err.append(self.__class__, 161, id, name=self.objectPath())
        try:
            self.addObject(object=container, name='guiAdmin')
        except:
            err.append(self.__class__, 162, name=self.objectPath())
        return err

    def loadDataFromXml(self, fileName=None, guiAdmin=False, check=True, function='PARAMS', loadHeader=False):
        from core import CCP4File
        from core import CCP4TaskManager
        f = CCP4File.CI2XmlDataFile(fullPath=fileName)
        #f.loadFile()
        #print 'CContainer.loadDataFromXml',fileName,f.header.function,self.pluginName()
        if f.header.function != function:
            raise CException(self.__class__, 107, 'Filename: ' + str(fileName), name=self.objectPath())
        if self.pluginName() is None or len(self.pluginName()) == 0:
            # Try loading for the appropriate plugin def file
            pluginName = str(f.header.pluginName)
            pluginVersion = str(f.header.pluginVersion)
            try:
                defFile = CCP4TaskManager.TASKMANAGER().lookupDefFile(pluginName, pluginVersion)
            except:
                raise CException(self.__class__, 150, pluginName, name=self.objectPath())
            #print 'CContainer.loadDataFromXml loading contents from defFile',defFile
            self.loadContentsFromXml(defFile, guiAdmin=guiAdmin)
        else:
            if check and f.header.pluginName != self.pluginName():
                #print 'loadDataFromXml', f.header.pluginName, self.pluginName()
                raise CException(self.__class__, 111, 'Filename: ' + fileName, name=self.objectPath())
        if loadHeader:
            self.__dict__['header'] = CCP4File.CI2XmlHeader()
            self.__dict__['header'].set(f.header)
        else:
            self.__dict__['paramsHeader'] = CCP4File.CI2XmlHeader()
            self.__dict__['paramsHeader'].set(f.header)
        errorReport = self.setEtree(f.getBodyEtree())
        return errorReport

    def saveContentsToXml(self, fileName=None, header=None):
        from core import CCP4File
        f = CCP4File.CI2XmlDataFile(fullPath=fileName)
        if header is not None:
            f.header.set(header)
        elif self.header is not None:
            f.header.set(self.header)
        f.header.setCurrent()
        bodyEtree, errorReport = self.saveContentsToEtree()
        #print 'CContainer.saveContentsToXml', fileName, errorReport.report()
        try:
            f.saveFile(bodyEtree=bodyEtree)
        except CException as e:
            errorReport.append(e)
        except:
            errorReport.append(self.__class__, 130, fileName, name=self.objectPath())
        return errorReport

    def getEtree(self, excludeUnset=True, useLXML=True):
        name = self.objectName()
        if name.count(' ') > 0:
            name = re.sub(' ', '_', name)
        return CCP4Data.CData.getEtree(self, excludeUnset=excludeUnset, name=name)

    def saveDataToXml(self, fileName=None, subContainer=None, function='PARAMS'):
        '''Save the data to xml PARAMS file with option to save just one
        specified subContainer (used for interruptStatus)'''
        from core import CCP4File
        errorReport = CErrorReport()
        f = CCP4File.CI2XmlDataFile(fullPath=fileName)
        if self.header is not None:
            f.header.set(self.header)
        f.header.function.set(function)
        f.header.setCurrent()
        bodyEtree = self.getEtree(excludeUnset=True)
        if subContainer is not None:
            subConList = bodyEtree.xpath('./' + subContainer)
            if len(subConList) > 0:
                #from lxml import etree
                #bodyEtree = etree.Element('dummy')
                #bodyEtree.append(subConList[0])
                bodyEtree = subConList[0]
        try:
            f.saveFile(bodyEtree=bodyEtree)
        except CException as e:
            errorReport.append(e)
        except:
            errorReport.append(self.__class__, 130, fileName, name=self.objectPath())
        return errorReport

    def loadFromEtree(self, element):
        myException = CErrorReport()
        self.__dict__['_dataOrder'] = []
        name = element.get('id', None)
        if name is not None:
            self.setObjectName(name)
        for ele in element.iterchildren():
            className = ele.tag
            cls = DATAMANAGER().getClass(className)
            if cls is None:
                myException.append(self.__class__, 142, className, name=self.objectPath())
            else:
                name = ele.get('id', None)
                if name is None:
                    myException.append(self.__class__, 143, className, name=self.objectPath())
                elif self._dataOrder.count(name):
                    myException.append(self.__class__, 144, name, name=self.objectPath())
                else:
                    try:
                        obj = cls(parent=self, name=name)
                    except CException as e:
                        myException.extend(e)
                    except Exception as e:
                        myException.append(self.__class__, 104, name, name=self.objectPath())
                    if obj is not None:
                        try:
                            self.addObject(object=obj)
                        except CException as e:
                            myException.extend(e)
                        except Exception as e:
                            myException.append(self.__class__, 116, name, name=self.objectPath())
                        rv = self._value[name].setEtree(ele)
                        if len(rv) > 0:
                            myException.extend(rv)
        return myException

    def loadContentsFromEtree(self, element, readHeader=True, overwrite=False):
        myException = CErrorReport()
        myOverwrite = overwrite
        if not myOverwrite:
            self.__dict__['_dataOrder'] = []
        name = element.get('id', None)
        if name is not None and readHeader:
            # Beware reading the ccp4i2_body id from a 'base' def file
            # readHeader will be False in this case
            self.setObjectName(name)
        for ele in element.iterchildren():
            obj = None
            if ele.tag == 'content':
                #try:
                e,obj = self.makeDataObjectFromEtree(ele)
                #except:
                #  myException.append(self.__class__, 125, 'Content:' + str(ele.get('id')), name=self.objectPath())
                myException.extend(e)
            elif ele.tag == 'container':
                name = ele.get('id', None)
                if name is None:
                    CErrorReport.append(self.__class__, 126, name=self.objectPath())
                    name = 'noName'
                # Allow that container may already exist and just load to it
                if name in self.__dict__['_value']:
                    obj = self.__dict__['_value'][name]
                else:
                    try:
                        obj = CContainer(parent=self,name=name)
                    except:
                        myException.append(self.__class__,110,name=self.objectPath())
                if  obj is not None:
                    if DEVELOPER():
                        obj.loadContentsFromEtree(ele,overwrite=myOverwrite)
                    else:
                        try:
                            obj.loadContentsFromEtree(ele,overwrite=myOverwrite)
                        except  CException as e:
                            myException.extend(e)
                            obj = None
                        except:
                            myException.append(self.__class__, 108, 'Tag: ' + name, name=self.objectPath())
                            obj = None
            elif ele.tag == 'header':
                if readHeader:
                    try:
                        if not hasattr(self.__dict__,'header'):
                            from core import CCP4File
                            self.__dict__['header']=CCP4File.CI2XmlHeader(parent=self, name='header')
                        self.__dict__['header'].setEtree(ele)
                        #print 'CContainer.loadContentsFromEtree header',self.__dict__['header'],self.__dict__['header'].pluginName
                        if self.__dict__['header'].pluginName is not None:
                            self.setObjectName(str(self.__dict__['header'].pluginName))
                    except CException as e:
                        myException.extend(e)
                    except:
                        myException.append(self.__class__, 123, name=self.objectPath())
            elif ele.tag == 'qualifiers':
                self.setQualifiersEtree(ele)
            elif ele.tag == 'file':
                fileNameEle = ele.find('CI2XmlDataFile')
                if fileNameEle is None:
                    myException.append(self.__class__, 147, name=self.objectPath())
                else:
                    from core import CCP4File
                    fileObj = CCP4File.CI2XmlDataFile()
                    err = fileObj.setEtree(fileNameEle)
                    if len(err) > 0:
                        myException.append(self.__class__, 147, name=self.objectPath())
                    else:
                        if not fileObj.exists():
                            myException.append(self.__class__, 148,str(fileObj), name=self.objectPath())
                        else:
                            try:
                                fileContentEle = fileObj.getBodyEtree()
                            except:
                                myException.append(self.__class__, 149,str(fileObj), name=self.objectPath())
                            else:
                                err = self.loadContentsFromEtree(fileContentEle, readHeader=False)
                                if len(err) > 0:
                                    myException.extend(err)
                                # We have read in a 'base' definition from another file and may now modify it by overwriting
                                # with definitions from this file
                                myOverwrite = True
            if obj is not None:
                if obj.objectName() in self.__dict__['CONTENTS'] and myOverwrite:
                    self.replaceObject(object=obj, name=obj.objectName())
                else:
                    self.addObject(object=obj)
        # Not sure if we treat this as exception or return
        #if len(myException)>0: raise myException
        #print 'CContainer.loadContentsFromEtree myException',len(myException),myException.report()
        #print self.__dict__, myException
        return myException

    def saveContentsToEtree(self):
        if XMLPARSER() == 'lxml':
            from lxml import etree
        errors = CErrorReport()
        # Create element
        element = etree.Element('container')
        name = self.objectName()
        if name is not None and len(name) > 0:
            element.set('id', name)
        # Loop over contents
        for name in self._dataOrder:
            if name in self.CONTENTS:
                cls = self.CONTENTS[name].get('class', None)
                obj = self.__dict__['_value'].get(name, None)
                #print 'saveContentsToEtree',name,cls,obj
                if cls is None or obj is None:
                    errors.append(self.__class__, 129, name, name=self.objectPath())
                elif cls == CContainer:
                    #try:
                    ele,errs = self.__dict__['_value'][name].saveContentsToEtree()
                    errors.extend(errs)
                    # also keep qualifiers
                    qualiEle, errs = obj.qualifiersEtree(customOnly=True, tag='qualifiers', recurse=False)
                    #print 'saveContentsToEtree qualiEle'
                    #print etree.tostring(qualiEle, pretty_print=True)
                    errors.extend(errs)
                    ele.append(qualiEle)
                    element.append(ele)
                    #except:
                    #  errors.append(self.__class__, 128, name)
                else:
                    #try:
                        ele = etree.Element('content')
                        ele.set('id',name)
                        classEle = etree.Element('className')
                        cls =  self.CONTENTS[name].get('class', None)
                        if cls is not None:
                            classEle.text = cls.__name__
                        else:
                            classEle.text = 'UNKNOWN'
                        ele.append(classEle)
                        qualiEle,errs = obj.qualifiersEtree(customOnly=True, tag='qualifiers', recurse=False)
                        errors.extend(errs)
                        ele.append(qualiEle)
                        if CCP4Data.isCollectionClass(cls):
                            listEle,errs = obj.subItemEtree(customOnly=True)
                            ele.append(listEle)
                            errors.extend(errs)
                        element.append(ele)
                    #except:
                    #  errors.append(self.__class__, 130, name)
        return element, errors

    def makeDataObjectFromEtree(self, element):
        rv = CErrorReport()
        if element.get('id') is not None:
            name = str(element.get('id'))
        else:
            name = None
        className = None
        qEle = None
        cEleList = []
        for subEle in element.iterchildren():
            if subEle.tag == 'className':
                className = str(subEle.text)
            elif subEle.tag == 'qualifiers':
                # Don't try to parse the qualifiers just
                # hand them over to the data class object
                qEle = subEle
            elif ['content', 'subItem', 'rowHeader', 'columnHeader'].count(str(subEle.tag)):
                # There is definition of some item of content
                cEleList.append(subEle)
        #print 'makeDataObjectFromEtree', name
        if len(cEleList) > 0:
            e,subContents = self.subContentsFromEtree(cEleList)
            CErrorReport.extend(e)
        else:
            subContents = {}
        #if len(cEleList) > 0: print 'makeDataObjectFromEtree',name,subContents
        if className is None:
            rv.append(self.__class__, 102, 'No className for tag: ' + str(name), name=self.objectPath())
            return rv,None
        else:
            cls = DATAMANAGER().getClass(className)
        if cls == None:
            rv.append(self.__class__, 103, 'Class: ' + className + ' specified for: ' + str(name), name=self.objectPath())
            return rv, None
        else:
            buildInInit = True
            if DEVELOPER():
                if len(subContents) == 0:
                    obj = cls(parent=self, name=name)
                elif CCP4Data.isCollectionClass(cls):
                    buildInInit = False
                    obj = cls(parent=self, name=name, subItem=subContents.get('subItem', {}) ,build=False)
                else:
                    obj = cls(parent=self, name=name, contents=subContents)
            else:
                try:
                    if len(subContents) == 0:
                        obj = cls(parent=self, name=name, build=False)
                    elif CCP4Data.isCollectionClass(cls):
                        buildInInit = False
                        obj = cls(parent=self, name=name, subItem=subContents.get('subItem', {}), build=False)
                    else:
                        obj = cls(parent=self, name=name, contents=subContents)
                except CException as e:
                    rv.extend(e)
                except:
                    rv.append(self.__class__, 104, 'Class: ' + className + ' specified for: ' + name, name=self.objectPath())
            #print 'makeDataObjectFromEtree',name,type(obj),qEle
            if obj is not None:
                if qEle is not None:
                    rv.extend(obj.setQualifiersEtree(qEle))
                    if not buildInInit:
                        obj.build(qualifiers=obj.__dict__['_qualifiers'])
                    #print 'CContainer.makeDataObjectFromEtree',name,obj.__dict__['_value'],qEle.find('default')
                    if qEle.find('default') is not None:
                        obj.setDefault()
        return rv, obj

    def subContentsFromEtree(self, elementList=[]):
        rv = CErrorReport()
        tmpContainer = CContainer(name='temporary', parent=self)
        contents = {}
        for ele in elementList:
            e,subObj = tmpContainer.makeDataObjectFromEtree(ele)
            if subObj is not None:
                name = subObj.objectName()
                if name is None or len(name) == 0:
                    name = str(ele.tag)
                contents[name] = {'class' : subObj.__class__, 'qualifiers' : subObj.qualifiers(default=False)}
        tmpContainer.deleteLater()
        rv.extend(e)
        #print 'subContentsFromEtree', contents
        return rv,contents


    def addContent(self, name=None, cls=None, qualifiers={}, subItem={}, value=None):
        # Add an item of content that is defined in Python code
        if name is None:
            raise CException(self.__class__, 114, 'Content name: ' + name, name=self.objectPath())
        elif name in self.__dict__['_value']:
            raise CException(self.__class__, 115, 'Content name: ' + name, name=self.objectPath())
        if cls is None:
            raise CException(self.__class__, 120, name=self.objectPath())
        if isinstance(cls,str):
            cls = DATAMANAGER().getClass(cls)
        if cls is None or not issubclass(cls, CCP4Data.CData):
            raise CException(self.__class__, 120, 'Content name: ' + name, name=self.objectPath())
        try:
            self.__dict__['CONTENTS'][name] = {'class' : cls , 'qualifiers' : qualifiers}
        except:
            raise CException(self.__class__, 121, 'Content name: ' + name, name=self.objectPath())
        if cls == CContainer:
            self.__dict__['_value'][name] = CContainer(name=name)
        else:
            try:
                self.buildItem(key=name, cls=cls, qualifiers=qualifiers, subItem=subItem)
                self.__dict__['_dataOrder'].append(name)
            except CException as e:
                del self.__dict__['CONTENTS'][name]
                e.append(self.__class__, 122, 'Content name: ' + name, name=self.objectPath())
                raise e
            except:
                del self.__dict__['CONTENTS'][name]
                raise CException(self.__class__, 122, 'Content name: ' + name, name=self.objectPath())
        if value is not None:
            self.__dict__['_value'][name].set(value)
        return self.__dict__['_value'][name]

    def addObject(self, object=None, name=None, reparent=True, afterObject=None):
        #print 'CContainer.addObject',name, object.__class__,object.objectName()
        import inspect
        derivesFromCData = CCP4Data.CData.__name__ in [candidate.__name__ for candidate in inspect.getmro(type(object))[1:]]
        if object is None or not derivesFromCData:#isinstance(object, CCP4Data.CData):
            print("Object is",object, type(object))
            raise CException(self.__class__, 113, name=self.objectPath())
        if name is None:
            name = object.objectName()
            if len(name) <= 0:
                raise CException(self.__class__, 114, name=self.objectPath())
        if name in self.__dict__['_value']:
            raise CException(self.__class__, 115, 'Object name: ' + name, name=self.objectPath())
        if reparent:
            try:
                object.setParent(self)
            except:
                raise CException(self.__class__, 119, 'Object name: ' + str(name), name=self.objectPath())
        try:
            self.__dict__['CONTENTS'][name] = {'class' : object.__class__, 'qualifiers' : object.qualifiers(default=False)}
            self.__dict__['_value'][name] = object
        except Exception as e:
            raise CException(self.__class__, 116, 'Object name: ' + name + ' ' + str(e), name=self.objectPath(), exc_info=sys.exc_info())
        try:
            if afterObject is not None and self.__dict__['_dataOrder'].count(afterObject):
                ipos = self.__dict__['_dataOrder'].index(afterObject)
            else:
                ipos = -1
            if ipos <0:
                self.__dict__['_dataOrder'].append(name)
            else:
                self.__dict__['_dataOrder'].insert(ipos+1,name)
        except Exception as e:
            raise CException(self.__class__,116,'Object name: ' + name + ' ' + str(e), name=self.objectPath(), exc_info=sys.exc_info())

    def replaceObject(self, object=None, name=None, reparent=True):
        if name is None or name not in self.__dict__['CONTENTS']:
            raise CException(self.__class__, 117, 'Object name: ' + str(name), name=self.objectPath())
        import inspect
        derivesFromCData = CCP4Data.CData.__name__ in [candidate.__name__ for candidate in inspect.getmro(type(object))[1:]]
        if object is None or not derivesFromCData:#isinstance(object, CCP4Data.CData):
            print("Object is",object, type(object))
            raise CException(self.__class__, 113,name=self.objectPath())
        if reparent:
            try:
                object.setParent(self)
            except:
                raise CException(self.__class__, 119, 'Object name: ' + str(name), name=self.objectPath())
        try:
            self.__dict__['CONTENTS'][name] = {'class' : object.__class__, 'qualifiers' : object.qualifiers(default=False)}
            self.__dict__['_value'][name] = object
        except:
            raise CException(self.__class__, 116, 'Object name: ' + name, name=self.objectPath())

    def deleteObject(self, name=None):
        if name is None or name not in self.__dict__['CONTENTS']:
            raise CException(self.__class__, 117, 'Object name: ' + str(name), name=self.objectPath())
        try:
            del self.__dict__['CONTENTS'][name]
            del self.__dict__['_value'][name]
        except:
            raise CException(self.__class__, 118, 'Object name: ' + str(name), name=self.objectPath())
        if self.__dict__['_dataOrder'].count(name):
            self.__dict__['_dataOrder'].remove(name)

    def renameObject(self,oldName=None,newName=None):
        if oldName is None or newName is None or len(newName) == 0:
            raise CException(self.__class__, 132, 'Old name: ' + str(oldName) + ' new name:' + str(newName), name=self.objectPath())
        elif newName in self.__dict__['_value']:
            raise CException(self.__class__, 133, 'Old name: ' + str(oldName) + ' new name:' + str(newName), name=self.objectPath())
        if oldName not in self.__dict__['_value']:
            raise CException(self.__class__, 134, 'Old name: ' + str(oldName) + ' new name:' + str(newName), name=self.objectPath())
        try:
            self.__dict__['CONTENTS'][newName] = self.__dict__['CONTENTS'].pop(oldName)
            self.__dict__['_value'][newName] = self.__dict__['_value'].pop(oldName)
            ii = self.__dict__['_dataOrder'].index(oldName)
            if ii>=0: self.__dict__['_dataOrder'][ii] = newName
            self.__dict__['_value'][newName].setObjectName(newName)
        except:
            raise CException(self.__class__, 135, 'Old name: ' + str(oldName) + ' new name:' + str(newName), name=self.objectPath())

    def clear(self):
        self.__dict__['CONTENTS'] = {}
        self.__dict__['_value'] = {}
        self.__dict__['_dataOrder'] = []

    def addHeader(self):
        #print 'CContainer.addHeader',repr(self), self.__dict__.has_key('header'),
        if 'header' not in self.__dict__:
            from core import CCP4File
            self.__dict__['header'] = CCP4File.CI2XmlHeader(parent=self)
        #print repr(self.__dict__['header']),self.__dict__['header']
        return self.__dict__['header']

    def getHeader(self):
        if 'header' in self.__dict__:
            return self.__dict__['header']
        else:
            return None

    def getHeaderEtree(self,jobId=None, project=None, function='PARAMS'):
        # Return an eTree for a data file header
        # This method could be member of CTaskViewer
        from core import CCP4File
        header = CCP4File.CI2XmlHeader(parent=self, name='header')
        if 'header' in self.__dict__:
            header.set(self.__dict__['header'])
        else:
            header.pluginName = self.objectName()
        header.function.set(function)
        header.userId.setCurrentUser()
        header.creationTime.setCurrentTime()
        header.ccp4iVersion.set(CCP4Config.VERSION().ccp4iVersion)
        if jobId is not None:
            header.jobId = jobId
        if project is not None:
            header.project = project
        return header.getEtree()

    def nonexistantFiles(self):
        #Return a list of files that mustExist but do not!    
        from core import CCP4File
        import inspect
        nonexistantList = []
        for key,obj0 in list(self.__dict__['_value'].items()):
            objList, xmlText, keyValues = obj0.saveToDb()
            for obj in objList:
                if isinstance(obj,CCP4File.CDataFile):
                    #print 'CContainer.nonexistantFiles', key, obj.qualifiers('mustExist'), obj.exists(), obj.isSet(), obj.qualifiers('allowUndefined'), obj.fullPath
                    if not obj.isSet():
                        if not obj.qualifiers('allowUndefined'):
                            nonexistantList.append(key)
                    else:
                        if obj.qualifiers('mustExist') and not obj.exists():
                            nonexistantList.append(key)
        return nonexistantList

    def inputFilesFileIds(self):
        from core import CCP4File
        ret= {}
        for key,obj0 in list(self.__dict__['_value'].items()):
            objList, xmlText, keyValues = obj0.saveToDb()
            for obj in objList:
                if isinstance(obj, CCP4File.CDataFile) and obj.dbFileId.isSet():
                    ret[str(obj.dbFileId)] = key
        return ret

    def dataOrder(self):
        return self.__dict__['_dataOrder']

    def validity(self):
        err = CErrorReport()
        for key in self.dataOrder():
            instance = getattr(self, key)
            if isinstance(instance, CContainer):
                err.extend(instance.validity())
            else:            
                err.extend(instance.validity(instance.get()))
        return err  

    def firstDataObject(self):
        #Return the first data (not container) object -
        # useful for initialising guis
        ii = -1
        while ii+1 < len(self.__dict__['_dataOrder']):
            ii = ii + 1
            name = self.__dict__['_dataOrder'][ii]
            #print 'firstDataItem', repr(self), ii, name
            if self.__dict__['CONTENTS'][name]['class'] == CContainer:
                firstItem = self.__dict__['_value'][name].firstDataObject()
                if firstItem is not None:
                    return firstItem
            else:
                return self.__dict__['_value'][name]
        return None

    def pluginName(self):
        if 'header' in self.__dict__ and self.__dict__['header'].pluginName is not None:
            return self.__dict__['header'].pluginName
        else:
            return self.objectName()

    def find(self, name=None):
        # Find an item in this or any child container
        if name in self.__dict__['_value']:
            return self.__dict__['_value'].get(name)
        else:
            splitName = name.split('.')
            if len(splitName)>1:
                if splitName[0] in self.__dict__['_value'] and isinstance(self.__dict__['_value'][splitName[0]], CContainer):
                    newName = splitName[1]
                    for item in splitName[2:]:
                        newName = newName + '.' + item
                    return self.__dict__['_value'][splitName[0]].find(newName)
                return None
            else:
                for key, obj in list(self.__dict__['_value'].items()):
                    if isinstance(obj,CContainer):
                        rv = obj.find(name)
                        if rv is not None:
                            return rv
                return None

    def copyData(self, otherContainer=None, dataList=None):
        report = CErrorReport()
        if otherContainer is None or not isinstance(otherContainer, CContainer):
            report.append(self.__class__, 138, name=self.objectPath())
            return report
        if dataList is None:
            dataList = self.dataOrder()
        for name in dataList:
            if isinstance(name,(list, tuple)):
                if len(name) >= 2:
                    otherName = name[0]
                    name = name[1]
                else:
                    name = name[0]
                    otherName = name[0]
            else:
                otherName = name
            myData = self.find(name)
            otherData = otherContainer.find(otherName)
            #print 'CContainer.copyData', name, myData,myData
            if myData is None:
                report.append(self.__class__, 139, details=name, name=self.objectPath())
            elif otherData is None:
                report.append(self.__class__, 140, details=otherName, name=self.objectPath())
            else:
                try:
                    myData.set(otherData)
                except CException as e:
                    report.extend(e)
                except Exception as e:
                    report.append(self.__class__, 141, details=name)
        return report

    def parseCommandLine(self, commandLine=[], template=[]):
        errorReport = CErrorReport()
        indx = -1
        while indx+1 < len(commandLine):
            indx = indx + 1
            # Try if command line keyword matches name of data object
            # If so then assume it is followed by one data value
            # This should work for the usual CCP4 com line
            dataObj = self.find(commandLine[indx])
            if dataObj is not None:
                if isinstance(dataObj,CCP4Data.CBoolean):
                    dataObj.set(True)
                else:
                    if indx+1 < len(commandLine):
                        indx = indx + 1
                        try:
                            dataObj.set(commandLine[indx])
                        except CException as e:
                            errorReport.extend(e)
                    else:
                        errorReport.append(self.__class__, 146, commandLine[indx], name=self.objectPath())
            # Alternatively try looking for keyword in the commandLineLookup
            else:
                for defn in template:
                    found = 0
                    if isinstance(defn[0], list):
                        found = defn[0].count(commandLine[indx])
                    elif commandLine[indx] == defn[0]:
                        found = 1
                    if found:
                        dataObj = self.find(defn[1])
                        if dataObj is not None:
                            if len(defn) == 2:
                                if indx + 1 < len(commandLine):
                                    indx = indx + 1
                                    #print 'parseCommandLine found dataObj', dataObj.objectName(), commandLine[indx]
                                    try:
                                        dataObj.set(commandLine[indx])
                                    except CException as e:
                                        errorReport.extend(e)
                                    except Exception as e:
                                        errorReport.append(self.__class__, 145, dataObj.objectName(), name=self.objectPath())
                                else:
                                    errorReport.append(self.__class__, 146, commandLine[indx], name=self.objectPath())
                            else:
                                dataDict = {}
                                if indx + len(defn[2]) < len(commandLine):
                                    for item in defn[2]:
                                        indx = indx + 1
                                        # There should be no neeed to convert type?
                                        dataDict[item] = commandLine[indx]
                                    #print 'CContainer.parseCommandLine dataDict',dataDict
                                    try:
                                        dataObj.set(dataDict)
                                    except CException as e:
                                        errorReport.extend(e)
                                    except Exception as e:
                                        errorReport.append(self.__class__, 145, dataObj.objectName(), name=self.objectPath())
                                else:
                                    errorReport.append(self.__class__, 146, commandLine[indx], name=self.objectPath())
        return errorReport


#===============================================================================================
import unittest

class testCContainer(unittest.TestCase):

    #Removing broken tests
    def broken_test1(self):
        from core import CCP4Utils
        self.summat = CContainer(definitionFile=os.path.join(CCP4Utils.getCCP4I2Dir(),'test/data/summat.contents.xml'))
        #print summat.name,summat.CONTENTS
        self.assertEqual(self.summat.header.pluginName,'Summat','Wrong header.pluginName')
        #print 'summat.CONTENTS',self.summat.CONTENTS.keys()
        self.assertEqual(len(list(self.summat.CONTENTS.keys())),4,'Wrong content length')
        self.assertEqual(self.summat.CONTENTS['range']['class'],CCP4Data.CIntRange,'range parameter has wrong class')
        #print "summat.CONTENTS['range']['qualifiers']",self.summat.CONTENTS['range']['qualifiers']
        self.assertEqual(self.summat.CONTENTS['range']['qualifiers']['compare'],1,'range has wrong qualifiers')
        #print 'summat.gubbins.startValue',self.summat.gubbins.startValue
        self.assertEqual(self.summat.gubbins.startValue,12.0,'sub-container CFloat wrong initial value')

    def broken_test2(self):
        contents = {'nCycles': {'class': CCP4Data.CInt, 'qualifiers': { 'default' : None, 'max': 20, 'allowUndefined': True, 'min': 1}}, 'cutoff': {'class': CCP4Data.CInt, 'qualifiers': { 'default' :3,'max': 100, 'min': 1}}, 'range': {'class': CCP4Data.CIntRange, 'qualifiers': { 'start' : { 'default' : 0 } ,  'end' : { 'default' : 10 } , 'compare': 1 }}}
        self.summat = CContainer(name='Summat',contents=contents)
        #self.summat.range.set(start=10,end=20)
        self.summat.range.start.set(10)
        self.summat.range.end.set(20)
        self.assertEqual(self.summat.get('range'), {'start':10,'end':20},'Failed setting range parameter')
    
    def broken_test3(self):
        from core import CCP4Utils
        defFile = os.path.join(CCP4Utils.getCCP4I2Dir(),'tasks','demo','demo.def.xml')
        demo =  CContainer(definitionFile=defFile)
        self.assertEqual(demo.controlParameters.QUICKNDIRTY,'quick','Loading demo.def.xml does not give correct value for controlParameters.QUICKNDIRTY')
      
  
    def broken_test4(self):
        from core import CCP4XtalData
        c = CContainer()
        c.addObject(CCP4XtalData.CMtzDataFile(name='MTZIN'))
        c.addObject(CCP4Data.CFloatRange(name='RESOLUTION_RANGE',qualifiers = {'compare':-1,'start' :{'min':0.0,'default':999.9},'end':{'min':0.0,'default':0.0}}))
        c.MTZIN.fullPath = '/foo/bar/wotsit.tmp'
        c.RESOLUTION_RANGE.start=67.9
        c.RESOLUTION_RANGE.end = 2.1
        self.assertEqual(c.RESOLUTION_RANGE.start.qualifiers('min'),0.0,'Error adding CFloatRange to CContainer - wrong qualifier')
        self.assertEqual(c.RESOLUTION_RANGE.end,2.1,'Error adding  CFloatRange to CContainer - wrong value')
        self.assertEqual(c.MTZIN.baseName,'wotsit.tmp','Error adding CMtzDataFile to CContainer')
        c.deleteObject('MTZIN')
        self.assertEqual(len(c),1,'Error deleteing object from CContainer')

    def test5(self):
        from core import CCP4XtalData
        c = CContainer()
        c.addContent(name='MTZIN',cls=CCP4XtalData.CMtzDataFile)
        c.addContent(name='RESOLUTION_RANGE',cls='CFloatRange',qualifiers={'compare':-1,'start' : { 'min':0.0,'default': 1.0}, 'end': {'min':0.0,'default':0.0}})
        self.assertEqual(c.RESOLUTION_RANGE.start.qualifiers('min'),0.0,'Error adding CFloatRange to CContainer - wrong qualifier')

    '''
    def test6(self):
        from core import CCP4Utils
        defFile = os.path.join(CCP4Utils.getCCP4I2Dir(),'tasks','demo','demo.def.xml')
        demo =  CContainer(definitionFile=defFile)
        p = demo.getPersistent()
        other = CContainer()
        e = other.setPersistent(p)
        self.assertEqual(other.controlParameters.QUICKNDIRTY,'quick','Loading demo.def.xml does not give correct value for controlParameters.QUICKNDIRTY')
    '''

    def broken_test7(self):
        from core import CCP4Utils
        defFile = os.path.join(CCP4Utils.getCCP4I2Dir(),'tasks','demo','demo.def.xml')
        demo =  CContainer(definitionFile=defFile)
        #print 'CContainer.test7',demo.inputData.dataOrder()
        demo.inputData.PDBIN = os.path.join(CCP4Utils.getCCP4I2Dir(),'test','data','1df7.pdb')
        #print 'CContainer.test7', demo.inputData.PDBIN
        demo.controlParameters.CHOOSEONE = 5
        demo.controlParameters.CHOOSEOTHER = 2
        other = CContainer(definitionFile=defFile)
        rv = other.copyData(demo,['PDBIN','CHOOSEONE','CHOOSEOTHER','FOO'])
        #print 'CContainer.test7',rv,other.inputData.PDBIN
        self.assertEqual(other.inputData.PDBIN.baseName,'1df7.pdb','Copying between containers failed to give correct filename')
        self.assertEqual(other.controlParameters.CHOOSEONE,5,'Copying between containers failed to give correct value for CHOOSEONE')
        self.assertEqual(len(rv),1,'Copying between containers failed to give correct length error report')
        self.assertEqual(rv[0]['code'],139,'Copying between containers failed to give correct error report code')

    def test8(self):
        from core import CCP4XtalData
        from core import CCP4ModelData 
        c = CContainer()
        i = CContainer(name='inputData')
        i.addContent(name='HKLIN',cls=CCP4XtalData.CMtzDataFile)
        i.addContent(name='XYZIN',cls=CCP4ModelData.CPdbDataFile)
        c.addObject(i)
        err = c.parseCommandLine(['HKLIN','/mydir/myfile.mtz','XYZIN','/mydir/myotherfile.pdb'])
        #print 'parseCommandLine err',err.report()
        self.assertEqual(len(err),0,'Parsing command line returns an error')
        self.assertEqual(str(c.inputData.XYZIN.fullPath),'/mydir/myotherfile.pdb','Parsing command line XYZIN misset')

    def test9(self):
        from core import CCP4MathsData
        from core import CCP4ModelData
        c = CContainer()
        c.addContent(name='XYZIN',cls=CCP4ModelData.CPdbDataFile)
        c.addContent(name='COORD',cls = CCP4MathsData.CXyz)
        err = c.parseCommandLine(['-p', '/mydir/myotherfile.pdb', '-x', '1.0', '2.0', '3.0'],
                                 [['-p' , 'XYZIN'], ['-x' , 'COORD', ['x', 'y', 'z']]])
        #print 'test9 err',err.report()
        self.assertEqual(len(err), 0, 'Parsing command line returns an error')
        self.assertEqual(str(c.XYZIN.fullPath), '/mydir/myotherfile.pdb', 'Parsing command line XYZIN misset')
        self.assertEqual(c.COORD.z, 3.0, 'Parsing command line COORDS.z missed')

    def broken_test10(self):
        from core import CCP4Utils
        defFile = os.path.join(CCP4Utils.getCCP4I2Dir(), 'test', 'data', 'summat_1.def.xml')
        c = CContainer(definitionFile=defFile)
        self.assertEqual(c.PDBIN.baseName.isSet(), False, 'Weirdness in loading from second def file')

    def broken_test11(self):
        from core import CCP4Utils
        c = CContainer()
        c.loadDataFromXml(os.path.join(CCP4Utils.getCCP4I2Dir(), 'pipelines', 'demo_copycell', 'test_data', 'test_1.params.xml'))
        self.assertTrue(c.inputData.HKLIN.isSet(), 'Error using loadDataFromXml with auto load of def file')
    
def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testCContainer)
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)
