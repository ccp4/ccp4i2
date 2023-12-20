from __future__ import print_function

"""
     CCP4ErrorHandling.py: CCP4 GUI Project
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
   Liz Potterton Aug 2010 - Exception class and definitions of severity
"""
SEVERITY_OK = 0
SEVERITY_UNDEFINED = 1
SEVERITY_WARNING = 2
SEVERITY_UNDEFINED_ERROR = 3
SEVERITY_ERROR = 4
SEVERITY_CRITICAL = 5
SEVERITY_TEXT = ['OK', 'WARNING DATA UNDEFINED', 'WARNING', 'ERROR DATA UNDEFINED', 'ERROR', 'CRITICAL']

STACK_LIMIT = 5

import sys
import time
import traceback
from core import CCP4Config
# from core import CCP4Data  # - KJS : Causes issues. Circular dependency with globals.
from xml.etree import ElementTree as ET

def formatExceptionInfo(maxTBlevel=5):
    cla, exc, trbk = sys.exc_info()
    excName = cla.__name__
    try:
        excArgs = exc.__dict__["args"]
    except KeyError:
        excArgs = "<no args>"
    excTb = traceback.format_tb(trbk, maxTBlevel)
    return (excName, excArgs, excTb)


class CErrorReport():
    '''Holds list of errors and warnings'''
    def __init__(self, cls=None, code=0, details=None, name=None, label=None, recordTime=False, stack=True, exc_info=None):
        self._reports = []
        if cls is not None:
            self.append(cls=cls, code=code, details=details, name=name, label=label, recordTime=recordTime, stack=stack, exc_info=exc_info)

    def append(self, cls=None, code=0, details=None, name=None, label=None, recordTime=False, stack=True, exc_info=None):
        report = {'class' : cls, 'name' : name, 'label' : label, 'code' : code, 'details' : details}
        if recordTime:
            report['time'] = time.mktime(time.localtime())
        if stack:
            report['stack'] = self.getStack(exc_info=exc_info)
        self._reports.append(report)
        try:
            from core import CCP4Data         # KJS : Did this work before without this ?
            severity = CCP4Data.errorCodeSeverity(report['class'], report['code'])
            if severity == SEVERITY_CRITICAL:
                print(self.report())
        except:
            pass

    def extend(self, other=None, label=None, recordTime=False, stack=True):
        if other is None or len(other) == 0 or not hasattr(other,"_reports"):
            return
        for item in other._reports:
            #if label is not None: item['name'] = label + '_' + item['name']
            if recordTime and item.get('time', None) is None:
                item['time'] = time.mktime(time.localtime())
            if stack and item.get('stack', None) is None :
                item['stack'] = self.getStack()
        self._reports.extend(other._reports)

    def appendPythonException(self, cls=None, exception=None):
        details = str(exception)
        self._reports.append({'class' : cls, 'name' : 'Python error', 'code' : -2, 'details' : details})

    def count(self, cls=None, code=None):
        n = 0
        for report in self._reports:
            #print 'CErrorReport.count',cls,report['class']
            if (cls is None or report['class'] == cls) and (code is None or report['code'] == code):
                n = n + 1
        return n

    def appendDetails(self, details=''):
        for report in self._reports:
            report['details'] = report['details'] + details

    def setName(self, name=''):
        for indx in range(len(self._reports)):
            self._reports[indx]['name'] = name

    def prependName(self, name=''):
        #print 'prependName',name,self._reports
        for item in self._reports:
            if item['name'] is None:
                item['name'] = name
            if len(item['name']) > 0 and not name == item['name']:
                item['name'] = name + '_' + item['name']
            else:
                item['name'] = name

    def maxSeverity(self):
        from core import CCP4Data
        maxSev = 0
        for report in self._reports:
            severity = -1
            try:
                code = report['code']
                if issubclass(report['class'], CCP4Data.CData):
                    severity = CCP4Data.errorCodeSeverity(report['class'], report['code'])
                elif  report['code'] in report['class'].ERROR_CODES:
                    severity = report['class'].ERROR_CODES[report['code']].get('severity', SEVERITY_ERROR)
                else:
                    for baseClass in report['class'].__bases__:
                        if 'ERROR_CODES' in baseClass.__dict__ and report['code'] in baseClass.ERROR_CODES:
                            severity = baseClass.ERROR_CODES[report['code']].get('severity', SEVERITY_ERROR)
                if severity < 0:
                    print('ERROR in CErrorHandling - error code not found', report['code'], report['class'])
                maxSev = max(maxSev, severity)
            except:
                print('maxSeverity error', report)
        return maxSev

    def __len__(self):
        return len(self._reports)

    def __getitem__(self, arg):
        return self._reports[arg]

    def __str__(self):
        return self.report()

    def description(self, report=None, inclName=True, user=False):
        from core import CCP4Data
        if report is None:
            report = self._reports[0]
        if 'description' in report and report['description'] is not None:
            desc = report['description']
            severity = report.get('severity', SEVERITY_ERROR)
        elif issubclass(report['class'], CCP4Data.CData):
            severity = CCP4Data.errorCodeSeverity(report['class'], report['code'])
            desc = CCP4Data.errorCodeDescription(report['class'], report['code'])
        else:
            if report['code'] == -2:
                severity = 4
                desc = 'Python error'
            elif report['code'] in report['class'].ERROR_CODES:
                severity = report['class'].ERROR_CODES[report['code']].get('severity',SEVERITY_ERROR)
                desc = report['class'].ERROR_CODES[report['code']].get('description','NO DESCRIPTION')
            else:
                severity = 4
                desc = 'No description available'
                for baseClass in report['class'].__bases__:
                    if 'ERROR_CODES' in baseClass.__dict__ and report['code'] in baseClass.ERROR_CODES:
                        severity = baseClass.ERROR_CODES[report['code']].get('severity',SEVERITY_ERROR)
                        desc = baseClass.ERROR_CODES[report['code']].get('description','NO DESCRIPTION')
                        break
        if inclName:
            if user and report.get('label', None) is not None and report['label'] is not NotImplemented:
                #print 'CErrorReport description', report['label'], '*', desc
                desc = str(report['label']) + ': ' +  desc 
                #desc = report['class'].QUALIFIERS['guiLabel'] + ': ' + desc 
            elif 'name' in report and report['name'] is not None:
                desc = str(report['name']) + ': ' +  str(desc)
        return desc, severity

    def report(self, user=False, ifStack=True, mode=0, minSeverity=SEVERITY_UNDEFINED):
        text = ''
        #print 'CErrorReport.report len',len(self._reports)
        if len(self._reports) > 0:
            for report in self._reports:
                #print 'CErrorReport.report', report['class'],report['code'], user
                desc, severity = self.description(report, inclName=True, user=user)
                if severity == SEVERITY_CRITICAL:
                    text = text + "\nCRITICAL ERROR PLEASE REPORT TO CCP4:"
                if severity >= minSeverity:
                    try:
                        className = report['class'].__name__
                    except:
                        className = str(report['class'])
                    name = str(report.get('name', ''))
                    #if user:
                    #  text = text + "\n{3}\n{0} : {2}".format(className, SEVERITY_TEXT[severity], report['code'], desc)
                    #else:
                    if mode == 0:
                        text = text + "\n{0:20} -{1}- {2}:{3} {4}".format(name, SEVERITY_TEXT[severity], className, report['code'], desc)
                    elif mode == 1:
                        text = text + "\n {0:20} -{1}- \n{2}:{3} {4}".format(name, SEVERITY_TEXT[severity], className, report['code'], desc)
                    else:
                        if user:
                            if severity == SEVERITY_WARNING:
                                text = text + "\nWarning: " + desc
                            else:
                                if (len(text)+len(desc)) < 60:
                                    text = text + " " + desc
                                else:
                                    text = text + "\n" + desc
                        else:
                            if (len(text)+len(desc)) < 60:
                                text = text + " " + desc
                            else:
                                text = text + "\n" + desc
                        #print 'report', className, name, '*', text, '*', desc, mode, user
                    if report['details'] is not None and len(str(report['details'])) > 0 and report['details'] != 'None':
                        if mode == 1:
                            text = text + ' ' + str(report['details'])
                        else:
                            if (len(text) + len(str(report['details']))) < 60:
                                text = text + ' ' + str(report['details']) + '\n'
                            else:
                                text = text + '\n' + str(report['details']) + '\n'
                    if ifStack and report.get('stack', None) is not None:
                        text = text + '\n'
                        for line in report['stack']:
                            text = text + line
            return text[1:]
        else:
            return ''

    def getStack(self, exc_info=None):
        formatted_stack = []
        if exc_info is not None:
            try:
                formatted_stack = traceback.format_exception(exc_info[0], exc_info[1], exc_info[2])
            except:
                pass
        else:
            try:
                # Try getting a traceback from sys.exc_info()
                stack = traceback.format_stack(sys.exc_info()[2])
            except:
                pass
                 
            # Just try getting a traceback from here
            if len(formatted_stack) == 0:
                try:
                    formatted_stack = traceback.format_stack()[0:-2]
                except:
                    pass
        #print 'getStack',formatted_stack,len(formatted_stack)
        if len(formatted_stack) > 0:
            return formatted_stack
        else:
            return None

    def warningMessage(self, windowTitle='', message='', jobId=None, parent=None, ifStack=True, minSeverity=SEVERITY_UNDEFINED):
        if len(message) > 0 and message[-1] !='\n':
            message = message + '\n'
        if CCP4Config.GRAPHICAL() and parent is not None:
            print('CException.warningMessage GRAPHICAL', CCP4Config.GRAPHICAL())
            #from PySide2 import QtGui, QtWidgets
            #QtWidgets.QMessageBox.warning(None, windowTitle, message + self.report())
            from qtgui import CCP4MessageBox
            m = CCP4MessageBox.CMessageBox(parent, title=windowTitle, message=message,
                                           details=self.report(ifStack=ifStack, minSeverity=minSeverity), jobId=jobId)
            m.show()
            #print 'CErrorReport.warningMessage', m
        else:
            print(self.report(ifStack=ifStack, minSeverity=minSeverity))

    '''
    def getEtree(self):
        from core import CCP4Config
        if CCP4Config.XMLPARSER() == 'lxml':
           from lxml import etree
        element = ET.Element('report')
        element.text = self.report()
        return element
    '''

    def getEtree(self):
        if CCP4Config.XMLPARSER() == 'lxml':
            from lxml import etree
        element = ET.Element('errorReportList')
        for item in self._reports:
            try:
                ele = ET.Element('errorReport')
                e = ET.Element('className')
                e.text = item['class'].__name__
                ele.append(e)
                e = ET.Element('code')
                e.text = str(item['code'])
                ele.append(e)
                e = ET.Element('description')
                desc,severity = self.description(item)
                e.text = desc
                ele.append(e)
                e = ET.Element('severity')
                e.text = SEVERITY_TEXT[severity]
                ele.append(e)
                if item['details'] is not None:
                    e = ET.Element('details')
                    e.text = str(item['details'])
                    ele.append(e)
                if item.get('time',None) is not None:
                    e = ET.Element('time')
                    e.text = str(item['time'])
                    ele.append(e)
                if item.get('stack',None) is not None:
                    e = ET.Element('stack')
                    #print 'CErrorReport.getEtree stack',item['stack'],type(item['stack']),type(item['stack'][0])
                    text = ''
                    for line in item['stack']:
                        text = text + line
                    e.text = text
                    ele.append(e)
                element.append(ele)
            except Exception as e:
                print('CErrorReport.getEtree error', e)
                traceback.print_exc()
        return element

    def setEtree(self, element=None, checkValidity=True):
        from core import CCP4DataManager
        from core import CCP4TaskManager
        DM = CCP4DataManager.DATAMANAGER()
        TM = CCP4TaskManager.TASKMANAGER()
        if CCP4Config.XMLPARSER() == 'lxml':
            from lxml import etree
        body = element.find('ccp4i2_body')
        if body is None:
            body = element
        for ele in body:
            if str(ele.tag) == 'errorReport':
                report = {}
                for e in ele:
                    name = str(e.tag)
                    if name == 'className':
                        clsName = str(e.text)
                        if clsName == 'CPluginScript':
                            from core import CCP4PluginScript
                            report['class'] = CCP4PluginScript.CPluginScript
                        else:
                            report['class'] = DM.getClass(clsName)
                            if report['class'] is None:
                                report['class'] = TM.getClass(clsName)
                    elif name == 'code':
                        try:
                            report['code'] = int(e.text)
                        except:
                            report['code'] = 0
                    elif name == 'stack':
                        try:
                            report['stack'] = str(e.text)
                        except:
                            report['stack'] = ''
                    elif name == 'details':
                        try:
                            report['details'] = str(e.text)
                        except:
                            report['details'] = None
                    elif name == 'description':
                        try:
                            report['description'] = str(e.text)
                        except:
                            pass
                    elif name == 'severity':
                        if SEVERITY_TEXT.count(str(e.text)):
                            report['severity'] = SEVERITY_TEXT.index(str(e.text))
                self._reports.append(report)

    def classesInReport(self):
        classList = []
        for report in self._reports:
            if report['class'] is not None and classList.count(report['class']) == 0:
                classList.append(report['class'])
        classNameList = []
        for item in classList:
            classNameList.append(item.__name__)
        classNameList.sort()
        return classNameList


class CException(CErrorReport, Exception):
    def __init__(self, cls=None, code=0, details='', name=None, recordTime=False, stack=True, report=None, exc_info=None):
        CErrorReport.__init__(self)
        Exception.__init__(self)
        if cls is not None:
            self.append(cls, code, details, name, recordTime, stack, exc_info)
        if report is not None:
            self.extend(report)



#===========================================================================================
import unittest
def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testError)
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)

class testError(unittest.TestCase):

    def test1(self):
        if CCP4Config.XMLPARSER() == 'lxml':
            from lxml import etree
        from core import CCP4Data
        e = CException(CCP4Data.CData, 1, 'foo')
        tree = e.getEtree()
        #text = ET.tostring(tree, xml_declaration=True)
        #print text
        f = CException()
        f.setEtree(tree)
        self.assertEqual(f[0]['code'],1,'Error save/restore CException to etree = wrong code')

