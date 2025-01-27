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

   Liz Potterton Aug 2010 - Exception class and definitions of severity
"""

from enum import Enum
from inspect import isclass
import time
import traceback
import unittest

from lxml import etree

from . import CCP4Config
from . import CCP4Data
from . import CCP4DataManager
from . import CCP4PluginScript
from . import CCP4TaskManager
from ..qtgui import CCP4MessageBox


class Severity(Enum):
    OK = 0, 'OK'
    UNDEFINED = 1, 'WARNING DATA UNDEFINED'
    WARNING = 2, 'WARNING'
    UNDEFINED_ERROR = 3, 'ERROR DATA UNDEFINED'
    ERROR = 4, 'ERROR'
    CRITICAL = 5, 'CRITICAL'

    def __lt__(self, other):
        if hasattr(other, "value"):
            return self.value < other.value
        return self.value < int(other)

    def __str__(self):
        return self.value[1]

    @classmethod
    def __getitem__(cls, text: str):
        return {severity.value[1]: severity for severity in cls}[text]


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
            report['stack'] = getStack(exc_info)
        self._reports.append(report)
        severity = errorCodeSeverity(cls, code)
        if severity == Severity.CRITICAL:
            print(self.report())

    def extend(self, other=None, recordTime=False, stack=True):
        if other is None or len(other) == 0 or not hasattr(other,"_reports"):
            return
        for item in other._reports:
            if recordTime and item.get('time', None) is None:
                item['time'] = time.mktime(time.localtime())
            if stack and item.get('stack', None) is None :
                item['stack'] = getStack()
        self._reports.extend(other._reports)

    def appendPythonException(self, cls=None, exception=None):
        details = str(exception)
        self._reports.append({'class' : cls, 'name' : 'Python error', 'code' : -2, 'details' : details})

    def count(self, cls=None, code=None):
        n = 0
        for report in self._reports:
            if (cls is None or report['class'] == cls) and (code is None or report['code'] == code):
                n = n + 1
        return n

    def appendDetails(self, details=''):
        for report in self._reports:
            report['details'] += details

    def setName(self, name=''):
        for report in self._reports:
            report['name'] = name

    def maxSeverity(self):
        maxSev = 0
        for report in self._reports:
            severity = errorCodeSeverity(report['class'], report['code'])
            if severity < 0:
                print('ERROR in CErrorHandling - error code not found', report['code'], report['class'])
            maxSev = max(maxSev, severity)
        return maxSev

    def __len__(self):
        return len(self._reports)

    def __getitem__(self, arg):
        return self._reports[arg]

    def __str__(self):
        return self.report()

    def description(self, report=None, user=False):
        if report is None:
            report = self._reports[0]
        if 'description' in report and report['description'] is not None:
            desc = report['description']
            severity = report.get('severity', Severity.ERROR)
        elif issubclass(report['class'], CCP4Data.CData):
            severity = errorCodeSeverity(report['class'], report['code'])
            desc = errorCodeDescription(report['class'], report['code'])
        elif report['code'] == -2:
            severity = Severity.ERROR
            desc = 'Python error'
        else:
            severity = errorCodeSeverity(report['class'], report['code'], Severity.ERROR)
            desc = errorCodeDescription(report['class'], report['code'], 'No description available')
        if user and report.get('label', None) is not None and report['label'] is not NotImplemented:
            desc = str(report['label']) + ': ' + desc
        elif 'name' in report and report['name'] is not None:
            desc = str(report['name']) + ': ' + str(desc)
        return desc, severity

    def report(self, user=False, ifStack=True, mode=0, minSeverity=Severity.UNDEFINED):
        text = ''
        if len(self._reports) > 0:
            for report in self._reports:
                desc, severity = self.description(report, user=user)
                if severity == Severity.CRITICAL:
                    text = text + "\nCRITICAL ERROR PLEASE REPORT TO CCP4:"
                if severity >= minSeverity:
                    try:
                        className = report['class'].__name__
                    except:
                        className = str(report['class'])
                    name = str(report.get('name', ''))
                    if mode == 0:
                        text = text + "\n{0:20} -{1}- {2}:{3} {4}".format(name, severity, className, report['code'], desc)
                    elif mode == 1:
                        text = text + "\n {0:20} -{1}- \n{2}:{3} {4}".format(name, severity, className, report['code'], desc)
                    else:
                        if user:
                            if severity == Severity.WARNING:
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
                    if 'details' in report and report['details'] is not None and len(str(report['details'])) > 0 and report['details'] != 'None':
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
        return ''

    def warningMessage(self, windowTitle='', message='', jobId=None, parent=None, ifStack=True, minSeverity=Severity.UNDEFINED):
        if len(message) > 0 and message[-1] !='\n':
            message = message + '\n'
        if CCP4Config.GRAPHICAL() and parent is not None:
            print('CException.warningMessage GRAPHICAL', CCP4Config.GRAPHICAL())
            m = CCP4MessageBox.CMessageBox(parent, title=windowTitle, message=message,
                                           details=self.report(ifStack=ifStack, minSeverity=minSeverity), jobId=jobId)
            m.show()
        else:
            print(self.report(ifStack=ifStack, minSeverity=minSeverity))

    def getEtree(self):
        element = etree.Element('errorReportList')
        for item in self._reports:
            try:
                ele = etree.Element('errorReport')
                e = etree.Element('className')
                e.text = item['class'].__name__
                ele.append(e)
                e = etree.Element('code')
                e.text = str(item['code'])
                ele.append(e)
                e = etree.Element('description')
                desc,severity = self.description(item)
                e.text = desc
                ele.append(e)
                e = etree.Element('severity')
                e.text = str(severity)
                ele.append(e)
                if item['details'] is not None:
                    e = etree.Element('details')
                    e.text = str(item['details'])
                    ele.append(e)
                if item.get('time',None) is not None:
                    e = etree.Element('time')
                    e.text = str(item['time'])
                    ele.append(e)
                if item.get('stack',None) is not None:
                    e = etree.Element('stack')
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

    def setEtree(self, element=None):
        DM = CCP4DataManager.DATAMANAGER()
        TM = CCP4TaskManager.TASKMANAGER()
        body = element.find('ccp4i2_body') or element
        for ele in body:
            if str(ele.tag) == 'errorReport':
                report = {}
                for e in ele.iterchildren():
                    name = str(e.tag)
                    if name == 'className':
                        clsName = str(e.text)
                        if clsName == 'CPluginScript':
                            report['class'] = CCP4PluginScript.CPluginScript
                        else:
                            report['class'] = DM.getClass(clsName)
                            if report['class'] is None:
                                report['class'] = TM.getClass(clsName)
                    elif name == 'code':
                        report['code'] = int(e.text) if e.text and e.text.isdecimal() else 0
                    elif name in {'stack', 'details', 'description'}:
                        report[name] = str(e.text)
                    elif name == 'severity':
                        report['severity'] = Severity[str(e.text)]
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
    def __init__(self, cls=None, code=0, details='', name=None, label=None, stack=True, exc_info=None):
        CErrorReport.__init__(self)
        Exception.__init__(self)
        if cls is not None:
            self.append(cls, code, details, name, label, False, stack, exc_info)


def getStack(exc_info=None):
    if exc_info is None:
        return traceback.format_stack()[0:-2] or None
    try:
        return traceback.format_exception(exc_info[0], exc_info[1], exc_info[2])
    except:
        return None


def errorCodeSeverity(class_, code, default=-1):  # KJS - Revise
    if isclass(class_):
        for cls in class_.__mro__:
            if hasattr(cls, 'ERROR_CODES') and code in cls.ERROR_CODES:
                return cls.ERROR_CODES[code].get('severity', Severity.ERROR)
    return default


def errorCodeDescription(class_, code, default=-1):  # KJS : Needs revision.
    if isclass(class_):
        for cls in class_.__mro__:
            if hasattr(cls, 'ERROR_CODES') and code in cls.ERROR_CODES:
                return cls.ERROR_CODES[code].get('description', 'NO DESCRIPTION')
    return default


#===========================================================================================
def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testError)
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)

class testError(unittest.TestCase):

    def test1(self):
        e = CException(CCP4Data.CData, 1, 'foo')
        tree = e.getEtree()
        f = CException()
        f.setEtree(tree)
        self.assertEqual(f[0]['code'],1,'Error save/restore CException to etree = wrong code')
