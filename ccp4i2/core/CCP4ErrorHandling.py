"Exception class and definitions of severity"

from enum import Enum
from inspect import isclass
from time import localtime, mktime
import traceback
import xml.etree.ElementTree as ET


class Severity(Enum):
    OK = 0, "OK"
    UNDEFINED = 1, "WARNING DATA UNDEFINED"
    WARNING = 2, "WARNING"
    UNDEFINED_ERROR = 3, "ERROR DATA UNDEFINED"
    ERROR = 4, "ERROR"
    CRITICAL = 5, "CRITICAL"

    def __str__(self):
        return self.value[1]

    def __gt__(self, other):
        if self.__class__ == other.__class__:
            return self.value > other.value
        return self.value[0] > other

    def __ge__(self, other):
        if self.__class__ == other.__class__:
            return self.value >= other.value
        return self.value[0] >= other

    @classmethod
    def __getitem__(cls, text: str):
        return {s.value[1]: s for s in cls}[text]


def getStack(exc_info=None):
    if exc_info is None:
        return traceback.format_stack()[0:-2] or None
    try:
        return traceback.format_exception(exc_info[0], exc_info[1], exc_info[2])
    except:
        return None


def getReport(
    cls: object = None,
    className: str = None,
    code: int = 0,
    name: str = None,
    label=None,
    details: str = None,
    description: str = None,
    severity: Severity = None,
    time: float = None,
    stack: list[str] = None,
):
    report = {
        'className': className or getattr(cls, "__name__", str(cls)),
        'code': code,
        'name': name,
        'label': label,
        'details': details,
        'description': description or 'No description available',
        'severity': severity or Severity.ERROR,
        'time': time,
        'stack': stack,
    }
    if isclass(cls):
        for class_ in cls.__mro__:
            if hasattr(class_, 'ERROR_CODES') and code in class_.ERROR_CODES:
                error = class_.ERROR_CODES[code]
                if description is None and 'description' in error:
                    report['description'] = error['description']
                if severity is None and 'severity' in error:
                    report['severity'] = error['severity']
                break
        else:
            print('ERROR in CErrorHandling - error code not found', code, cls)
    return report


class CErrorReport():
    '''Holds list of errors and warnings'''
    def __init__(self, cls=None, code=0, details=None, name=None, label=None, recordTime=False, stack=True, exc_info=None):
        self._reports = []
        if cls is not None:
            self.append(cls=cls, code=code, details=details, name=name, label=label, recordTime=recordTime, stack=stack, exc_info=exc_info)

    def append(self, cls=None, code=0, details=None, name=None, label=None, recordTime=False, stack=True, exc_info=None):
        report = getReport(
            cls=cls,
            code=code,
            details=details,
            name=name,
            label=label,
            time=mktime(localtime()) if recordTime else None,
            stack=getStack(exc_info) if stack else None,
        )
        self._reports.append(report)
        if report['severity'] == Severity.CRITICAL:
            print(self.report())

    def extend(self, other=None, recordTime=False, stack=True):
        if other is None or len(other) == 0 or not hasattr(other,"_reports"):
            return
        for item in other._reports:
            if recordTime and item.get('time', None) is None:
                item['time'] = mktime(localtime())
            if stack and item.get('stack', None) is None:
                item['stack'] = getStack()
        self._reports.extend(other._reports)

    def appendPythonException(self, cls=None, exception=None):
        report = getReport(cls=cls, name='Python error', code=-2, details=str(exception))
        self._reports.append(report)

    def count(self, cls=None, code=None):
        className = getattr(cls, "__name__", str(cls))
        n = 0
        for report in self._reports:
            if (cls is None or report['className'] == className) and (code is None or report['code'] == code):
                n = n + 1
        return n

    def appendDetails(self, details=''):
        for report in self._reports:
            report['details'] += details

    def setName(self, name=''):
        for report in self._reports:
            report['name'] = name

    def maxSeverity(self):
        return max((r['severity'] for r in self._reports), default=Severity.OK)

    def __len__(self):
        return len(self._reports)

    def __getitem__(self, arg: int):
        return self._reports[arg]

    def __str__(self):
        return self.report()

    def description(self, report=None, user=False):
        if report is None:
            report = self._reports[0]
        desc = report['description']
        if user and report.get('label', None) is not None and report['label'] is not NotImplemented:
            desc = str(report['label']) + ': ' + desc
        elif 'name' in report and report['name'] is not None:
            desc = str(report['name']) + ': ' + str(desc)
        return desc, report['severity']

    def report(self, user=False, ifStack=True, mode=0, minSeverity=Severity.UNDEFINED):
        text = ''
        for report in self._reports:
            desc, severity = self.description(report, user=user)
            className = report['className']
            code = report['code']
            if severity == Severity.CRITICAL:
                text = text + "\nCRITICAL ERROR PLEASE REPORT TO CCP4:"
            if severity >= minSeverity:
                name = str(report['name'])
                if mode == 0:
                    text += f"\n{name:20} -{severity}- {className}:{code} {desc}"
                elif mode == 1:
                    text += f"\n{name:20} -{severity}- \n{className}:{code} {desc}"
                else:
                    if user and severity == Severity.WARNING:
                        text += f"\nWarning: {desc}"
                    if (len(text) + len(desc)) < 60:
                        text += f" {desc}"
                    else:
                        text += f"\n{desc}"
                details = str(report['details'])
                if len(details) > 0 and details != 'None':
                    if mode == 1:
                        text += f' {details}'
                    elif (len(text) + len(details)) < 60:
                        text += f' {details}\n'
                    else:
                        text += f'\n{details}\n'
                if ifStack and report['stack'] is not None:
                    text += '\n'
                    text += ''.join(report['stack'])
        return text[1:]

    def getEtree(self):
        element = ET.Element('errorReportList')
        for item in self._reports:
            ele = ET.Element('errorReport')
            e = ET.Element('className')
            e.text = item['className']
            ele.append(e)
            e = ET.Element('code')
            e.text = str(item['code'])
            ele.append(e)
            e = ET.Element('description')
            desc, severity = self.description(item)
            e.text = desc
            ele.append(e)
            e = ET.Element('severity')
            e.text = str(severity)
            ele.append(e)
            if item['details'] is not None:
                e = ET.Element('details')
                e.text = str(item['details'])
                ele.append(e)
            if item['time'] is not None:
                e = ET.Element('time')
                e.text = str(item['time'])
                ele.append(e)
            if item['stack'] is not None:
                e = ET.Element('stack')
                e.text = ''.join(item['stack'])
                ele.append(e)
            element.append(ele)
        return element

    def setEtree(self, element=None):
        body = element.find('ccp4i2_body') or element
        for ele in body:
            if str(ele.tag) == 'errorReport':
                report = {}
                for e in ele.iterchildren():
                    name = str(e.tag)
                    if name == 'code':
                        report['code'] = int(e.text) if e.text and e.text.isdecimal() else 0
                    elif name == 'severity':
                        report['severity'] = Severity[str(e.text)]
                    elif name in {'className', 'stack', 'details', 'description'}:
                        report[name] = str(e.text)
                self._reports.append(report)


class CException(CErrorReport, Exception):
    def __init__(self, cls=None, code=0, details='', name=None, label=None, stack=True, exc_info=None):
        CErrorReport.__init__(self, cls, code, details, name, label, stack=stack, exc_info=exc_info)
        Exception.__init__(self)
