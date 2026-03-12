"""
Report metadata elements.

Title, JobDetails, JobLogFiles, Reference, ReferenceGroup, GenericReport.
"""

import re

from ccp4i2.core.CCP4ErrorHandling import CException
from ccp4i2 import I2_TOP
from ccp4i2.report.core import ReportClass, Container, Report


class Title(ReportClass):
    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super().__init__()
        import time

        self.title0 = 'Job ' + \
            str(jobInfo['jobnumber']) + ': ' + jobInfo['tasktitle']
        if jobInfo.get(
                'jobtitle',
                None) is not None and len(
                jobInfo['jobtitle']) > 0:
            self.title1 = jobInfo['jobtitle']
        else:
            self.title1 = None

        self.title2 = time.strftime(
            '%H:%M %d-%b-%Y',
            time.localtime(
                jobInfo['creationtime']))

    def as_data_etree(self):
        root = super().as_data_etree()
        title1 = getattr(self, 'title1', '')
        if title1 is None:
            title1 = ''
        root.set('title1', title1)
        title2 = getattr(self, 'title2', '')
        if title2 is None:
            title2 = ''
        root.set('title2', title2)
        return root


class JobDetails(ReportClass):
    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super().__init__()
        self.id = kw.get('id', None)
        self.class_ = kw.get('class_', None)
        self.jobInfo = {}
        self.jobInfo.update(jobInfo)

    def as_data_etree(self):
        import time
        root = super().as_data_etree()
        root.set(
            'creationtime',
            time.strftime(
                '%H:%M %d-%b-%Y',
                time.localtime(
                    self.jobInfo['creationtime'])))
        root.set(
            'finishtime',
            time.strftime(
                '%H:%M %d-%b-%Y',
                time.localtime(
                    self.jobInfo['finishtime'])))
        root.set('status', self.jobInfo.get('status', 'Unknown'))
        return root



class JobLogFiles(ReportClass):
    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super().__init__()
        self.id = kw.get('id', None)
        self.class_ = kw.get('class_', None)
        self.jobInfo = {}
        self.jobInfo.update(jobInfo)

    def as_data_etree(self):
        import time
        root = super().as_data_etree()
        root.set(
            'creationtime',
            time.strftime(
                '%H:%M %d-%b-%Y',
                time.localtime(
                    self.jobInfo['creationtime'])))
        root.set(
            'finishtime',
            time.strftime(
                '%H:%M %d-%b-%Y',
                time.localtime(
                    self.jobInfo['finishtime'])))
        root.set('status', self.jobInfo.get('status', 'Unknown'))
        return root



class GenericReport(Report):
    def __init__(self, xmlnode=None, jobInfo={}, **kw):
        Report. __init__(self, xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        title = jobInfo.get('tasktitle', '')
        self.addText(text=title)


class Reference(ReportClass):
    ERROR_CODES = {}

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super().__init__()
        self.id = kw.get('id', None)
        if xrtnode is not None:
            data = xrtnode
        else:
            data = kw
        self.href = data.get('href', None)
        self.authorList = data.get('authorList', [])
        if data.get('author', None) is not None:
            self.authorList.append(data.get('author', None))
        self.source = data.get('source', None)
        self.articleTitle = data.get('articleTitle', None)
        self.articleLink = data.get('articleLink', None)

    def as_data_etree(self):
        root = super().as_data_etree()
        if self.articleTitle is not None:
            root.set('articleTitle', self.articleTitle)
        if self.articleLink is not None:
            root.set('articleLink', self.articleLink)
        if self.source is not None:
            root.set('source', self.source)
        if len(self.authorList) > 0:
            root.set('authorList', str(self.authorList))
        return root


class ReferenceGroup(Container):
    ERROR_CODES = {
        100: {
            'description': 'Failed attempting to load MedLine file - file not found'}}

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        Container.__init__(
            self,
            xrtnode=xrtnode,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.label = 'References'
        self.tag = 'div'
        self._class = 'bibreference_group'
        self.taskName = kw.get('taskName', None)

    def loadFromMedLine(self, taskName):
        path = I2_TOP / "references" / f"{taskName}.medline.txt"
        if not path.exists():
            self.errReport.append(
                self.__class__,
                100,
                f'Taskname: {taskName} Filename: {path}')
            return
        self.taskName = taskName

        from ccp4i2.core import CCP4Utils
        try:
            text = CCP4Utils.readFile(fileName=path)
        except CException as e:
            self.errReport.extend(e)
            return
        textList = text.split('\nPMID- ')
        for text in textList:
            ref = Reference()
            m = re.search(r'TI  -(.*)', text)
            if m is not None:
                ref.articleTitle = m.groups()[0].strip()
            m = re.search(r'SO  -(.*)', text)
            if m is not None:
                ref.source = m.groups()[0].strip()
            m = re.findall(r'AU  -(.*)', text)
            for item in m:
                ref.authorList.append(item.strip())
            m = re.search(r'URL -(.*)', text)
            if m is not None:
                ref.articleLink = m.groups()[0].strip()
            if ref.source is not None:
                self.append(ref)
