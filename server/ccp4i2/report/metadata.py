"""
Report metadata elements.

Title, JobDetails, JobLogFiles, Reference, ReferenceGroup, GenericReport.

These elements display job metadata (times, status, references) in reports.
"""

from __future__ import annotations

import logging
import re
import xml.etree.ElementTree as etree
from typing import Any

from ccp4i2.core.CCP4ErrorHandling import CException
from ccp4i2 import I2_TOP
from ccp4i2.report.core import ReportClass, Container, Report

logger = logging.getLogger(__name__)


class Title(ReportClass):
    """Job title bar showing job number, task name, user title, and timestamp."""

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        super().__init__()
        if jobInfo is None:
            jobInfo = {}
        import time

        self.title0: str = 'Job ' + \
            str(jobInfo['jobnumber']) + ': ' + jobInfo['tasktitle']
        if jobInfo.get(
                'jobtitle',
                None) is not None and len(
                jobInfo['jobtitle']) > 0:
            self.title1: str | None = jobInfo['jobtitle']
        else:
            self.title1 = None

        self.title2: str = time.strftime(
            '%H:%M %d-%b-%Y',
            time.localtime(
                jobInfo['creationtime']))

    def as_data_etree(self) -> etree.Element:
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
    """Job details section showing creation/finish times and status."""

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        super().__init__()
        if jobInfo is None:
            jobInfo = {}
        self.id: str | None = kw.get('id', None)
        self.class_: str | None = kw.get('class_', None)
        self.jobInfo: dict[str, Any] = {}
        self.jobInfo.update(jobInfo)

    def as_data_etree(self) -> etree.Element:
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
    """Job log files section showing creation/finish times and status."""

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        super().__init__()
        if jobInfo is None:
            jobInfo = {}
        self.id: str | None = kw.get('id', None)
        self.class_: str | None = kw.get('class_', None)
        self.jobInfo: dict[str, Any] = {}
        self.jobInfo.update(jobInfo)

    def as_data_etree(self) -> etree.Element:
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
    """Fallback report for tasks that ship no report class of their own.

    Many tasks — the whole ``Import*`` family, ``mtzheader``, ``pisa`` and
    friends — have no ``*_report.py``.  Qt-i2 still showed them a plain report
    listing what went in and what came out; the Django port instead showed
    "No report because: No report class found for task", which reads as a
    failure rather than as "this task has nothing to plot".

    There is no program XML to parse for these tasks, hence
    ``USEPROGRAMXML = False``.  Everything the report shows — title, input and
    output files, job details, log files, references — is added by
    :meth:`Report.standardisePythonReport`, so this class only needs to name
    the task and, while the job is still running or has failed, say so.
    """

    USEPROGRAMXML = False

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        # TASKNAME is a class attribute on hand-written reports; set it per
        # instance so ReferenceGroup.loadFromMedLine() finds this task's
        # bibliography, and so the class stays usable for any task.
        self.TASKNAME = jobInfo.get('taskname', '') or ''
        # A standardised report gets a Title bar, but that bar only carries the
        # job's own title - it renders nothing when the job was never given one.
        # Name the task here in that case, and whenever there is no bar at all.
        title = jobInfo.get('tasktitle', '') or self.TASKNAME
        if title and not (self.standardise and jobInfo.get('jobtitle')):
            self.addText(text=title, style='font-size:1.1em;')
        status = kw.get('jobStatus') or jobInfo.get('status')
        if status and status != 'Finished':
            self.addText(text=f'Job status: {status}')


class Reference(ReportClass):
    """A single bibliographic reference (article title, authors, source, link)."""

    ERROR_CODES: dict = {}

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        super().__init__()
        if jobInfo is None:
            jobInfo = {}
        self.id: str | None = kw.get('id', None)
        data = kw
        self.href: str | None = data.get('href', None)
        self.authorList: list[str] = data.get('authorList', [])
        if data.get('author', None) is not None:
            self.authorList.append(data.get('author', None))
        self.source: str | None = data.get('source', None)
        self.articleTitle: str | None = data.get('articleTitle', None)
        self.articleLink: str | None = data.get('articleLink', None)

    def as_data_etree(self) -> etree.Element:
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
    """Group of bibliographic references, loadable from MedLine files."""

    ERROR_CODES: dict = {
        100: {
            'description': 'Failed attempting to load MedLine file - file not found'}}

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        Container.__init__(
            self,
            xmlnode=xmlnode,
            jobInfo=jobInfo,
            **kw)
        self.label: str = 'References'
        self.tag: str = 'div'
        self._class: str = 'bibreference_group'
        self.taskName: str | None = kw.get('taskName', None)

    def loadFromMedLine(self, taskName: str) -> None:
        """Parse a MedLine-format file and populate references."""
        from ccp4i2.core import CCP4Utils as _utils
        # Resolved case-insensitively: report classes name the file by TASKNAME
        # while the bibliography builder names it by citation key, and the two
        # differ in case for AceDRG. See CCP4Utils.findReferenceFile.
        path = _utils.findReferenceFile(taskName)
        if path is None:
            path = I2_TOP / "references" / f"{taskName}.medline.txt"
            self.errReport.append(
                self.__class__,
                100,
                f'Taskname: {taskName} Filename: {path}')
            # errReport alone is not enough: nothing reads it, so this used to
            # render as a silently empty <CCP4i2ReportReferenceGroup/> and the
            # AceDRG-on-Linux bug sat undetected for as long as the file existed.
            # Tasks with no citable upstream program are declared in
            # core.citations.NON_CITABLE and stay quiet; anything else is a real
            # gap and says so once, at WARNING, naming the task.
            from ccp4i2.core.citations import NON_CITABLE
            if taskName not in NON_CITABLE:
                logger.warning(
                    "No reference file for report task %r (looked for %s) - "
                    "this report's bibliography will be empty. Add the file, "
                    "add the task to bibliography._NON_CITABLE, or - if the "
                    "citation lives under another key - see "
                    "docs/error-handling-remediation.md#bibliography",
                    taskName, path.name)
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
