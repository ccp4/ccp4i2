# Copyright (C) 2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
"""
Report metadata elements.

Title, JobDetails, JobLogFiles, Reference, ReferenceGroup, GenericReport.

These elements display job metadata (times, status, references) in reports.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as etree
from typing import Any

from ccp4i2.core.CCP4ErrorHandling import CException
from ccp4i2 import I2_TOP
from ccp4i2.report.core import ReportClass, Container, Report


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
    """Fallback report for tasks without a custom report class."""

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        title = jobInfo.get('tasktitle', '')
        self.addText(text=title)


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
