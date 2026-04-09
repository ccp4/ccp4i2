"""
Report action elements.

Help, Launch, CopyToClipboard, CopyUrlToClipboard, Download, LaunchTask, FileLink.

These elements represent interactive actions that appear in reports —
buttons for launching viewers, downloading data, copying information,
or linking to files in the job directory.
"""

from __future__ import annotations

import os
import re
import sys
import xml.etree.ElementTree as etree
from typing import Any


class Help:
    """Help button linking to documentation.

    Resolves ``$CCP4I2/`` paths to the actual CCP4I2 directory.
    Used via ``Container.addHelp()`` or inline in reports.
    """

    def __init__(
        self,
        config: etree.Element | None = None,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        self.id: str | None = kw.get('id', None)
        if config is not None:
            self.ref: str | None = config.get('ref', None)
        else:
            self.ref = kw.get('ref', None)
        if self.ref is not None and self.ref[0] == '$':

            from ccp4i2.core import CCP4Utils
            if sys.platform == "win32":
                # This had better be sane.
                tweak = CCP4Utils.getCCP4I2Dir().replace('\\', '/')
                tweakref = self.ref.replace('\\', '/')              # Ditto
                self.ref = re.sub(r'\$CCP4I2', tweak, tweakref)
                self.ref = os.path.normpath(self.ref)
            else:
                self.ref = re.sub(
                    r'\$CCP4I2', CCP4Utils.getCCP4I2Dir(), self.ref)
        if config is not None:
            self.label: str = config.get(
                'label',
                'About this ' +
                kw.get(
                    'mode',
                    ''))
        else:
            self.label = kw.get('label', 'About this ' + kw.get('mode', ''))


class Launch:
    """Launch button for external viewers (CCP4mg, Coot, loggraph).

    Tracks scene files and CCP4 data IDs for the viewer to open.
    """

    counter: int = 0

    def __init__(
        self,
        config: etree.Element | None = None,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}

        Launch.counter += 1
        self.id: str | None = kw.get('id', None)
        self.jobId: str | None = jobInfo.get('jobid', None)
        self.exe: str | None = None
        self.label: str | None = None
        # This is a list - could be more than one graph
        self.ccp4_data_id: list[str] = []
        self.sceneFile: str | None = None

        if config is not None:
            self.exe = config.get('exe', None)
            self.label = config.get('label', None)
            if config.get('ccp4_data_id', None) is not None:
                self.ccp4_data_id.append(config.get('ccp4_data_id'))
            self.sceneFile = config.get('sceneFile', None)
        self.exe = kw.get('exe', self.exe)
        self.label = kw.get('label', self.label)
        if kw.get('ccp4_data_id', None) is not None:
            self.ccp4_data_id.append(kw.get('ccp4_data_id'))
        # Use relative path in case project moved - Launcher widget will know
        # jobId and be able to find file
        self.sceneFile = kw.get('sceneFile', self.sceneFile)
        if self.sceneFile is not None:
            self.sceneFile = './' + os.path.split(self.sceneFile)[-1]

    def appendDataId(self, ccp4_data_id: str | None = None) -> None:
        """Append a data ID if not already present."""
        if self.ccp4_data_id.count(ccp4_data_id) == 0:
            self.ccp4_data_id.append(ccp4_data_id)


class CopyToClipboard:
    """Button to copy text to the clipboard."""

    def __init__(self, text: str = "", label: str = "Copy to clipboard", **kw: Any) -> None:
        self.text = text
        self.label = label


class CopyUrlToClipboard:
    """Button to copy a job URL to the clipboard."""

    def __init__(self, text: str = "", label: str = "Copy to clipboard", **kw: Any) -> None:
        self.text = text
        self.label = label
        self.projectId: str | None = kw.get("projectId")
        self.jobnumber: int | None = kw.get("jobnumber")


class Download:
    """Download button for report data files."""

    counter: int = 0

    def __init__(
        self,
        config: etree.Element | None = None,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        Download.counter += 1
        self.id: str | None = kw.get('id', None)
        self.jobId: str | None = jobInfo.get('jobid', None)
        self.dataName: str | None = None
        self.label: str | None = None

        if config is not None:
            self.dataName = config.get('dataName', None)
            self.label = config.get('label', None)
        self.dataName = kw.get('dataName', self.dataName)
        self.label = kw.get('label', self.label)


class LaunchTask:
    """Button to launch a follow-on task from the report.

    Note: Currently unused — no callers exist in the codebase.
    """

    counter: int = 0

    def __init__(
        self,
        config: etree.Element | None = None,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        LaunchTask.counter += 1
        self.id: str | None = kw.get('id', None)
        self.jobId: str | None = jobInfo.get('jobid', None)
        self.taskName: str | None = None
        self.label: str | None = None
        self.ccp4_data_id: str | None = None
        if config is not None:
            self.taskName = config.get('taskName', None)
            self.label = config.get('label', None)
            self.ccp4_data_id = config.get('ccp4_data_id', None)

        self.taskName = kw.get('taskName', self.taskName)
        self.label = kw.get('label', self.label)
        self.ccp4_data_id = kw.get('ccp4_data_id', self.ccp4_data_id)


class FileLink:
    """Clickable link to a file in the job directory hierarchy.

    Used by reports to offer inspection of log files, HTML reports, etc.
    The frontend renders this as a button that either opens the file in
    the Monaco preview dialog (text/log files) or in a new browser tab
    (HTML reports).

    Attributes:
        label: Display text for the link (e.g. "Show Pointless logfile")
        relativePath: Path relative to the job directory (e.g. "job_1/log.txt")
        fileType: "text" for Monaco preview, "html" for browser tab
        projectId: Project integer PK (for constructing project_file URL)
        jobId: Job UUID string
    """

    counter: int = 0

    def __init__(
        self,
        config: etree.Element | None = None,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        FileLink.counter += 1
        self.id: str | None = kw.get('id', None)
        self.internalId: str = kw.get(
            'internalId', f'FileLink_{FileLink.counter}')
        self.label: str | None = kw.get('label', None)
        self.relativePath: str | None = kw.get('relativePath', None)
        self.fileType: str = kw.get('fileType', 'text')
        self.projectId: str | None = kw.get('projectId', None)
        self.jobId: str | None = jobInfo.get('jobid', None)

        if config is not None:
            self.label = config.get('label', self.label)
            self.relativePath = config.get('relativePath', self.relativePath)
            self.fileType = config.get('fileType', self.fileType)
            self.projectId = config.get('projectId', self.projectId)

    def as_data_etree(self) -> etree.Element:
        """Serialise to XML for the React frontend."""
        root = etree.Element(
            'CCP4i2ReportFileLink',
            key=self.internalId,
        )
        if self.label is not None:
            root.set('label', self.label)
        if self.relativePath is not None:
            root.set('relativePath', self.relativePath)
        root.set('fileType', self.fileType)
        if self.projectId is not None:
            root.set('projectId', str(self.projectId))
        if self.jobId is not None:
            root.set('jobId', str(self.jobId))
        return root
