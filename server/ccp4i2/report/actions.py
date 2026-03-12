"""
Report action elements.

Help, Launch, CopyToClipboard, CopyUrlToClipboard, Download, LaunchTask.

These elements represent interactive actions that appear in reports —
buttons for launching viewers, downloading data, or copying information.
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
    Used via ``Container.addHelp()`` or inline in XRT templates.
    """

    def __init__(
        self,
        xrtnode: etree.Element | None = None,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        self.id: str | None = kw.get('id', None)
        if xrtnode is not None:
            self.ref: str | None = xrtnode.get('ref', None)
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
        if xrtnode is not None:
            self.label: str = xrtnode.get(
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
        xrtnode: etree.Element | None = None,
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

        if xrtnode is not None:
            self.exe = xrtnode.get('exe', None)
            self.label = xrtnode.get('label', None)
            if xrtnode.get('ccp4_data_id', None) is not None:
                self.ccp4_data_id.append(xrtnode.get('ccp4_data_id'))
            self.sceneFile = xrtnode.get('sceneFile', None)
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
        xrtnode: etree.Element | None = None,
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

        if xrtnode is not None:
            self.dataName = xrtnode.get('dataName', None)
            self.label = xrtnode.get('label', None)
        self.dataName = kw.get('dataName', self.dataName)
        self.label = kw.get('label', self.label)


class LaunchTask:
    """Button to launch a follow-on task from the report.

    Note: Currently unused — no callers exist in the codebase.
    """

    counter: int = 0

    def __init__(
        self,
        xrtnode: etree.Element | None = None,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        LaunchTask.counter += 1
        self.id: str | None = kw.get('id', None)
        self.jobId: str | None = jobInfo.get('jobid', None)
        if xrtnode is not None:
            self.taskName: str | None = xrtnode.get('taskName', None)
            self.label: str | None = xrtnode.get('label', None)
            self.ccp4_data_id: str | None = xrtnode.get('ccp4_data_id', ccp4_data_id)

        self.taskName = kw.get('taskName', self.taskName)
        self.label = kw.get('label', self.label)
        self.ccp4_data_id = kw.get('ccp4_data_id', self.ccp4_data_id)
