"""
Report action elements.

Help, Launch, CopyToClipboard, CopyUrlToClipboard, Download, LaunchTask.
"""

import os
import re
import sys
import xml.etree.ElementTree as etree


class Help:
    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        self.id = kw.get('id', None)
        if xrtnode is not None:
            self.ref = xrtnode.get('ref', None)
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
            self.label = xrtnode.get(
                'label',
                'About this ' +
                kw.get(
                    'mode',
                    ''))
        else:
            self.label = kw.get('label', 'About this ' + kw.get('mode', ''))


class Launch:

    counter = 0

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):

        Launch.counter += 1
        self.id = kw.get('id', None)
        self.jobId = jobInfo.get('jobid', None)
        self.exe = None
        self.label = None
        # This is a list - could be more than one graph
        self.ccp4_data_id = []
        self.sceneFile = None

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

    def appendDataId(self, ccp4_data_id=None):
        if self.ccp4_data_id.count(ccp4_data_id) == 0:
            self.ccp4_data_id.append(ccp4_data_id)


class CopyToClipboard:
    def __init__(self, text="", label="Copy to clipboard", **kw):
        self.text = text
        self.label = label


class CopyUrlToClipboard:
    def __init__(self, text="", label="Copy to clipboard", **kw):
        self.text = text
        self.label = label
        self.projectId = kw.get("projectId")
        self.jobnumber = kw.get("jobnumber")


class Download:

    counter = 0

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        Launch.counter += 1
        self.id = kw.get('id', None)
        self.jobId = jobInfo.get('jobid', None)
        self.dataName = None
        self.label = None

        if xrtnode is not None:
            self.dataName = xrtnode.get('dataName', None)
            self.label = xrtnode.get('label', None)
        self.dataName = kw.get('dataName', self.dataName)
        self.label = kw.get('label', self.label)


class LaunchTask:

    counter = 0

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        Launch.counter += 1
        self.id = kw.get('id', None)
        self.jobId = jobInfo.get('jobid', None)
        if xrtnode is not None:
            self.taskName = xrtnode.get('taskName', None)
            self.label = xrtnode.get('label', None)
            self.ccp4_data_id = xrtnode.get('ccp4_data_id', ccp4_data_id)

        self.taskName = kw.get('taskName', self.taskName)
        self.label = kw.get('label', self.label)
        self.ccp4_data_id = kw.get('ccp4_data_id', self.ccp4_data_id)
