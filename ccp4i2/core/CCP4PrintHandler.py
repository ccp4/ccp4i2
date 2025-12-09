from __future__ import print_function

import os
import sys
import glob
import time
import shutil
import threading

from PySide2 import QtCore

from ccp4i2.core import CCP4Config
from ccp4i2.core import CCP4Utils

class CPrintHandler:
    insts = None

    def __init__(self, parent=None, tmpDir=None):
        CPrintHandler.insts = self
        self.openFileObjects = {}
        self.threadNames = {}
        self.nThreads = 0
        #if tmpDir is None: tmpDir= tempfile.mkdtemp()
        self.tmpDir = os.path.join(CCP4Utils.getDotDirectory(), 'logs', 'started_' + str(int(time.time())))
        os.mkdir(self.tmpDir)
        sys.__stdout__.write('CPrintHandler saving print output to directory: ' + self.tmpDir + '\n')
        sys.__stdout__.flush()

    def cleanupLogs(self):
        # KJS : Fixed. Now properly removes logfiles over one week old in /<area>/.CCP4/logs/
        oneweek = 604800.0  # seconds
        logDirs = glob.glob(os.path.join(CCP4Utils.getDotDirectory(), 'logs', 'started_*'))
        for logDir in logDirs:
            time_since_sdmod = time.time() - os.stat(logDir).st_mtime
            if time_since_sdmod > oneweek:
                shutil.rmtree(logDir)

    def getFileObject(self, thread=None, name=None):
        if thread is None:
            thread = threading.currentThread()
        n = self.threadNames.get(thread, None)
        if n is None:
            self.nThreads += 1
            if name is None:
                name = 'thread_' + str(self.nThreads)
            f = open(os.path.join(self.tmpDir, name), "a+")
            self.openFileObjects[name] = f
            self.threadNames[thread] = name
        else:
            return self.openFileObjects[n]
        return f

    def write(self, value):
        f = self.getFileObject()
        f.write(str(value))
        if CCP4Config.CConfig.insts.developer:      # KJS : There seems to be a problem here.
            #sys.__stdout__.write(str(threading.currentThread())+value)
            sys.__stdout__.write(str(value))

    @QtCore.Slot()
    def exit(self):
        for key, f in list(self.openFileObjects.items()):
            f.close()

    def flush(self):
        pass

    def getContent(self, name=None):
        if not os.path.exists(os.path.join(self.tmpDir, str(name))):
            return ''
        f = self.openFileObjects.get(name)
        try:
            f.seek(0)
            text = f.read()
        except:
            pass
        finally:
            f.seek(0, os.SEEK_END)
        return text
