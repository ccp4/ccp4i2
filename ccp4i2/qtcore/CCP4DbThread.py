from __future__ import print_function
import sys

from PySide2 import QtCore
class CDbThread(QtCore.QThread):
    insts = None
    databaseCalls = [
        'getJobsWithOutputFiles',
        'getProjectJobListInfo',
        'getProjectDirectory',
        'getFileInfo',
        'getProjectFiles',
        'getFullPath',
        'getAllInfo',
        'getProjectJobFile',
        'listProjects',
        'updateProject'
    ]
    def __init__(self, fileName=None, *arg, **kw):
        super().__init__(*arg, **kw)
        #super().__init__(self)
        from utils.startup import startDb
        self.db = startDb(fileName=fileName)
        if sys.version_info >= (3,0):
            import queue
            self.queue = queue.Queue()
        else:
            import Queue
            self.queue = Queue.Queue()
        self.setObjectName('DbServer')
        CDbThread.insts = self

    def run(self):
        while True:
            callback = self.queue.get(True) #does block
            if callback == "ShutdownSignal":
                self.queue.task_done()
                break
            else:
                response = self.handleRequest(callback['path'])
                callback['responseQueue'].put(response)
                self.queue.task_done()
        print('CDbThread shut down')

    def handleRequest(self, requestPath):
        tokens = requestPath.split('?')
        tokensDict = {}
        for token in tokens:
            splitToken = token.strip().split('=')
            if len(splitToken) == 1: tokensDict[splitToken[0]] = True
            else: tokensDict[splitToken[0]] = '='.join(splitToken[1:])

        kwargs = {}
        
        for token in tokens[2:]:
            subtokens=token.split("=")
            if len(subtokens) == 2:
                if subtokens[1] == 'True':
                    kwargs[subtokens[0]]=True
                elif subtokens[1] == 'False':
                    kwargs[subtokens[0]]=False
                else:
                    kwargs[subtokens[0]]=subtokens[1]
            else:
                print('\n\n** parsing dodgy token:' + token)
        response = (getattr(self.db, tokens[1]))(**kwargs)
        return response
