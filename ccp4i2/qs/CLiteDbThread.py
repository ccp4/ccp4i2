from __future__ import print_function

import sys
import os
import time
import sqlite3
from PyQt4 import QtCore
import queue
if sys.version_info >= (3,0):
    import urllib.parse
else:
    import urlparse
import signal
from dbapi.CCP4DbApi import CDbApi

class CLiteDbThread(QtCore.QThread):
    insts = None
    #Here a dictionary to avoid execution of arbitrary code
    def __init__(self, *arg, **kw):
        #Queue to which to send signal when thread is running
        self.hasStartedQueue = kw.get('hasStartedQueue')
        del kw['hasStartedQueue']
        self.dbfile = kw.get('dbFile', None)
        del kw['dbFile']
        super(CLiteDbThread, self).__init__(*arg, **kw)
        self.databaseCalls = {
                'getProjectsAndJobs':self.getProjectsAndJobs,
                'getFileListForJob':self.getFileListForJob,
                'getProjectFiles':self.getProjectFiles,
                'getJobDirectory':self.getJobDirectory,
        }
        self._queue = queue.Queue()
        self.setObjectName('DbServer')
        self.jobDirectoryCache = {}
        CLiteDbThread.insts = self
    
    def getQueue(self):
        return self._queue
    queue = property(getQueue)
    
    def setDbfile(self, filename):
        self._dbfile = filename
    def getDbfile(self):
        return self._dbfile
    dbfile = property(getDbfile, setDbfile)

    def run(self):
        self.hasStartedQueue.put("Db thread successfully started")
        while True:
            callback = self.queue.get(True) #does block
            try:
                if callback == "ShutdownSignal":
                    self.queue.task_done()
                    break
                else:
                    response = self.handleRequest(callback['path'])
                    callback['responseQueue'].put(response)
                    self.queue.task_done()
            except:
                callback['responseQueue'].put('FailedDbInteraction')

    def shutdown(self):
        self.queue.put("ShutdownSignal")
            
    def performAction(self, action=None, responseQueue=None, timeout=None):
        self.queue.put({'path':action, 'responseQueue': responseQueue})
        response = None
        if responseQueue is not None:
            response = responseQueue.get(True, timeout)
        return response

    def handleRequest(self, requestPath):
        if sys.version_info >= (3,0):
            parsed = urllib.parse.urlparse(requestPath)
            if parsed.path in self.databaseCalls:
                response = (self.databaseCalls[parsed.path])(urllib.parse.parse_qs(parsed.query))
        else:
            parsed = urlparse.urlparse(requestPath)
            if parsed.path in self.databaseCalls:
                response = (self.databaseCalls[parsed.path])(urlparse.parse_qs(parsed.query))
        return response

    def getProjectsAndJobs(self, requestDict):
        try:
            t1 = time.time()
            f = self.dbfile
            conn = sqlite3.connect(f)
            cursor = conn.cursor()
            cursor.execute("SELECT {} FROM Projects ORDER BY ProjectName ASC".format(",".join(CDbApi.PROJECTITEMS)))
            rows = cursor.fetchall()
            projSums = []
            projIds = {}
            projJobs = {}
            projDirs = {}
            totalJobs = 0
            for row in rows:
                t = row[1]
                projIds[row[1]] = row[0]
                cursor.execute("SELECT {} FROM Jobs WHERE ProjectID=?".format(",".join(CDbApi.JOBITEMS)),(row[0],))
                jobs = cursor.fetchall()
                #t += " ("+str(len(jobs))+" jobs)"
                totalJobs += len(jobs)
                projSums.append(t)
                projJobs[row[1]] = jobs
                projDirs[row[1]] = row[6]
            conn.close()
            t2 = time.time()
            print("Time to readDB",t2-t1)
            return {'projSums':projSums, 'projJobs':projJobs,'projDirs':projDirs,'projIds':projIds,}
        except Exception as err:
            print("Exception in getProjectsAndJobs", err)
            raise

    def getFileListForJob(self, requestDict):
        try:
            #urlparse.parse_qs returns an array of matches
            jobId = requestDict['jobId'][0]
            f = self.dbfile
            conn = sqlite3.connect(f)
            cursor = conn.cursor()
            cursor.execute("SELECT * FROM Files WHERE JobID='"+jobId+"'"+";")
            response = cursor.fetchall()
            conn.close()
            return response
        except Exception as err:
            print("Exception in getFileListForJob", err)
            raise

    def getProjectFiles(self, requestDict):
        #print 'getProjectFiles',requestDict
        try:
            projId = requestDict['projId'][0]
            com, args = "SELECT Jobs.TaskName,Jobs.JobNumber, Files.FileId, Files.Annotation,ImportFiles.SourceFilename,Files.FiletypeID,Files.FileSubType,Files.JobParamName, Jobs.CreationTime FROM Jobs INNER JOIN Files ON Files.JobId=Jobs.JobId LEFT OUTER JOIN ImportFiles ON Files.FileId=ImportFiles.FileId WHERE Jobs.ProjectId=?  AND ImportFiles.ImportId IS NULL AND Jobs.ParentJobId IS NULL ORDER BY Jobs.CreationTime DESC;", [projId]
            f = self.dbfile
            conn = sqlite3.connect(f)
            cursor = conn.cursor()
            cursor.execute(com,args)
            files = cursor.fetchall()
            conn.close()
            #print "----------------------------"
            #print "Got",len(files),"files"
            #print "----------------------------"
            return files
        except Exception as err:
            print("Exception in getProjectFiles", err)
            raise

    def getJobDirectory(self, requestDict):
        #print "getJobDirectory", requestDict
        try:
            jobId = requestDict['jobId'][0]
            #Caching jobDirectories....Look out for this if we start to move projects around
            if jobId in self.jobDirectoryCache: return self.jobDirectoryCache[jobId]
            com, args = "SELECT Jobs.JobNumber, Projects.ProjectDirectory FROM Jobs INNER JOIN Projects ON Jobs.ProjectID=Projects.ProjectID WHERE Jobs.JobID=?", [jobId]
            f = self.dbfile
            conn = sqlite3.connect(f)
            cursor = conn.cursor()
            cursor.execute(com,args)
            results = cursor.fetchall()
            conn.close()
            answer = os.path.join(results[0][1], "CCP4_JOBS", "job_{}".format(results[0][0]))
            self.jobDirectoryCache[jobId] = answer
            return answer
        except Exception as err:
            print("Exception in getJobDirectory", err)
            raise

if __name__ == "__main__":
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    application = QtCore.QCoreApplication(sys.argv)
    hasStartedQueue = queue.Queue()
    a=CLiteDbThread(hasStartedQueue=hasStartedQueue)
    a.dbfile = "/Volumes/Home/martin/.CCP4I2/db/database.sqlite"
    a.start()
    callback = {
        'path': 'getProjectsAndJobs',
        'responseQueue':queue.Queue()
    }
    a.queue.put(callback)
    response = callback['responseQueue'].get(True)
    print(response)
