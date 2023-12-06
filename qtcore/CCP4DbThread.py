from __future__ import print_function
from urllib.parse import parse_qs
from urllib.parse import urlparse
import sys
import mimetypes
import traceback

from PySide2 import QtCore
from dbapi import CCP4DjangoApi
from core import CCP4Modules

class CDbThread(QtCore.QThread):
    insts = None
    databaseCalls = [
        'modelValues'
    ]

    def __init__(self, fileName=None, *arg, **kw):
        super().__init__(*arg, **kw)
        # super().__init__(self)
        from utils.startup import startDb
        self.db = startDb(fileName=fileName)
        if sys.version_info >= (3, 0):
            import queue
            self.queue = queue.Queue()
        else:
            import Queue
            self.queue = Queue.Queue()
        self.setObjectName('DbServer')
        CDbThread.insts = self

    def run(self):
        while True:
            request = self.queue.get(True)  # does block
            if request == "ShutdownSignal":
                self.queue.task_done()
                break
            elif request["method"] == "GET":
                response = self.handleGetRequest(request['apiEndpoint'], request['query'])
                request['responseQueue'].put(response)
                self.queue.task_done()
            elif request["method"] == "POST":
                response = self.handlePostRequest(request)
                request['responseQueue'].put(response)
                self.queue.task_done()

        print('CDbThread shut down')

    def handleGetRequest(self, apiEndpoint, queryDict):
        #print('In dbThread handleGetRequest with params', apiEndpoint, queryDict)
        try:
            if apiEndpoint == "ModelValues":
                resultDict = self.handleModelValues(queryDict)
            elif apiEndpoint == "getProjectJobFile":
                resultDict = self.handleGetProjectJobFile(queryDict)
            elif apiEndpoint == "getJobFile":
                resultDict = self.handleGetJobFile(queryDict)
            elif apiEndpoint == "getProjectJobFileName":
                resultDict = self.handleGetProjectJobFileName(queryDict)
            elif apiEndpoint == "getFileWithPredicate":
                resultDict = self.handleGetFileWithPredicate(queryDict)
            elif apiEndpoint == "makeTerminateFile":
                resultDict = self.handleMakeTerminateFile(queryDict)
            elif apiEndpoint == "getProjectFileData":
                resultDict = self.handleGetProjectFileData(queryDict)
            elif apiEndpoint == "getReportXML":
                resultDict = self.handleGetReportXML(queryDict)
            return resultDict
        except Exception as err:
            resultDict = {"status":f"Exception {err}"}
            traceback.print_exc()
            return resultDict
    
    def handlePostRequest(self, request):
        #print('In dbThread handlePostRequest with endPoint', request["apiEndpoint"])
        try:
            if request["apiEndpoint"] == "uploadFileToJob":
                resultDict = self.handleUploadFileToJob(request['fields'])
            if request["apiEndpoint"] == "uploadFileForJobObject":
                resultDict = self.handleUploadFileForJobObject(request['fields'])
            elif request["apiEndpoint"] == "cloneJob":
                resultDict = self.handleCloneJob(request['fields'])
            elif request["apiEndpoint"] == "createJob":
                resultDict = self.handleCreateJob(request['fields'])
            elif request["apiEndpoint"] == "runJob":
                resultDict = self.handleRunJob(request['fields'])
            elif request["apiEndpoint"] == "setJobParameterByXML":
                resultDict = self.handleSetJobParameterByXML(request['fields'])
            return resultDict
        except Exception as err:
            resultDict = {"status":f"Exception {err}"}
            return resultDict
        
    def handleCloneJob(self, queryDict):
        disArrayedQueryDict = {key: queryDict[key][0] for key in queryDict}
        newJobId, newJobNumber, projectId, parentJobId = CCP4DjangoApi.cloneJob(disArrayedQueryDict['jobId'])
        resultDict = {"status":"Success", "jobId":newJobId, "jobNumber":newJobNumber, "projectId":projectId, "parentJobId":parentJobId}
        CCP4Modules.PROJECTSMANAGER().db().jobCreated.emit(
            {'jobId': resultDict['jobId'], 
                'projectId': resultDict['projectId'], 
                'parentJobId': resultDict['parentJobId']})
        return resultDict
    
    def handleCreateJob(self, queryDict):
        disArrayedQueryDict = {key: queryDict[key][0] for key in queryDict}
        newJobId, newJobNumber, projectId, parentJobId = CCP4DjangoApi.createJob(**disArrayedQueryDict)
        resultDict = {"status":"Success", "jobId":newJobId, "jobNumber":newJobNumber, "projectId":projectId, "parentJobId":parentJobId}
        CCP4Modules.PROJECTSMANAGER().db().jobCreated.emit(
            {'jobId': resultDict['jobId'], 
                'projectId': resultDict['projectId'], 
                'parentJobId': resultDict['parentJobId']})
        return resultDict
    
    def handleRunJob(self, queryDict):
        disArrayedQueryDict = {key: queryDict[key][0] for key in queryDict}
        jobId, projectId, status = CCP4DjangoApi.updateJobStatus(jobId=disArrayedQueryDict['jobId'], statusId=2)
        resultDict = {'jobId': jobId, 'projectId': projectId, 'key': 'status', 'value': status}
        #Am a bit surprised I don't need to emit this signal !
        #CCP4Modules.PROJECTSMANAGER().db().jobUpdated.emit(resultDict)
        return resultDict

    def handleSetJobParameterByXML(self, queryDict):
        disArrayedQueryDict = {key: queryDict[key][0] for key in queryDict}
        newValue = CCP4DjangoApi.setJobParameterByXML(**disArrayedQueryDict)
        return {"status":"Success", "newValue":newValue}

    def handleUploadFileToJob(self, queryDict):
        disArrayedQueryDict = {key: queryDict[key][0] for key in queryDict}
        fileDict = CCP4DjangoApi.uploadFileToJob(**disArrayedQueryDict)
        return {"status":"Success", "fileDict":fileDict}

    def handleUploadFileForJobObject(self, queryDict):
        disArrayedQueryDict = {key: queryDict[key][0] for key in queryDict}
        fileDict = CCP4DjangoApi.uploadFileForJobObject(**disArrayedQueryDict)
        return {"status":"Success", "fileDict":fileDict}

    def handleGetProjectJobFile(self, queryDict):
        disArrayedQueryDict = {key: queryDict[key][0] for key in queryDict}
        data = CCP4DjangoApi.getProjectJobFile(**disArrayedQueryDict)
        contentType = mimetypes.guess_type(disArrayedQueryDict['fileName'])[0]
        if queryDict['fileName'] in ['report.html', 'report_tmp.html']:
            contentType = 'application/xhtml+xml'
        isBinary = False
        return {'contentType': contentType, 'data': data, 'isBinary': isBinary}

    def handleGetProjectFileData(self, queryDict):
        disArrayedQueryDict = {key: queryDict[key][0] for key in queryDict}
        data, fileName = CCP4DjangoApi.getProjectFileData(**disArrayedQueryDict)
        contentType = mimetypes.guess_type(fileName)[0]
        if queryDict['baseName'] in ['report.html', 'report_tmp.html']:
            contentType = 'application/xhtml+xml'
        isBinary = True
        return {'contentType': contentType, 'data': data, 'isBinary': isBinary}

    def handleGetFileWithPredicate(self, queryDict):
        disArrayedQueryDict = {key: queryDict[key][0] for key in queryDict}
        data, filePath, contentType = CCP4DjangoApi.getFileWithPredicate(disArrayedQueryDict)
        isBinary = True
        return {'contentType': contentType, 'data': data, 'isBinary': isBinary}

    def handleMakeTerminateFile(self, queryDict):
        disArrayedQueryDict = {key: queryDict[key][0] for key in queryDict}
        data, contentType = CCP4DjangoApi.makeTerminateFile(disArrayedQueryDict)
        return {'contentType': contentType, 'data': data, 'isBinary': False}
    
    def handleGetJobFile(self, queryDict):
        disArrayedQueryDict = {key: queryDict[key][0] for key in queryDict}
        data = CCP4DjangoApi.getJobFile(**disArrayedQueryDict)
        contentType = mimetypes.guess_type(disArrayedQueryDict['fileName'])[0]
        if queryDict['fileName'] in ['report.html', 'report_tmp.html']:
            contentType = 'application/xhtml+xml'
        isBinary = False
        return {'contentType': contentType, 'data': data, 'isBinary': isBinary}
    
    def handleGetProjectJobFileName(self, queryDict):
        disArrayedQueryDict = {key: queryDict[key][0] for key in queryDict}
        fullPath = CCP4DjangoApi.getProjectJobFileName(**disArrayedQueryDict)
        contentType = mimetypes.guess_type(disArrayedQueryDict['fileName'])[0]
        if queryDict['fileName'] in ['report.html', 'report_tmp.html']:
            contentType = 'application/xhtml+xml'
        isBinary = False
        return {'contentType': contentType, 'fullPath': fullPath}

    def handleModelValues(self, queryDict):
        order_byArray = []
        valuesArray = []
        predicate = {}

        for key in queryDict:
            # print(key)
            if key == '__type__':
                objectType = queryDict.get(key)[0]
            elif key == '__order_by__':
                order_byArray.append(queryDict.get(key))
            elif key == '__order_by__[]':
                order_byArray += queryDict.get(key)
            elif key == '__values__':
                valuesArray.append(queryDict.get(key))
            elif key == '__values__[]':
                valuesArray += queryDict.get(key)
            else:
                if len(queryDict.get(key)) > 1 or key.endswith('[]'):
                    nonArrayKey = key[:-2]
                    if nonArrayKey not in predicate:
                        predicate[nonArrayKey] = []
                    for value in queryDict.get(key):
                        if value == "__True__":
                            value = True
                        elif value == "__False__":
                            value = False
                        predicate[nonArrayKey].append(value)
                else:
                    value = queryDict.get(key)[0]
                    if value == "__True__":
                        value = True
                    elif value == "__False__":
                        value = False
                    predicate[key] = value
        resultDict = CCP4DjangoApi.modelValues(
            objectType=objectType, predicate=predicate, order_byArray=order_byArray, valuesArray=valuesArray)
        return resultDict

    def handleGetReportXML(self, queryDict):
        disArrayedQueryDict = {key: queryDict[key][0] for key in queryDict}
        reportXml = CCP4DjangoApi.getReportXML(**disArrayedQueryDict)
        contentType = 'application/xhtml+xml'
        isBinary = False
        return {'contentType': contentType, 'data': reportXml, 'isBinary': isBinary}
    
