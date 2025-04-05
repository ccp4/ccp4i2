from http.server import SimpleHTTPRequestHandler
from http.server import ThreadingHTTPServer
import glob
import json
import mimetypes
import os
import queue
import re

from PySide2 import QtCore
import dials
import mrparse

from ..core import CCP4Utils
from ..core.CCP4Modules import DBSERVER, PRINTHANDLER
from ..utils.QApp import QTAPPLICATION


DEFAULT_PORT = 43434


def HTTPSERVER(fileName=None):
    if CHTTPServerThread.insts is None:
        CHTTPServerThread(fileName=fileName)
    return CHTTPServerThread.insts


def makeServer(port: int):
      CHTTPRequestHandler.protocol_version = "HTTP/1.0"
      try:
        return ThreadingHTTPServer(('127.0.0.1', port), CHTTPRequestHandler)
      except:
        return None


class CHTTPServerThread(QtCore.QThread):

    insts = None

    def __init__(self, parent=None, parentDir=None, fileName=None):
        if parent is None:
            parent = QTAPPLICATION()
        QtCore.QThread.__init__(self, parent)
        self.setObjectName('HTTPServer')
        CHTTPServerThread.insts = self
        self.port = None
        self.parentDir = parentDir
        self.dbThread = DBSERVER(fileName)
      
    def quitServer(self):
        self.dbThread.queue.put("ShutdownSignal")
        if hasattr(self, "httpd"):
            self.httpd.shutdown()

    def run(self):
      # Code from http://www.linuxjournal.com/content/tech-tip-really-simple-http-server-python
      # Using this short script so that it will only serve localhost

      if self.parentDir is not None:
          os.chdir(self.parentDir)
      self.httpd = makeServer(DEFAULT_PORT)
      if self.httpd is None:
          self.httpd = makeServer(0)

      if self.httpd is not None:
           sa = self.httpd.socket.getsockname()
           message = f"CCP4i2 starting HTTP server on {sa[0]} port {sa[1]}\n"
           f = PRINTHANDLER().getFileObject(thread='HTTPServer',name='HTTPServer')
           f.write(message)
           print(message)
           self.port = sa[1]
           self.httpd.serve_forever()


class CHTTPRequestHandler(SimpleHTTPRequestHandler):
    def end_headers (self):
        self.send_header('Cross-Origin-Opener-Policy', 'same-origin')
        self.send_header('Cross-Origin-Embedder-Policy', 'require-corp')
        SimpleHTTPRequestHandler.end_headers(self)

    #Here I am going to do some hackery to allow the HTTP server to return information about the
    #database
    def do_GET(self):
        #print('CHTTPRequestHandler.do_GET',self.path)
        if (
            '/database/projectid/' in self.path
            and (projectIdMatch := re.search(r'/projectid/([^/]+)', self.path))
            and (jobNumberMatch := re.search(r'/jobnumber/([^/]+)', self.path))
            and (fileMatch := re.search(r'/file/(.+)/?', self.path))
        ):
            oldUrl = (
                "/database/?getProjectJobFile"
                f"?projectId={projectIdMatch.group(1)}"
                f"?fileName={fileMatch.group(1)}"
                f"?jobNumber={jobNumberMatch.group(1)}"
            )
            #print("Old style url",oldUrl)
            return self.do_GET_main(oldUrl)
        return self.do_GET_main(self.path)

    def do_GET_main(self,self_path):
        #print("do_GET_main",self_path)
        if (i := self_path.find("site-packages/dials/static")) > -1:
            # This mangling is done to help with imported projects
            # which might have links to different file locations
            subPath = self_path[i:].removeprefix("site-packages/dials/")
            # Leads to a mixture of path separators on Windows
            newPath = os.path.join(os.path.dirname(dials.__file__), subPath)
            try:
                fileType = mimetypes.guess_type(newPath.split("?")[0])[0]
                self.returnFileAsData(contentType=fileType, fullPath=newPath.split("?")[0])
            except:
                self.send_response(404)
            return
        if (i := self_path.find("site-packages/mrparse/html")) > -1:
            # This mangling is done to help with imported projects
            # which might have links to different file locations
            subPath = self_path[i:].removeprefix("site-packages/mrparse/")
            # Leads to a mixture of path separators on Windows
            newPath = os.path.join(os.path.dirname(mrparse.__file__), subPath)
            try:
                fileType = mimetypes.guess_type(newPath.split("?")[0])[0]
                self.returnFileAsData(contentType=fileType, fullPath=newPath.split("?")[0])
            except:
                self.send_response(404)
            return
        if not self_path.startswith('/database'):
            #serve files normally
            #print("Normal service")
            SimpleHTTPRequestHandler.do_GET(self)
        else:
            dbQueue = HTTPSERVER().dbThread.queue

            #Create my own Queue for reading response
            dbRequest = {'responseQueue':queue.Queue()}

            newPath = None
            if self_path.startswith("/database/projectId/"):
                if("Referer" in self.headers):
                    a = self_path
                    projectId = a[len("/database/projectId/"):a.find("/jobNumber")]
                    jobNumber = a[a.find("/jobNumber")+len("/jobNumber/"):a.find("/fileName")]
                    fileName = a[a.find("/fileName")+len("/fileName/"):]
                    if "?" in fileName:
                        fileName = fileName[:fileName.find("?")]
                    path = self.headers["Referer"]
                    path = path[:path.find('/database')]
                    dbRequest['path'] = path + "/database/?getProjectJobFileName?projectId="+projectId+"?fileName="+fileName+"?jobNumber="+jobNumber
                    dbQueue.put(dbRequest)
                    response = dbRequest['responseQueue'].get(True, 10)
                    newPath = response

            elif self_path.startswith("/database/") and not self_path.startswith("/database/?"):
                if("Referer" in self.headers):
                    refTokens = self.headers["Referer"].split('?')
                    if refTokens[1] == "getProjectJobFile":
                        dbRequest['path'] = self.headers["Referer"].replace("getProjectJobFile","getProjectJobFileName")
                        dbQueue.put(dbRequest)
                        response = dbRequest['responseQueue'].get(True, 10)
                        newPath = os.path.join(os.path.dirname(response),self_path.replace("/database/",""))

            #print("newPath",newPath)
            if newPath:
                        fileType = 'text/plain'
                        if newPath.lower().endswith(".pdb") or newPath.lower().endswith(".ent"):
                            fileType = 'chemical/x-pdb'
                            self.returnFileAsData(contentType=fileType, fullPath=newPath)
                            dbRequest['responseQueue'].task_done()
                            return
                        elif newPath.lower().endswith(".pdf"):
                            fileType = 'application/pdf'
                            self.returnFileAsDownload(contentType=fileType, fullPath=newPath)
                            dbRequest['responseQueue'].task_done()
                            return
                        elif newPath.lower().endswith(".jpg") or newPath.lower().endswith(".jpeg"):
                            fileType = 'image/jpeg'
                            self.returnFileAsData(contentType=fileType, fullPath=fullPath)
                            dbRequest['responseQueue'].task_done()
                            return
                        elif newPath.lower().endswith(".png"):
                            fileType = 'image/png'
                            self.returnFileAsData(contentType=fileType, fullPath=newPath,isBinary=True)
                            dbRequest['responseQueue'].task_done()
                            return
                        elif newPath.lower().endswith(".svg"):
                            fileType = 'image/svg+xml'
                            self.returnFileAsData(contentType=fileType, fullPath=fullPath)
                            dbRequest['responseQueue'].task_done()
                            return
                        elif newPath.lower().endswith(".html") or newPath.lower().endswith(".htm"):
                            fileType = 'text/html; charset=UTF-8'
                            self.returnFileAsData(contentType=fileType, fullPath=newPath)
                        elif newPath.lower().endswith(".css"):
                            fileType = 'text/css'
                            self.returnFileAsData(contentType=fileType, fullPath=newPath)
                        self.returnFileAsDownload(contentType=fileType, fullPath=newPath)
                        dbRequest['responseQueue'].task_done()
                        return

            tokens = self_path.split('?')
            tokensDict = {}
            for token in tokens:
                splitToken = token.strip().split('=')
                if len(splitToken) == 1: tokensDict[splitToken[0]] = True
                else: tokensDict[splitToken[0]] = '='.join(splitToken[1:])
            print(tokensDict)

            if len(tokens) == 1:
                docTemplate = '''
                    <!DOCTYPE html>
                    <html xmlns="http://www.w3.org/1999/xhtml" lang="en" xml:lang="en">
                    <head>
                    <meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
                    <script>
                    var serverName = '__serverName__';
                    var serverAddress = '127.0.0.1:__serverPort__';
                    </script>
                    
                    <!--CCP4i2 javscript and css-->
                    <link rel="stylesheet" type="text/css" href="/report_files/0.1.0/xreport.css"></link>
                    
                    <link rel="stylesheet" type="text/css" href="/report_files/0.1.0/dist/themes/default/style.min.css"></link>
                    
                    <!-- requirejs javascript -->
                    <script src="/report_files/0.1.0/require.min.js" data-main="/report_files/0.1.0/mainApp"></script>
                    
                    <!-- jspimple javascript -->
                    <link rel="stylesheet" type="text/css" href="/report_files/0.1.0/jspimple.css"></link>

                    <!-- jquery-ui javascript and css -->
                    <script src="/report_files/0.1.0/jquery-ui-1.10.4.custom/js/jquery-ui-1.10.4.custom.min.js"></script>
                    <link rel="stylesheet" href="/report_files/0.1.0/jquery-ui-1.10.4.custom/css/ui-lightness/jquery-ui-1.10.4.custom.css"></link>
                    
                    <style>body {font-size:10px;}</style>
                    <title>CCP4i2 Database view</title>
                    </head>
                    <body></body>
                    </html>
                    '''
                docString1 = re.sub('__serverName__',str(self.server.server_name),docTemplate)
                docString = re.sub('__serverPort__', str(self.server.server_port),docString1)
                self.returnData('text/html; charset=UTF-8', docString)
                # SJM 18/07/2018 - Commented out this, no other reference to db anywhere...
                # db.close()
                return

            from ..qtcore.CCP4DbThread import CDbThread
            if tokens[1] in CDbThread.databaseCalls:
                dbRequest['path'] = self_path
                dbQueue.put(dbRequest)
                response = dbRequest['responseQueue'].get(True, 10)
                isLog = False
                isPng = False
                isSvg = False
                isCss = False
                isJs = False
                isJson = False
                theFileName = ""
                for q in tokens[1:]:
                    if q.startswith("fileName="):
                        if q.split("=")[1] != "report.html" and q.split("=")[1] !=  "report_tmp.html":
                            theFileName = q.split("=")[1]
                            if theFileName.endswith(".txt") or theFileName.endswith(".log"):
                                isLog = True
                            if theFileName.endswith(".png") or theFileName.endswith(".PNG"):
                                isPng = True
                            if theFileName.endswith(".svg") or theFileName.endswith(".SVG"):
                                isSvg = True
                            if theFileName.endswith(".js") or theFileName.endswith(".JS"):
                                isJs = True
                            if theFileName.endswith(".json") or theFileName.endswith(".JSON"):
                                isJson = True
                            if theFileName.endswith(".css") or theFileName.endswith(".CSS"):
                                isCss = True
                if response is None:
                    self.send_response(404)
                    dbRequest['responseQueue'].task_done()
                    return
                else:
                    #print("Doing something!!!!!!!!",theFileName,isPng,isSvg,isLog,isCss,isJs,isJson)
                    if isPng:
                        self.returnData('image/png',response,isBinary=True)
                    elif isSvg:
                        self.returnData('image/svg+xml',response,isBinary=True)
                    elif hasattr(response,"startswith") and (response.strip().startswith(b'<!DOCTYPE html') or response.strip().startswith(b'<HTML') or response.strip().startswith(b'<html')):
                        #Hopefully this is html
                        self.returnData('text/html; charset=UTF-8',response,isBinary=True)
                    elif isLog:
                        self.returnData('text/plain',response,isBinary=True)
                    elif isCss:
                        self.returnData('text/css',response,isBinary=True)
                    elif isJs:
                        self.returnData('text/javascript',response,isBinary=True)
                    elif isJson:
                        self.returnData('application/json',response,isBinary=True)
                    else:
                        self.returnData('application/javascript',json.dumps(response),isBinary=True)
                    dbRequest['responseQueue'].task_done()
                    return
    
            if 'File' in tokensDict:
                fileId = tokensDict.get('fileId', None)
                jobId = tokensDict.get('jobId', None)
                filePath = tokensDict.get('filePath', None)
                iconType = tokensDict.get('icon', None)
                #print 'Digested:',fileId, jobId, filePath
                if fileId is not None:
                    #Get info about file
                    dbRequest['path'] = '/database?getFileInfo?mode=filetype?fileId='+str(fileId)
                    dbQueue.put(dbRequest)
                    fileType = dbRequest['responseQueue'].get(True, 10)
                    if fileType is None:
                        self.send_response(404)
                        dbRequest['responseQueue'].task_done()
                        return
                    dbRequest['responseQueue'].task_done()
                    dbRequest['path'] = '/database?getFullPath?fileId='+str(fileId)
                    dbQueue.put(dbRequest)
                    fullPath = dbRequest['responseQueue'].get(True, 10)
                    #if 'Download' in tokensDict: self.returnFileAsFile(fileType, )
                    if 'Download' in tokensDict: self.returnFileAsDownload(contentType=fileType, fullPath=fullPath)
                    else: self.returnFileAsData(contentType=fileType, fullPath=fullPath)
                    dbRequest['responseQueue'].task_done()
                    return
            
                elif jobId is not None and filePath is not None:
                    dbRequest['path'] = '/database?jobDirectory?jobId='+str(jobId)
                    dbQueue.put(dbRequest)
                    jobDirectory = dbRequest['responseQueue'].get(True, 10)
                    fullPath = None
                    if filePath.startswith(jobDirectory):
                        fullPath = filePath
                    else:
                        fullPath = os.path.join(jobDirectory,filePath)
                    fileType = mimetypes.guess_type(fullPath)[0]
                    #Here apply patches where content type is explicitly known
                    if filePath == 'report.html': fileType = 'application/xhtml+xml'
                    if 'Download' in tokensDict: self.returnFileAsDownload(contentType=fileType, fullPath=fullPath)
                    else: self.returnFileAsData(contentType=fileType, fullPath=fullPath)
                    dbRequest['responseQueue'].task_done()
                    return
                        
                elif iconType is not None:
                    from ..dbapi import CCP4DbApi
                    itype = 0
                    try:
                        itype = CCP4DbApi.FILETYPES_TEXT.index(iconType)
                    except:
                        pass
                    pattern = os.path.join(CCP4Utils.getCCP4I2Dir(), 'qticons','')+CCP4DbApi.FILETYPES_CLASS[itype]+'.*'
                    possibleFiles = glob.glob(pattern)
                    if len(possibleFiles) > 0:
                        fullPath = possibleFiles[0]
                        fileType = mimetypes.guess_type(fullPath)[0]
                        self.returnFileAsDownload(contentType=fileType, fullPath=fullPath)
                    else:
                        self.send_response(404)
                    return
                else:
                    self.send_response(404)
                    return

    def returnFileAsDownload(self, contentType='application/xhtml+xml', fullPath=None):
        if not os.path.isfile(fullPath):
            self.send_response(404)
            return
        with open(fullPath, 'rb') as content_file:
            content = content_file.read()
            self.send_response(200)
            self.send_header('Content-type',contentType)
            self.send_header('Content-disposition','attachment;filename='+os.path.split(fullPath)[1])
            self.end_headers()
            self.wfile.write(content)

    def returnFileAsData(self, contentType='application/xhtml+xml', fullPath=None, isBinary=False):
        if not os.path.isfile(fullPath):
            self.send_response(404)
            return
        if isBinary:
            with open(fullPath, 'rb') as content_file:
                content = content_file.read()
                self.returnData(contentType, content, isBinary=True)
        else:
            with open(fullPath, 'r') as content_file:
                content = content_file.read()
                self.returnData(contentType, content)

    def returnData(self, contentType=None, data=None, isBinary=False):
        if contentType is not None and data is not None:
            self.send_response(200)
            self.send_header('Content-type',contentType)
            self.end_headers()
            if not isBinary:
                self.wfile.write(bytes(data,"utf-8"))
            else:
                self.wfile.write(data)
