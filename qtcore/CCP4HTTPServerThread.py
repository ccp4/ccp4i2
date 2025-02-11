from __future__ import print_function

import sys,os
import re

if sys.version_info >= (3,0):
    import http.server
    from http.server import SimpleHTTPRequestHandler
else:
    import BaseHTTPServer
    from SimpleHTTPServer import SimpleHTTPRequestHandler
    
from PySide2 import QtCore
from core import CCP4Modules

DEFAULT_PORT = 43434

if sys.version_info >= (3,0):
    from socketserver import ThreadingMixIn
    from http.server import HTTPServer
else:
    from SocketServer import ThreadingMixIn
    from BaseHTTPServer import HTTPServer

class MultiThreadedHTTPServer(ThreadingMixIn, HTTPServer):
    pass

def makeServer(port):
      HandlerClass = CHTTPRequestHandler
      ServerClass  = MultiThreadedHTTPServer
      Protocol     = "HTTP/1.0"
      server_address = ('127.0.0.1',port )
      HandlerClass.protocol_version = Protocol
      try:
        httpd = ServerClass(server_address, HandlerClass)
      except:
        return None
      else:
        return httpd

def testServer(port):
      if sys.version_info >= (3,0):
          import http.client
      else:
          import httplib

      ret = True
      try:
        if sys.version_info >= (3,0):
            h = http.client.HTTPConnection('127.0.0.1',port,timeout=3)
        else:
            h = httplib.HTTPConnection('127.0.0.1',port,timeout=3)
        h.connect()
      except:
        ret=False
      print('testServer',port,h,ret)
      try:
        h.close()
      except:
        pass
      return ret

class CHTTPServerThread(QtCore.QThread):

    insts = None

    def __init__(self,parent=None,parentDir=None,port=None,fileName=None):
      if parent is None: parent = CCP4Modules.QTAPPLICATION()
      QtCore.QThread.__init__(self,parent)
      self.setObjectName('HTTPServer')
      CHTTPServerThread.insts = self
      self.port = None
      if port is not None:
        self.defaultPort = port
      else:
        self.defaultPort = '43434'
      self.parentDir = parentDir
      self.dbThread = CCP4Modules.DBSERVER(fileName)
      
    def quitServer(self):
      self.dbThread.queue.put("ShutdownSignal")
      #sys.__stdout__.write('Try to shut down server cleanly\n');sys.__stdout__.flush()
      if hasattr(self,"httpd"): self.httpd.shutdown()

    def run(self):

      # Code from http://www.linuxjournal.com/content/tech-tip-really-simple-http-server-python
      # Using this short script so that it will only serve localhost

      if self.parentDir is not None: os.chdir(self.parentDir)
      self.httpd =  makeServer(self.defaultPort)  
      if self.httpd is None: self.httpd = makeServer(0)
      
      if self.httpd is not None:
           sa = self.httpd.socket.getsockname()
           f = CCP4Modules.PRINTHANDLER().getFileObject(thread='HTTPServer',name='HTTPServer')
           f.write( "CCP4i2 starting HTTP server on "+str( sa[0])+ " port "+str(sa[1])+'\n' )
           print("CCP4i2 starting HTTP server on "+str( sa[0])+ " port "+str(sa[1])+'\n')
           self.port = sa[1]
           self.httpd.serve_forever()
      return

from dbapi import CCP4DbApi

class CHTTPRequestHandler(SimpleHTTPRequestHandler):
    def end_headers (self):
        self.send_header('Cross-Origin-Opener-Policy', 'same-origin')
        self.send_header('Cross-Origin-Embedder-Policy', 'require-corp')
        SimpleHTTPRequestHandler.end_headers(self)
        
    def log_message(self,format,*args):
      # More programming by stackoverflow...
      # http://stackoverflow.com/questions/10651052/how-to-quiet-simplehttpserver/10651257#10651257
      # The base class code
      '''
      sys.stderr.write("%s - - [%s] %s\n" %
                     (self.address_string(),
                      self.log_date_time_string(),
                      format%args))
      '''
                    
      f = CCP4Modules.PRINTHANDLER().getFileObject(thread='HTTPServer',name='HTTPServer')
      f.write("%s - - [%s] %s\n" %
                     (self.address_string(),
                      self.log_date_time_string(),
                      format%args))
    
    """
    def do_POST(self):
        import cgi
        form = cgi.FieldStorage(fp=self.rfile, headers=self.headers, \
            environ={'REQUEST_METHOD':'POST', 'CONTENT_TYPE':self.headers['Content-Type'], })
        upfile = form['upfile']
        project = form['project']

        # extract basename of input filename, remove non-alphanumeric characters
        if '\\' in upfile.filename:
            filename = upfile.filename.split('\\')[-1]
        else:
            filename = upfile.filename.split('/')[-1]
        filename = re.sub('[ \t]','-', filename)
        filename = re.sub('[^a-zA-Z0-9_.:-]','', filename)

        data = ''
        while True:
            chunk = upfile.file.read(8192)
            if len(chunk) == 0:
                break
            else:
                data += chunk
        self.send_response(200)
        self.end_headers()

        CHTTPServerThread.insts.emit(QtCore.SIGNAL('insertFileInDBRequest'),(os.path.basename(upfile.filename),data,project.value))

        self.wfile.write('<html><head><title>Upload</title></head><body>\
            <h1>Requested uploaded</h1><br>From: %s<br>To project:%s</body></html>' % \
            (upfile.filename, project.value) )
    """

    #Here I am going to do some hackery to allow the HTTP server to return information about the
    #database
    def do_GET(self):
        #print('CHTTPRequestHandler.do_GET',self.path)
        if "/database/projectid/" in self.path and "/jobnumber/" in self.path and "/file/" in self.path:
            #print("Hello!!!!")
            u = self.path
            projidfind = re.compile(r'/projectid/[a-z0-9\-]*/')
            jobnumberfind = re.compile(r'/jobnumber/[0-9\.]*/')
            filefind = re.compile(r'/file/.*')
            if projidfind.search(u) is not None and jobnumberfind.search(u) is not None and filefind.search(u) is not None:
                theProjectId = projidfind.search(u).group()[11:].rstrip('/')
                theJobNumber = jobnumberfind.search(u).group()[11:].rstrip('/')
                theFile = filefind.search(u).group()[6:]
                #print("Found",theProjectId,theJobNumber,theFile)
                old_url =  "/database/?getProjectJobFile?projectId="+theProjectId+"?fileName="+theFile+"?jobNumber="+theJobNumber
                #print("Old style url",old_url)
                return self.do_GET_main(old_url)
        return self.do_GET_main(self.path)

    def do_GET_main(self,self_path):
        if "site-packages/dials/static" in self_path:
            import dials
            f = self_path
            #This mangling is done to help with imported projects which might have links to different file locations
            newPath = os.path.join(os.path.dirname(dials.__file__),f[f.find("site-packages/dials/static"):][len("site-packages/dials")+1:])
            try:
                import mimetypes
                fileType = mimetypes.guess_type(newPath.split("?")[0])[0]
                self.returnFileAsData(contentType=fileType, fullPath=newPath.split("?")[0])
                return
            except:
                self.send_response(404)
                return
        if "site-packages/mrparse/html" in self_path:
            import mrparse
            f = self_path
            #This mangling is done to help with imported projects which might have links to different file locations
            newPath = os.path.join(os.path.dirname(mrparse.__file__),f[f.find("site-packages/mrparse/html"):][len("site-packages/mrparse")+1:])
            try:
                import mimetypes
                fileType = mimetypes.guess_type(newPath.split("?")[0])[0]
                self.returnFileAsData(contentType=fileType, fullPath=newPath.split("?")[0])
                return
            except:
                self.send_response(404)
                return
        if not self_path.startswith('/database'):
            #serve files normally
            #print("Normal service")
            SimpleHTTPRequestHandler.do_GET(self)
        else:
            from core.CCP4Modules import HTTPSERVER
            dbQueue = HTTPSERVER().dbThread.queue
            
            #Create my own Queue for reading response
            if sys.version_info >= (3,0):
                import queue
                dbRequest = {'responseQueue':queue.Queue()}
            else:
                import Queue
                dbRequest = {'responseQueue':Queue.Queue()}

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
                            fileType = 'text/html'
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
            #print(tokensDict)
            import json
            
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
                self.returnData('text/html', docString)
#SJM 18/07/2018 - Commented out this, no other reference to db anywhere...
                #db.close()
                return
            
            from qtcore.CCP4DbThread import CDbThread
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
                        self.returnData('image/svg+xml',response)
                    elif hasattr(response,"startswith") and (response.strip().startswith('<!DOCTYPE html') or response.strip().startswith('<HTML') or response.strip().startswith('<html')):
                        #Hopefully this is html
                        self.returnData('text/html',response)
                    elif isLog:
                        self.returnData('text/plain',response)
                    elif isCss:
                        self.returnData('text/css',response)
                    elif isJs:
                        self.returnData('text/javascript',response)
                    elif isJson:
                        self.returnData('application/json',response)
                    else:
                        print("else .....")
                        self.returnData('application/javascript',json.dumps(response))
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
                    import mimetypes
                    fileType = mimetypes.guess_type(fullPath)[0]
                    #Here apply patches where content type is explicitly known
                    if filePath == 'report.html': fileType = 'application/xhtml+xml'
                    if 'Download' in tokensDict: self.returnFileAsDownload(contentType=fileType, fullPath=fullPath)
                    else: self.returnFileAsData(contentType=fileType, fullPath=fullPath)
                    dbRequest['responseQueue'].task_done()
                    return
                        
                elif iconType is not None:
                    from core import CCP4Utils
                    import glob
                    from dbapi import CCP4DbApi
                    itype = 0
                    try:
                        itype = CCP4DbApi.FILETYPES_TEXT.index(iconType)
                    except:
                        pass
                    pattern = os.path.join(CCP4Utils.getCCP4I2Dir(), 'qticons','')+CCP4DbApi.FILETYPES_CLASS[itype]+'.*'
                    possibleFiles = glob.glob(pattern)
                    if len(possibleFiles) > 0:
                        import mimetypes
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

