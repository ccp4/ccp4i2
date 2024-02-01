from __future__ import print_function
from dbapi import CCP4DbApi
from urllib.parse import parse_qs
from urllib.parse import urlparse
from urllib.parse import unquote

import sys
import os
import queue
import json
import mimetypes
import posixpath


if sys.version_info >= (3, 0):
    import http.server
    from http.server import SimpleHTTPRequestHandler
else:
    import BaseHTTPServer
    from SimpleHTTPServer import SimpleHTTPRequestHandler

from PySide2 import QtCore
from core import CCP4Modules
from core import CCP4Utils

DEFAULT_PORT = 43434

if sys.version_info >= (3, 0):
    from socketserver import ThreadingMixIn
    from http.server import HTTPServer
else:
    from SocketServer import ThreadingMixIn
    from BaseHTTPServer import HTTPServer


class MultiThreadedHTTPServer(ThreadingMixIn, HTTPServer):
    pass


def makeServer(port):
    HandlerClass = CHTTPRequestHandler
    ServerClass = MultiThreadedHTTPServer
    Protocol = "HTTP/1.0"
    server_address = ('127.0.0.1', port)
    HandlerClass.protocol_version = Protocol
    try:
        httpd = ServerClass(server_address, HandlerClass)
    except:
        return None
    else:
        return httpd


def testServer(port):
    if sys.version_info >= (3, 0):
        import http.client
    else:
        import httplib

    ret = True
    try:
        if sys.version_info >= (3, 0):
            h = http.client.HTTPConnection('127.0.0.1', port, timeout=3)
        else:
            h = httplib.HTTPConnection('127.0.0.1', port, timeout=3)
        h.connect()
    except:
        ret = False
    print('testServer', port, h, ret)
    try:
        h.close()
    except:
        pass
    return ret


class CHTTPServerThread(QtCore.QThread):

    insts = None

    def __init__(self, parent=None, parentDir=None, port=None, fileName=None):
        if parent is None:
            parent = CCP4Modules.QTAPPLICATION()
        QtCore.QThread.__init__(self, parent)
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
        # sys.__stdout__.write('Try to shut down server cleanly\n');sys.__stdout__.flush()
        if hasattr(self, "httpd"):
            self.httpd.shutdown()

    def run(self):

        # Code from http://www.linuxjournal.com/content/tech-tip-really-simple-http-server-python
        # Using this short script so that it will only serve localhost

        if self.parentDir is not None:
            os.chdir(self.parentDir)
        self.httpd = makeServer(self.defaultPort)
        if self.httpd is None:
            self.httpd = makeServer(0)

        if self.httpd is not None:
            sa = self.httpd.socket.getsockname()
            f = CCP4Modules.PRINTHANDLER().getFileObject(
                thread='HTTPServer', name='HTTPServer')
            f.write("CCP4i2 starting HTTP server on " +
                    str(sa[0]) + " port "+str(sa[1])+'\n')
            print("CCP4i2 starting HTTP server on " +
                  str(sa[0]) + " port "+str(sa[1])+'\n')
            self.port = sa[1]
            self.httpd.serve_forever()
        return


class CHTTPRequestHandler(SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header('Cross-Origin-Opener-Policy', 'same-origin')
        self.send_header('Cross-Origin-Embedder-Policy', 'require-corp')
        self.send_header('Access-Control-Allow-Origin', '*')

        SimpleHTTPRequestHandler.end_headers(self)

    def translate_path(self, path):
        parsedRequest = urlparse(self.path)
        pathElements = parsedRequest.path[1:].split("/")

        # Factor out some plugin-specific paths
        if "site-packages/dials/static" in path:
            import dials
            fsPath = os.path.join(os.path.dirname(dials.__file__), path[path.find(
                "site-packages/dials/static"):][len("site-packages/dials")+1:])
            return fsPath

        elif "site-packages/mrparse/html" in path:
            import mrparse
            fsPath = os.path.join(os.path.dirname(mrparse.__file__), path[path.find(
                "site-packages/mrparse/html"):][len("site-packages/mrparse")+1:])
            return fsPath

        elif len(pathElements) > 1 and pathElements[-2] == 'svgicons':
            newPath = os.path.join(CCP4Utils.getCCP4I2Dir(
            ), pathElements[-2], f"{pathElements[-1]}.svg")
            if not os.path.exists(newPath):
                newPath = os.path.join(
                    CCP4Utils.getCCP4I2Dir(), pathElements[-2], "ccp4.svg")
            return newPath

        elif len(pathElements) > 1 and pathElements[-2] == 'qticons':
            newPath = os.path.join(CCP4Utils.getCCP4I2Dir(
            ), pathElements[-2], f"{pathElements[-1]}.png")
            if not os.path.exists(newPath):
                newPath = os.path.join(
                    CCP4Utils.getCCP4I2Dir(), pathElements[-2], "ccp4.png")
            return newPath

        return SimpleHTTPRequestHandler.translate_path(self, path)

    def log_message(self, format, *args):
        f = CCP4Modules.PRINTHANDLER().getFileObject(
            thread='HTTPServer', name='HTTPServer')
        f.write("%s - - [%s] %s\n" %
                (self.address_string(),
                 self.log_date_time_string(),
                 format % args))

    def do_POST(self):
        #print('CHTTPRequestHandler.do_POST', self.path)
        parsedRequest = urlparse(self.path)
        print('parsedRequest', parsedRequest)
        pathElements = parsedRequest.path[1:].split("/")
        print('pathElements', pathElements)

        if not (parsedRequest.path.startswith("/api") or parsedRequest.path.startswith("/database")):
            self.send_response(404)

        elif pathElements[-1] in ["uploadFileToJob", "uploadFileForJobObject",
                                  "cloneJob", "createJob", "runJob", 
                                  "setJobParameterByXML"]:
            apiEndpoint = pathElements[-1]
            import cgi
            from core.CCP4Modules import HTTPSERVER
            ctype, pdict = cgi.parse_header(self.headers['Content-Type'])
            print(f"ctype {ctype} pdict {pdict}")
            if ctype == 'multipart/form-data':
                #print("Handling multipart/form-data")
                pdict['boundary'] = bytes(pdict['boundary'], 'utf-8')
                dbQueue = HTTPSERVER().dbThread.queue
                dbRequest = {'responseQueue': queue.Queue()}
                dbRequest['apiEndpoint'] = apiEndpoint
                dbRequest['path'] = self.path
                dbRequest['method'] = "POST"
                dbRequest['fields'] = cgi.parse_multipart(self.rfile, pdict)
                dbQueue.put(dbRequest)
                resultDict = dbRequest['responseQueue'].get(True, 10)
            elif ctype == "application/x-www-form-urlencoded":
                # Parse the form data
                form = cgi.FieldStorage(
                    fp=self.rfile,
                    headers=self.headers,
                    environ={'REQUEST_METHOD': 'POST',
                             'CONTENT_TYPE': self.headers['Content-Type'],
                             }
                    )
                fields = {}
                for key in form.keys(  ):
                    fields[key] = [form[key].value]
                dbQueue = HTTPSERVER().dbThread.queue
                dbRequest = {'responseQueue': queue.Queue()}
                dbRequest['apiEndpoint'] = apiEndpoint
                dbRequest['path'] = self.path
                dbRequest['method'] = "POST"
                dbRequest['fields'] = fields
                dbQueue.put(dbRequest)
                resultDict = dbRequest['responseQueue'].get(True, 10)
        else:
            self.send_response(404)

        self.returnData(contentType="application/json",
                        data=json.dumps(resultDict), isBinary=False)

    # Here I am going to do some hackery to allow the HTTP server to return information about the
    # database
    def do_GET(self):
        #print('CHTTPRequestHandler.do_GET', self.path)
        path = self.path
        if("Referer" in self.headers):
            #This deals with the css/images in the embedded report from pairef.
            try:
                b = self.headers["Referer"].split('?')[1]
                a = self.path
                if a.startswith("/database/") and not "getProjectJobFile?projectId=" in a:
                    old = b[b.find("fileName=")+len("fileName="):].split("&")[0].split("/")[-1]
                    new = a[10:].split("?")[0]
                    newpath = "/database/getProjectJobFile?"+b.replace(old,new)
                    path = newpath
            except:
                pass

        parsedRequest = urlparse(path)
        #print('parsedRequest', parsedRequest)
        pathElements = parsedRequest.path[1:].split("/")
        #print('pathElements', pathElements)
        if pathElements[0] in ["api", "database"] and \
                pathElements[-1] in ["ModelValues", "getProjectJobFile", "getProjectJobFileName", "getJobFile",
                                     "getFileWithPredicate", "makeTerminateFile", "getProjectFileData", "getReportXML"]:
            apiEndpoint = pathElements[-1]
            #print('This is a dbRequest API call')
            from core.CCP4Modules import HTTPSERVER
            dbQueue = HTTPSERVER().dbThread.queue
            dbRequest = {'responseQueue': queue.Queue()}
            dbRequest['apiEndpoint'] = apiEndpoint
            dbRequest['query'] = parse_qs(parsedRequest.query)
            dbRequest['method'] = "GET"
            dbQueue.put(dbRequest)
            resultDict = dbRequest['responseQueue'].get(True, 10)

            if "status" in resultDict and resultDict["status"] == "Exception":
                dbRequest['responseQueue'].task_done()
                self.send_response(404)
                return

            if apiEndpoint == "ModelValues":
                self.returnData(contentType="application/json",
                                data=json.dumps(resultDict), isBinary=False)
            elif apiEndpoint in ["getProjectJobFile", "getJobFile"]:
                self.returnData(**resultDict)
            elif apiEndpoint == "getProjectJobFileName":
                self.returnFileAsData(**resultDict)
            elif apiEndpoint in ["getFileWithPredicate", "makeTerminateFile", "getProjectFileData", "getReportXML"]:
                self.returnData(**resultDict)
            else:
                print("Unknown request", apiEndpoint)

            dbRequest['responseQueue'].task_done()

        else:
            return SimpleHTTPRequestHandler.do_GET(self)

    def returnFileAsDownload(self, contentType='application/xhtml+xml', fullPath=None):
        if not os.path.isfile(fullPath):
            self.send_response(404)
            return
        with open(fullPath, 'rb') as content_file:
            content = content_file.read()
            self.send_response(200)
            self.send_header('Content-type', contentType)
            self.send_header('Content-disposition',
                             'attachment;filename='+os.path.split(fullPath)[1])
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

    def returnData(self, contentType=None, data=None, isBinary=False, status="Success"):
        if contentType is not None and data is not None:
            self.send_response(200)
            self.send_header('Content-type', contentType)
            self.send_header('Content-length', len(data))
            self.end_headers()
            if not isBinary and not type(data) == bytes:
                self.wfile.write(bytes(data, "utf-8"))
            else:
                self.wfile.write(data)
