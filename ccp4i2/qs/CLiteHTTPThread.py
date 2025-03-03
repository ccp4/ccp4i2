from __future__ import print_function

import sys
import os
import signal
import http.server
from http.server import SimpleHTTPRequestHandler
from PyQt4 import QtCore
import queue

from socketserver import ThreadingMixIn
class MultiThreadedHTTPServer(ThreadingMixIn, http.server.HTTPServer):
    #Let there be peace !
    def log_message(*argv, **kw):
        pass

import CLiteDbThread

class CLiteHTTPThread(QtCore.QThread):
    insts= None
    def __init__(self,*args, **kw):
        #Queue to which to send signal when thread is running
        self.hasStartedQueue = kw.get('hasStartedQueue')
        del kw['hasStartedQueue']
        #Queue to which to send requests for db Informatin
        self.dbQueue = None
        self.dbQueue = kw.get('dbQueue')
        del kw['dbQueue']
        #Directory to form principal route for file searches
        self.parentDir = kw.get('parentDir')
        del kw['parentDir']
        #Default port on which to serve
        self.defaultPort = 43434
        if 'defaultPort' in kw:
            self.defaultPort = kw.get('defaultPort')
            del kw['defaultPort']
        super(CLiteHTTPThread,self).__init__(*args, **kw)
        self.setObjectName('HTTPServer')
        self.port = None
        CLiteHTTPThread.insts = self

    def shutdown(self):
        self.httpd.shutdown()

    def run(self):
        self.hasStartedQueue.put("HTTP thread successfully started")
        # Code from http://www.linuxjournal.com/content/tech-tip-really-simple-http-server-python
        # Using this short script so that it will only serve localhost
        os.chdir(self.parentDir)
        self.httpd =  self.makeServer()
        if self.httpd is not None:
            sa = self.httpd.socket.getsockname()
            self.port = sa[1]
            self.httpd.serve_forever()
        return

    def makeServer(self):
        HandlerClass = CLiteHTTPRequestHandler
        ServerClass  = MultiThreadedHTTPServer
        Protocol     = "HTTP/1.0"
        server_address = ('127.0.0.1', self.defaultPort )
        #print server_address
        HandlerClass.protocol_version = Protocol
        try:
            httpd = ServerClass(server_address, HandlerClass)
            return httpd
        except:
            server_address = ('127.0.0.1', 0 )
            #print server_address
            httpd = ServerClass(server_address, HandlerClass)
            return httpd

class CLiteHTTPRequestHandler(SimpleHTTPRequestHandler):
    fileTypesMap= {"Unknown":"DataFile","application/CCP4-seq":"SeqDataFile","chemical/x-pdb":"PdbDataFile","MultiPDB":"","application/CCP4-mtz":"MtzDataFile","application/CCP4-unmerged-mtz":"MtzDataFile","application/CCP4-unmerged-experimental":"UnmergedDataFile","application/CCP4-map":"MapDataFile","application/refmac-dictionary":"DictDataFile","application/refmac-TLS":"TLSDataFile","application/CCP4-mtz-freerflag":"FreeRDataFile","application/CCP4-mtz-observed":"ObsDataFile","application/CCP4-mtz-phases":"PhsDataFile","application/CCP4-mtz-map":"MapCoeffsDataFile","":"","application/CCP4-seqalign":"SeqAlignDataFile","application/CCP4-mtz-mini":"MiniMtzDataFile","application/coot-script":"CootHistoryDataFile","application/refmac-external-restraints":"RefmacRestraintsDataFile","application/CCP4-scene":"SceneDataFile","application/CCP4-shelx-FA":"ShelxFADataFile","application/phaser-sol":"PhaserSolDataFile","chemical/x-mdl-molfile":"MDLMolDataFile","application/iMosflm-xml":"ImosflmXmlDataFile","application/CCP4-image":"ImageFile","application/CCP4-generic-reflections":"GenericReflDataFile","application/HHPred-alignments":"HhpredDataFile","application/Blast-alignments":"BlastDataFile","chemical/x-pdb-ensemble":"EnsemblePdbDataFile","application/CCP4-asu-content":"AsuDataFile"}
    
    #Let there be peace
    def log_message(*args, **kw):
        pass
    def translate_path(self,path):
        #print "Translating path ", path
        if path.startswith("/jobId/"):
            pathElements = path.split("/")
            self.responseQueue = queue.Queue()
            callback = {'responseQueue':self.responseQueue, 'path':'getJobDirectory?jobId={}'.format(pathElements[2])}
            CLiteDbThread.CLiteDbThread.insts.queue.put(callback)
            result = self.responseQueue.get(True)
            rerooted = os.path.join(result,os.sep.join(pathElements[3:]))
            return rerooted
        elif path.startswith("/icon/"):
            pathElements = path.split("/")
            mimeName = "/".join(pathElements[-2:])
            rerooted = path
            if mimeName in CLiteHTTPRequestHandler.fileTypesMap:
                oneDirUp = os.path.split(os.path.abspath(__file__))[0]
                twoDirsUp = os.path.split(oneDirUp)[0]
                className = CLiteHTTPRequestHandler.fileTypesMap[mimeName]
                rerooted = os.path.join(twoDirsUp,"qticons","{}.svg".format(className))
                if os.path.isfile(rerooted): return rerooted
                rerooted = os.path.join(twoDirsUp,"qticons","{}.png".format(className))
                if os.path.isfile(rerooted): return rerooted
            return rerooted
        else:
            result = SimpleHTTPRequestHandler.translate_path(self, path)
        return result

def sigint_handler(*args):
    QtCore.QCoreApplication.quit()

if __name__ == "__main__":
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    application = QtCore.QCoreApplication(sys.argv)
    serverThread = CLiteHTTPThread(parent=application, defaultPort=43434, parentDir=os.path.split(os.path.abspath(__file__))[0])
    serverThread.start()
    sys.exit(application.exec_())
