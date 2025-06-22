import os
import socketserver
import threading
import time
import xml.etree.ElementTree as ET

from ....core import CCP4Utils
from ....core.CCP4PluginScript import CPluginScript


class MosflmRequestHandler(socketserver.StreamRequestHandler):
    commandLines = []
    xmlRoot = None
    xmlFilepath = None
    
    def handle(self):
        responseText = "Squazgar"
        currentXmlLength = ''
        
        while responseText != '<done>':
            responseText = self.rfile.readline().strip()
        
        #Enter into command issuing dialog
        for commandLine in MosflmRequestHandler.commandLines:
            if not ( commandLine.lower().startswith('hklo') or commandLine.lower().startswith('exit') or commandLine.lower().startswith('mtzd')):
                self.wfile.write(commandLine.strip()+'\n')
                self.wfile.flush()
                responseText = self.rfile.readline().strip()
        
        #print '\n\nInto listen loop\n\n'
        inBlockIntegrate = False
        while responseText != "<done>":
            #Block until readable
            responseText = self.rfile.readline().strip()
            if responseText != "<done>":
                try:
                    aDom = ET.fromstring(responseText)
                    MosflmRequestHandler.xmlRoot.append(aDom)
                    if aDom.tag == 'integration_response':
                        inBlockIntegrate = True
                    elif inBlockIntegrate:
                        inBlockIntegrate = False
                        newXmlLength = len(ET.tostring(MosflmRequestHandler.xmlRoot))
                        if newXmlLength > currentXmlLength:
                            currentXmlLength = newXmlLength
                            CCP4Utils.writeXml(MosflmRequestHandler.xmlRoot, MosflmRequestHandler.xmlFilepath+'.tmp')
                            os.rename(MosflmRequestHandler.xmlFilepath+'.tmp',MosflmRequestHandler.xmlFilepath)
                except:
                    print('Couldnt interpret ['+responseText+']')
            self.wfile.write('continue\n')
            self.wfile.flush()
        CCP4Utils.writeXml(MosflmRequestHandler.xmlRoot, MosflmRequestHandler.xmlFilepath+'.tmp')
        os.rename(MosflmRequestHandler.xmlFilepath+'.tmp',MosflmRequestHandler.xmlFilepath)

class MosflmServer(socketserver.ThreadingMixIn, socketserver.TCPServer):
    allow_reuse_address = True
    daemon_threads = True

class mosflm(CPluginScript):

    TASKTITLE = 'Run mosflm and capture output'     # A short title for gui menu
    TASKNAME = 'mosflm'                                  # Task name - should be same as class name
    TASKCOMMAND = 'ipmosflm'                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template
    MAINTAINER = 'martin.noble@newcastle.ac.uk'

    '''
    def __init__(self,parent=None,name=None,workDirectory=''):
      CPluginScript. __init__(self,parent=parent,name=name)
    '''
    
    def makeCommandAndScript(self):
        #self.simpleServer = MyServer()
        #self.simpleServer.start()
        
        address = ('127.0.0.1', 0)

        MosflmRequestHandler.xmlRoot = ET.Element('MosflmXML')
        MosflmRequestHandler.xmlFilepath = self.makeFileName( 'PROGRAMXML' )
        MosflmRequestHandler.commandLines = self.container.controlParameters.SCRIPT.split('\n')

        myMosflmServer = MosflmServer(address, MosflmRequestHandler)
        t = threading.Thread(target=myMosflmServer.serve_forever)
        t.setDaemon(True) # don't hang on exit
        t.start()
        print('Server loop running in thread:', t.getName())

        time.sleep(0.1)
        address,port = myMosflmServer.server_address
        self.appendCommandLine(['MOSFLMSOCKET', str(port)])
        self.appendCommandLine(['HKLOUT',self.container.outputData.UNMERGEDMTZ.fullPath])
        return CPluginScript.SUCCEEDED
