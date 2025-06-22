import argparse
import functools
import glob
import os
import sys
import time
import xml.etree.ElementTree as ET

from PySide2 import QtCore

from ..core import CCP4Utils
from ..core.CCP4Config import CONFIG
from ..utils import startup
from ..utils.QApp import QTAPPLICATION


def quitThread(thread):
    print('quitThread',thread, flush=True)
    QTAPPLICATION(graphical=False).quit()
    sys.exit()


def onTestRunnerComplete(logXmlPath, logXmlRoot, app):
    threads = app.findChildren(QtCore.QThread)
    for t in threads:
        if hasattr(t,"quitServer"):
            t.quitServer()
        t.wait()
    if logXmlRoot is not None:
        CCP4Utils.writeXml(logXmlRoot, logXmlPath)
    sys.exit()


def main():
    # Redirect stderr to stdout
    sys.stderr = sys.stdout
    # print 'testi2sys sys.argv',sys.argv
    parser = argparse.ArgumentParser(description='Run CCP4i2 test project')
    parser.add_argument('-x', '--xmlOut',action='store_true', default=False)
    parser.add_argument('-c', '--configFile',help='Path to config file', default=None)
    parser.add_argument('-o', '--outputDirectory',help='Path to directory where job test will be unpacked and run', default=None)
    parser.add_argument('-p', '--project',help='Path to compressed project archive')
    parser.add_argument('-j','--jobs', help = "Run only seleced jobs, input in form '1,2,7-8' with no spaces", default=None)
    #Following arguments are ones that might be passed on to the CProjectBasedTesting instance
    parser.add_argument('-r', '--runMode', help = "Runmode - 0 : 'no testing' , 1 : 'run jobs', 2: 'run tests', 3: 'run jobs and tests'", type=int, default=3, choices=[0,1,2,3])
    parser.add_argument('-v', '--verbosity', help = "Verbosity - 0 : 'fails only' , 1 : 'jobs and parameters tested', 3 : 'jobs, parameters and parameter values'", type=int, default=1, choices=[0,1,3])
    parser.add_argument('-t', '--testSubJobs', help = "Tset sub jobs - 0 : 'no' , 1 : 'for failed jobs', 2: 'for all jobs'", type=int, default=0,choices=[0,1,2])
    parser.add_argument('-f', '--copyFiles', help = "Copy files 0 : 'no' , 1 : 'yes'", type=int, default=0, choices=[0,1])
    parser.add_argument('-i', '--copyCCP4ImportedFiles', help = "Copy CCP4_IMPORTED_FILES files 0 : 'no' , 1 : 'yes'", type=int, default=0, choices=[0,1])
    parser.add_argument('-l', '--resetBaseline', help = "Reset baseline - 0 : 'no' , 1 : 'yes", type=int, default=0, choices=[0,1])
    pns = parser.parse_args()

    configFile = pns.configFile
    compressedProjectFile = pns.project
    outputDirectory = pns.outputDirectory
    selectedJobs = pns.jobs

    from ..core import CCP4ProjectBasedTesting
    kw = {}
    for name in CCP4ProjectBasedTesting.MODES:
        kw[name] = getattr(pns, name)

    # truncating sys.argv to stop the annoying error messages about 'usage'
    # but need to leave something for the CBazaar.initialise to pop()
    del sys.argv[2:]

    source = CCP4Utils.getCCP4I2Dir()
    t = time.localtime()
    startTime = time.strftime('%y-%m-%d-%H-%M',t)
    startTime0 = time.strftime('%H-%M %d %b %y',t)

    if compressedProjectFile is None:
        print("ERROR: no compressed project file")
        sys.exit()
    compressedProjectFile = os.path.abspath(compressedProjectFile)
    if not os.path.exists(compressedProjectFile):
        print("ERROR: compressed project file not found: "+compressedProjectFile)
        sys.exit()
    elif os.path.isdir(compressedProjectFile):
        from ..qtcore import CCP4Export
        compressedProjectFileList = glob.glob(os.path.join(compressedProjectFile,"*."+CCP4Export.COMPRESSED_SUFFIX))
        if len(compressedProjectFileList)<1:
            print("ERROR: no compressed project files found in directory: "+compressedProjectFile)
            sys.exit()
    elif os.path.isfile(compressedProjectFile):
        compressedProjectFileList = [compressedProjectFile]
    else:
        print("ERROR: compressed project file recognised as file or directory: "+compressedProjectFile)
        sys.exit()

    # Use the specified config file or dbFile
    if configFile is not None:
        print('Running tests using config file:', configFile)
    config = CONFIG(configFile, graphical=False)

    if outputDirectory is None:
        outputDirectory = os.getcwd()
    else:
        outputDirectory = os.path.abspath(outputDirectory)
    
    if not os.path.exists(outputDirectory):
        try:
            os.makedirs(outputDirectory)
        except:
            print("ERROR: failed to create test directory:",outputDirectory)
            sys.exit()
        else:
            print('Directory for output files: '+outputDirectory)

    dbFile = os.path.join(outputDirectory,'db-'+startTime+'.sqlite')
    if os.path.exists(dbFile):
        print("ERROR: database file exists:",dbFile)
        sys.exit()

    config.dbFile = dbFile

    #sys.argv =  [sys.argv[0]]

    app = QTAPPLICATION(graphical=False)
    pm = startup.startProjectsManager(dbFileName=dbFile,loadDiagnostic=False)
    pm.startCheckForFinishedJobs()
    jc = startup.startJobController()
    jc.setDiagnostic(False)
    if dbFile is not None: jc.setDbFile(dbFile)
    pm.doCheckForFinishedJobs.connect(pm.checkForFinishedJobs)
    from ..core import CCP4DataManager
    CCP4DataManager.DATAMANAGER().buildClassLookup()

    print('Test results will be saved to: '+os.path.join(outputDirectory,'test-'+startTime+'.log')+' and database '+dbFile)
    log = CCP4ProjectBasedTesting.Logger(os.path.join(outputDirectory,'test-'+startTime+'.log'))
    log.write('Running CCP4i2 from: '+str(source)+'\n')
    log.write('Started: '+startTime0+'\n\n\n')
    if pns.xmlOut:
        logXmlPath = os.path.join(outputDirectory,'test-'+startTime+'.xml')
        logXmlRoot = ET.Element('ProjectTesting',startTime=startTime)
    else:
        logXmlPath = None
        logXmlRoot = None
    testRunner = CCP4ProjectBasedTesting.CProjectBasedTesting(sourceProjectList=compressedProjectFileList,outputDirectory=outputDirectory,parent=app,log=log,logXmlRoot=logXmlRoot,selectedJobs=selectedJobs,useCurrentDb=True,**kw)
    testRunner.finished.connect(functools.partial(onTestRunnerComplete, logXmlPath, logXmlRoot,app))
    testRunner.runTests()
    #print 'testi2sys from testRunner.runTests()'
    sys.exit(app.exec_())
