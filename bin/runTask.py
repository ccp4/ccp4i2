from __future__ import print_function


import os
import sys
import argparse
import functools
from PySide2 import QtCore

def getCCP4I2Dir(up=1):
    target = os.path.join(os.path.realpath(sys.argv[0]),"..")
    abstarget = os.path.abspath(target)
    splittarget = abstarget.split()
    if splittarget.count('ccp4i2'):
        splittarget.reverse()
        up = splittarget.index('ccp4i2')
    while up>0:
        abstarget = os.path.dirname(abstarget)
        up = up -1
    return abstarget

@QtCore.Slot('CRunPlugin')
def quitThread(thread):
    print('quitThread',thread); sys.stdout.flush()
    from core import CCP4Modules
    CCP4Modules.QTAPPLICATION(graphical=False).quit()
    sys.exit()

#--------------------------------------------------------------------------

if __name__ == '__main__':
    sys.stderr = sys.stdout     # Redirect stderr to stdout
    parser = argparse.ArgumentParser()     # KJS : Switch to argparse.
    parser.add_argument("xmlIn", help="Default Input File")
    parser.add_argument("-dbxml", "--dbXmlFile", help="Read in a i2 database xml file")
    parser.add_argument("-c", "--configFile", help="Read in an i2 Configuration file")
    parser.add_argument("-db", "--dbFile", help="Read in an i2 database file")
    args = parser.parse_args()
    top_path = getCCP4I2Dir()
    print('runTask Script raw Input Arguments', sys.argv)
    print("runTask Script top_path(CCPI2dir) and __file__ variables", top_path, __file__)
    exec(compile(open(os.path.join(top_path, 'utils', 'startup.py')).read(), os.path.join(top_path, 'utils', 'startup.py'), 'exec'))  # This seems to be a hack to expose setupPythonpath() ?
    graphical = False
    if graphical:
        setupPythonpath(mode='qtgui')
    else:
        setupPythonpath(mode='qtcore')
    # Use the specified config file or dbFile
    from core import CCP4Config
    if args.configFile is not None:
        config = CCP4Config.CONFIG(args.configFile)
        print('Running plugin using config file:', args.configFile)
    else:
        print("what is going on ? loadConfig")
        config = loadConfig()
    config.set('graphical', graphical)
    if args.dbFile is not None and os.path.exists(args.dbFile):
        print('Running plugin using database file:', args.dbFile)
        config.set('dbFile', args.dbFile)
    # Get name of plugin and whether it is asynchronous
    sXmlIn = str(args.xmlIn)
    if os.path.splitext(sXmlIn)[1] == '.xml':
        from core import CCP4File
        comFilePath = sXmlIn
        compressedFile = None
        xmlHeader = CCP4File.CI2XmlHeader()
        xmlHeader.loadFromXml(sXmlIn, checkValidity=False)
        #print 'runTask pluginName', xmlHeader, xmlHeader.pluginName
    else:
        comFilePath = None
        compressedFile = sXmlIn
    from core import CCP4PluginScript
    from core import CCP4Modules
    if 0:
        print('Run non-asynchronous task',xmlHeader.pluginName)
        runPlugin = CCP4PluginScript.CRunPlugin(None, top_path, comFilePath=comFilePath, compressedFile=compressedFile, dbXmlFile=args.dbXmlFile)
        runPlugin.run()
    else:
        #print 'Run ASYNCHRONOUS task',xmlHeader.pluginName
        print('Running QApplication to support asyncronous sub-processes')
        app = CCP4Modules.QTAPPLICATION(graphical=graphical)
        runPlugin = CCP4PluginScript.CRunPlugin(app, top_path, comFilePath=comFilePath, compressedFile=compressedFile, dbXmlFile=args.dbXmlFile)
        app.aboutToQuit.connect(functools.partial(quitThread, runPlugin))
        runPlugin.finished.connect(functools.partial(quitThread, runPlugin))
        runPlugin.run()
        sys.exit(app.exec_())
