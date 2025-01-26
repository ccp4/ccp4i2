import argparse
import functools
import os
import sys

from PySide2 import QtCore

from ..core import CCP4Config
from ..core import CCP4File
from ..core import CCP4Modules
from ..core import CCP4PluginScript
from ..core.CCP4Utils import getCCP4I2Dir


@QtCore.Slot('CRunPlugin')
def quitThread(thread):
    print('quitThread',thread); sys.stdout.flush()
    CCP4Modules.QTAPPLICATION(graphical=False).quit()
    sys.exit()


def main():
    sys.stderr = sys.stdout     # Redirect stderr to stdout
    parser = argparse.ArgumentParser()
    parser.add_argument("xmlIn", help="Default Input File")
    parser.add_argument("-dbxml", "--dbXmlFile", help="Read in a i2 database xml file")
    parser.add_argument("-c", "--configFile", help="Read in an i2 Configuration file")
    parser.add_argument("-db", "--dbFile", help="Read in an i2 database file")
    args = parser.parse_args()
    top_path = getCCP4I2Dir()
    print('runTask Script raw Input Arguments', sys.argv)
    print("runTask Script top_path(CCPI2dir) and __file__ variables", top_path, __file__)
    graphical = False
    # Use the specified config file or dbFile
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
        comFilePath = sXmlIn
        compressedFile = None
        xmlHeader = CCP4File.CI2XmlHeader()
        xmlHeader.loadFromXml(sXmlIn, checkValidity=False)
        #print 'runTask pluginName', xmlHeader, xmlHeader.pluginName
    else:
        comFilePath = None
        compressedFile = sXmlIn
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
