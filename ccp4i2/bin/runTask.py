import argparse
import functools
import os
import sys

from PySide2 import QtCore

from ..core.CCP4Config import CONFIG
from ..core.CCP4Utils import getCCP4I2Dir
from ..utils.QApp import QTAPPLICATION


@QtCore.Slot('CRunPlugin')
def quitThread(thread):
    print('quitThread',thread); sys.stdout.flush()
    QTAPPLICATION(graphical=False).quit()
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
    # Use the specified config file or dbFile
    if args.configFile is not None:
        print('Running plugin using config file:', args.configFile)
    config = CONFIG(args.configFile, graphical=False)
    if args.dbFile is not None and os.path.exists(args.dbFile):
        print('Running plugin using database file:', args.dbFile)
        config.dbFile = args.dbFile
    # Get name of plugin and whether it is asynchronous
    sXmlIn = str(args.xmlIn)
    if os.path.splitext(sXmlIn)[1] == '.xml':
        from ..core import CCP4File
        comFilePath = sXmlIn
        compressedFile = None
        xmlHeader = CCP4File.CI2XmlHeader()
        xmlHeader.loadFromXml(sXmlIn, checkValidity=False)
        #print 'runTask pluginName', xmlHeader, xmlHeader.pluginName
    else:
        comFilePath = None
        compressedFile = sXmlIn
    from ..core import CCP4PluginScript
    if 0:
        print('Run non-asynchronous task',xmlHeader.pluginName)
        runPlugin = CCP4PluginScript.CRunPlugin(None, top_path, comFilePath=comFilePath, compressedFile=compressedFile, dbXmlFile=args.dbXmlFile)
        runPlugin.run()
    else:
        #print 'Run ASYNCHRONOUS task',xmlHeader.pluginName
        print('Running QApplication to support asyncronous sub-processes')
        app = QTAPPLICATION(graphical=False)
        runPlugin = CCP4PluginScript.CRunPlugin(app, top_path, comFilePath=comFilePath, compressedFile=compressedFile, dbXmlFile=args.dbXmlFile)
        app.aboutToQuit.connect(functools.partial(quitThread, runPlugin))
        runPlugin.finished.connect(functools.partial(quitThread, runPlugin))
        runPlugin.run()
        sys.exit(app.exec_())
