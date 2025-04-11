import argparse
import functools
import os
import sys

from PySide2 import QtCore

from .. import I2_TOP
from ..core.CCP4Config import CONFIG
from ..utils.QApp import QTAPPLICATION


@QtCore.Slot("CRunPlugin")
def quitThread(thread):
    print("quitThread", thread, flush=True)
    QTAPPLICATION(graphical=False).quit()
    sys.exit()


def main():
    sys.stderr = sys.stdout  # Redirect stderr to stdout
    parser = argparse.ArgumentParser()
    parser.add_argument("xml", help="Default Input File")
    parser.add_argument("-db", help="Read in an i2 database file")
    parser.add_argument("-dbxml", help="Read in a i2 database xml file")
    args = parser.parse_args()
    print(f"runTask Script {sys.argv = }")
    print(f"runTask Script {I2_TOP = }")
    print(f"runTask Script {__file__ = }")
    config = CONFIG(graphical=False)
    if args.db is not None and os.path.exists(args.db):
        print("Running plugin using database file:", args.db)
        config.dbFile = args.db
    if args.xml.endswith(".xml"):
        comFilePath = args.xml
        compressedFile = None
    else:
        comFilePath = None
        compressedFile = args.xml
    print("Running QApplication to support asyncronous sub-processes")
    app = QTAPPLICATION(graphical=False)
    from ..core.CCP4PluginScript import CRunPlugin
    runPlugin = CRunPlugin(
        parent=app,
        ccp4i2Path=str(I2_TOP),
        comFilePath=comFilePath,
        compressedFile=compressedFile,
        dbXmlFile=args.dbxml,
    )
    app.aboutToQuit.connect(functools.partial(quitThread, runPlugin))
    runPlugin.finished.connect(functools.partial(quitThread, runPlugin))
    runPlugin.run()
    sys.exit(app.exec_())
