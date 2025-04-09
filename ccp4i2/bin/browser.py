import atexit
import cProfile
import faulthandler
import io
import os
import pstats
import sys
import time

from PySide2 import QtCore, QtGui, QtWidgets, QtWebEngineWidgets
import ccp4mg # Sets sys.path so import of MG modules will work from here onwards

from .. import I2_TOP
from ..utils.QApp import QTAPPLICATION
from ..utils.startup import startBrowser


def main():
    # faulthandler.enable()  # Help with debugging segfaults
    print("CCP4", os.environ["CCP4"])
    print('Running CCP4i2 browser from:', I2_TOP)
    print('Python', sys.version)
    try:
        print('Qt version', QtCore.qVersion())
    except:
        print('Failed finding Qt verion')
    print(' ')
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_ShareOpenGLContexts)
    app = QTAPPLICATION(graphical=True)
    QtWebEngineWidgets.QWebEngineProfile.defaultProfile().clearHttpCache()
    splash = None
    # This splash screen is not used properly as it does not wait for any window.
    # However it does seem to set whatever is necessary to get the app to switch to the correct menu bar on OS X.
    pixmap = QtGui.QPixmap(str(I2_TOP / "qticons" / "ccp4i2.png"))
    pix2 = QtGui.QPixmap(pixmap)
    pix2.fill(QtCore.Qt.white)
    painter = QtGui.QPainter()
    painter.begin(pix2)
    painter.drawPixmap(0, 0, pixmap)
    painter.end()
    splash = QtWidgets.QSplashScreen(pix2)
    version_string = '1.0'
    splash.showMessage(splash.tr("Starting ") + version_string + "...   ", QtCore.Qt.AlignRight | QtCore.Qt.AlignBottom, QtCore.Qt.black)
    splash.show()
    time.sleep(0.05)  # KJS : This+repaint fixes the Qt problem with the Splash Screen on Linux.
    splash.repaint()
    splash.raise_()
    splash.update()
    app.processEvents()
    if 'CCP4' not in os.environ:
        box = QtWidgets.QMessageBox()
        box.setText('CCP4 not setup - programs will not run.\nPlease exit and setup CCP4 before restarting.')
        box.setWindowTitle('No CCP4 setup')
        box.addButton('Continue', QtWidgets.QMessageBox.YesRole)
        box.addButton('Exit', QtWidgets.QMessageBox.NoRole)
        box.show()
        box.raise_()
        reply = box.exec_()
        if reply == 1:
            sys.exit()
    if "--profile" in sys.argv:
        pr = cProfile.Profile()
        pr.enable()
        while "--profile" in sys.argv:
            sys.argv.remove("--profile")
        print("Profiling")
        startBrowser(sys.argv[1:], app=app, splash=splash)
        #-----------------------------------------------------------------------
        def printStats():
            pr.disable()
            s = io.StringIO()
            sortby = 'cumulative'
            ps = pstats.Stats(pr, stream=s).sort_stats(sortby).reverse_order()
            ps.print_stats()
            print(s.getvalue())
        #-----------------------------------------------------------------------
        atexit.register(printStats)
        sys.exit(app.exec_())
    else:
        startBrowser(sys.argv[1:], app=app, splash=splash)
        sys.exit(app.exec_())
