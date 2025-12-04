from __future__ import print_function

import sys

from baselayer import QtCore, QtGui, QtWidgets

import ccp4mg
import PdbView

if __name__ == "__main__":

  from sys import platform

  # Check if we're on OS X, first.
  if platform == 'darwin':
      try:
          from Foundation import NSBundle
          bundle = NSBundle.mainBundle()
          if bundle:
              info = bundle.localizedInfoDictionary() or bundle.infoDictionary()
              if info and info['CFBundleName'] == 'Python':
                  info['CFBundleName'] = "PdbView"
      except:
          # We cannot do this without the PyObjC libraries
          print("PyObjC libraries not found - cannot change application menubar name.")

  app = QtWidgets.QApplication(sys.argv)

  widget = PdbView.PDBVIEWWidget()

  win = PdbView.PDBVIEWMainWindow(None,widget)

  splash = QtWidgets.QSplashScreen()
  version_string = '1.0'
  splash.showMessage(splash.tr("Starting ")+version_string+"...   ",QtCore.Qt.AlignRight|QtCore.Qt.AlignBottom,QtCore.Qt.black)
  splash.show()
  splash.raise_()

  splash.finish(win)

  win.show()
  win.raise_()

  if len(sys.argv)>1:
    win.openFile(sys.argv[1])
    print("Opened",sys.argv[1])

  for arg in sys.argv[2:]:
    widget2 = PdbView.PDBVIEWWidget()
    win2 = PdbView.PDBVIEWMainWindow(None,widget2)
    win2.openFile(arg)
    win2.show()
    win2.raise_()

  sys.exit(app.exec_())
