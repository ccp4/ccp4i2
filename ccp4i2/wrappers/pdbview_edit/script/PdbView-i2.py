import functools
import glob
import os
import sys

from ccp4i2.baselayer import QtCore, QtWidgets
from ccp4i2.core.mgimports import PdbView

#workDirectory = "XXXXX_WORK_DIR_XXXXX"

class PdbViewMainWindowI2(PdbView.PDBVIEWMainWindow):

  @QtCore.Slot(str)
  def WriteFileToI2(self,workDirectory):
      print("Triggered!!")
      dropDir = os.path.join(workDirectory,"CCP4MG_FILE_DROP")
      outList = glob.glob(os.path.join(dropDir,'output*.pdb'))

      print("dropDir",dropDir)
      print("outList",outList)

      maxIndx = 0
      for f in outList:
        fpath,fname = os.path.split(f)
        maxIndx =  max(maxIndx,int(fname[6:-4]))

      print("maxIndx",maxIndx)

      fname = os.path.join(dropDir,'output'+str(maxIndx+1)+'.pdb')
      print("Save",fname)

      self.save(fname)

  def __init__(self,parent=None,widget=None,workDir=None):
    PdbView.PDBVIEWMainWindow.__init__(self,parent,widget)
    menus = self.menuBar().findChildren(QtWidgets.QMenu);
    for m in menus:
      if m.title() == "File":
        act = m.addAction("Save all atoms to i2 database")
        act.setShortcut("Ctrl+2")
        act.triggered.connect(functools.partial(self.WriteFileToI2,workDir))
        print("Connected signal")

if __name__ == "__main__":
  app = QtWidgets.QApplication(sys.argv)

  widget = PdbView.PDBVIEWWidget()

  workDir = sys.argv[1]
  win = PdbViewMainWindowI2(None,widget,workDir=workDir)

  splash = QtWidgets.QSplashScreen()
  version_string = '1.0'
  splash.showMessage(splash.tr("Starting ")+version_string+"...   ",QtCore.Qt.AlignRight|QtCore.Qt.AlignBottom,QtCore.Qt.black)
  splash.show()
  splash.raise_()

  splash.finish(win)

  win.show()
  win.raise_()

  if len(sys.argv)>1:
    win.openFile(sys.argv[2])
    print("Opened",sys.argv[2])

  for arg in sys.argv[3:]:
    widget2 = PdbView.PDBVIEWWidget()
    win2 = PdbViewMainWindowI2(None,widget2,workDir=workDir)
    win2.openFile(arg)
    win2.show()
    win2.raise_()

  sys.exit(app.exec_())
