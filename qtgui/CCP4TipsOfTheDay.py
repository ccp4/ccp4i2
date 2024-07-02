import os
import random
import glob

from PySide6 import QtGui, QtWidgets, QtCore, QtWebEngine, QtWebEngineWidgets

from core import CCP4Utils, CCP4Modules

class CTipsOfTheDay(QtWidgets.QDialog):

    def __init__(self,parent=None):
        QtWidgets.QDialog.__init__(self,parent)
        self.setWindowTitle("Tip of the Day")
        layout = QtWidgets.QVBoxLayout()
        headLayout = QtWidgets.QHBoxLayout()
        self.setLayout(layout)
        lightBulb = QtWidgets.QLabel()
        lightBulbPixmap = QtGui.QPixmap(os.path.join(CCP4Utils.getCCP4I2Dir(), 'qticons', "lightbulb.svg"))
        lightBulb.setPixmap(lightBulbPixmap.scaledToHeight(50))
        headLayout.addWidget(lightBulb)
        didYouKnow = QtWidgets.QLabel("Did you know ...?")
        didYouKnow.setStyleSheet("QLabel { font-size: 32pt;}")
        headLayout.addWidget(didYouKnow)
        layout.addLayout(headLayout)
        tipWindow = QtWebEngineWidgets.QWebEngineView()
        tipWindow.setMaximumWidth(530)
        layout.addWidget(tipWindow)
        tipWindow.setContextMenuPolicy(QtCore.Qt.NoContextMenu) 

        showCheck = QtWidgets.QCheckBox("Show Tips on Startup")
        layout.addWidget(showCheck)
        showTipsOfTheDay = CCP4Modules.PREFERENCES().SHOW_TIPS_OF_THE_DAY
        if showTipsOfTheDay:
            showCheck.setChecked(True)
        @QtCore.Slot(bool)
        def showTipsOfTheDayHandler(v):
            showTipsOfTheDay.set(v)
        showCheck.clicked.connect(showTipsOfTheDayHandler)

        dialogButtonBox = QtWidgets.QDialogButtonBox()
        layout.addWidget(dialogButtonBox)
        cb = dialogButtonBox.addButton(QtWidgets.QDialogButtonBox.Close)
        pb = dialogButtonBox.addButton("Previous tip",QtWidgets.QDialogButtonBox.ActionRole)
        nb = dialogButtonBox.addButton("Next tip",QtWidgets.QDialogButtonBox.ActionRole)
        nb.setDefault(True)

        def numList(elem):
            return int(os.path.basename(elem).split(".")[0])
        numTips = int(os.path.basename(sorted(glob.glob(os.path.join(CCP4Utils.getCCP4I2Dir(),'tipsOfTheDay','*.html')),key=numList)[-1]).split(".")[0])

        #In Python2 we need to deal with a list so that we can edit the variable in an inner-scope. (Can use nonlocal in Python3.)
        day = [int(random.random()*(numTips-1))]

        @QtCore.Slot()
        def setNextTipText():
           day[0] += 1
           if day[0] > numTips:
               day[0] = 1
           dayTextFile = os.path.join(CCP4Utils.getCCP4I2Dir(), 'tipsOfTheDay', str(day[0])+".html")
           self.setWindowTitle(str(day[0]))
           tipWindow.load(QtCore.QUrl.fromLocalFile(dayTextFile))

        @QtCore.Slot()
        def setPreviousTipText():
           day[0] -= 1
           if day[0] == 0:
               day[0] = numTips
           dayTextFile = os.path.join(CCP4Utils.getCCP4I2Dir(), 'tipsOfTheDay', str(day[0])+".html")
           self.setWindowTitle(str(day[0]))
           tipWindow.load(QtCore.QUrl.fromLocalFile(dayTextFile))

        setNextTipText()
 
        cb.clicked.connect(self.accept)
        nb.clicked.connect(setNextTipText)
        pb.clicked.connect(setPreviousTipText)
        
        @QtCore.Slot()
        def zoomIn():
            tipWindow.setZoomFactor(tipWindow.zoomFactor()*1.2)

        @QtCore.Slot()
        def zoomOut():
            tipWindow.setZoomFactor(tipWindow.zoomFactor()/1.2)

        @QtCore.Slot()
        def resetZoom():
            tipWindow.setZoomFactor(1.0)

        zoomOutShortcut = QtGui.QShortcut(QtGui.QKeySequence.ZoomOut,self)
        zoomInShortcut = QtGui.QShortcut(QtGui.QKeySequence.ZoomIn,self)
        zoomInShortcut2 = QtGui.QShortcut(QtGui.QKeySequence(self.tr("Ctrl+=", "Zoom in")),self)
        resetZoomShortcut = QtGui.QShortcut(QtGui.QKeySequence(self.tr("Ctrl+0", "Reset Zoom")),self)

        zoomInShortcut.activated.connect(zoomIn)
        zoomInShortcut2.activated.connect(zoomIn)
        zoomOutShortcut.activated.connect(zoomOut)
        resetZoomShortcut.activated.connect(resetZoom)
