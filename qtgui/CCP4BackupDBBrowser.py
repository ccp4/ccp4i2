from __future__ import print_function

import sys
import os
import sqlite3
import glob
from datetime import datetime
import functools

from PySide2 import QtCore, QtGui, QtWidgets
import sqlite3

def datetimesort(k1,k2):
    t1 = datetime.strptime(k1.lstrip("database_sqlite_backup-"),"%d%m%Y-%H%M%S")
    t2 = datetime.strptime(k2.lstrip("database_sqlite_backup-"),"%d%m%Y-%H%M%S")
    if t1 < t2:
        return -1
    if t1 > t2:
        return 1
    if t1 == t2:
        return 0
    return None

class SimpleTextDelegate(QtWidgets.QStyledItemDelegate):
    def __init__(self,parent=None):
        QtWidgets.QStyledItemDelegate.__init__(self,parent)

    def paint(self,painter,option,index):

        """
        Allows us to store the full file name as the data but display only the date part in the list.
        """

        painter.save();

        if option.state & QtWidgets.QStyle.State_Selected and option.state & QtWidgets.QStyle.State_Active:
            painter.setPen(option.palette.highlightedText().color());        
            painter.fillRect(option.rect, option.palette.highlight());

        style = QtWidgets.QApplication.style();
        textRect = style.subElementRect(QtWidgets.QStyle.SE_ItemViewItemText, option);

        t1 = datetime.strptime(str(index.data()).lstrip("database_sqlite_backup-"),"%d%m%Y-%H%M%S")
        dispStr = datetime.strftime(t1,"%a %d %b %Y %H:%M:%S")
        painter.drawText(textRect, option.displayAlignment, dispStr);

        painter.restore();

class CBackupDBBrowser(QtWidgets.QDialog):

    databaseRecoveryFile = QtCore.Signal(str)

    def backupFromXML(self):
        import tempfile
        from lxml import etree

        from dbapi.CCP4DbApi import CDbXml, CDbApi
        from core import CCP4Utils

        projectDirectories = []

        dbListBackupName = os.path.join(CCP4Utils.getDotDirectory(),'projectList-backup.xml')
        dbListBackupFile = open(dbListBackupName)
        dbListBackup = dbListBackupFile.read()
        dbListBackupFile.close()

        parser = ET.XMLParser()
        backupListTree = ET.fromstring(dbListBackup, parser)

        for p in backupListTree.findall("project"):
            projectDirectories.append(p.text)

        tFile = tempfile.NamedTemporaryFile(delete=True)
        dbFileName = tFile.name
        tFile.close()

        db = CDbApi(mode='sqlite',fileName=dbFileName,createDb=True)

        #Hmm, occasional foreign key problem. Do I have some dodgy projects?

        for projectDirectory in projectDirectories[:]:
            dbXML = os.path.join(projectDirectory,"backup_db.xml")
            if os.path.exists(dbXML):
                print("Trying",dbXML)
                try:
                    parser = ET.XMLParser()
                    f = open(dbXML)
                    s = f.read()
                    f.close()
                    tree = ET.fromstring(s, parser)
                    projectId = tree.findall('ccp4i2_header/projectId')[0].text
                    projectName = tree.findall('ccp4i2_header/projectName')[0].text
                    db.createProject(projectName,projectId=projectId,projectDirectory=projectDirectory)
                    dbImport = CDbXml(db=db,xmlFile=dbXML)
                    err = dbImport.loadTable()
                except:
                    pass

        db.close()
        return dbFileName


    def __init__(self,parent=None):
        QtWidgets.QDialog.__init__(self,parent)

        self.setWindowTitle("CCP4i2 database file corruption recovery")

        self.dbDirectory = None
        self.dbFiles = []
        self.mainLayout = QtWidgets.QVBoxLayout()
        self.setLayout(self.mainLayout)
        self.listWidget = QtWidgets.QListView()
        self.listWidget.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.listWidget.setItemDelegate(SimpleTextDelegate())
        self.model = QtCore.QStringListModel()
        self.infoModel = QtCore.QStringListModel()
        self.listWidget.setModel(self.model)

        warningWidget = QtWidgets.QWidget()
        warningLayout = QtWidgets.QHBoxLayout()
        warningWidget.setLayout(warningLayout)

        warningPixmapLabel = QtWidgets.QLabel()

        icon = QtWidgets.QApplication.style().standardIcon(QtWidgets.QStyle.SP_MessageBoxCritical)
        warningPixmap = icon.pixmap(QtCore.QSize(64, 64));
        warningPixmapLabel.setPixmap(warningPixmap)

        #warningLabel = QtWidgets.QLabel("<b>The default CCP4i2 database file seems to have become corrupted.<br>You can try to recover the database either from one of the backup files listed below or from a database file specified by you.<br/>Alternatively  you can attempt to reconstruct the database by scanning the project directories.  <br/><br/>You will need to restart CCP4i2 after this.")
        warningLabel = QtWidgets.QLabel("<b>The default CCP4i2 database file seems to have become corrupted.<br>You can try to recover the database either from one of the backup files listed below or from a database file specified by you.  <br/><br/>You will need to restart CCP4i2 after this.")

        warningLayout.addWidget(warningPixmapLabel)
        warningLayout.addWidget(warningLabel)
        self.mainLayout.addWidget(warningWidget)

        self.listLayout = QtWidgets.QGridLayout()
        self.knownRadio = QtWidgets.QRadioButton("Default backup files")
        self.mainLayout.addWidget(self.knownRadio)
        self.mainLayout.addLayout(self.listLayout)

        self.label = QtWidgets.QListView()
        self.label.setModel(self.infoModel)
        self.label.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.label.setSelectionMode(QtWidgets.QAbstractItemView.NoSelection)

        self.listLayout.addWidget(QtWidgets.QLabel("Backup file date"),0,0)
        self.listLayout.addWidget(self.listWidget,1,0)
        self.summaryLabel = QtWidgets.QLabel("Database File Project Summary")
        self.listLayout.addWidget(self.summaryLabel,0,1)
        self.listLayout.addWidget(self.label,1,1)

        self.specificRadio = QtWidgets.QRadioButton("Specific file")
        self.mainLayout.addWidget(self.specificRadio)

        specifyFileWidget = QtWidgets.QWidget()
        specifyFileLayout = QtWidgets.QHBoxLayout()
        specifyFileWidget.setLayout(specifyFileLayout)
        self.mainLayout.addWidget(specifyFileWidget)
        getFileButton = QtWidgets.QPushButton("Browse for database file")
        specifyFileLayout.addWidget(getFileButton)
        specificFileLabel = QtWidgets.QLabel("<None>")
        specifyFileLayout.addWidget(specificFileLabel)
        specifyFileLayout.addStretch(2)
        @QtCore.Slot()
        def getSpecificFile():
            specF,selectedFilter = QtWidgets.QFileDialog.getOpenFileName()
            if specF:
                specificFileLabel.setText(specF)
                self.setInfoList(str(specF))
                self.specificSelectedFile = str(specF)
                self.myStatusBar.showMessage("Selected file: "+self.specificSelectedFile)
        getFileButton.clicked.connect(getSpecificFile)

        self.specificRadio.clicked.connect(functools.partial(self.listWidget.setEnabled,False))
        self.knownRadio.clicked.connect(functools.partial(self.listWidget.setEnabled,True))
        self.specificRadio.clicked.connect(functools.partial(specifyFileWidget.setEnabled,True))
        self.knownRadio.clicked.connect(functools.partial(specifyFileWidget.setEnabled,False))

        specifyFileWidget.setEnabled(False)
        self.knownRadio.setChecked(True)

        self.xmlRadio = QtWidgets.QRadioButton("Attemept to recreate database by scanning project directories")
        #self.mainLayout.addWidget(self.xmlRadio)

        self.selectedFile = None
        self.specificSelectedFile = None
        self.defaultSelectedFile = None

        @QtCore.Slot()
        def setStatusBarToSpecific():
            if self.specificSelectedFile is not None:
                self.myStatusBar.showMessage("Selected file: "+self.specificSelectedFile)
                self.setInfoList(str(self.specificSelectedFile))
            else:
                self.myStatusBar.showMessage("Selected file: <None>")
                self.valid = False
                self.infoModel.setStringList([])

        @QtCore.Slot()
        def setStatusBarToDefault():
            if self.defaultSelectedFile is not None and len(self.listWidget.selectionModel().selectedIndexes()) > 0:
                self.myStatusBar.showMessage("Selected file: "+self.defaultSelectedFile)
                self.setInfoList(self.defaultSelectedFile)
            else:
                self.myStatusBar.showMessage("Selected file: <None>")
                self.valid = False
                self.infoModel.setStringList([])

        @QtCore.Slot()
        def setStatusBarToXML():
            self.infoModel.setStringList([])
            self.valid = True
            self.myStatusBar.showMessage("Selected file: <None> - Attempting to reconstruct database from XML files in project direcories")

        self.specificRadio.clicked.connect(setStatusBarToSpecific)
        self.knownRadio.clicked.connect(setStatusBarToDefault)
        self.xmlRadio.clicked.connect(setStatusBarToXML)

        self.xmlRadio.clicked.connect(functools.partial(self.listWidget.setEnabled,False))
        self.xmlRadio.clicked.connect(functools.partial(specifyFileWidget.setEnabled,False))

        @QtCore.Slot('QModelIndex')
        def changeLabel(index):
            f = os.path.join(self.dbOrigDir,str(self.model.data(index,QtCore.Qt.DisplayRole)))
            self.setInfoList(f)
            self.defaultSelectedFile = f
            self.myStatusBar.showMessage("Selected file: "+f)

        @QtCore.Slot('QItemSelection','QItemSelection')
        def changeLabelItem(sel,desel):
            if len(sel.indexes())>0:
                changeLabel(sel.indexes()[0])

        #self.listWidget.clicked.connect(changeLabel)
        self.listWidget.selectionModel().selectionChanged.connect(changeLabelItem)

        monFont = QtGui.QFont("Monospace");
        monFont.setStyleHint(QtGui.QFont.TypeWriter);
        if sys.platform == "darwin":
            #monFont = QtGui.QFont("Monaco");
            monFont = QtGui.QFont("Menlo");
        self.listWidget.setFont(monFont)
        self.label.setFont(monFont)

        buttons = QtWidgets.QDialogButtonBox()
        cancel = buttons.addButton(QtWidgets.QDialogButtonBox.Cancel)
        load1 = buttons.addButton("Load", QtWidgets.QDialogButtonBox.ActionRole)
        cancel.clicked.connect(self.reject)
        @QtCore.Slot()
        def askLoad():
            if self.knownRadio.isChecked():
                if len(self.listWidget.selectionModel().selectedIndexes()) == 0:
                    res = QtWidgets.QMessageBox.warning(None,"No backup database file specified","No backup database file has been specified. There are no default database backup files.")
                    return
                datafile = os.path.join(self.dbOrigDir,str(self.listWidget.selectionModel().selectedIndexes()[0].data()))
                if not self.valid:
                    res = QtWidgets.QMessageBox.warning(None,"Backup database invalid","The selected backup database file "+datafile+" appears to be invalid. You cannot recover the CCP4i2 database from this file. ")
                else:
                    if datafile is not None:
                        t1 = datetime.strptime(str(self.listWidget.selectionModel().selectedIndexes()[0].data()).lstrip("database_sqlite_backup-"),"%d%m%Y-%H%M%S")
                        dispStr = datetime.strftime(t1,"%a %d %b %Y %H:%M:%S")
                        res = QtWidgets.QMessageBox.question(None,"Load backup database","Copy backup database file "+datafile+" ("+dispStr+" ) to the default database location and use this as the CCP4i2 database? (This cannot be undone.)<br/><br/>You will need to restart CCP4i2 after this.",QtWidgets.QMessageBox.Cancel|QtWidgets.QMessageBox.Ok)
                        if res == QtWidgets.QMessageBox.Ok:
                            self.databaseRecoveryFile.emit(datafile)
                            self.selectedFile = datafile
                            self.accept()
            elif self.xmlRadio.isChecked():
                datafile = self.backupFromXML()
                print(datafile)
                res = QtWidgets.QMessageBox.question(None,"Load backup database","Copy reconstructed database file "+datafile+" to the default database location and use this as the CCP4i2 database? (This cannot be undone.)<br/><br/>You will need to restart CCP4i2 after this.",QtWidgets.QMessageBox.Cancel|QtWidgets.QMessageBox.Ok)
                if res == QtWidgets.QMessageBox.Ok:
                    self.databaseRecoveryFile.emit(datafile)
                    self.selectedFile = datafile
                    self.accept()
            else:
                datafile = self.specificSelectedFile
                if not self.valid:
                    if datafile is not None:
                        res = QtWidgets.QMessageBox.warning(None,"Backup database invalid","The selected backup database file "+datafile+" appears to be invalid. You cannot recover the CCP4i2 database from this file. ")
                    else:
                        res = QtWidgets.QMessageBox.warning(None,"No backup database specified","No backup database file has been specified")
                else:
                    if datafile is not None:
                        res = QtWidgets.QMessageBox.question(None,"Load backup database","Copy backup database file "+datafile+" ) to the default database location and use this as the CCP4i2 database? (This cannot be undone.)<br/><br/>You will need to restart CCP4i2 after this.",QtWidgets.QMessageBox.Cancel|QtWidgets.QMessageBox.Ok)
                        if res == QtWidgets.QMessageBox.Ok:
                            self.databaseRecoveryFile.emit(datafile)
                            self.selectedFile = datafile
                            self.accept()
                    else:
                        res = QtWidgets.QMessageBox.warning(None,"No backup database file specified","No backup database file has been specified.")
        load1.clicked.connect(askLoad)
        self.mainLayout.addWidget(buttons)
        self.myStatusBar = QtWidgets.QStatusBar()
        self.mainLayout.addWidget(self.myStatusBar)
        self.valid = False

    def setInfoList(self,f):
        try:
            conn = sqlite3.connect(f)
            cursor = conn.cursor()
            cursor.execute("SELECT * FROM Projects ORDER BY ProjectName ASC")
            rows = cursor.fetchall()
            projSums = []
            for row in rows:
                t = row[1]
                cursor.execute("SELECT * FROM Jobs WHERE ProjectID='"+row[0]+"'"+";")
                jobs = cursor.fetchall()
                t += " ("+str(len(jobs))+" jobs)"
                projSums.append(t)
            conn.close()
            self.infoModel.setStringList(projSums)
            self.summaryLabel.setText("Database File Project Summary ("+str(len(rows))+" projects)")
            self.valid = True
        except:
            self.infoModel.setStringList([])
            self.summaryLabel.setText("Database File Project Summary (0 projects, invalid file)")
            self.valid = False


    def setDBDirectory(self,dbDirectory):
        self.dbDirectory = dbDirectory
        self.draw()

    def draw(self):
        self.dbFiles = []
        if self.dbDirectory is None:
            return
        dbOrigName = os.path.join(self.dbDirectory,'database.sqlite')
        self.dbOrigDir = self.dbDirectory
        backups = glob.glob(os.path.join(self.dbOrigDir,"database_sqlite_backup-*[0-9]"))
        baseBackups = [ os.path.basename(f) for f in backups ]
        sortedBaseBackups = sorted(baseBackups,key=functools.cmp_to_key(datetimesort),reverse=True)
        print(sortedBaseBackups)
        sortedFullBaseBackups = [ os.path.join(self.dbOrigDir,f) for f in sortedBaseBackups ]
        for f in sortedFullBaseBackups:
            try:
                c = sqlite3.connect(f,5.0,1, check_same_thread = False)
                #FIXME - This needs to be more thorough. We need to check that all the required tables exist.
                c.execute('PRAGMA foreign_keys = ON')
                c.execute('PRAGMA synchronous=OFF')
                c.close()
                self.dbFiles.append(os.path.basename(f))
            except:
                continue
        self.model.setStringList(self.dbFiles)
        self.listWidget.setCurrentIndex(self.model.index(0, 0));

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    win = CBackupDBBrowser()
    win.setDBDirectory("/Users/stuart/.CCP4I2/db")

    @QtCore.Slot(str)
    def recoverTest(f):
        print(f)
    win.databaseRecoveryFile.connect(recoverTest)
    print(win.exec_())
