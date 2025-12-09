from __future__ import print_function
import sys
import os
import glob
import tempfile
import shutil
import sqlite3

from os.path import expanduser

from ccp4i2.baselayer import QtCore, QtGui, QtWidgets

import reconstructDBFromXML
if __name__ == "__main__":
    sys.path.append(os.path.join(os.path.dirname(__file__),".."))
import importDir

class ReconstructBrowserDialog(QtWidgets.QDialog):

    @QtCore.Slot()
    def addDirectory(self):
        root = str(QtWidgets.QFileDialog.getExistingDirectory(self, self.tr("Open Directory"), expanduser("~"), QtWidgets.QFileDialog.ShowDirsOnly | QtWidgets.QFileDialog.DontResolveSymlinks))
        print(root)
        if root:
            if os.path.exists(os.path.join(root,"CCP4_JOBS")) and  os.path.exists(os.path.join(root,"CCP4_PROJECT_FILES")):
                item = QtWidgets.QListWidgetItem(root,self.listWidget)
                item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
                item.setCheckState(QtCore.Qt.Checked) 
            else:
                QtWidgets.QMessageBox.critical(self,self.tr("Not a project directory"),self.tr("Invalid directory: does not contain CCP4_JOBS and CCP4_PROJECT_FILES"))

    def __init__(self,parent=None):
        QtWidgets.QDialog.__init__(self,parent)
        layout = QtWidgets.QVBoxLayout()
        self.setLayout(layout)
        self.listWidget = QtWidgets.QListWidget()
        layout.addWidget(self.listWidget)

        ccp4i2_root_dir = os.path.join(expanduser("~"),"CCP4I2_PROJECTS")
        for root in glob.glob(os.path.join(ccp4i2_root_dir,"*")):
            if os.path.exists(os.path.join(root,"CCP4_JOBS")) and  os.path.exists(os.path.join(root,"CCP4_PROJECT_FILES")):
                item = QtWidgets.QListWidgetItem(root,self.listWidget)
                item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
                item.setCheckState(QtCore.Qt.Checked) 

        backup_file_layout = QtWidgets.QHBoxLayout()
        self.backup_file_name = os.path.join(expanduser("~"),"ccp4i2_backup_db.sqlite")
        backup_file_label = QtWidgets.QLabel("Output database file: "+self.backup_file_name)
        backup_file_layout.addWidget(backup_file_label)
        change_backup_file_but = QtWidgets.QPushButton("Change")
        backup_file_layout.addWidget(change_backup_file_but)
        layout.addLayout(backup_file_layout)
        @QtCore.Slot()
        def change_backup_location():
            save_browser = QtWidgets.QFileDialog()
            save_browser.selectFile(self.backup_file_name)
            save_browser.setAcceptMode(QtWidgets.QFileDialog.AcceptSave);
            save_browser.setOption(QtWidgets.QFileDialog.DontUseNativeDialog);
            save_browser.setLabelText(QtWidgets.QFileDialog.Accept,"Set file name");
            fn = save_browser.exec_()
            if fn:
                print("Got a new filename",save_browser.selectedFiles()[0])
                self.backup_file_name = str(save_browser.selectedFiles()[0])
                backup_file_label.setText("Output database file: "+self.backup_file_name)
        change_backup_file_but.clicked.connect(change_backup_location)

        dbb = QtWidgets.QDialogButtonBox()
        layout.addWidget(dbb)

        okButton = dbb.addButton(QtWidgets.QDialogButtonBox.Ok)
        cancelButton = dbb.addButton(QtWidgets.QDialogButtonBox.Cancel)
        clearAllButton = dbb.addButton("Clear",QtWidgets.QDialogButtonBox.ActionRole)
        selectAllButton = dbb.addButton("Select all",QtWidgets.QDialogButtonBox.ActionRole)
        addButton = dbb.addButton("Add",QtWidgets.QDialogButtonBox.ActionRole)

        @QtCore.Slot()
        def clearAll():
            for i in range(self.listWidget.count()):
                self.listWidget.item(i).setCheckState(QtCore.Qt.Unchecked)

        @QtCore.Slot()
        def selectAll():
            for i in range(self.listWidget.count()):
                self.listWidget.item(i).setCheckState(QtCore.Qt.Checked)

        clearAllButton.clicked.connect(clearAll)
        selectAllButton.clicked.connect(selectAll)
        cancelButton.clicked.connect(self.reject)
        okButton.clicked.connect(self.accept)
        addButton.clicked.connect(self.addDirectory)

if __name__ == "__main__":
    from lxml import etree

    app = QtWidgets.QApplication(sys.argv)
    win = ReconstructBrowserDialog()
    win.setWindowTitle("Select project directories to reconstruct from")
    dirs = []
    if win.exec_():
        mb = QtWidgets.QMessageBox(QtWidgets.QMessageBox.Information,"Generating database","Generating database from project files. Please wait.")
        mb.setStandardButtons(QtWidgets.QMessageBox.NoButton)
        mb.show()
        mb.raise_()
        app.processEvents()
        fileName = win.backup_file_name
        for i in range(win.listWidget.count()):
            if win.listWidget.item(i).checkState() == QtCore.Qt.Checked:
                t = win.listWidget.item(i).text()
                if hasattr(t,"toString"):
                    dirs.append(str(t.toString()))
                else:
                    dirs.append(t)

        tfile = tempfile.NamedTemporaryFile(delete=False)
        tfn = tfile.name+".sqlite"
        tfile2 = tempfile.NamedTemporaryFile(delete=False)
        tfn2 = tfile2.name
        tfile.close()
        tfile2.close()

        for d in dirs:
            mb.setText("Generating database from project files. Please wait.\nCurrently processing project directory: "+str(d))
            app.processEvents()
            print("------------------------------------------------------------")
            print(d)
            project_tree = reconstructDBFromXML.generate_xml_from_project_directory(str(d))
            if len(project_tree.xpath("//ccp4i2_body/jobTable")) == 0:
                continue
            if len(project_tree.xpath("//ccp4i2_body/jobTable/job")) == 0:
                continue
            if sys.version_info < (3,0):
                outl = etree.tostring(project_tree,pretty_print=True)
            else:
                outl = etree.tostring(project_tree,pretty_print=True).decode()
            dbxmlout = os.path.join(str(d),"DATABASE.db.xml")
            with open(dbxmlout,"w+") as outfd:
                outfd.write(outl)

            importDir.importFilesFromDirXML(tfn,str(d))

        #Now we use sqlite api to make fresh, "clean" hopefully copy of new db...
        con = sqlite3.connect(tfn)
        sql = "".join([s+"\n" for s in con.iterdump()])

        f2trunc = open(tfn2,"w")
        f2trunc.close()
        print("Writing to",tfn2)

        conbak = sqlite3.connect(tfn2)
        cur = conbak.cursor()

        try:
            cur.executescript(sql)
        except:
            print("Fail",com)
            conbak.close()
            raise

        conbak.commit()
        
        #... and then copy to desited location.
        shutil.copy(tfn2,fileName)

        print()
        print("############################################################")
        print("Saved backup database to",fileName)
        print("############################################################")
        conbak.close()

