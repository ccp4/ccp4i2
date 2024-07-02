"""
     CCP4ManageDbGui.py: CCP4 GUI Project
     Copyright (C) 2014 STFC

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
"""

"""
     Liz Potterton Sep 2014 - Tools to manage database
"""

from PySide6 import QtGui,QtWidgets,QtCore
  
class CChallengeUnknownUser(QtWidgets.QDialog):
    # This class is being used to create a pop-up (in the initialiser) as part of a db user/passwd system
    # I've changed this to make it a bit easier for users, & less misleading, should they ever encounter a db username mix-up.

    def __init__(self, parent=None, currentUserName='UNKNOWN', origUserName='UNKNOWN'):
        QtWidgets.QDialog.__init__(self, parent=None)
        self.origUserName= origUserName
        self.setModal(True)
        self.setWindowTitle('Unknown user')
        layout = QtWidgets.QVBoxLayout()
        self.setLayout(layout)
        self.layout().addWidget(QtWidgets.QLabel('Your current username \'' + currentUserName 
                                             + '\' is not in the ccp4i2 database', self))
        self.layout().addWidget(QtWidgets.QLabel('This can happen if the i2 database was created under a different username (eg. \'' 
                                             + self.origUserName + '\' )', self))
        font = QtGui.QFont()
        font.setBold(True)
        lc = QtWidgets.QLabel('\nYou can add \'' + currentUserName + '\' to the database by clicking apply and i2 can be then used with this account.', self)
        lc.setFont(font)
        self.layout().addWidget(lc)
        self.layout().addWidget(QtWidgets.QLabel('If you do not wish to do this, you will need to log into the computer account used to create the i2 database.\n', self))
        line = QtWidgets.QHBoxLayout()
        self.layout().addLayout(line)
        # This is clearly not right, currently a game of guess the prior owner. Need to just get this info from the db.
        #line.addWidget(QtWidgets.QLabel('Enter the name of the database owner', self))
        #self.userNameWidget = QtWidgets.QLineEdit(self)
        #line.addWidget(self.userNameWidget)
        # The following two lines are def not ok, checking that tick-box will cause i2 to fail. 
        # Will make a database look corrupt to users, when it's actually not corrupt at all.
        #self.resetName = QtWidgets.QCheckBox('Reset name of database owner to your username', self)
        #self.layout().addWidget(self.resetName)
        box = QtWidgets.QDialogButtonBox(self)
        button_y = box.addButton(QtWidgets.QDialogButtonBox.Apply)
        self.connect(button_y, QtCore.SIGNAL('clicked()'), self.accept)
        button_n = box.addButton(QtWidgets.QDialogButtonBox.Cancel)
        self.connect(button_n, QtCore.SIGNAL('clicked()'), self.reject)
        self.layout().addWidget(box)
        self.raise_()

    def getOwner(self):
        return self.origUserName
        #return self.userNameWidget.text().__str__()

    def getReset(self):
        return False #self.resetName.isChecked()

class CDbManagerDialog(QtWidgets.QDialog):

    def __init__(self,parent=None):
        QtWidgets.QDialog.__init__(self,parent)
        self.setWindowTitle('Database manager')
        layout = QtWidgets.QVBoxLayout()
        self.setLayout(layout)
