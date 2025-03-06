from __future__ import print_function


import os, sys, time
import subprocess
from PySide2 import QtGui, QtWidgets, QtCore
from core.CCP4Bazaar import CUpdateUser

class CUpdateThread(QtCore.QThread, CUpdateUser):

  @staticmethod
  def _is_exec(name):
    proc = subprocess.Popen(('which', name), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdo, stde = proc.communicate()
    return os.path.basename(stdo.strip()) == name

  @classmethod
  def _set_sudo_cmd(cls):
    cls._sudo_cmd = None
    try:
      if sys.platform.startswith('darwin') and cls._is_exec('osascript'):
        command = '\\"%s\\" \\"%s\\"' %(sys.executable, cls._script)
        script ='do shell script "%s" with administrator privileges without altering line endings'
        cls._sudo_cmd = ['osascript', '-e', script %command]

      elif sys.platform.startswith('linux') and cls._is_exec('gksudo'):
        message = 'CCP4I2 Update needs your sudo password'
        cls._sudo_cmd = ['gksudo', '-g', '-m', message, sys.executable, cls._script]

    except Exception as e:
      if cls._verbose:
        print('ERRORi in_set_sudo_cmd:', e.args[0])
         

  @classmethod
  def _init_starting(cls):
    cls._verbose = False
    cls._force_compile = False

  @classmethod
  def _init_finished(cls):
    pass

  msg_info = QtCore.Signal(str)
  msg_action = QtCore.Signal(str)

  def __init__(self, updating):
    QtCore.QThread.__init__(self)
    self._updating = updating
    unknown = 'Unknown CCP4I2 version is running'
    running = 'CCP4I2 version %d is running'
    self._msg1_running = running %self._original if self._original else unknown
    self._msg1_current = self._msg1_running + ', %d is installed, %d is the latest.'
    self._msg1_running = self._msg1_running + '.'

  def run(self):
    self.msg_info.emit(self._msg1_running)
    self._bzr_run(self._updating)

  def _bzr_check_starting(self):
    self.msg_action.emit('Checking for updates, please wait...')

  def _bzr_apply_starting(self):
    self.msg_action.emit('Updating CCP4I2, please wait...')

  def _subprocess_starting(self):
    self.msg_action.emit('Updating CCP4I2, please wait...')

  def _bzr_check_finished(self):
    self._bzr_finished()

  def _bzr_apply_finished(self):
    self._bzr_finished()

  def _bzr_run_failed(self):
    self.msg_info.emit(self._msg1_running)
    self.msg_action.emit('Sorry, updater is off.')
    self._continue = False

  def _bzr_finished(self):
    ready = False
    if self._current < self._latest:
      if self._writable or self._sudo_cmd:
        ready = True
        msg2 = 'Please update and then restart CCP4I2.'

      else:
        msg2 = 'CCP4I4 is out of date. Please ask sysadmin to update.'

    elif self._original == 0:
      msg2 = 'The running version number will be set during the next update.'

    elif self._original < self._current:
      msg2 = 'The running version is out of date. Please restart CCP4I2.'

    else:
      msg2 = 'The running version is up to date.'

    msg1 = self._msg1_current %(self._current, self._latest)
    self.msg_info.emit(msg1)
    self.msg_action.emit(msg2)
    self._continue = ready

  def canContinue(self):
    return self._continue

CUpdateThread.initialise()

class CUpdateDialog(QtWidgets.QDialog):
  def __init__(self, parent=None):
    QtWidgets.QDialog.__init__(self, parent=parent)
    self.setWindowTitle('GUI Update')
    self.setModal(True)
    if parent:
      self.move(parent.x(), parent.y() + 99)

    self.label1 = QtWidgets.QLabel(108* ' ')
    self.label2 = QtWidgets.QLabel()

    self.buttonStart = QtWidgets.QPushButton('Update', self)
    self.buttonFinish = QtWidgets.QPushButton('Close', self)

    buttons = QtWidgets.QHBoxLayout()
    buttons.addWidget(self.buttonStart)
    buttons.addWidget(self.buttonFinish)

    self.setLayout(QtWidgets.QVBoxLayout())
    self.layout().addWidget(self.label1)
    self.layout().addWidget(self.label2)
    self.layout().addLayout(buttons)

    self.buttonStart.released.connect(self.checkStarted)
    self.buttonFinish.released.connect(self.close)

    self.show()
    self.checkStarted(False)

  @QtCore.Slot(bool)
  def checkStarted(self, updating=True):
    self.buttonStart.setEnabled(False)
    self.buttonFinish.setEnabled(False)
    self.check_thread = CUpdateThread(updating)
    self.check_thread.finished.connect(self.checkFinished)
    self.check_thread.msg_info.connect(self.label1.setText)
    self.check_thread.msg_action.connect(self.label2.setText)
    self.check_thread.start()

  @QtCore.Slot()
  def checkFinished(self):
    self.buttonFinish.setEnabled(True)
    if self.check_thread.canContinue():
      self.buttonStart.setEnabled(True)
      self.buttonStart.setFocus()

    else:
      self.buttonStart.setEnabled(False)
      self.buttonFinish.setFocus()

class CUpdateLauncher(QtWidgets.QDialog):
  def __init__(self):
    QtWidgets.QDialog.__init__(self)
    self.setWindowTitle('Main')

    self.buttonStart = QtWidgets.QPushButton('Start', self)
    self.buttonFinish = QtWidgets.QPushButton('Finish', self)
    self.buttonFinish.activateWindow()
    buttons = QtWidgets.QHBoxLayout()
    buttons.addWidget(self.buttonStart)
    buttons.addWidget(self.buttonFinish)
    self.setLayout(QtWidgets.QVBoxLayout())
    self.layout().addLayout(buttons)

    self.buttonStart.released.connect(self.start)
    self.buttonFinish.released.connect(self.close)

    self.activateWindow()
    self.show()

  @QtCore.Slot()
  def start(self):
    self.q = CUpdateDialog(self)


if __name__ == '__main__':
  app = QtWidgets.QApplication(sys.argv)
  q = CUpdateLauncher()
  app.exec_()

