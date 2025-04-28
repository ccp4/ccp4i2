
import subprocess as SP
import sys

from PySide2 import QtCore, QtWidgets

from ..qtgui import CCP4WebBrowser


class CUpdateMessageBox(QtWidgets.QMessageBox):

  exit_requested = QtCore.Signal()

  def __init__(self, available, w_access):
    QtWidgets.QMessageBox.__init__(self)
    self.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
    self.setWindowTitle('CCP4 Updates')
    if w_access > 0:
      self.setStandardButtons(QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
      self.setDefaultButton(QtWidgets.QMessageBox.Yes)
      self.setIcon(QtWidgets.QMessageBox.Question)

    else:
      self.setStandardButtons(QtWidgets.QMessageBox.Ok)
      self.setDefaultButton(QtWidgets.QMessageBox.Ok)
      self.setIcon(QtWidgets.QMessageBox.Information)

    t_new = 'New CCP4 updates are available.'
    t_proceed = 'Proceed with installation?'

    t_um = 'CCP4 Update Manager allows to examine,'
    t_um += ' apply or remove updates. Proceed?'

    t_admin = 'You will be asked for administrative password.'

    t_restart = 'The CCP4 interface will be turned off'
    t_restart += ' while you are managing updates.'

    t_run_cmd = 'Run Update Manager from administrative account'
    t_run_cmd += ' or from the command line (ccp4um).'

    t_release = 'NEW CCP4 VERSION HAS BEEN RELEASED'
    t_no_new = 'The running version of the CCP4 software'
    t_no_new += ' is no longer supported and will not be updated.'
    t_no_new += ' However, you can manage existing updates.'
    t_continue = 'Continue?'

    if w_access > 0:
      if available == 254:
        self.summary = t_release
        self.next_summary = t_release
        if w_access == 1:
          info = (t_no_new, t_admin, t_continue)

        else:
          info = (t_no_new, t_continue)

      else:
        self.summary = ' '.join((t_new, t_proceed))
        self.next_summary = t_um
        if w_access == 1:
          info = (t_admin, t_restart)

        else:
          info = (t_restart,)

    else:
      if available == 254:
        self.summary = t_release
        self.next_summary = t_release
        info = (t_no_new, t_run_cmd)

      else:
        self.summary = ' '.join((t_new, t_run_cmd))
        self.next_summary = t_run_cmd

        info = ()

    self.setInformativeText('\n\n'.join(info))
    self.first = True
    self.first_show = 0 < available and available < 255

  def execute(self):
    text = self.summary if self.first else self.next_summary
    self.setText('<b><font size="+1">' + text + '</font></b>')
    if not self.first or self.first_show:
      if self.exec_() == self.Yes:
        self.exit_requested.emit()

    self.first = False

class CUpdateManager(QtCore.QThread):

  cmd = sys.executable, '-m', 'um.wrapper', '-stamp', '2'
  arg = '-check-silent', '-timeout', '33000'

  def __init__(self):
    QtCore.QThread.__init__(self)
    self.finished.connect(self.next)
    self.lock = True
    self.mbox = None
    self.info = None
    QtCore.QTimer.singleShot(1000, self.start)

  def run(self):
    cmd_check = self.cmd + self.arg
    try:
      proc = SP.Popen(cmd_check, stdout=SP.PIPE, stderr=SP.PIPE)
      stdo, stde = proc.communicate()
      if not proc.wait():
        items = stdo.split()
        self.info = int(items[0]), int(items[1])

    except:
      pass

  @QtCore.Slot()
  def next(self):
    if self.info:
      self.mbox = CUpdateMessageBox(*self.info)
      self.mbox.exit_requested.connect(self.exit)
      self.mbox.execute()
      self.lock = False

  def is_unlocked(self):
    return not self.lock

  def manage(self):
    if not self.lock:
      self.mbox.execute()

  @QtCore.Slot()
  def exit(self):
    app = QtWidgets.QApplication.instance()
    app.setQuitOnLastWindowClosed(False)
    CCP4WebBrowser.exitBrowser()
    app.exit(0)
    SP.Popen(self.cmd)
    # bug on windows: does not accept paths with spaces 
    # os.execv(self.cmd[0], self.cmd)

um = CUpdateManager()

