from __future__ import print_function


# This file should only be used in the script $CCP4/bin/ccp4i2-update
# Currently it is source in several palces in order to call get_revno

import os, sys
from core.CCP4Bazaar import CUpdateUser
from core.CCP4Bazaar import bzrlib_exists

class CUpdate(CUpdateUser):

    # reusing from super:
    # _init_starting
    # _init_finished
    # _bzr_apply_starting
    # _bzr_unlock_message

    @classmethod
    def _set_sudo_cmd(cls):
        cls._sudo_cmd = ['sudo', sys.executable, os.path.relpath(cls._script, os.getcwd())]
        cls._sudo_cmd.extend(sys.argv[1:])
        cls._sudo_cmd.append(cls._command)

    def _bzr_check_starting(self):
        print('Checking for new updates, please wait...')

    def _subprocess_starting(self):
        fmt = 'Updater needs your sudo password to make changes within "%s"'
        print(fmt %self._checkout)

    def _bzr_check_finished(self):
        fmt = 'Currently installed CCP4I2 version %s; the latest available version %s'
        print(fmt %(self._current, self._latest))

    def _bzr_apply_finished(self):
        fmt = 'Your CCP4I2 installation has been updated to version %s.'
        print(fmt %self._current)
        print('Please Restart CCP4I2 if it is currently running.')

    def _bzr_run_failed(self):
        print('ERROR:', self._errmsg)
        if self._verbose:
            print('See bzr messages above for more details.')

        else:
            print('For more detailes, run "%s -v"' %self._command)

        sys.exit(1)

    def _bzr_unlock_message(cls, line, retcode):
        print('bzr break-lock: retcode=%s, message=%s' %(str(retcode), line))

    def main(self):
        self._bzr_run(updating=False)
        if self._current < self._latest:
            self._bzr_run(updating=True)

        else:
            print("Nothing to do.")

  
class CUpdateForTestSys(CUpdate):
  @classmethod
  def initialise(cls):
    CUpdate.initialise()
    from core import CCP4Utils
    cls._command =CCP4Utils.getCCP4I2Dir()     

if bzrlib_exists:
    CUpdate.initialise()
    if CUpdate._original:
        running_version = str(CUpdate._original)

    else:
        running_version = '0'

else:
    try:
        istr = open(os.path.join(os.environ['CCP4'], 'lib', 'ccp4', 'MAJOR_MINOR'))
        running_version = istr.read().strip()
        istr.close()
    except:
        running_version = '9.0.000'

def get_revno():
    numeric = ''.join(running_version.split('.'))
    return int(numeric) if numeric.isdigit() else 0

def get_ccp4_str():
    if bzrlib_exists:
        return 'CCP4i2 alpha-' + running_version

    else:
        return 'CCP4-' + running_version

if __name__ == '__main__':
    if bzrlib_exists:
        CUpdate().main()

