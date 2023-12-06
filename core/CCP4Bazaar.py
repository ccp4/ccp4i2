from __future__ import print_function


import os, sys
import compileall
import subprocess
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import shutil

try:
  import bzrlib
  bzrlib_exists = True

except ImportError:
  bzrlib_exists = False


class BazaarContext(list):

    # stdout and stderr cannot be redirected using
    # bzrlib.initialize because of bugs in bzrlib

    verbose = False

    @classmethod
    def set_verbose(cls, value):
        cls.verbose = value

    def __init__(self):
        self.ostream = StringIO()

    def __enter__(self):
        sys.stderr = sys.stdout = self.ostream
        os.environ["BZR_PLUGIN_PATH"]=''
        self.bzrlib_state = bzrlib.initialize()
        self.bzrlib_state.__enter__()
        return self

    def _dump(self):
        if self.verbose:
            sys.__stdout__.write(self.ostream.getvalue())
            sys.__stdout__.flush()

        self.ostream.truncate(0)

    def bzr(self, *args):
        cmd = ('bzr',) + args
        print('cmd:')
        print(' '.join(['"%s"' %w if ' ' in w else w for w in cmd]))
        self._dump()
        retval = bzrlib.commands.main(cmd)
        self.ostream.seek(0)
        self[:] = self.ostream.readlines()
        print("returned:", retval)
        print()
        return retval

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.bzrlib_state.__exit__(None, None, None)
        del os.environ["BZR_PLUGIN_PATH"]
        self._dump()
        self.ostream.close()
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        return False


class CUpdateUser(object):

    @classmethod
    def _diff_uid(cls, path, uid):
        if os.stat(path).st_uid != uid:
            return True

        if os.path.isdir(path):
            for name in os.listdir(path):
                if cls._diff_uid(os.path.join(path, name), uid):
                    return True

    @classmethod
    def _set_writable(cls):
        cls._writable = False
        cls._writable = not (cls._uid and cls._diff_uid(cls._checkout, cls._uid))

    def _write_original(self):
        with open(self._version_file, 'w') as ostream:
            print(self._current, file=ostream)

    @classmethod
    def _set_original(cls):
        cls._original = 0
        try:
            with open(cls._version_file) as istream:
                cls._original = int(istream.read().strip())

        except:
            pass

    @classmethod
    def _copy_uid(cls, dir):
        uid = os.stat(dir).st_uid
        for name in os.listdir(dir):
            path = os.path.join(dir, name)
            if os.stat(path).st_uid != uid:
                os.chown(path, uid, -1)

            if os.path.isdir(path):
                cls._copy_uid(path)

    @classmethod
    def _set_path(cls, path):
        return os.path.relpath(path, os.getcwd()) if cls._verbose else os.path.realpath(path)

    @classmethod
    def _set_uid(cls):
        try:
            cls._uid = os.getuid()

        except:
            cls._uid = None

    @classmethod
    def initialise(cls):
        cls._init_starting()
        cls._script = cls._set_path(__file__)
        cls._checkout = cls._set_path(os.path.join(os.path.dirname(__file__), '..'))
        cls._version_file = os.path.join(cls._checkout, '.bzr', 'running_version')
        BazaarContext.set_verbose(cls._verbose)
        cls._set_uid()
        cls._set_sudo_cmd()
        cls._set_original()
        cls._init_finished()

    def _bzr_run(cls, updating):
        cls._errmsg = None
        try:
            cls._set_writable()
            if updating:
                if cls._writable:
                    cls._bzr_apply_starting()
                    for bname in ('limbo', 'pending-deletion'):
                        dname = os.path.join(cls._checkout, '.bzr', 'checkout', bname)
                        if os.path.exists(dname):
                            shutil.rmtree(dname)

                    retcode = 99
                    line = ''
                    with BazaarContext() as rc:
                        retcode = rc.bzr('break-lock', '--force', cls._checkout)
                        line = rc[0] if rc else ''

                    if retcode or line:
                        cls._bzr_unlock_message(line, retcode)

                    if retcode:
#                       cls._original = 0
#                       cls._write_original()
                        raise Exception('Repository is buisy, please try later.')

                    with BazaarContext() as rc:
                        rc.bzr('update', cls._checkout)
                        rc.bzr('resolve', '--all', '--take-other', '-d', cls._checkout)
                        rc.bzr('revert', '--no-backup', cls._checkout)
                        rc.bzr('diff', cls._checkout)

                elif cls._sudo_cmd:
                    cls._subprocess_starting()
                    subprocess.Popen(cls._sudo_cmd).wait()

            else:
                cls._bzr_check_starting()

            cls._latest = 0
            cls._current = 0
            with BazaarContext() as rc:
                rc.bzr('revno', cls._checkout)
                word = rc[0].strip()
                if word.isdigit():
                    cls._latest = int(word)

                rc.bzr('revno', cls._checkout, '--tree')
                word = rc[0].strip()
                if word.isdigit():
                    cls._current = int(word)

            if not (cls._latest and cls._current and cls._current <= cls._latest):
                raise Exception('Update check failed')

            if updating:
                if not (cls._latest and cls._current and cls._current == cls._latest):
                    raise Exception('Update failed')

                if cls._writable:
                    cls._write_original()
                    compileall.compile_dir(cls._checkout, force=cls._force_compile, quiet=not cls._verbose)
                    if cls._verbose:
                        print()

                    if cls._uid == 0:
                        cls._copy_uid(cls._checkout)

                else:
                    cls._set_writable()

                return cls._bzr_apply_finished()

            else:
                return cls._bzr_check_finished()

        except Exception as e:
            cls._errmsg = e.args[0]
            cls._bzr_run_failed()

    @classmethod
    def _init_starting(cls):
        cls._command = sys.argv.pop()
        import argparse
        parser = argparse.ArgumentParser(
           prog=cls._command,
           description=' '.join((
              'This is an Updater for CCP4I2.',
              'NB: sudo is invoked autoamtically if needed.',
           ))
        )
        parser.add_argument(
           '-v', '--verbose',
           help='print stdout and stderr from bzr',
           action='store_true'
        )
        parser.add_argument(
           '-f', '--force-compile',
           help='force (re)compilation of all py-files',
           action='store_true'
        )
        opt = parser.parse_args()
        cls._verbose = opt.verbose
        cls._force_compile = opt.force_compile

    @classmethod
    def _init_finished(cls):
        if cls._verbose:
            print('print_info:')
            print()
            print('checkout', cls._checkout)
            print('version_file', cls._version_file)
            print('command', cls._command)
            print('verbose', cls._verbose)
            print('force_compile', cls._force_compile)
            print('uid', cls._uid)
            print('sudo_cmd', cls._sudo_cmd)
            print('original', cls._original)
            print()

    @classmethod
    def _set_sudo_cmd(cls):
        cls._sudo_cmd = None

    def _bzr_check_starting(cls):
        sys.exit(1)

    def _bzr_apply_starting(cls):
        print('Updating, please wait...')

    def _subprocess_starting(cls):
        sys.exit(1)

    def _bzr_check_finished(cls):
        sys.exit(1)

    def _bzr_apply_finished(cls):
        pass

    def _bzr_run_failed(cls):
        sys.exit(1)

    def _bzr_unlock_message(cls, *args):
        pass

if __name__ == '__main__':

    if bzrlib_exists:
        CUpdateUser.initialise()
        CUpdateUser()._bzr_run(updating=True)

