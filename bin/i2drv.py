from __future__ import print_function

import sys
import os

def getCCP4MG_DIR():
    return 'C:\Program Files (x86)\CCP4MG'

def getCCP4I2_DIR():
    target = os.path.join(os.path.realpath(sys.argv[0]),"..")
    abstarget = os.path.abspath(target)
    targetdir = os.path.normpath(os.path.dirname(abstarget))
    return targetdir

def setup_pythonpath():
    CCP4I2 = getCCP4I2_DIR()
    sys.path.append(CCP4I2)
    sys.path.append(os.path.join(CCP4I2,"lib"))

def setup_environment():
    ROOT = getCCP4I2_DIR()
    PYTHONHOME=os.path.join(ROOT,"pythondist")

    path_sep = ';'

    DYLD_FALLBACK_LIBRARY_PATH=ROOT + path_sep + os.path.join(ROOT,"lib")


    os.environ["DYLD_FALLBACK_LIBRARY_PATH"] = DYLD_FALLBACK_LIBRARY_PATH
    os.environ["LD_LIBRARY_PATH"] = DYLD_FALLBACK_LIBRARY_PATH
    os.environ["LIBPATH"] = DYLD_FALLBACK_LIBRARY_PATH

    OLDPATH=os.environ["PATH"]
    PATH=os.path.join(ROOT,"bin") + path_sep + os.path.join(ROOT,"pythondist","bin") + path_sep + DYLD_FALLBACK_LIBRARY_PATH + path_sep + OLDPATH
    os.environ["PATH"] = PATH
    os.environ["PYTHONPATH"] = PATH
    if "HOME" not in os.environ or os.environ["HOME"] == "":
        try:
            HOME=os.environ["USERPROFILE"]
            os.environ["HOME"] = HOME
        except:
            pass

    os.environ["CCP4MG"] = getCCP4MG_DIR()
    os.environ["CCP4I2"] = ROOT
    os.environ["CCP4I2_TOP"] = ROOT
    os.environ["PYTHONHOME"] = PYTHONHOME


if __name__ == '__main__':
    setup_environment()
    ccp4i2 = getCCP4I2_DIR()
    mainpy = os.path.join("bin","browser.py")
    import subprocess
    args = ["python",mainpy]
    for sysarg in sys.argv[1:]:
        args.append(sysarg)
    retval = subprocess.call(args)

    print("Finished")
