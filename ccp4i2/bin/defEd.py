from __future__ import print_function


import os,sys
from PySide2 import QtCore

def getCCP4I2Dir(up=1):
    target = os.path.join(os.path.realpath(sys.argv[0]),"..")
    abstarget = os.path.abspath(target)
    splittarget = abstarget.split()
    if splittarget.count('ccp4i2'):
        splittarget.reverse()
        up = splittarget.index('ccp4i2')
    while up>0:
        abstarget = os.path.dirname(abstarget)
        up = up -1
    return abstarget

if __name__ == '__main__':

    top_path = getCCP4I2Dir()
    print('Running CCP4i2 browser from: '+top_path)
    print('Python '+sys.version)
    print(' ')
    exec(compile(open(os.path.join(top_path,'utils','startup.py')).read(), os.path.join(top_path,'utils','startup.py'), 'exec'))
    setupEnvironment()
    setupPythonpath(top=top_path,mode='qtgui')
    setupGuiPluginsPath(top=top_path)
    from ccp4i2.core.CCP4Modules import QTAPPLICATION
    app = QTAPPLICATION(graphical=True)

    defEd = startDefEd()

    sys.exit(app.exec_())
