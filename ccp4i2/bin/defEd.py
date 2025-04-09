import sys

from ..core.CCP4Config import CONFIG
from ..core.CCP4Utils import getCCP4I2Dir
from ..utils.QApp import QTAPPLICATION


def main():
    print('Running CCP4i2 browser from:', getCCP4I2Dir())
    print('Python', sys.version)
    print(' ')
    app = QTAPPLICATION(graphical=True)
    from ..qtgui.CCP4DefEd import CDefEd
    CONFIG(graphical=True, developer=True)
    defEd = CDefEd()
    defEd.raise_()
    sys.exit(app.exec_())
