import sys

from ..core.CCP4Config import CConfig
from ..core.CCP4Modules import QTAPPLICATION
from ..core.CCP4Utils import getCCP4I2Dir
from ..qtgui.CCP4DefEd import CDefEd
from ..utils.startup import setupEnvironment


def main():
    print('Running CCP4i2 browser from:', getCCP4I2Dir())
    print('Python', sys.version)
    print(' ')
    setupEnvironment()
    app = QTAPPLICATION(graphical=True)
    CConfig(qt=True, graphical=True, developer=True)
    defEd = CDefEd()
    defEd.raise_()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
