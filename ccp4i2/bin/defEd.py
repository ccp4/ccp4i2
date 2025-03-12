import sys

from ..core.CCP4Utils import getCCP4I2Dir


def main():
    print('Running CCP4i2 browser from:', getCCP4I2Dir())
    print('Python', sys.version)
    print(' ')
    from ..utils.startup import setupEnvironment()
    setupEnvironment()
    from ..core.CCP4Modules import QTAPPLICATION
    app = QTAPPLICATION(graphical=True)
    from ..core.CCP4Config import CConfig
    from ..qtgui.CCP4DefEd import CDefEd
    CConfig(qt=True, graphical=True, developer=True)
    defEd = CDefEd()
    defEd.raise_()
    sys.exit(app.exec_())
