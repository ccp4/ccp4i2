import os
import sys

from ..core.CCP4Modules import QTAPPLICATION
from ..core.CCP4Utils import getCCP4I2Dir


if __name__ == '__main__':
    top_path = getCCP4I2Dir()
    print('Running CCP4i2 browser from: '+top_path)
    print('Python', sys.version)
    print(' ')
    path = os.path.join(top_path,'utils','startup.py')
    exec(compile(open(path).read(), path, 'exec'))
    setupEnvironment()
    app = QTAPPLICATION(graphical=True)
    defEd = startDefEd()
    sys.exit(app.exec_())
