import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core import CCP4TaskManager
from core import CCP4Config
a=CCP4Config.CConfig()
a.insts.graphical=True
taskManager=CCP4TaskManager.CTaskManager()
taskManager.buildLookupFromScratch()
