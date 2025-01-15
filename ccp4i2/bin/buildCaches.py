from ..core import CCP4TaskManager
from ..core import CCP4Config

a=CCP4Config.CConfig()
a.insts.graphical=True
taskManager=CCP4TaskManager.CTaskManager()
taskManager.buildLookupFromScratch()
