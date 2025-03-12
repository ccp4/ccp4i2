from ..core import CCP4TaskManager
from ..core import CCP4Config


def main():
    config = CCP4Config.CConfig()
    config.insts.graphical = True
    taskManager=CCP4TaskManager.CTaskManager()
    taskManager.buildLookupFromScratch()
