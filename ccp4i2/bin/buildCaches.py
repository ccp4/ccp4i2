from ..core import CCP4Config
from ..core import CCP4TaskManager


def main():
    CCP4Config.CConfig(graphical=True)
    taskManager=CCP4TaskManager.CTaskManager()
    taskManager.buildLookupFromScratch()
