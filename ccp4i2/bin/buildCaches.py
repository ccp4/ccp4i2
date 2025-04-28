from ..core import CCP4TaskManager
from ..core.CCP4Config import CONFIG


def main():
    CONFIG(graphical=True)
    taskManager=CCP4TaskManager.CTaskManager()
    taskManager.buildLookupFromScratch()
