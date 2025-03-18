from argparse import ArgumentParser
from datetime import datetime
import getpass
import os
import sys
import time

from ..utils.startup import setupEnvironment


def checkForPythonNameClash(nameRoot):
    from ..core import CCP4TaskManager
    CCP4TaskManager.TASKMANAGER()
    scriptClass = nameRoot
    guiClass = nameRoot+"_gui"
    reportClass = nameRoot +"_report"
    for className in [scriptClass, guiClass, reportClass]:
        print("Trying to import ",className,"...")
        try:
            __import__(className)
            print("Proposed class name clashes with something else in the namespace (i.e. import "+className+" already succeeds)")
            return True
        except:
            print("... No clash for className "+className)
    print("Proposed class "+nameRoot+" does not clash with anything in the ccp4i2 python namespace")
    return False


def makeScriptDirectory(wrapperScriptDirectory):
    try:
        os.makedirs(wrapperScriptDirectory)
        return 0
    except Exception as err:
        print("makeScriptDirectory generated error ", err)
        return 1


def main():
    CCP4I2_TOP= os.path.abspath(os.environ["CCP4I2"])
    setupEnvironment(path=CCP4I2_TOP)
    from ..core import CCP4TaskManager
    destinations = CCP4TaskManager.MODULE_ORDER
    parser = ArgumentParser(description='Initiate CCP4i2 plugin from boiler plate')
    parser.add_argument('-n','--name', required=True, help='name of the plugin...this will end up as a class name and so should have no spaces')
    parser.add_argument('-m','--maintainer',required=False, default="person@server.com",help='e-mail address of user to register as the maintainer of this plugin')
    parser.add_argument('-u','--user', required=False, help='user to be credited with initiating the def.xml of this task')
    parser.add_argument('-s','--shortTitle', nargs='+', required=True, help='short title of task')
    parser.add_argument('-l','--longTitle', nargs='+', required=True, help='more descriptive title of task to be presented in task selection gui')
    parser.add_argument('-d','--description', nargs='+', required=True, help='longer description of task to be presented in task selection gui')
    parser.add_argument('-o','--overwrite', action='store_const', const=True, help='Override namespace clashes and pre-existing files')
    parser.add_argument('-t', '--type', choices=['pipeline','plugin'], default='plugin', help='Use plugin or pipeline boilerplate template')
    parser.add_argument('-w', '--where', choices=destinations, default='refinement', help='Section of the task selection GUI where task should show up')
    parser.add_argument('-f','--firstPlugin', type=str, default='firstPlugin', help='For pipeline code, optionally provide name of first plugin to be called')
    parameterNamespace = parser.parse_args()
    
    user = getpass.getuser()
    nameRoot = parameterNamespace.name
    shortTitle = " ".join(parameterNamespace.shortTitle)
    longTitle = " ".join(parameterNamespace.longTitle)
    description = " ".join(parameterNamespace.description)
    
    timestampString=datetime.fromtimestamp(time.time()).isoformat()

    if parameterNamespace.user is not None: user = parameterNamespace.user

    nameClashed = checkForPythonNameClash(nameRoot)
    if nameClashed and not parameterNamespace.overwrite: sys.exit(1)
    
    boilerplateRoot = os.path.join(CCP4I2_TOP, "wrappers", "boilerplate", "script")
    pluginRoot = os.path.join(CCP4I2_TOP, "wrappers", nameRoot, "script")
    nameToReplace="ZZPluginNameZZ"
    if parameterNamespace.type == "pipeline":
        boilerplateRoot = os.path.join(CCP4I2_TOP, "pipelines", "boilerplate", "script")
        pluginRoot = os.path.join(CCP4I2_TOP, "pipelines", nameRoot, "script")
        nameToReplace="ZZPipelineNameZZ"

    directoryClashed = makeScriptDirectory(pluginRoot)
    if directoryClashed and not parameterNamespace.overwrite: sys.exit(1)

    for nameModifier in [".def.xml",".py","_gui.py","_report.py"]:
        srcFilePath = os.path.join(boilerplateRoot,"boilerplate"+nameModifier)
        dstFilePath = os.path.join(pluginRoot,nameRoot+nameModifier)
        with open(dstFilePath,"w") as dstFile:
            with open(srcFilePath,"r") as srcFile:
                dstFile.write(srcFile.read().replace(nameToReplace, nameRoot)\
                              .replace ("ZZPluginMaintainerZZ", parameterNamespace.maintainer)\
                              .replace ("ZZPluginUserZZ", user)\
                              .replace ("ZZPluginCreatedZZ", timestampString)\
                              .replace ("ZZFirstPluginNameZZ", parameterNamespace.firstPlugin)\
                              .replace ("ZZShortTitleZZ", shortTitle)\
                              .replace ("ZZLongTitleZZ", longTitle)\
                              .replace ("ZZFolderZZ", parameterNamespace.where)\
                              .replace ("ZZDescriptionZZ", description)\
                              )
    with open(os.path.join(pluginRoot,"__init__.py"),"w") as pluginFile:
        pluginFile.write("")
    pluginRoot = os.path.split(pluginRoot)[0]
    with open(os.path.join(pluginRoot,"__init__.py"),"w") as pluginFile:
        pluginFile.write("")
